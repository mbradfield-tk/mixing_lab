"""
Page 4 – Mixing Sensitivity Calculation Workflow
=================================================
Select a reactor, reaction, and fluid system → compute hydrodynamic parameters,
Damköhler numbers, and assess mixing sensitivity.
"""

import streamlit as st

import pandas as pd
import numpy as np
import re
import pathlib
import plotly.graph_objects as go

from utils.calculations import (
    compute_reactor_hydro,
    compute_damkohler_numbers,
    kolmogorov_length,
    batchelor_length,
    micromixing_time_engulfment,
    blend_time_turbulent,
    epsilon_max_estimate,
    impeller_power,
    reynolds_number,
    power_number_correlation,
    mixing_sensitivity_assessment,
    settling_velocity,
    particle_reynolds,
    zwietering_njs,
    solid_liquid_mass_transfer,
    particle_suspension_criterion,
    liquid_height_from_volume,
    reaction_rate_mol_per_s,
    heat_generation_rate,
    estimate_jacket_area,
    estimate_U,
    estimate_U_detailed,
    heat_removal_capacity,
    heat_balance_assessment,
)
from utils.rom_registry import (
    compute_reactor_hydro_with_mode,
    get_correlations,
    has_any_alt_correlations,
    PARAM_DISPLAY,
    SUPPORTED_PARAMS,
)
from utils.corr_widgets import (
    render_correlation_matrix_multi,
    priority_mode_dict,
    active_modes_set,
    build_mode_dict_for,
    MODE_COLORS,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── helpers to load databases (reuse session state if available) ─────────

def _load_csv(key, filename, columns):
    if key not in st.session_state:
        p = DATA_DIR / filename
        if p.exists():
            st.session_state[key] = pd.read_csv(p)
        else:
            st.session_state[key] = pd.DataFrame(columns=columns)
    return st.session_state[key]


reactors = _load_csv("reactor_db", "reactors.csv", ["reactor_name"])
reactions = _load_csv("reaction_db", "reactions.csv", ["reaction_name"])
fluids = _load_csv("fluid_db", "fluids.csv", ["fluid_name"])
particles = _load_csv("particle_db", "particles.csv", ["particle_name"])

st.title("⚙️ Mixing Sensitivity Calculation Workflow")

if reactors.empty or reactions.empty or fluids.empty:
    st.warning("Please populate the Reactor, Reaction, and Fluid databases first.")
    st.stop()

# ── Step 1: Select system ────────────────────────────────────────────────
st.header("1 · Select System")

# Persist selections across page navigations
_reactor_list = reactors["reactor_name"].tolist()
_reaction_list = reactions["reaction_name"].tolist()
_fluid_list = fluids["fluid_name"].tolist()

def _idx(lst, key, default=0):
    """Return list index of session_state[key] if present, else default."""
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col1, col2, col3 = st.columns(3)

with col1:
    reactor_name = st.selectbox("Reactor", _reactor_list, index=_idx(_reactor_list, "_sel_reactor"), key="ms_reactor")
    st.session_state["_sel_reactor"] = reactor_name
with col2:
    reaction_name = st.selectbox("Reaction", _reaction_list, index=_idx(_reaction_list, "_sel_reaction"), key="ms_reaction")
    st.session_state["_sel_reaction"] = reaction_name
with col3:
    fluid_name = st.selectbox("Fluid", _fluid_list, index=_idx(_fluid_list, "_sel_fluid"), key="ms_fluid")
    st.session_state["_sel_fluid"] = fluid_name

reactor = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
reaction = reactions[reactions["reaction_name"] == reaction_name].iloc[0]
fluid = fluids[fluids["fluid_name"] == fluid_name].iloc[0]

# ── Helper: safely read a numeric field with a fallback ──────────────────
def _safe(series_val, default):
    """Return float(series_val) if non-NaN, else default."""
    try:
        v = float(series_val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default

# Use selection names in widget keys so they reset automatically on change
_rk = reactor_name   # reactor key fragment
_fk = fluid_name     # fluid key fragment
_xk = reaction_name  # reaction key fragment

# ── Step 2: Allow overrides ──────────────────────────────────────────────
st.divider()
st.header("2 · Review / Override Parameters")

with st.expander("Reactor geometry & agitation", expanded=False):
    oc1, oc2, oc3, oc4 = st.columns(4)
    D_tank = oc1.number_input("D_tank (m)", value=_safe(reactor.get("D_tank_m"), 0.10), format="%.4f", key=f"ov_Dt_{_rk}")
    D_imp = oc2.number_input("D_imp (m)", value=_safe(reactor.get("D_imp_m"), 0.05), format="%.4f", key=f"ov_Di_{_rk}")
    _rpm_lo = _safe(reactor.get("N_rpm_min"), 0.0)
    _rpm_hi = _safe(reactor.get("N_rpm_max"), 0.0)
    _rpm_default = (_rpm_lo + _rpm_hi) / 2 if _rpm_lo > 0 and _rpm_hi > 0 else (_rpm_hi if _rpm_hi > 0 else _safe(reactor.get("N_rps"), 5.0) * 60)
    _N_rpm_input = oc3.number_input("N (RPM)", value=_rpm_default, step=2.0, format="%.0f", key=f"ov_N_{_rk}")
    N_rps = _N_rpm_input / 60.0
    oc5, oc6 = st.columns(2)
    Np_in = oc5.number_input("Np", value=_safe(reactor.get("Np"), 1.27), format="%.2f", key=f"ov_Np_{_rk}")
    Nq_in = oc6.number_input("Nq", value=_safe(reactor.get("Nq"), 0.79), format="%.2f", key=f"ov_Nq_{_rk}")

# ── Volume selection ─────────────────────────────────────────────────────
V_L_min = _safe(reactor.get("V_L_min"), 0.0)
V_L_max = _safe(reactor.get("V_L_max"), 0.0)
V_L_nom = _safe(reactor.get("V_L"), 0.0)
H_max = _safe(reactor.get("H_max_m"), _safe(reactor.get("H_m"), 0.13))
_bottom_dish = str(reactor.get("bottom_dish", "")) if pd.notna(reactor.get("bottom_dish")) else ""

# Compute a sensible default volume
if V_L_min > 0 and V_L_max > 0:
    V_L_default = (V_L_min + V_L_max) / 2.0
elif V_L_nom > 0:
    V_L_default = V_L_nom
else:
    # Fall back to geometric estimate
    V_L_default = np.pi / 4 * D_tank**2 * _safe(reactor.get("H_m"), 0.13) * 1000

V_L = st.number_input(
    "Liquid volume (L)", value=V_L_default, min_value=0.01,
    format="%.2f", key=f"ov_VL_{_rk}",
    help=f"Reactor range: {V_L_min:.1f} – {V_L_max:.1f} L" if V_L_max > 0 else None,
)

# Validate volume is within the reactor's operating range
_vol_ok = True
if V_L_min > 0 and V_L_max > 0:
    if V_L < V_L_min or V_L > V_L_max:
        st.error(f"Volume {V_L:.1f} L is outside the reactor operating range "
                 f"({V_L_min:.1f} – {V_L_max:.1f} L).")
        _vol_ok = False

# Derive liquid height from volume, accounting for bottom dish geometry
H = liquid_height_from_volume(V_L, D_tank, H_max, _bottom_dish)

with st.expander("Fluid properties", expanded=False):
    fc1, fc2, fc3 = st.columns(3)
    rho = fc1.number_input("ρ (kg/m³)", value=float(fluid["rho_kg_m3"]), format="%.1f", key=f"ov_rho_{_fk}")
    mu = fc2.number_input("μ (Pa·s)", value=float(fluid["mu_Pa_s"]), format="%.6f", key=f"ov_mu_{_fk}")
    D_mol = fc3.number_input("D_mol (m²/s)", value=float(fluid["D_mol_m2_s"]), format="%.2e", key=f"ov_Dmol_{_fk}")

# ── Particle (solid) phase ────────────────────────────────────────────────
include_particles = st.checkbox("Include solid particles", value=False, key="include_particles")
if include_particles and not particles.empty:
    with st.expander("Particle properties", expanded=True):
        def _on_particle_change():
            """Reset override widgets when a new particle is selected."""
            pname = st.session_state["sel_particle"]
            p = particles[particles["particle_name"] == pname].iloc[0]
            st.session_state["ov_rhop"] = float(p["rho_p_kg_m3"])
            st.session_state["ov_d50"] = float(p["d50_um"])
            st.session_state["ov_phi"] = float(p["shape_factor"])

        particle_name = st.selectbox(
            "Particle", particles["particle_name"].tolist(),
            key="sel_particle", on_change=_on_particle_change,
        )
        particle = particles[particles["particle_name"] == particle_name].iloc[0]
        # Initialise session state on very first render
        if "ov_rhop" not in st.session_state:
            st.session_state["ov_rhop"] = float(particle["rho_p_kg_m3"])
        if "ov_d50" not in st.session_state:
            st.session_state["ov_d50"] = float(particle["d50_um"])
        if "ov_phi" not in st.session_state:
            st.session_state["ov_phi"] = float(particle["shape_factor"])

        pc1, pc2, pc3 = st.columns(3)
        rho_p = pc1.number_input("ρ_p (kg/m³)", min_value=100.0, max_value=25000.0,
                                 step=10.0, format="%.0f", key="ov_rhop")
        d50_um = pc2.number_input("D50 (µm)", min_value=0.1, max_value=10000.0,
                                  step=1.0, format="%.1f", key="ov_d50")
        phi_p = pc3.number_input("Shape factor φ",
                                 min_value=0.01, max_value=1.0, step=0.05,
                                 format="%.2f", key="ov_phi")
        pc4, pc5 = st.columns(2)
        X_wt = pc4.number_input("Solids loading X (wt-%)", value=5.0, min_value=0.01,
                                format="%.2f", key="ov_Xwt",
                                help="Mass fraction of solids as weight-%")
        S_zw = pc5.number_input("Zwietering S constant", value=5.5, min_value=0.5,
                                max_value=20.0, format="%.1f", key="ov_Szw",
                                help="Geometry-dependent constant (typ. 1–10, PBT ≈ 5.5)")
elif include_particles and particles.empty:
    st.warning("Particle database is empty. Add particles on the Particle Database page.")

with st.expander("Gas sparging (for kLa)", expanded=False):
    gc1, gc2 = st.columns(2)
    v_s = gc1.number_input("Superficial gas velocity v_s (m/s)", value=0.0, min_value=0.0, format="%.4f",
                           key=f"ov_vs_{_rk}", help="Set > 0 to compute kLa. Typical: 0.001 – 0.05 m/s.")
    coalescing = gc2.selectbox("Liquid type", ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                               key=f"ov_coal_{_fk}")
    is_coalescing = coalescing.startswith("Coalescing")

with st.expander("Reaction parameters", expanded=False):
    rc1, rc2, rc3 = st.columns(3)
    k_val = rc1.number_input("k", value=float(reaction["k_value"]), format="%.6g", key=f"ov_k_{_xk}")
    C0 = rc2.number_input("C₀ (mol/L)", value=float(reaction["C0_mol_L"]), format="%.4g", key=f"ov_C0_{_xk}")
    t_rxn_input = rc3.number_input("t_rxn (s)", value=float(reaction["t_rxn_s"]), format="%.4g", key=f"ov_trxn_{_xk}",
                                    help="Characteristic reaction time. 0 = auto-compute.")

# ── Step 3: Correlations ─────────────────────────────────────────────
st.divider()
st.header("3 · Correlation Source Selection")

if has_any_alt_correlations(reactor_name):
    with st.expander("Correlation source selection", expanded=False):
        st.caption(
            "Select one or more of **Literature**, **ROM**, **Experimental** per "
            "parameter.  Each selection generates its own operating envelope.  "
            "Metrics show the highest-priority value "
            "(Experimental → ROM → Literature)."
        )
        corr_selections = render_correlation_matrix_multi(
            reactor_name, key_prefix="ms_corr",
        )
else:
    st.info(f"Only **Literature** correlations available for **{reactor_name}**. "
            "Register ROM or Experimental correlations to enable source comparison.")
    corr_selections = {p: ["Literature"] for p in SUPPORTED_PARAMS}

# Derive a single priority mode dict for the main metrics
_priority_mode = priority_mode_dict(corr_selections)
# Collect all active modes for envelope computation
_active_modes = active_modes_set(corr_selections)


# ── Heat balance options ──────────────────────────────────────────────────────
st.divider()
st.header("4 · Heat Balance Options")
def _parse_fluid_temp(name: str, fallback: float = 25.0) -> float:
    m = re.search(r'\(([-\d.]+)\s*°?C\)', name)
    return float(m.group(1)) if m else fallback

fluid_T_C = _parse_fluid_temp(fluid_name)
rxn_T_C = _safe(reaction.get("T_C"), 25.0)
rxn_delta_H = _safe(reaction.get("delta_H_kJ_mol"), 0.0)
rxn_order = str(reaction.get("order", "1"))

_ms_heat_context = f"{fluid_name}|{reaction_name}"
if st.session_state.get("_ms_heat_context") != _ms_heat_context:
    if st.session_state.get("_ms_heat_active", False):
        st.session_state["_ms_heat_active"] = False
        st.session_state["_ms_heat_stale"] = True
    st.session_state["_ms_heat_context"] = _ms_heat_context

if st.button("🔥 Run Heat Balance"):
    st.session_state["_ms_heat_active"] = True
    st.session_state.pop("_ms_heat_stale", None)
    _default_T = rxn_T_C if rxn_T_C > 0 else fluid_T_C
    st.session_state["ms_T_process"] = _default_T
    st.session_state["ms_T_cool"] = _default_T - 20.0
    st.rerun()

if st.session_state.get("_ms_heat_stale"):
    st.info("Fluid or reaction changed — click **🔥 Run Heat Balance** to update.")

include_heat = st.session_state.get("_ms_heat_active", False)
ms_T_process = fluid_T_C
ms_T_coolant = fluid_T_C - 20.0
if include_heat:
    hcol1, hcol2, hcol3 = st.columns(3)
    with hcol1:
        ms_T_process = st.number_input(
            "Process temperature (°C)",
            format="%.1f", key="ms_T_process",
            help="Defaults to reaction temperature; adjust as needed",
            step=1.0)
    with hcol2:
        ms_T_coolant = st.number_input(
            "Coolant temperature (°C)",
            format="%.1f", key="ms_T_cool",
            help="Jacket coolant inlet temperature",
            step=1.0)
    with hcol3:
        st.markdown(f"**ΔH** = {rxn_delta_H:.0f} kJ/mol  ")
        st.markdown(f"**ΔT** = {ms_T_process - ms_T_coolant:.1f} °C")
    if rxn_delta_H == 0:
        st.info("Selected reaction has no enthalpy value – add ΔH in the Reaction Database.")

# ── Step 5: Compute ──────────────────────────────────────────────────────
st.divider()
st.header("5 · Centerpoint Results")

# Compute reaction time if needed
t_rxn = t_rxn_input
if t_rxn <= 0 and k_val > 0:
    order = str(reaction.get("order", "1"))
    if order in ("1", "pseudo-1"):
        t_rxn = np.log(2) / k_val
    elif order in ("2", "pseudo-2") and C0 > 0:
        t_rxn = 1.0 / (k_val * C0)
    else:
        t_rxn = 1.0  # fallback

hydro, rom_sources = compute_reactor_hydro_with_mode(
    _priority_mode, reactor_name,
    N=N_rps, D_imp=D_imp, D_tank=D_tank, H=H,
    rho=rho, mu=mu, Np=Np_in, Nq=Nq_in,
    v_s=v_s, coalescing=is_coalescing,
    D_mol=D_mol,
)

da = compute_damkohler_numbers(
    t_blend=hydro["Blend time 95% (s)"],
    t_micro=hydro["Micromix time t_E (s)"],
    t_rxn=t_rxn,
    kLa=hydro["kLa (1/s)"],
    kLa_surface=hydro["kLa_surface (1/s)"],
)

# Show which parameters are using ROM / Experimental correlations
if rom_sources:
    _src_lines = [f"- **{PARAM_DISPLAY.get(k, k)}**: {v}" for k, v in rom_sources.items()]
    _modes_label = ", ".join(m for m in _active_modes if m != "Literature") or "Literature"
    st.info(
        f"**Metrics show highest-priority source** ({_modes_label}) for {reactor_name}:\n"
        + "\n".join(_src_lines)
    )

# ── Display ──────────────────────────────────────────────────────────────

# Key metrics row
m1, m2, m3, m4 = st.columns(4)
m1.metric("Re", f"{hydro['Re']:.0f}")
m2.metric("P/V (W/L)", f"{hydro['P/V (W/L)']:.3g}")
m3.metric("Blend time (s)", f"{hydro['Blend time 95% (s)']:.2f}")
m4.metric("Micromix t_E (s)", f"{hydro['Micromix time t_E (s)']:.4g}")

m5, m6, m7, m8 = st.columns(4)
m5.metric("Tip speed (m/s)", f"{hydro['Tip speed (m/s)']:.2f}")
m6.metric("Kolmogorov η (µm)", f"{hydro['Kolmogorov η (µm)']:.1f}")
m7.metric("Da_macro", f"{da['Da_macro']:.3g}")
m8.metric("Da_micro", f"{da['Da_micro']:.3g}")

# New parameters row
m9, m10, m11, m12 = st.columns(4)
m9.metric("Avg shear rate (1/s)", f"{hydro['Avg shear rate (1/s)']:.1f}")
m10.metric("Max shear rate (1/s)", f"{hydro['Max shear rate (1/s)']:.0f}")
m11.metric("Avg shear stress (Pa)", f"{hydro['Avg shear stress (Pa)']:.3g}")
m12.metric("t_E local (s)", f"{hydro['Micromix time t_E_local (s)']:.4g}")

if v_s > 0:
    m13, m14, m15, _ = st.columns(4)
    m13.metric("kLa sparged (1/s)", f"{hydro['kLa (1/s)']:.4g}")
    m14.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
    m15.metric("Da_GL", f"{da['Da_GL']:.3g}")
else:
    m13, m14, _, _ = st.columns(4)
    m13.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
    da_gl_val = da['Da_GL']
    if da_gl_val > 0:
        m14.metric("Da_GL (surface)", f"{da_gl_val:.3g}")
    st.caption("Sparged kLa not shown (set v_s > 0 in *Gas sparging* section). "
               "Surface kLa uses the Lamont-Scott free-surface model.")

# Assessment banners – one per sensitivity type
def _da_banner(label: str, Da: float, thresholds: dict[str, tuple[float, float]]) -> None:
    """Display a coloured banner for a single Damköhler-based assessment."""
    if Da < 0.01:
        st.success(f"🟢 **{label}:** Not sensitive (Da = {Da:.3g})")
    elif Da < 0.1:
        st.success(f"🟢 **{label}:** Likely not sensitive (Da = {Da:.3g})")
    elif Da < 1:
        st.warning(f"🟡 **{label}:** Potentially sensitive (Da = {Da:.3g})")
    elif Da < 10:
        st.error(f"🔴 **{label}:** Likely sensitive (Da = {Da:.3g})")
    else:
        st.error(f"🔴 **{label}:** Highly sensitive (Da = {Da:.3g})")

_da_banner("Macromixing", da["Da_macro"], {})
_da_banner("Micromixing", da["Da_micro"], {})
if da["Da_GL"] > 0:
    _da_banner("Gas-liquid mass transfer", da["Da_GL"], {})

# Batchelor length
lam_B = batchelor_length(mu / rho, hydro["P/V (W/kg)"], D_mol)
st.info(f"Batchelor microscale λ_B = {lam_B * 1e6:.2f} µm  |  Reaction time t_rxn = {t_rxn:.4g} s")

# ── Particle (solid-phase) results ──────────────────────────────────────
particle_results = {}
particle_meta = {}
if include_particles and not particles.empty:
    st.subheader("Solid-Phase Parameters")
    d_p_m = d50_um * 1e-6  # convert µm → m
    nu = mu / rho
    delta_rho = abs(rho_p - rho)

    v_t = settling_velocity(d_p_m, rho_p, rho, mu, phi_p)
    Re_p = particle_reynolds(d_p_m, v_t, rho, mu)
    N_js = zwietering_njs(S_zw, nu, d_p_m, delta_rho, rho, X_wt, D_imp)
    k_SL = solid_liquid_mass_transfer(d_p_m, v_t, rho, mu, D_mol)
    susp = particle_suspension_criterion(N_rps, N_js)

    particle_results = {
        "d50 (µm)": d50_um,
        "ρ_p (kg/m³)": rho_p,
        "Shape factor φ": phi_p,
        "X (wt-%)": X_wt,
        "v_t (m/s)": v_t,
        "Re_p": Re_p,
        "N_js (rev/s)": N_js,
        "N_js (RPM)": N_js * 60,
        "k_SL (m/s)": k_SL,
    }
    particle_meta = {
        "Particle": particle_name,
        "Suspension": susp,
    }

    sp1, sp2, sp3, sp4 = st.columns(4)
    sp1.metric("Settling velocity (m/s)", f"{v_t:.3e}")
    sp2.metric("Re_p", f"{Re_p:.3g}")
    sp3.metric("N_js (RPM)", f"{N_js * 60:.1f}")
    sp4.metric("k_SL (m/s)", f"{k_SL:.3e}")

    # Suspension assessment
    if "Poorly" in susp:
        st.error(f"🔴 **{susp}** — current speed is below just-suspended speed")
    elif "Partially" in susp:
        st.warning(f"🟡 **{susp}** — increase speed to achieve full suspension")
    elif "Just" in susp:
        st.info(f"🟢 **{susp}** — operating near just-suspended conditions")
    else:
        st.success(f"🟢 **{susp}** — particles are well suspended")

# ── Heat Balance Results ─────────────────────────────────────────────────
heat_results = {}
if include_heat and rxn_delta_H != 0:
    st.subheader("🌡️ Heat Balance")
    _r_material = str(reactor.get("material", ""))
    _r_U_override = _safe(reactor.get("U_W_m2K"), 0.0)
    _r_A_override = _safe(reactor.get("A_ht_m2"), 0.0)
    _r_wall_mm = _safe(reactor.get("wall_thickness_mm"), 0.0)

    _r_mol_s = reaction_rate_mol_per_s(rxn_order, k_val, C0, V_L)
    _Q_gen = heat_generation_rate(rxn_delta_H, _r_mol_s)
    _A_ht = _r_A_override if _r_A_override > 0 else estimate_jacket_area(D_tank, H, _bottom_dish)
    if _r_U_override > 0:
        _U_ht = _r_U_override
        _u_warns = []
    else:
        _U_ht, _u_warns = estimate_U_detailed(
            N_rps=N_rps, D_imp=D_imp, D_tank=D_tank,
            rho=rho, mu=mu,
            material=_r_material,
            wall_thickness_mm=_r_wall_mm,
            fluid_name=fluid_name,
        )
    _dT = ms_T_process - ms_T_coolant
    _Q_cool = heat_removal_capacity(_U_ht, _A_ht, _dT)
    _ratio_pct = _Q_gen / _Q_cool * 100 if _Q_cool > 0 else np.inf

    hm1, hm2, hm3, hm4 = st.columns(4)
    hm1.metric("Q_gen (W)", f"{_Q_gen:.1f}")
    hm2.metric("Q_cool (W)", f"{_Q_cool:.1f}")
    hm3.metric("U (W/m²·K)", f"{_U_ht:.0f}")
    hm4.metric("A_ht (m²)", f"{_A_ht:.3f}")

    _assessment = heat_balance_assessment(_Q_gen, _Q_cool)
    if _ratio_pct < 100:
        st.success(f"🟢 Q_gen/Q_cool = {_ratio_pct:.1f}% — **{_assessment}**")
    elif _ratio_pct < 10000:
        st.error(f"🔴 Q_gen/Q_cool = {_ratio_pct:.1f}% — **{_assessment}**")
    else:
        st.error(f"🔴 Q_gen/Q_cool = ∞ — **{_assessment}**")

    if _u_warns:
        with st.expander("ℹ️ U estimation notes", expanded=False):
            for _w in _u_warns:
                st.markdown(f"- {_w}")

    heat_results = {
        "Q_gen (W)": _Q_gen,
        "Q_cool (W)": _Q_cool,
        "U (W/m²·K)": _U_ht,
        "A_ht (m²)": _A_ht,
        "Q_gen/Q_cool (%)": _ratio_pct,
    }

# ── Operating Envelope Charts ────────────────────────────────────────────
st.divider()
st.header("6. Operating Envelopes")
st.caption("Parameter variation across the reactor's full RPM and volume range.")

# Read RPM range
_rpm_min = _safe(reactor.get("N_rpm_min"), 0.0)
_rpm_max = _safe(reactor.get("N_rpm_max"), 0.0)
_n_rps_default = _safe(reactor.get("N_rps"), 0.0)
if _rpm_max == 0 and _n_rps_default > 0:
    _rpm_max = _n_rps_default * 60
    _rpm_min = _rpm_max
if _rpm_min == 0:
    _rpm_min = _rpm_max

_N_lo = _rpm_min / 60.0 if _rpm_min > 0 else N_rps
_N_hi = _rpm_max / 60.0 if _rpm_max > 0 else N_rps

_env_V_max = V_L_max if V_L_max > 0 else V_L
_env_V_min = V_L_min if V_L_min > 0 else V_L

_can_envelope = _N_lo > 0 and _N_hi > 0 and _env_V_max > 0

if _can_envelope:
    # Reactor heat-transfer metadata
    _r_material_env = str(reactor.get("material", ""))
    _r_U_override_env = _safe(reactor.get("U_W_m2K"), 0.0)
    _r_A_override_env = _safe(reactor.get("A_ht_m2"), 0.0)
    _r_wall_mm_env = _safe(reactor.get("wall_thickness_mm"), 0.0)

    PLOT_PARAMS = [
        "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
        "Micromix time t_E (s)", "Micromix time t_E_local (s)",
        "Kolmogorov η (µm)", "Re",
        "Avg shear rate (1/s)", "Max shear rate (1/s)", "Avg shear stress (Pa)",
        "Da_macro", "Da_micro", "Da_GL", "ε_max (W/kg)",
        "kLa (1/s)", "kLa_surface (1/s)",
    ]
    HEAT_PARAMS = ["Q_gen (W)", "Q_cool (W)", "U (W/m²·K)", "A_ht (m²)", "Q_gen/Q_cool (%)"]
    if include_heat and rxn_delta_H != 0:
        PLOT_PARAMS = PLOT_PARAMS + HEAT_PARAMS
    PARTICLE_PARAMS = ["N_js (RPM)", "N/N_js", "v_t (m/s)", "Re_p", "k_SL (m/s)"]
    if include_particles and not particles.empty:
        PLOT_PARAMS = PLOT_PARAMS + PARTICLE_PARAMS

    _N_INTERP = 50
    N_arr = np.linspace(_N_lo, _N_hi, _N_INTERP)
    pct_arr = N_arr / _N_hi * 100 if _N_hi > 0 else np.zeros(_N_INTERP)

    # ── Pre-compute RPM-independent particle quantities (hoisted) ────────
    _part_static: dict | None = None
    if include_particles and not particles.empty:
        _dp = d50_um * 1e-6
        _nu_p = mu / rho
        _drho = abs(rho_p - rho)
        _vt = settling_velocity(_dp, rho_p, rho, mu, phi_p)
        _rep = particle_reynolds(_dp, _vt, rho, mu)
        _njs = zwietering_njs(S_zw, _nu_p, _dp, _drho, rho, X_wt, D_imp)
        _ksl = solid_liquid_mass_transfer(_dp, _vt, rho, mu, D_mol)
        _part_static = {
            "N_js_rps": _njs,
            "v_t (m/s)": _vt,
            "Re_p": _rep,
            "k_SL (m/s)": _ksl,
            "N_js (RPM)": _njs * 60,
        }

    # ── Deduplicate mode dicts so identical sweeps run only once ──────────
    _mode_dicts: dict[str, dict[str, str]] = {}
    _mode_label_to_key: dict[str, str] = {}  # mode_label → canonical key
    _seen_tuples: dict[tuple, str] = {}       # frozen mode_dict → first label
    for _mode_label in _active_modes:
        _md = build_mode_dict_for(_mode_label, corr_selections)
        _frozen = tuple(sorted(_md.items()))
        if _frozen in _seen_tuples:
            _mode_label_to_key[_mode_label] = _seen_tuples[_frozen]
        else:
            _seen_tuples[_frozen] = _mode_label
            _mode_label_to_key[_mode_label] = _mode_label
            _mode_dicts[_mode_label] = _md

    # ── Sweep only unique mode dicts ─────────────────────────────────────
    # curve_data: mode_label → vol_key → {param: np.array}
    _unique_curve_data: dict[str, dict[str, dict[str, np.ndarray]]] = {}

    for _mode_label, _mode_dict in _mode_dicts.items():
        mode_curves: dict[str, dict[str, np.ndarray]] = {}
        for vol_key, _vl in [("maxV", _env_V_max), ("minV", _env_V_min)]:
            H_v = liquid_height_from_volume(_vl, D_tank, H_max, _bottom_dish)
            param_arrs = {p: np.empty(_N_INTERP) for p in PLOT_PARAMS}

            # Pre-compute RPM-independent heat quantities for this volume
            _Q_gen_v = _A_ht_v = 0.0
            if include_heat and rxn_delta_H != 0:
                _r_mol_s_e = reaction_rate_mol_per_s(rxn_order, k_val, C0, _vl)
                _Q_gen_v = heat_generation_rate(rxn_delta_H, _r_mol_s_e)
                _A_ht_v = _r_A_override_env if _r_A_override_env > 0 else estimate_jacket_area(D_tank, H_v, _bottom_dish)

            for j, _N in enumerate(N_arr):
                _h, _ = compute_reactor_hydro_with_mode(
                    _mode_dict, reactor_name,
                    N=_N, D_imp=D_imp, D_tank=D_tank, H=H_v,
                    rho=rho, mu=mu, Np=Np_in, Nq=Nq_in,
                    v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                )
                _da = compute_damkohler_numbers(
                    _h["Blend time 95% (s)"], _h["Micromix time t_E (s)"], t_rxn,
                    kLa=_h["kLa (1/s)"], kLa_surface=_h["kLa_surface (1/s)"],
                )
                _vals = {**_h, **_da}
                if _part_static is not None:
                    _vals["N_js (RPM)"] = _part_static["N_js (RPM)"]
                    _vals["N/N_js"] = _N / _part_static["N_js_rps"] if _part_static["N_js_rps"] > 0 else 0.0
                    _vals["v_t (m/s)"] = _part_static["v_t (m/s)"]
                    _vals["Re_p"] = _part_static["Re_p"]
                    _vals["k_SL (m/s)"] = _part_static["k_SL (m/s)"]
                if include_heat and rxn_delta_H != 0:
                    if _r_U_override_env > 0:
                        _U_ht_e = _r_U_override_env
                    else:
                        _U_ht_e, _ = estimate_U_detailed(
                            N_rps=_N, D_imp=D_imp, D_tank=D_tank,
                            rho=rho, mu=mu,
                            material=_r_material_env,
                            wall_thickness_mm=_r_wall_mm_env,
                            fluid_name=fluid_name,
                        )
                    _dT_e = ms_T_process - ms_T_coolant
                    _Q_cool_e = heat_removal_capacity(_U_ht_e, _A_ht_v, _dT_e)
                    _vals["Q_gen (W)"] = _Q_gen_v
                    _vals["Q_cool (W)"] = _Q_cool_e
                    _vals["U (W/m²·K)"] = _U_ht_e
                    _vals["A_ht (m²)"] = _A_ht_v
                    _vals["Q_gen/Q_cool (%)"] = _Q_gen_v / _Q_cool_e * 100 if _Q_cool_e > 0 else np.inf
                for p in PLOT_PARAMS:
                    param_arrs[p][j] = _vals.get(p, np.nan)
            mode_curves[vol_key] = param_arrs
        _unique_curve_data[_mode_label] = mode_curves

    # Map all mode labels (including duplicates) to their computed data
    curve_data: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for _mode_label in _active_modes:
        curve_data[_mode_label] = _unique_curve_data[_mode_label_to_key[_mode_label]]

    # Mark current operating point as RPM %
    _current_pct = N_rps / _N_hi * 100 if _N_hi > 0 else 50.0

    # Priority mode label for the ★ operating-point marker
    # (Experimental > ROM > Literature among active modes)
    _priority_mode_label = "Literature"
    for _ml in reversed(_active_modes):  # _active_modes is ordered Lit, ROM, Exp
        if _ml in curve_data:
            _priority_mode_label = _ml
            break

    _DISPLAY_NAMES: dict[str, str] = {
        "Da_macro": "Macromixing (Da_macro)",
        "Da_micro": "Micromixing (Da_micro)",
        "Da_GL": "Gas-Liquid Mass Transfer (Da_GL)",
        "Q_gen/Q_cool (%)": "Heat Transfer Capacity (Q_gen/Q_cool (%))",
    }
    _display = lambda p: _DISPLAY_NAMES.get(p, p)

    _DEFAULT_PARAMS = ["Da_micro", "Da_macro", "Da_GL", "Q_gen/Q_cool (%)", "Blend time 95% (s)", "P/V (W/L)"]
    _defaults = [p for p in _DEFAULT_PARAMS if p in PLOT_PARAMS]

    with st.expander("Show / hide envelope charts", expanded=True):
        params_to_plot = st.multiselect(
            "Parameters to plot",
            PLOT_PARAMS,
            default=_defaults,
            key="ms_env_params",
            format_func=_display,
        )

        for param in params_to_plot:
            fig = go.Figure()

            # One envelope per active mode
            for _mode_label in _active_modes:
                _COLOR = MODE_COLORS.get(_mode_label, "#999999")
                _mc = curve_data[_mode_label]
                y_maxV = _mc["maxV"][param]
                y_minV = _mc["minV"][param]

                # Filled polygon
                poly_x = np.concatenate([pct_arr, pct_arr[::-1], [pct_arr[0]]])
                poly_y = np.concatenate([y_maxV, y_minV[::-1], [y_maxV[0]]])
                fig.add_trace(go.Scatter(
                    x=poly_x, y=poly_y,
                    fill="toself", fillcolor=_COLOR, opacity=0.15,
                    line=dict(color=_COLOR, width=1), mode="lines",
                    name=f"{_mode_label} envelope",
                    legendgroup=_mode_label, showlegend=True,
                    hoverinfo="skip",
                ))
                # Max-volume boundary (solid)
                fig.add_trace(go.Scatter(
                    x=pct_arr, y=y_maxV,
                    mode="lines", line=dict(color=_COLOR, width=2),
                    name=f"{_mode_label} max V ({_env_V_max:.1f} L)",
                    legendgroup=_mode_label, showlegend=True,
                    hovertemplate=f"{_mode_label} | %% max RPM: %{{x:.1f}}%%<br>%{{y:.3g}}<extra>max V</extra>",
                ))
                # Min-volume boundary (dotted)
                fig.add_trace(go.Scatter(
                    x=pct_arr, y=y_minV,
                    mode="lines", line=dict(color=_COLOR, width=2, dash="dot"),
                    name=f"{_mode_label} min V ({_env_V_min:.1f} L)",
                    legendgroup=_mode_label, showlegend=True,
                    hovertemplate=f"{_mode_label} | %% max RPM: %{{x:.1f}}%%<br>%{{y:.3g}}<extra>min V</extra>",
                ))

            # Current operating point (uses priority mode)
            _priority_curves = curve_data[_priority_mode_label]
            _y_maxV_p = _priority_curves["maxV"][param]
            _y_minV_p = _priority_curves["minV"][param]
            if abs(_env_V_max - _env_V_min) > 1e-6:
                _frac = (V_L - _env_V_min) / (_env_V_max - _env_V_min)
                _frac = max(0.0, min(1.0, _frac))
                _y_interp = (np.interp(_current_pct, pct_arr, _y_minV_p) * (1 - _frac)
                             + np.interp(_current_pct, pct_arr, _y_maxV_p) * _frac)
            else:
                _y_interp = np.interp(_current_pct, pct_arr, _y_maxV_p)
            fig.add_trace(go.Scatter(
                x=[_current_pct], y=[_y_interp],
                mode="markers", marker=dict(size=12, color="red", symbol="star",
                                             line=dict(width=1, color="white")),
                name="Current", showlegend=True,
                hovertemplate=f"Current: {N_rps*60:.0f} RPM, {V_L:.1f} L<br>" + param + " = %{y:.3g}<extra></extra>",
            ))

            # Reference lines for Da parameters
            if param in ("Da_macro", "Da_micro", "Da_GL"):
                import math
                for da_val, da_color, label in [
                    (0.1, "orange", "Da=0.1 (onset of sensitivity)"),
                    (1.0, "red", "Da=1 (limited)"),
                ]:
                    fig.add_shape(
                        type="line", x0=0, x1=1, y0=da_val, y1=da_val,
                        xref="paper", yref="y",
                        line=dict(color=da_color, width=1.5, dash="dash"),
                    )
                    fig.add_annotation(
                        x=1.0, xref="paper", xanchor="right",
                        y=math.log10(da_val), yref="y",
                        yanchor="top", yshift=-2,
                        text=label, font=dict(size=11, color=da_color),
                        showarrow=False,
                    )
                fig.update_yaxes(type="log")

            if param == "N/N_js":
                fig.add_shape(
                    type="line", x0=0, x1=1, y0=1.0, y1=1.0,
                    xref="paper", yref="y",
                    line=dict(color="red", width=1.5, dash="dash"),
                )
                fig.add_annotation(
                    x=1.0, xref="paper", xanchor="right",
                    y=1.0, yref="y", yanchor="bottom", yshift=2,
                    text="N/N_js = 1 (just suspended)",
                    font=dict(size=11, color="red"), showarrow=False,
                )

            if param == "Q_gen/Q_cool (%)":
                fig.add_shape(
                    type="line", x0=0, x1=1, y0=100.0, y1=100.0,
                    xref="paper", yref="y",
                    line=dict(color="red", width=1.5, dash="dash"),
                )
                fig.add_annotation(
                    x=1.0, xref="paper", xanchor="right",
                    y=100.0, yref="y", yanchor="bottom", yshift=2,
                    text="100% (cooling limit)",
                    font=dict(size=11, color="red"), showarrow=False,
                )

            _param_label = _display(param)
            _yaxis_opts: dict = dict(
                showspikes=True, spikemode="across",
                spikethickness=1, spikecolor="grey", spikedash="dot",
            )
            if param == "Q_gen/Q_cool (%)":
                _yaxis_opts["ticksuffix"] = "%"
            fig.update_layout(
                title=_param_label,
                xaxis_title=f"Stir speed (% of max RPM = {_rpm_max:.0f})",
                yaxis_title=_param_label,
                xaxis=dict(range=[0, 105], dtick=10,
                           showspikes=True, spikemode="across",
                           spikethickness=1, spikecolor="grey", spikedash="dot"),
                yaxis=_yaxis_opts,
                height=500,
                margin=dict(t=50, b=50),
                hovermode="closest",
            )
            st.plotly_chart(fig, width="content")

        with st.expander("Chart legend"):
            st.markdown("""
**Shaded region:** Operational envelope spanning min-to-max fill volume.

**Boundary lines:**
- **Solid line** = max fill volume edge
- **Dotted line** = min fill volume edge

**★ Red star** = Current operating point (selected RPM & volume)

**Colours (when multiple sources selected):**
- 🔵 Blue = Literature
- 🟢 Green = ROM
- 🟠 Orange = Experimental
""")
else:
    st.info("Reactor has no RPM range defined — cannot compute operating envelope.")

# Full table
st.subheader("Full Hydrodynamic Parameter Table")
hydro_df = pd.DataFrame([hydro]).T
hydro_df.columns = ["Value"]
if particle_results:
    part_df = pd.DataFrame([particle_results]).T
    part_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, part_df])
if heat_results:
    heat_df = pd.DataFrame([heat_results]).T
    heat_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, heat_df])
st.dataframe(hydro_df, width="stretch")
if particle_results:
    st.caption(f"Particle: **{particle_meta['Particle']}**  ·  {particle_meta['Suspension']}")

# ── Step 4: Save to recorded results ────────────────────────────────────
st.header("4 · Save Result")

if st.button("📌 Save this result to Recorded Results"):
    result_row = {
        "reactor": reactor_name,
        "reaction": reaction_name,
        "fluid": fluid_name,
        "Re": hydro["Re"],
        "P/V (W/L)": hydro["P/V (W/L)"],
        "Tip speed (m/s)": hydro["Tip speed (m/s)"],
        "Blend time (s)": hydro["Blend time 95% (s)"],
        "Micromix t_E (s)": hydro["Micromix time t_E (s)"],
        "Micromix t_E_local (s)": hydro["Micromix time t_E_local (s)"],
        "Kolmogorov η (µm)": hydro["Kolmogorov η (µm)"],
        "Batchelor λ_B (µm)": lam_B * 1e6,
        "Avg shear rate (1/s)": hydro["Avg shear rate (1/s)"],
        "Max shear rate (1/s)": hydro["Max shear rate (1/s)"],
        "Avg shear stress (Pa)": hydro["Avg shear stress (Pa)"],
        "kLa (1/s)": hydro["kLa (1/s)"],
        "kLa_surface (1/s)": hydro["kLa_surface (1/s)"],
        "t_rxn (s)": t_rxn,
        "Da_macro": da["Da_macro"],
        "Da_micro": da["Da_micro"],
        "Da_GL": da["Da_GL"],
        "Assessment": da["Assessment"],
    }
    if particle_results:
        result_row.update({
            "Particle": particle_meta.get("Particle", ""),
            "d50 (µm)": particle_results.get("d50 (µm)", ""),
            "ρ_p (kg/m³)": particle_results.get("ρ_p (kg/m³)", ""),
            "v_t (m/s)": particle_results.get("v_t (m/s)", ""),
            "Re_p": particle_results.get("Re_p", ""),
            "N_js (RPM)": particle_results.get("N_js (RPM)", ""),
            "k_SL (m/s)": particle_results.get("k_SL (m/s)", ""),
            "Suspension": particle_meta.get("Suspension", ""),
        })
    if "recorded_results" not in st.session_state:
        st.session_state.recorded_results = pd.DataFrame()
    st.session_state.recorded_results = pd.concat(
        [st.session_state.recorded_results, pd.DataFrame([result_row])],
        ignore_index=True,
    )
    # Also persist to CSV
    results_csv = DATA_DIR / "recorded_results.csv"
    st.session_state.recorded_results.to_csv(results_csv, index=False)
    st.success("Result saved! View it on the **Recorded Results** page.")
