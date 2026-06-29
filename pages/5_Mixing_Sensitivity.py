"""
Page 4 – Mixing Sensitivity Calculation Workflow
=================================================
Select a reactor, reaction, and fluid system → compute hydrodynamic parameters,
Damköhler numbers, and assess mixing sensitivity.
"""

import streamlit as st

import pandas as pd
import numpy as np
import pathlib
import plotly.graph_objects as go

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, list_solvents, is_known_solvent,
    hansen_distance, miscibility_assessment, get_hsp, solvent_miscibility,
)
from utils.calculations import (
    compute_reactor_hydro,
    compute_damkohler_numbers,
    damkohler_sl,
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
    solid_liquid_kla,
    archimedes_number,
    particle_suspension_criterion,
    liquid_height_from_volume,
    reaction_rate_mol_per_s,
    heat_generation_rate,
    estimate_jacket_area,
    estimate_U,
    estimate_U_detailed,
    heat_removal_capacity,
    heat_balance_assessment,
    mesomixing_time,
    cooling_rate,
    time_to_cool_or_heat,
    gmb_njs,
    gas_holdup_hughmark,
    sauter_bubble_diameter,
    gas_flooding_speed,
    gas_flow_rate_from_vs,
    weber_number,
    sauter_drop_diameter,
    phase_separation_check,
    minimum_dispersion_speed,
    liquid_liquid_mass_transfer,
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
from utils.data_helpers import load_db, safe_float as _safe, all_fluid_names, safe_iloc, find_reactor_image, reactor_search_name, find_reactor_model_3d, render_reactor_3d

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

reactors = load_db("reactor_db", "reactors.csv", ["reactor_name"])
reactions = load_db("reaction_db", "reactions.csv", ["reaction_name"])
custom_fluids = load_db("fluid_db", "fluids.csv", ["fluid_name"])
particles = load_db("particle_db", "particles.csv", ["particle_name"])

_all_fluid_names = all_fluid_names(custom_fluids)

st.title("🌀 Mixing Assessment")

if reactors.empty or reactions.empty:
    st.warning("Please populate the Reactor and Reaction databases first.")
    st.stop()

# ── Step gating: track confirmed input sections ──────────────────────────
if "_p5_step" not in st.session_state:
    st.session_state["_p5_step"] = 0

def _reset_p5():
    st.session_state["_p5_step"] = 0

# ── Step 1: Select system ────────────────────────────────────────────────
st.header("1 · Define System")

# Persist selections across page navigations
_reactor_list = reactors["reactor_name"].tolist()
_reaction_list = reactions["reaction_name"].tolist()

def _idx(lst, key, default=0):
    """Return list index of session_state[key] if present, else default."""
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col_r, col_rx, col_f = st.columns(3)

with col_r:
    reactor_name = st.selectbox("Reactor", _reactor_list, index=_idx(_reactor_list, "_sel_reactor"), key="ms_reactor", on_change=_reset_p5,
                                format_func=lambda n: reactor_search_name(reactors, n))
    st.session_state["_sel_reactor"] = reactor_name
with col_rx:
    reaction_name = st.selectbox("Reaction", _reaction_list, index=_idx(_reaction_list, "_sel_reaction"), key="ms_reaction", on_change=_reset_p5)
    st.session_state["_sel_reaction"] = reaction_name

col_T, col_P, col_cw = st.columns(3)
with col_f:
    fluid_name = st.selectbox("Fluid (continuous phase)", _all_fluid_names, index=_idx(_all_fluid_names, "_sel_fluid"), key="ms_fluid", on_change=_reset_p5)
    st.session_state["_sel_fluid"] = fluid_name
with col_T:
    _is_solvent = is_known_solvent(fluid_name)
    if _is_solvent:
        fluid_T_C = st.number_input("Temperature (°C)", value=25.0, step=1.0, format="%.1f", key="ms_fluid_T")
    else:
        st.caption("Custom fluid — fixed properties (no T dependence)")
        fluid_T_C = 25.0
with col_P:
    if _is_solvent:
        fluid_P_atm = st.number_input("Pressure (atm)", value=1.0, min_value=0.01, step=0.1, format="%.2f", key="ms_fluid_P")
    else:
        fluid_P_atm = 1.0
with col_cw:
    ms_T_coolant = st.number_input(
        "Coolant temperature (°C)", value=15.0, step=1.0,
        format="%.1f", key="ms_T_cool",
        help="Jacket coolant inlet temperature",
    )

# Compute fluid properties (continuous phase)
if _is_solvent:
    _fprops = get_properties(fluid_name, fluid_T_C, fluid_P_atm)
    _rho_default = _fprops["rho_kg_m3"]
    _mu_default = _fprops["mu_Pa_s"]
    _D_mol_default = _fprops["D_mol_m2_s"]
    _sigma_default = _fprops.get("surface_tension_N_m", 0.072)
    if not _fprops["in_range"]:
        st.warning(f"⚠️ {fluid_T_C:.1f} °C is outside the liquid range "
                   f"({_fprops['mp_C']:.0f} – {_fprops['bp_at_P_C']:.0f} °C) for {fluid_name}. "
                   f"Values are extrapolated.")
else:
    _cust = safe_iloc(custom_fluids, "fluid_name", fluid_name, "Custom fluid")
    _rho_default = float(_cust["rho_kg_m3"])
    _mu_default = float(_cust["mu_Pa_s"])
    _D_mol_default = float(_cust["D_mol_m2_s"])
    _sigma_default = float(_cust.get("surface_tension_N_m", 0.072))

reactor = safe_iloc(reactors, "reactor_name", reactor_name, "Reactor")
reaction = safe_iloc(reactions, "reaction_name", reaction_name, "Reaction")

# Use selection names in widget keys so they reset automatically on change
_rk = reactor_name   # reactor key fragment
_fk = f"{fluid_name}_{fluid_T_C:.1f}"  # fluid key fragment (includes T)
_xk = reaction_name  # reaction key fragment

# ── System type selection ────────────────────────────────────────────────
st.subheader("System Type")
st.caption(
    "Select multiphase types as applicable. Leave empty for single-phase liquid."
)
_SYSTEM_OPTIONS = ["Gas-Liquid (GL)", "Liquid-Liquid (LL)", "Solid-Liquid (SL)"]
system_types = st.multiselect(
    "Multiphase system types",
    _SYSTEM_OPTIONS,
    default=st.session_state.get("_sel_sys_types", []),
    key="ms_sys_types",
    on_change=_reset_p5,
    help="Select one or more. Leave empty for single-phase liquid.",
)
st.session_state["_sel_sys_types"] = system_types
include_GL = "Gas-Liquid (GL)" in system_types
include_LL = "Liquid-Liquid (LL)" in system_types
include_SL = "Solid-Liquid (SL)" in system_types

if not system_types:
    st.info("**Liquid-only** mode.")

# ── Process mode ─────────────────────────────────────────────────────────
is_semi_batch = st.checkbox(
    "Semi-batch (fed-batch) process",
    value=st.session_state.get("_sel_semi_batch", False),
    key="ms_semi_batch",
    on_change=_reset_p5,
    help="Reagent or anti-solvent is dosed during the reaction.",
)
st.session_state["_sel_semi_batch"] = is_semi_batch

# ── GL inputs ────────────────────────────────────────────────────────────
v_s = 0.0
is_coalescing = True
gl_bubble_d_mm = 3.0
gl_sparged = False
if include_GL:
    with st.expander("Gas-Liquid — sparging & gas phase", expanded=True):
        gl_mode = st.radio(
            "GL mass-transfer mode",
            ["Sparged (gas fed through sparger)", "Headspace only (free-surface transfer)"],
            horizontal=True, key="ov_gl_mode",
        )
        gl_sparged = gl_mode.startswith("Sparged")
        if gl_sparged:
            gc1, gc2, gc3 = st.columns(3)
            v_s = gc1.number_input(
                "Superficial gas velocity v_s (m/s)", value=0.005, min_value=0.0,
                format="%.4f", key=f"ov_vs_{_rk}",
                help="Typical: 0.001 – 0.05 m/s.",
            )
            coalescing = gc2.selectbox(
                "Liquid type",
                ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                key=f"ov_coal_{_fk}",
            )
            is_coalescing = coalescing.startswith("Coalescing")
            gl_bubble_d_mm = gc3.number_input(
                "Estimated bubble diameter (mm)", value=3.0, min_value=0.1,
                max_value=50.0, format="%.1f", key="ov_bubble_d",
                help="Initial estimate; d₃₂ computed via Calderbank.",
            )
        else:
            st.caption("Headspace only — kLa_surface via Lamont-Scott model.")

# ── LL inputs ────────────────────────────────────────────────────────────
ll_rho_d = 800.0
ll_mu_d = 0.001
ll_sigma_LL = 0.030
ll_phi_d = 0.10
ll_D_mol_LL = 1e-9
if include_LL:
    with st.expander("Liquid-Liquid — dispersed phase", expanded=True):
        st.caption("Define the dispersed phase (continuous phase = fluid above).")
        _ll_source = st.radio(
            "Dispersed-phase definition",
            ["Select from fluid database", "Enter properties manually"],
            horizontal=True, key="ov_ll_source",
        )
        if _ll_source.startswith("Select"):
            _ll_col1, _ll_col2 = st.columns(2)
            _ll_fluid_name = _ll_col1.selectbox(
                "Dispersed-phase fluid", _all_fluid_names,
                key="ov_ll_fluid",
            )
            _ll_is_solvent = is_known_solvent(_ll_fluid_name)
            if _ll_is_solvent:
                _ll_T = _ll_col2.number_input(
                    "Dispersed-phase T (°C)", value=fluid_T_C, step=1.0,
                    format="%.1f", key="ov_ll_T",
                )
                _ll_fprops = get_properties(_ll_fluid_name, _ll_T, fluid_P_atm)
                _ll_rho_val = _ll_fprops["rho_kg_m3"]
                _ll_mu_val = _ll_fprops["mu_Pa_s"]
            else:
                _ll_cust = custom_fluids[custom_fluids["fluid_name"] == _ll_fluid_name].iloc[0]
                _ll_rho_val = float(_ll_cust["rho_kg_m3"])
                _ll_mu_val = float(_ll_cust["mu_Pa_s"])
                _ll_col2.caption("Custom fluid — fixed properties")
            _ll_ov1, _ll_ov2 = st.columns(2)
            ll_rho_d = _ll_ov1.number_input(
                "ρ_d (kg/m³)", value=_ll_rho_val, min_value=100.0,
                format="%.1f", key="ov_ll_rho_d_db",
                help="From database — override if needed",
            )
            ll_mu_d = _ll_ov2.number_input(
                "μ_d (Pa·s)", value=_ll_mu_val, min_value=1e-6,
                format="%.6f", key="ov_ll_mu_d_db",
                help="From database — override if needed",
            )
        else:
            _ll_m1, _ll_m2 = st.columns(2)
            _ll_m1.text_input(
                "Dispersed-phase name", value="Dispersed liquid",
                key="ov_ll_name",
            )
            _ll_m2.caption("")  # spacer
            _ll_m3, _ll_m4 = st.columns(2)
            ll_rho_d = _ll_m3.number_input(
                "ρ_d (kg/m³)", value=800.0, min_value=100.0,
                format="%.1f", key="ov_ll_rho_d",
                help="Dispersed-phase density",
            )
            ll_mu_d = _ll_m4.number_input(
                "μ_d (Pa·s)", value=0.001, min_value=1e-6,
                format="%.6f", key="ov_ll_mu_d",
                help="Dispersed-phase viscosity",
            )
        ll_c1, ll_c2, ll_c3 = st.columns(3)
        ll_sigma_LL = ll_c1.number_input(
            "Interfacial tension σ (N/m)", value=0.030, min_value=1e-5,
            format="%.4f", key="ov_ll_sigma",
            help="Between continuous and dispersed phases",
        )
        ll_phi_d = ll_c2.number_input(
            "Dispersed-phase volume fraction φ_d", value=0.10,
            min_value=0.001, max_value=0.70, format="%.3f",
            key="ov_ll_phi_d",
            help="Range 0–0.7",
        )
        ll_D_mol_LL = ll_c3.number_input(
            "D_mol for LL transfer (m²/s)", value=1e-9,
            format="%.2e", key="ov_ll_Dmol",
            help="Solute diffusivity across L-L interface",
        )

        # ── Miscibility screening ─────────────────────────────────────────
        _cont_name = fluid_name
        _disp_name = _ll_fluid_name if _ll_source.startswith("Select") else None

        if _disp_name:
            _misc = solvent_miscibility(_cont_name, _disp_name)
            st.divider()
            st.markdown("**Miscibility Screening**")
            _hc1, _hc2 = st.columns(2)
            _hc1.metric("Assessment", _misc["assessment"])
            if _misc["Ra"] is not None:
                _hc2.metric("Hansen R_a (MPa½)", f"{_misc['Ra']:.1f}")
            else:
                _hc2.metric("Hansen R_a (MPa½)", "—")

            if _misc["hsp_1"] and _misc["hsp_2"]:
                _hd1, _hd2 = st.columns(2)
                _hd1.caption(f"δ ({_cont_name}): ({_misc['hsp_1'][0]:.1f}, {_misc['hsp_1'][1]:.1f}, {_misc['hsp_1'][2]:.1f}) MPa½")
                _hd2.caption(f"δ ({_disp_name}): ({_misc['hsp_2'][0]:.1f}, {_misc['hsp_2'][1]:.1f}, {_misc['hsp_2'][2]:.1f}) MPa½")

            if _misc["source"] == "lookup":
                _src_note = "Experimental data."
            elif _misc["source"] == "Hansen estimate":
                _src_note = "Hansen estimate — verify experimentally."
            else:
                _src_note = ""

            if _misc["miscible"] is True:
                if "Partially" in _misc["assessment"]:
                    st.warning(f"🟡 **{_misc['assessment']}** — limited mutual solubility; L-L may apply at higher φ.")
                else:
                    st.error(f"🔴 **{_misc['assessment']}** — single phase; L-L dispersion N/A.")
            elif _misc["miscible"] is False:
                st.success(f"🟢 **{_misc['assessment']}** — suitable for L-L dispersion.")
            else:
                st.info("Miscibility unknown for this pair.")

            if _src_note:
                st.caption(_src_note)
        else:
            st.caption("Miscibility screening requires dispersed phase from fluid database.")

# ── SL inputs ────────────────────────────────────────────────────────────
include_particles = include_SL
if include_SL and not particles.empty:
    with st.expander("Solid-Liquid — particle properties", expanded=True):
        st.caption("Solid-phase properties for suspension and S-L mass transfer.")

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
        # Solids loading
        _solids_basis = st.radio("Solids loading basis",
                                 ["Mass (wt-%)", "Volume (vol-%)"],
                                 horizontal=True, key="ov_solids_basis")
        pc4, pc5, pc6 = st.columns(3)
        if _solids_basis == "Mass (wt-%)":
            X_wt = pc4.number_input("Solids loading X (wt-%)", value=5.0,
                                    min_value=0.01, format="%.2f", key="ov_Xwt",
                                    help="(mass solids / mass liquid) × 100")
            X_vol = 100.0 * X_wt * rho_p / (X_wt * rho_p + 100.0 * _rho_default) if rho_p > 0 else 0.0
            pc5.metric("Xv (vol-%)", f"{X_vol:.2f}")
        else:
            X_vol = pc4.number_input("Solids loading Xv (vol-%)", value=2.0,
                                     min_value=0.01, max_value=99.0,
                                     format="%.2f", key="ov_Xv",
                                     help="(vol solids / vol slurry) × 100")
            X_wt = 100.0 * rho_p * X_vol / (_rho_default * (100.0 - X_vol)) if (100.0 - X_vol) > 0 and _rho_default > 0 else 0.0
            pc5.metric("X (wt-%)", f"{X_wt:.2f}")
        S_zw = pc6.number_input("Zwietering S constant", value=5.5, min_value=0.5,
                                max_value=20.0, format="%.1f", key="ov_Szw",
                                help="Typ. 1–10, PBT ≈ 5.5")

        st.markdown("**Grenville, Mak & Brown (GMB) parameters**")
        _gmb_z_default = _safe(reactor.get("GMB_z"), 0.0)
        _C_imp_default = _safe(reactor.get("imp1_clearance_m"), 0.0)
        _C_D_default = _C_imp_default / _safe(reactor.get("D_imp_m"), 0.05) if _safe(reactor.get("D_imp_m"), 0.05) > 0 and _C_imp_default > 0 else 0.33
        pc7, pc8 = st.columns(2)
        gmb_z = pc7.number_input("GMB z constant", value=_gmb_z_default if _gmb_z_default > 0 else 3.0,
                                  min_value=0.1, max_value=30.0, format="%.2f", key="ov_gmb_z",
                                  help="Impeller-type dependent")
        C_D_ratio = pc8.number_input("C/D (clearance / impeller dia)",
                                      value=_C_D_default, min_value=0.01,
                                      max_value=2.0, format="%.3f", key="ov_CD",
                                      help="Clearance / impeller diameter")
elif include_SL and particles.empty:
    st.warning("Particle database is empty. Add particles on the Particle Database page.")

# ── Gate 1: confirm system selection ─────────────────────────────────────
if st.session_state["_p5_step"] < 1:
    st.info("Confirm the system selection above to continue.")
    if st.button("✅ Confirm system selection", key="p5_gate1", type="primary"):
        st.session_state["_p5_step"] = 1
        st.rerun()
    st.stop()

# ── Show iso image of selected reactor ───────────────────────────────────

_model_path = find_reactor_model_3d(reactors, reactor_name)
_iso_path = find_reactor_image(reactors, reactor_name, "iso")
_side_path = find_reactor_image(reactors, reactor_name, "side")
if _model_path or _iso_path or _side_path:
    with st.container(border=True):
        st.markdown(f"**{reactor_name}**")
        if _model_path:
            render_reactor_3d(_model_path, height=360, auto_rotate=False)
            st.caption("3D model — drag to rotate · scroll to zoom · right-drag to pan")
        else:
            _imgs = [(p, lbl) for p, lbl in [(_iso_path, "Iso view"), (_side_path, "Side view")] if p]
            _cols = st.columns(len(_imgs))
            for _i, (_path, _lbl) in enumerate(_imgs):
                with _cols[_i]:
                    st.image(str(_path), caption=_lbl, width='stretch')

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
    V_L_default = np.pi / 4 * D_tank**2 * _safe(reactor.get("H_m"), 0.13) * 1000

V_L = st.number_input(
    "Liquid volume (L)", value=V_L_default, min_value=0.01,
    format="%.2f", key=f"ov_VL_{_rk}",
    help=f"Reactor range: {V_L_min:.1f} – {V_L_max:.1f} L" if V_L_max > 0 else None,
)

_vol_ok = True
if V_L_min > 0 and V_L_max > 0:
    if V_L < V_L_min or V_L > V_L_max:
        st.error(f"Volume {V_L:.1f} L is outside the reactor operating range "
                 f"({V_L_min:.1f} – {V_L_max:.1f} L).")
        _vol_ok = False

H = liquid_height_from_volume(V_L, D_tank, H_max, _bottom_dish)

with st.expander("Fluid properties (continuous phase)", expanded=False):
    if _is_solvent:
        st.caption(f"Computed from **{fluid_name}** correlations at **{fluid_T_C:.1f} °C**")
    fc1, fc2, fc3, fc4 = st.columns(4)
    rho = fc1.number_input("ρ (kg/m³)", value=_rho_default, min_value=0.1, format="%.1f", key=f"ov_rho_{_fk}")
    mu = fc2.number_input("μ (Pa·s)", value=_mu_default, min_value=1e-7, format="%.6f", key=f"ov_mu_{_fk}")
    D_mol = fc3.number_input("D_mol (m²/s)", value=_D_mol_default, format="%.2e", key=f"ov_Dmol_{_fk}")
    sigma_c = fc4.number_input("σ surface tension (N/m)", value=_sigma_default, format="%.4f", key=f"ov_sigma_{_fk}",
                               help="Used for GL/LL calculations")

with st.expander("Reaction parameters", expanded=False):
    rc1, rc2, rc3 = st.columns(3)
    k_val = rc1.number_input("k", value=float(reaction["k_value"]), format="%.6g", key=f"ov_k_{_xk}")
    C0 = rc2.number_input("C₀ (mol/L)", value=float(reaction["C0_mol_L"]), format="%.4g", key=f"ov_C0_{_xk}")
    t_rxn_input = rc3.number_input("t_rxn (s)", value=float(reaction["t_rxn_s"]), format="%.4g", key=f"ov_trxn_{_xk}",
                                    help="0 = auto-compute from k and C₀")

# ── Gate 2: confirm parameter overrides ──────────────────────────────────
if st.session_state["_p5_step"] < 2:
    st.info("Confirm parameters above to continue.")
    if st.button("✅ Confirm parameters", key="p5_gate2", type="primary"):
        st.session_state["_p5_step"] = 2
        st.rerun()
    st.stop()

# ── Step 3: Correlations ─────────────────────────────────────────────
st.divider()
st.header("3 · Correlation Source Selection")

if has_any_alt_correlations(reactor_name):
    with st.expander("Correlation source selection", expanded=False):
        st.caption(
            "Select sources per parameter. Priority: Experimental → ROM → Literature."
        )
        corr_selections = render_correlation_matrix_multi(
            reactor_name, key_prefix="ms_corr",
        )
else:
    st.info(f"Only **Literature** correlations for **{reactor_name}**. Register ROM/Experimental to compare sources.")
    corr_selections = {p: ["Literature"] for p in SUPPORTED_PARAMS}

# Derive a single priority mode dict for the main metrics
_priority_mode = priority_mode_dict(corr_selections)
# Collect all active modes for envelope computation
_active_modes = active_modes_set(corr_selections)


# ── Reaction thermodynamics (used later for auto heat balance) ────────────────
rxn_T_C = _safe(reaction.get("T_C"), 25.0)
rxn_delta_H = _safe(reaction.get("delta_H_kJ_mol"), 0.0)
rxn_order = str(reaction.get("order", "1"))
ms_T_process = fluid_T_C
include_heat = rxn_delta_H != 0

# ── Gate 3: confirm correlations & compute ───────────────────────────────
if st.session_state["_p5_step"] < 3:
    st.info("Confirm correlation sources to compute results.")
    if st.button("🔬 Confirm & Compute", key="p5_gate3", type="primary"):
        st.session_state["_p5_step"] = 3
        st.rerun()
    st.stop()

# ── Step 4: Compute ──────────────────────────────────────────────────────
st.divider()
st.header("4 · Centerpoint Results")

# Compute reaction time if needed
t_rxn = t_rxn_input
if t_rxn <= 0 and k_val > 0:
    order = str(reaction.get("order", "1"))
    if order in ("1", "pseudo-1"):
        # Characteristic time = reciprocal rate constant (time constant 1/k),
        # the e-folding time used in the Damkohler definition.
        t_rxn = 1.0 / k_val
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

# Assessment banners — helper
def _da_banner(label: str, Da: float) -> None:
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

# ── 4a · Mixing ──────────────────────────────────────────────────────────
st.subheader("Mixing")

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

m9, m10, m11, m12 = st.columns(4)
m9.metric("Avg shear rate (1/s)", f"{hydro['Avg shear rate (1/s)']:.1f}")
m10.metric("Max shear rate (1/s)", f"{hydro['Max shear rate (1/s)']:.0f}")
m11.metric("Avg shear stress (Pa)", f"{hydro['Avg shear stress (Pa)']:.3g}")
m12.metric("t_E local (s)", f"{hydro['Micromix time t_E_local (s)']:.4g}")

m9b, m10b, m11b, m12b = st.columns(4)
m9b.metric("Circulation time (s)", f"{hydro['Circulation time (s)']:.2f}")
m10b.metric("Froude number", f"{hydro['Froude number']:.4g}")
m11b.metric("EDCF (W/kg/s)", f"{hydro['EDCF (W/kg/s)']:.3g}")
m12b.metric("Torque (N·m)", f"{hydro['Torque (N·m)']:.3g}")

_da_banner("Macromixing", da["Da_macro"])
_da_banner("Micromixing", da["Da_micro"])

lam_B = batchelor_length(mu / rho, hydro["P/V (W/kg)"], D_mol)
st.info(f"Batchelor microscale λ_B = {lam_B * 1e6:.2f} µm  |  Reaction time t_rxn = {t_rxn:.4g} s")

if is_semi_batch:
    st.warning(
        "⚠️ **Semi-batch process** — three mixing scales act in series at the feed point:\n\n"
        "1. **Macromixing** (θ₉₅) — bulk blending of added reagent into the vessel\n"
        "2. **Mesomixing** ($t_{meso}$) — turbulent dispersion of the feed plume; "
        "controls local concentration before molecular homogenisation\n"
        "3. **Micromixing** ($t_E$) — engulfment at the Kolmogorov scale\n\n"
        "To distinguish these experimentally, run the full **🅱️ Bourne Protocol**:\n"
        "- Vary **impeller speed** → probes micromixing (ε changes)\n"
        "- Vary **feed rate / addition time** → probes mesomixing\n"
        "- Vary **feed location** → probes meso- and macromixing"
    )

# ── 4b · Gas-Liquid Mass Transfer ─────────────────────────────────────────
gl_results = {}
if include_GL:
    st.subheader("Gas-Liquid Mass Transfer")

    _P_V_Wm3 = hydro["P/V (W/m³)"]
    _eps_kg = hydro["P/V (W/kg)"]

    if gl_sparged and v_s > 0:
        mt1, mt2, mt3, mt4 = st.columns(4)
        mt1.metric("kLa sparged (1/s)", f"{hydro['kLa (1/s)']:.4g}")
        mt2.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
        mt3.metric("Da_GL", f"{da['Da_GL']:.3g}")

        # Gas holdup, bubble size, flooding
        _eps_G = gas_holdup_hughmark(v_s, _P_V_Wm3, mu, sigma_c, rho)
        _d32_bubble = sauter_bubble_diameter(_P_V_Wm3, v_s, sigma_c, rho)
        _Q_gas = gas_flow_rate_from_vs(v_s, D_tank)
        _N_flood = gas_flooding_speed(Nq_in, D_imp, _Q_gas)

        mt4.metric("Gas holdup ε_G", f"{_eps_G:.3f}")
        mt5, mt6, mt7, mt8 = st.columns(4)
        mt5.metric("d₃₂ bubble (mm)", f"{_d32_bubble * 1e3:.2f}")
        mt6.metric("Q_gas (m³/s)", f"{_Q_gas:.3e}")
        mt7.metric("N_flood (RPM)", f"{_N_flood * 60:.1f}")
        _flood_ratio = N_rps / _N_flood if _N_flood > 0 else np.inf
        mt8.metric("N/N_flood", f"{_flood_ratio:.2f}")

        if _flood_ratio < 1.0:
            st.error("🔴 **Flooded** — below flooding speed.")
        elif _flood_ratio < 1.3:
            st.warning("🟡 **Near flooding** — increase speed.")
        else:
            st.success(f"🟢 **Good dispersion** — N/N_flood = {_flood_ratio:.2f}")

        gl_results = {
            "GL mode": "Sparged",
            "kLa sparged (1/s)": hydro["kLa (1/s)"],
            "kLa surface (1/s)": hydro["kLa_surface (1/s)"],
            "Gas holdup ε_G": _eps_G,
            "d32 bubble (mm)": _d32_bubble * 1e3,
            "Q_gas (m³/s)": _Q_gas,
            "N_flood (RPM)": _N_flood * 60,
            "N/N_flood": _flood_ratio,
        }
    else:
        # Headspace-only mode (or sparged with v_s=0)
        mt1, mt2, _, _ = st.columns(4)
        mt1.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
        da_gl_val = da["Da_GL"]
        if da_gl_val > 0:
            mt2.metric("Da_GL (surface)", f"{da_gl_val:.3g}")
        if gl_sparged and v_s <= 0:
            st.caption("Sparged selected but v_s = 0 — set v_s > 0.")
        else:
            st.caption("Headspace only — kLa via Lamont-Scott.")
        gl_results = {
            "GL mode": "Headspace",
            "kLa surface (1/s)": hydro["kLa_surface (1/s)"],
        }

    if da["Da_GL"] > 0:
        _da_banner("Gas-liquid mass transfer", da["Da_GL"])
else:
    st.subheader("Mass Transfer (surface only)")
    mt1, mt2, _, _ = st.columns(4)
    mt1.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
    da_gl_val = da["Da_GL"]
    if da_gl_val > 0:
        mt2.metric("Da_GL (surface)", f"{da_gl_val:.3g}")
    st.caption("GL not selected — free-surface kLa only (Lamont-Scott).")

# ── 4c · Liquid-Liquid ───────────────────────────────────────────────────
ll_results = {}
if include_LL:
    st.subheader("Liquid-Liquid Dispersion")

    _eps_kg_ll = hydro["P/V (W/kg)"]
    _ll_check = phase_separation_check(
        N=N_rps, D_imp=D_imp, D_tank=D_tank,
        rho_c=rho, rho_d=ll_rho_d, mu_c=mu,
        sigma_LL=ll_sigma_LL, phi_d=ll_phi_d,
    )
    _We_ll = _ll_check["We"]
    _d32_drop = _ll_check["d32 (m)"]
    _N_min_disp = minimum_dispersion_speed(D_imp, ll_sigma_LL, rho, ll_phi_d)
    _k_LL = liquid_liquid_mass_transfer(_d32_drop, ll_D_mol_LL, rho, mu, _eps_kg_ll)
    _a_LL = 6.0 * ll_phi_d / _d32_drop if _d32_drop > 0 else 0.0
    _kLa_LL = _k_LL * _a_LL

    ll1, ll2, ll3, ll4 = st.columns(4)
    ll1.metric("Weber number", f"{_We_ll:.1f}")
    ll2.metric("d₃₂ drop (µm)", f"{_ll_check['d32 (µm)']:.1f}")
    ll3.metric("N_min dispersion (RPM)", f"{_N_min_disp * 60:.1f}")
    ll4.metric("N/N_min", f"{N_rps / _N_min_disp:.2f}" if _N_min_disp > 0 else "—")

    ll5, ll6, ll7, ll8 = st.columns(4)
    ll5.metric("k_LL (m/s)", f"{_k_LL:.3e}")
    ll6.metric("a_LL (1/m)", f"{_a_LL:.1f}")
    ll7.metric("kLa_LL (1/s)", f"{_kLa_LL:.4g}")
    ll8.metric("Drop settling (m/s)", f"{_ll_check['Drop settling velocity (m/s)']:.3e}")

    ll9, ll10, _, _ = st.columns(4)
    ll9.metric("Separation time (s)", f"{_ll_check['Separation time (s)']:.0f}" if _ll_check["Separation time (s)"] < 1e6 else "∞")
    ll10.metric("Δρ (kg/m³)", f"{abs(ll_rho_d - rho):.1f}")

    # Dispersion stability assessment
    _disp_ratio = N_rps / _N_min_disp if _N_min_disp > 0 else 0
    if _disp_ratio < 0.8:
        st.error(f"🔴 **Poor dispersion** — N/N_min = {_disp_ratio:.2f}")
    elif _disp_ratio < 1.0:
        st.warning(f"🟡 **Marginal** — near minimum dispersing speed.")
    else:
        st.success(f"🟢 **Good dispersion** — N/N_min = {_disp_ratio:.2f}")

    st.info(f"Phase separation estimate: **{_ll_check['Assessment']}**")

    ll_results = {
        "We": _We_ll,
        "d32 drop (µm)": _ll_check["d32 (µm)"],
        "Drop settling velocity (m/s)": _ll_check["Drop settling velocity (m/s)"],
        "Separation time (s)": _ll_check["Separation time (s)"],
        "N_min dispersion (RPM)": _N_min_disp * 60,
        "k_LL (m/s)": _k_LL,
        "a_LL (1/m)": _a_LL,
        "kLa_LL (1/s)": _kLa_LL,
        "LL Assessment": _ll_check["Assessment"],
    }

# ── 4d · Heat Transfer ───────────────────────────────────────────────────
st.subheader("Heat Transfer")

heat_results = {}
if include_heat:
    _r_material = str(reactor.get("shell_material", ""))
    _r_lining = str(reactor.get("lining_material", ""))
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
            lining_material=_r_lining,
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

    hm5, hm6, _, _ = st.columns(4)
    hm5.metric("ΔT process–coolant (°C)", f"{_dT:.1f}")
    hm6.metric("ΔH (kJ/mol)", f"{rxn_delta_H:.1f}")

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
else:
    st.caption("No ΔH data — heat balance skipped.")

# ── 4e · Solid-Liquid ────────────────────────────────────────────────────
particle_results = {}
particle_meta = {}
if include_SL and not particles.empty:
    st.subheader("Solid-Liquid")
    d_p_m = d50_um * 1e-6
    nu = mu / rho if rho > 0 else 0.0
    delta_rho = abs(rho_p - rho)

    v_t = settling_velocity(d_p_m, rho_p, rho, mu, phi_p)
    Re_p = particle_reynolds(d_p_m, v_t, rho, mu)
    N_js_zw = zwietering_njs(S_zw, nu, d_p_m, delta_rho, rho, X_wt, D_imp)
    N_js_gmb = gmb_njs(gmb_z, Np_in if Np_in else 1.27, D_imp, d_p_m,
                       delta_rho, rho, X_vol, C_D_ratio)
    N_js = max(N_js_zw, N_js_gmb)
    _eps_kg_cp = hydro["P/V (W/kg)"] if hydro["P/V (W/kg)"] > 0 else 0.0
    _v_slip_cp = max(v_t, (_eps_kg_cp * d_p_m) ** (1.0 / 3.0) if _eps_kg_cp > 0 else 0.0)
    k_SL = solid_liquid_mass_transfer(d_p_m, _v_slip_cp, rho, mu, D_mol)
    _Ar = archimedes_number(d_p_m, rho, delta_rho, mu)
    _phi_s = X_vol / 100.0  # volume fraction
    _kLa_SL = solid_liquid_kla(k_SL, d_p_m, _phi_s)
    susp = particle_suspension_criterion(N_rps, N_js)

    particle_results = {
        "d50 (µm)": d50_um,
        "ρ_p (kg/m³)": rho_p,
        "Shape factor φ": phi_p,
        "X (wt-%)": X_wt,
        "Xv (vol-%)": X_vol,
        "v_t (m/s)": v_t,
        "Re_p": Re_p,
        "Archimedes Ar": _Ar,
        "N_js Zwietering (rev/s)": N_js_zw,
        "N_js Zwietering (RPM)": N_js_zw * 60,
        "N_js GMB (rev/s)": N_js_gmb,
        "N_js GMB (RPM)": N_js_gmb * 60,
        "N_js (rev/s)": N_js,
        "N_js (RPM)": N_js * 60,
        "k_SL (m/s)": k_SL,
        "kLa_SL (1/s)": _kLa_SL,
    }
    particle_meta = {
        "Particle": particle_name,
        "Suspension": susp,
    }

    sp1, sp2, sp3, sp4 = st.columns(4)
    sp1.metric("Settling velocity (m/s)", f"{v_t:.3e}")
    sp2.metric("Re_p", f"{Re_p:.3g}")
    sp3.metric("N_js Zwietering (RPM)", f"{N_js_zw * 60:.1f}")
    sp4.metric("N_js GMB (RPM)", f"{N_js_gmb * 60:.1f}")

    sp5, sp6, sp7, sp8 = st.columns(4)
    sp5.metric("N_js design (RPM)", f"{N_js * 60:.1f}",
              help="Higher of Zwietering and GMB estimates")
    sp6.metric("k_SL (m/s)", f"{k_SL:.3e}")
    sp7.metric("kLa_SL (1/s)", f"{_kLa_SL:.4g}")
    sp8.metric("Archimedes Ar", f"{_Ar:.3g}")

    # Recompute Damköhler numbers including Da_SL
    _Da_SL = damkohler_sl(_kLa_SL, t_rxn)
    da["Da_SL"] = _Da_SL
    da["Assessment"] = mixing_sensitivity_assessment(
        da["Da_macro"], da["Da_micro"], da["Da_GL"], _Da_SL,
    )
    particle_results["Da_SL"] = _Da_SL

    sp9, _, _, _ = st.columns(4)
    sp9.metric("Da_SL", f"{_Da_SL:.3g}")

    _da_banner("Solid-liquid mass transfer", _Da_SL)

    if "Poorly" in susp:
        st.error(f"🔴 **{susp}** — below N_js")
    elif "Partially" in susp:
        st.warning(f"🟡 **{susp}** — increase speed")
    elif "Just" in susp:
        st.info(f"🟢 **{susp}** — near N_js")
    else:
        st.success(f"🟢 **{susp}**")

# ── Operating Envelope Charts ────────────────────────────────────────────
st.divider()
st.header("5 · Operating Envelopes")
st.caption("Sweep across full RPM and volume range.")

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
    _r_material_env = str(reactor.get("shell_material", ""))
    _r_lining_env = str(reactor.get("lining_material", ""))
    _r_U_override_env = _safe(reactor.get("U_W_m2K"), 0.0)
    _r_A_override_env = _safe(reactor.get("A_ht_m2"), 0.0)
    _r_wall_mm_env = _safe(reactor.get("wall_thickness_mm"), 0.0)

    PLOT_PARAMS = [
        "Power (W)", "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
        "Circulation time (s)",
        "Micromix time t_E (s)", "Micromix time t_E_local (s)",
        "Kolmogorov η (µm)", "Re",
        "Avg shear rate (1/s)", "Max shear rate (1/s)", "Avg shear stress (Pa)",
        "Da_macro", "Da_micro", "Da_GL", "Da_SL", "ε_max (W/kg)",
        "EDCF (W/kg/s)", "Torque (N·m)", "Torque/V (N·m/m³)", "Froude number",
        "kLa (1/s)", "kLa_surface (1/s)",
    ]
    HEAT_PARAMS = ["Q_gen (W)", "Q_cool (W)", "U (W/m²·K)", "A_ht (m²)", "Q_gen/Q_cool (%)"]
    if include_heat and rxn_delta_H != 0:
        PLOT_PARAMS = PLOT_PARAMS + HEAT_PARAMS
    PARTICLE_PARAMS = ["N_js Zwietering (RPM)", "N_js GMB (RPM)", "N_js (RPM)", "N/N_js", "v_t (m/s)", "Re_p", "k_SL (m/s)", "kLa_SL (1/s)", "Da_SL"]
    if include_SL and not particles.empty:
        PLOT_PARAMS = PLOT_PARAMS + PARTICLE_PARAMS
    GL_PARAMS = ["Gas holdup ε_G", "d32 bubble (mm)", "N_flood (RPM)", "N/N_flood"]
    if include_GL and gl_sparged and v_s > 0:
        PLOT_PARAMS = PLOT_PARAMS + GL_PARAMS
    LL_PARAMS = ["We", "d32 drop (µm)", "N_min disp (RPM)", "N/N_min disp", "kLa_LL (1/s)"]
    if include_LL:
        PLOT_PARAMS = PLOT_PARAMS + LL_PARAMS

    _N_INTERP = 50
    N_arr = np.linspace(_N_lo, _N_hi, _N_INTERP)
    pct_arr = N_arr / _N_hi * 100 if _N_hi > 0 else np.zeros(_N_INTERP)

    # ── Pre-compute RPM-independent particle quantities (hoisted) ────────
    _part_static: dict | None = None
    if include_SL and not particles.empty:
        _dp = d50_um * 1e-6
        _nu_p = mu / rho
        _drho = abs(rho_p - rho)
        _vt = settling_velocity(_dp, rho_p, rho, mu, phi_p)
        _rep = particle_reynolds(_dp, _vt, rho, mu)
        _njs_zw = zwietering_njs(S_zw, _nu_p, _dp, _drho, rho, X_wt, D_imp)
        _njs_gmb = gmb_njs(gmb_z, Np_in if Np_in else 1.27, D_imp, _dp,
                           _drho, rho, X_vol, C_D_ratio)
        _njs = max(_njs_zw, _njs_gmb)
        _phi_s_env = X_vol / 100.0
        _part_static = {
            "N_js_zw_rps": _njs_zw,
            "N_js_gmb_rps": _njs_gmb,
            "N_js_rps": _njs,
            "v_t (m/s)": _vt,
            "Re_p": _rep,
            "N_js Zwietering (RPM)": _njs_zw * 60,
            "N_js GMB (RPM)": _njs_gmb * 60,
            "N_js (RPM)": _njs * 60,
            "d_p": _dp,
            "phi_s": _phi_s_env,
        }

    # ── Pre-compute LL RPM-independent quantities ────────────────────────
    _ll_N_min_static: float = 0.0
    if include_LL:
        _ll_N_min_static = minimum_dispersion_speed(D_imp, ll_sigma_LL, rho, ll_phi_d)

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

    # ── Sweep only unique mode dicts (memoized) ──────────────────────────
    # curve_data: mode_label → vol_key → {param: np.array}
    # Build a signature of every input that affects the sweep so the result
    # can be reused across reruns when nothing relevant changed.
    _L = locals()

    def _sig_round(x, n=6):
        try:
            return round(float(x), n)
        except (TypeError, ValueError):
            return x

    _env_sig = (
        reactor_name, _sig_round(D_imp), _sig_round(D_tank), _sig_round(H_max),
        str(_bottom_dish), _sig_round(_env_V_min), _sig_round(_env_V_max),
        _sig_round(_N_lo), _sig_round(_N_hi), _N_INTERP,
        _sig_round(rho), _sig_round(mu), _sig_round(Np_in), _sig_round(Nq_in),
        _sig_round(v_s), bool(is_coalescing), _sig_round(D_mol), _sig_round(t_rxn),
        tuple(_active_modes),
        tuple(sorted((k, tuple(v)) for k, v in corr_selections.items())),
        tuple(PLOT_PARAMS),
        # Heat
        bool(include_heat), _sig_round(rxn_delta_H),
        str(_L.get("rxn_order", "")), _sig_round(_L.get("k_val", 0.0)),
        _sig_round(_L.get("C0", 0.0)), _sig_round(_r_U_override_env),
        _sig_round(_r_A_override_env), _sig_round(_r_wall_mm_env),
        str(_r_material_env), str(_r_lining_env), str(fluid_name),
        _sig_round(_L.get("ms_T_process", 0.0)), _sig_round(_L.get("ms_T_coolant", 0.0)),
        # Solid-liquid
        bool(include_SL and not particles.empty),
        _sig_round(_L.get("d50_um", 0.0)), _sig_round(_L.get("rho_p", 0.0)),
        _sig_round(_L.get("phi_p", 0.0)), _sig_round(_L.get("S_zw", 0.0)),
        _sig_round(_L.get("gmb_z", 0.0)), _sig_round(_L.get("X_wt", 0.0)),
        _sig_round(_L.get("X_vol", 0.0)), _sig_round(_L.get("C_D_ratio", 0.0)),
        # Gas-liquid
        bool(include_GL and gl_sparged and v_s > 0), _sig_round(sigma_c),
        # Liquid-liquid
        bool(include_LL), _sig_round(_L.get("ll_sigma_LL", 0.0)),
        _sig_round(_L.get("ll_phi_d", 0.0)), _sig_round(_L.get("ll_D_mol_LL", 0.0)),
    )

    _env_cached = (
        st.session_state.get("_ms_env_sig") == _env_sig
        and "_ms_unique_curve_data" in st.session_state
    )
    _unique_curve_data: dict[str, dict[str, dict[str, np.ndarray]]] = (
        st.session_state["_ms_unique_curve_data"] if _env_cached else {}
    )

    if not _env_cached:
      with st.spinner("Computing operating envelopes…"):
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
                _kLa_SL_env = 0.0
                if _part_static is not None:
                    _eps_kg = _h["P/V (W/kg)"] if _h["P/V (W/kg)"] > 0 else 0.0
                    _v_slip = max(_part_static["v_t (m/s)"],
                                  (_eps_kg * _part_static["d_p"]) ** (1.0 / 3.0) if _eps_kg > 0 else 0.0)
                    _ksl_env = solid_liquid_mass_transfer(_part_static["d_p"], _v_slip, rho, mu, D_mol)
                    _kLa_SL_env = solid_liquid_kla(_ksl_env, _part_static["d_p"], _part_static["phi_s"])
                _da = compute_damkohler_numbers(
                    _h["Blend time 95% (s)"], _h["Micromix time t_E (s)"], t_rxn,
                    kLa=_h["kLa (1/s)"], kLa_surface=_h["kLa_surface (1/s)"],
                    kLa_SL=_kLa_SL_env,
                )
                _vals = {**_h, **_da}
                if _part_static is not None:
                    _vals["N_js Zwietering (RPM)"] = _part_static["N_js Zwietering (RPM)"]
                    _vals["N_js GMB (RPM)"] = _part_static["N_js GMB (RPM)"]
                    _vals["N_js (RPM)"] = _part_static["N_js (RPM)"]
                    _vals["N/N_js"] = _N / _part_static["N_js_rps"] if _part_static["N_js_rps"] > 0 else 0.0
                    _vals["v_t (m/s)"] = _part_static["v_t (m/s)"]
                    _vals["Re_p"] = _part_static["Re_p"]
                    _vals["k_SL (m/s)"] = _ksl_env
                    _vals["kLa_SL (1/s)"] = _kLa_SL_env
                    _vals["Da_SL"] = _da["Da_SL"]
                if include_GL and gl_sparged and v_s > 0:
                    _P_V_e = _h["P/V (W/m³)"]
                    _vals["Gas holdup ε_G"] = gas_holdup_hughmark(v_s, _P_V_e, mu, sigma_c, rho)
                    _vals["d32 bubble (mm)"] = sauter_bubble_diameter(_P_V_e, v_s, sigma_c, rho) * 1e3
                    _Q_gas_e = gas_flow_rate_from_vs(v_s, D_tank)
                    _N_flood_e = gas_flooding_speed(Nq_in, D_imp, _Q_gas_e)
                    _vals["N_flood (RPM)"] = _N_flood_e * 60
                    _vals["N/N_flood"] = _N / _N_flood_e if _N_flood_e > 0 else np.inf
                if include_LL:
                    _We_e = weber_number(rho, _N, D_imp, ll_sigma_LL)
                    _d32_e = sauter_drop_diameter(_We_e, D_imp, ll_phi_d)
                    _vals["We"] = _We_e
                    _vals["d32 drop (µm)"] = _d32_e * 1e6
                    _vals["N_min disp (RPM)"] = _ll_N_min_static * 60
                    _vals["N/N_min disp"] = _N / _ll_N_min_static if _ll_N_min_static > 0 else 0.0
                    _eps_kg_e = _h["P/V (W/kg)"]
                    _k_LL_e = liquid_liquid_mass_transfer(_d32_e, ll_D_mol_LL, rho, mu, _eps_kg_e)
                    _a_LL_e = 6.0 * ll_phi_d / _d32_e if _d32_e > 0 else 0.0
                    _vals["kLa_LL (1/s)"] = _k_LL_e * _a_LL_e
                if include_heat and rxn_delta_H != 0:
                    if _r_U_override_env > 0:
                        _U_ht_e = _r_U_override_env
                    else:
                        _U_ht_e, _ = estimate_U_detailed(
                            N_rps=_N, D_imp=D_imp, D_tank=D_tank,
                            rho=rho, mu=mu,
                            material=_r_material_env,
                            lining_material=_r_lining_env,
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
      st.session_state["_ms_env_sig"] = _env_sig
      st.session_state["_ms_unique_curve_data"] = _unique_curve_data

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
        "Da_SL": "Solid-Liquid Mass Transfer (Da_SL)",
        "Q_gen/Q_cool (%)": "Heat Transfer Capacity (Q_gen/Q_cool (%))",
        "Gas holdup ε_G": "Gas Holdup ε_G",
        "d32 bubble (mm)": "Sauter Mean Bubble Diameter d₃₂ (mm)",
        "N/N_flood": "N/N_flood (Gas Flooding Ratio)",
        "N_flood (RPM)": "Gas Flooding Speed (RPM)",
        "We": "Weber Number (LL)",
        "d32 drop (µm)": "Sauter Mean Drop Diameter d₃₂ (µm)",
        "N_min disp (RPM)": "Min Dispersion Speed (RPM)",
        "N/N_min disp": "N/N_min (LL Dispersion Ratio)",
        "kLa_LL (1/s)": "LL Mass Transfer kLa_LL (1/s)",
        "kLa_SL (1/s)": "SL Mass Transfer kLa_SL (1/s)",
    }
    _display = lambda p: _DISPLAY_NAMES.get(p, p)

    _DEFAULT_PARAMS = ["Da_micro", "Da_macro", "Da_GL", "Q_gen/Q_cool (%)", "Blend time 95% (s)", "P/V (W/L)"]
    if include_GL and gl_sparged and v_s > 0:
        _DEFAULT_PARAMS.extend(["N/N_flood", "Gas holdup ε_G"])
    if include_LL:
        _DEFAULT_PARAMS.extend(["N/N_min disp", "d32 drop (µm)"])
    if include_SL and not particles.empty:
        _DEFAULT_PARAMS.extend(["N/N_js"])
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
            if param in ("Da_macro", "Da_micro", "Da_GL", "Da_SL"):
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

            if param == "N/N_flood":
                fig.add_shape(
                    type="line", x0=0, x1=1, y0=1.0, y1=1.0,
                    xref="paper", yref="y",
                    line=dict(color="red", width=1.5, dash="dash"),
                )
                fig.add_annotation(
                    x=1.0, xref="paper", xanchor="right",
                    y=1.0, yref="y", yanchor="bottom", yshift=2,
                    text="N/N_flood = 1 (flooding)",
                    font=dict(size=11, color="red"), showarrow=False,
                )

            if param == "N/N_min disp":
                fig.add_shape(
                    type="line", x0=0, x1=1, y0=1.0, y1=1.0,
                    xref="paper", yref="y",
                    line=dict(color="red", width=1.5, dash="dash"),
                )
                fig.add_annotation(
                    x=1.0, xref="paper", xanchor="right",
                    y=1.0, yref="y", yanchor="bottom", yshift=2,
                    text="N/N_min = 1 (min dispersion)",
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
            st.plotly_chart(fig, width='content')

        with st.expander("Chart legend"):
            st.markdown("""
- **Shaded region** = volume envelope (min → max fill)
- **Solid line** = max volume · **Dotted** = min volume
- **★ Red star** = current operating point
- 🔵 Literature · 🟢 ROM · 🟠 Experimental
""")
else:
    st.info("No RPM range defined — envelope unavailable.")

# Full table
st.subheader("Full Hydrodynamic Parameter Table")
hydro_df = pd.DataFrame([hydro]).T
hydro_df.columns = ["Value"]
if gl_results:
    gl_df = pd.DataFrame([gl_results]).T
    gl_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, gl_df])
if ll_results:
    ll_df = pd.DataFrame([ll_results]).T
    ll_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, ll_df])
if particle_results:
    part_df = pd.DataFrame([particle_results]).T
    part_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, part_df])
if heat_results:
    heat_df = pd.DataFrame([heat_results]).T
    heat_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, heat_df])
st.dataframe(hydro_df, width='content')
if particle_results:
    st.caption(f"Particle: **{particle_meta['Particle']}**  ·  {particle_meta['Suspension']}")

# ── Persist snapshot for report page ─────────────────────────────────────
def _serialise_envelope() -> dict | None:
    """Convert envelope data (with numpy arrays) to JSON-safe form."""
    if not _can_envelope:
        return None
    cd_serial = {}
    for ml, vols in curve_data.items():
        cd_serial[ml] = {}
        for vk, params in vols.items():
            cd_serial[ml][vk] = {p: arr.tolist() for p, arr in params.items()}
    return {
        "curve_data": cd_serial,
        "pct_arr": pct_arr.tolist(),
        "active_modes": list(_active_modes),
        "priority_mode_label": _priority_mode_label,
        "current_pct": _current_pct,
        "env_V_max": _env_V_max,
        "env_V_min": _env_V_min,
        "rpm_max": _rpm_max,
        "rpm_min": _rpm_min,
    }

st.session_state["_ms_report_snapshot"] = {
    "reactor": reactor_name,
    "reaction": reaction_name,
    "fluid": fluid_name,
    "fluid_T_C": fluid_T_C,
    "N_rpm": _N_rpm_input,
    "V_L": V_L,
    "system_types": system_types,
    "hydro": dict(hydro),
    "da": dict(da),
    "t_rxn": t_rxn,
    "heat_results": dict(heat_results) if heat_results else {},
    "gl_results": dict(gl_results) if gl_results else {},
    "ll_results": dict(ll_results) if ll_results else {},
    "particle_results": dict(particle_results) if particle_results else {},
    "particle_meta": dict(particle_meta) if particle_meta else {},
    "batchelor_um": lam_B * 1e6,
    "envelope": _serialise_envelope(),
}

# ── Step 4: Save to recorded results ────────────────────────────────────
st.header("6 · Save Result")

if st.button("📌 Save this result to Recorded Results"):
    result_row = {
        "reactor": reactor_name,
        "reaction": reaction_name,
        "fluid": fluid_name,
        "fluid_T_C": fluid_T_C,
        "RPM": _N_rpm_input,
        "Volume (L)": V_L,
        "Re": hydro["Re"],
        "P/V (W/L)": hydro["P/V (W/L)"],
        "Tip speed (m/s)": hydro["Tip speed (m/s)"],
        "Blend time (s)": hydro["Blend time 95% (s)"],
        "Circulation time (s)": hydro["Circulation time (s)"],
        "Micromix t_E (s)": hydro["Micromix time t_E (s)"],
        "Micromix t_E_local (s)": hydro["Micromix time t_E_local (s)"],
        "Kolmogorov η (µm)": hydro["Kolmogorov η (µm)"],
        "Batchelor λ_B (µm)": lam_B * 1e6,
        "EDCF (W/kg/s)": hydro["EDCF (W/kg/s)"],
        "Torque (N·m)": hydro["Torque (N·m)"],
        "Froude number": hydro["Froude number"],
        "Avg shear rate (1/s)": hydro["Avg shear rate (1/s)"],
        "Max shear rate (1/s)": hydro["Max shear rate (1/s)"],
        "Avg shear stress (Pa)": hydro["Avg shear stress (Pa)"],
        "kLa (1/s)": hydro["kLa (1/s)"],
        "kLa_surface (1/s)": hydro["kLa_surface (1/s)"],
        "t_rxn (s)": t_rxn,
        "Da_macro": da["Da_macro"],
        "Da_micro": da["Da_micro"],
        "Da_GL": da["Da_GL"],
        "Da_SL": da.get("Da_SL", 0.0),
        "Assessment": da["Assessment"],
    }
    if particle_results:
        result_row.update({
            "Particle": particle_meta.get("Particle", ""),
            "d50 (µm)": particle_results.get("d50 (µm)", ""),
            "ρ_p (kg/m³)": particle_results.get("ρ_p (kg/m³)", ""),
            "v_t (m/s)": particle_results.get("v_t (m/s)", ""),
            "Re_p": particle_results.get("Re_p", ""),
            "Archimedes Ar": particle_results.get("Archimedes Ar", ""),
            "N_js Zwietering (RPM)": particle_results.get("N_js Zwietering (RPM)", ""),
            "N_js GMB (RPM)": particle_results.get("N_js GMB (RPM)", ""),
            "N_js (RPM)": particle_results.get("N_js (RPM)", ""),
            "k_SL (m/s)": particle_results.get("k_SL (m/s)", ""),
            "kLa_SL (1/s)": particle_results.get("kLa_SL (1/s)", ""),
            "Suspension": particle_meta.get("Suspension", ""),
        })
    if gl_results:
        result_row.update({
            "Gas holdup ε_G": gl_results.get("Gas holdup ε_G", ""),
            "d32 bubble (mm)": gl_results.get("d32 bubble (mm)", ""),
            "Q_gas (m³/s)": gl_results.get("Q_gas (m³/s)", ""),
            "N_flood (RPM)": gl_results.get("N_flood (RPM)", ""),
            "N/N_flood": gl_results.get("N/N_flood", ""),
        })
    if ll_results:
        result_row.update({
            "We (LL)": ll_results.get("We", ""),
            "d32 drop (µm)": ll_results.get("d32 drop (µm)", ""),
            "N_min dispersion (RPM)": ll_results.get("N_min dispersion (RPM)", ""),
            "k_LL (m/s)": ll_results.get("k_LL (m/s)", ""),
            "kLa_LL (1/s)": ll_results.get("kLa_LL (1/s)", ""),
            "LL Assessment": ll_results.get("LL Assessment", ""),
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
    st.success("Saved — see **Recorded Results** page.")

# ── Generate PDF Report ──────────────────────────────────────────────────
st.divider()
st.header("7 · Export Report")

if st.button("📥 Export PDF Report", type="primary", key="p5_export_pdf"):
    with st.spinner("Generating PDF…"):
        try:
            from utils.report_builder import build_mixing_assessment_pdf, report_filename
            _pdf_bytes = build_mixing_assessment_pdf(st.session_state["_ms_report_snapshot"])
            st.session_state["_p5_pdf_bytes"] = _pdf_bytes
            st.session_state["_p5_pdf_name"] = report_filename(
                "Mixing_Assessment", reactor_name
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_p5_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_p5_pdf_bytes"],
        file_name=st.session_state["_p5_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")
