"""
Page 5 – Hydrodynamic & Mixing Comparison of Selected Reactors
===============================================================
Compare key hydrodynamic parameters side-by-side across multiple reactors
for a given fluid system.  Each reactor's full operating envelope is mapped
from its RPM and volume ranges (4 corner conditions).
"""

import streamlit as st

import pandas as pd
import numpy as np
import re
import pathlib
import plotly.graph_objects as go

from utils.calculations import (
    compute_reactor_hydro, compute_damkohler_numbers,
    settling_velocity, particle_reynolds, zwietering_njs,
    solid_liquid_mass_transfer, particle_suspension_criterion,
    reaction_rate_mol_per_s, heat_generation_rate,
    estimate_jacket_area, estimate_U, estimate_U_detailed,
    heat_removal_capacity, heat_balance_assessment,
    liquid_height_from_volume,
)
from utils.rom_registry import (
    compute_reactor_hydro_with_mode,
    get_correlations,
    has_any_alt_correlations,
    PARAM_DISPLAY,
)
from utils.corr_widgets import render_correlation_matrix

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── Colour palette by scale ──────────────────────────────────────────────
SCALE_COLORS = {"Lab": "#1f77b4", "Pilot": "#ff7f0e", "Manufacturing": "#2ca02c"}


def _load(key, fn):
    if key not in st.session_state:
        p = DATA_DIR / fn
        st.session_state[key] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    return st.session_state[key]


def _safe_float(val, default=0.0):
    try:
        v = float(val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


reactors = _load("reactor_db", "reactors.csv")
fluids = _load("fluid_db", "fluids.csv")
reactions = _load("reaction_db", "reactions.csv")
particles_db = _load("particle_db", "particles.csv")

st.title("📊 Reactor Comparison")

if reactors.empty:
    st.warning("Populate the Reactor Database first.")
    st.stop()

# ── Selection ────────────────────────────────────────────────────────────
st.header("1 · Select Reactors & Conditions")

_all_reactor_names = reactors["reactor_name"].tolist()
_preferred_defaults = ["Nalas – EasyMax 102", "Cambrex – R-802", "Nalas – 20-L"]
_defaults = [n for n in _preferred_defaults if n in _all_reactor_names] or _all_reactor_names[:4]

# Restore previous selection if available, otherwise use defaults
_saved = st.session_state.get("_sel_cmp_reactors")
if _saved is not None:
    _initial = [n for n in _saved if n in _all_reactor_names]
else:
    _initial = _defaults

selected_names = st.multiselect(
    "Reactors to compare",
    _all_reactor_names,
    default=_initial,
    key="cmp_reactors",
)
st.session_state["_sel_cmp_reactors"] = selected_names

if not selected_names:
    st.info("Select at least one reactor above.")
    st.stop()

def _parse_fluid_temp(name: str, fallback: float = 25.0) -> float:
    """Extract temperature from fluid name like 'Water (25 °C)'."""
    m = re.search(r'\(([-\d.]+)\s*°?C\)', name)
    return float(m.group(1)) if m else fallback


def _find_fluid_for_reaction(solvent: str, T_C: float,
                              fluid_names: list[str]) -> tuple[str | None, str | None]:
    """Find the best matching fluid name for a reaction's solvent and temperature.

    Returns (matched_fluid_name, warning_message_or_None).
    """
    if not solvent or not fluid_names:
        return None, None

    solvent_lower = solvent.strip().lower()

    # Map of common reaction-solvent names → fluid-DB base names
    _SOLVENT_ALIASES: dict[str, list[str]] = {
        "meoh":  ["methanol"],
        "etoh":  ["ethanol"],
        "water": ["water"],
        "h2o":   ["water"],
        "thf":   ["thf", "tetrahydrofuran"],
        "dcm":   ["dcm", "dichloromethane"],
        "dmf":   ["dmf", "dimethylformamide"],
        "dmso":  ["dmso", "dimethyl sulfoxide"],
        "toluene": ["toluene"],
        "heptane": ["heptane"],
        "acetonitrile": ["acetonitrile"],
        "ipa": ["ipa", "isopropanol"],
        "ethyl acetate": ["ethyl acetate"],
        "mek": ["mek"],
        "acetone": ["acetone"],
    }

    # Build a list of base-name candidates for matching
    base_candidates = [solvent_lower]
    for alias, targets in _SOLVENT_ALIASES.items():
        if alias in solvent_lower:
            base_candidates.extend(targets)

    # 1) Try exact temperature match
    for base in base_candidates:
        for fn in fluid_names:
            fn_lower = fn.lower()
            fn_base = re.sub(r"\s*\(.*?\)\s*$", "", fn_lower).strip()
            if fn_base in base_candidates or base in fn_base:
                fn_T = _parse_fluid_temp(fn, 25.0)
                if abs(fn_T - T_C) < 0.5:
                    return fn, None

    # 2) Fallback to 25 °C
    for base in base_candidates:
        for fn in fluid_names:
            fn_lower = fn.lower()
            fn_base = re.sub(r"\s*\(.*?\)\s*$", "", fn_lower).strip()
            if fn_base in base_candidates or base in fn_base:
                fn_T = _parse_fluid_temp(fn, 25.0)
                if abs(fn_T - 25.0) < 0.5:
                    warn = (f"Fluid \"{fn}\" at 25 °C selected – "
                            f"reaction specifies {solvent} at {T_C:.0f} °C "
                            f"which is not in the fluid database.")
                    return fn, warn

    return None, None


# -- Auto-select fluid when reaction changes --------------------------------
_fluid_names = fluids["fluid_name"].tolist() if not fluids.empty else []

if not reactions.empty and _fluid_names:
    # Detect reaction change
    _prev_rxn = st.session_state.get("_prev_cmp_rxn")
    _cur_rxn = st.session_state.get("cmp_rxn", None)
    if _cur_rxn is not None and _cur_rxn != _prev_rxn:
        st.session_state["_prev_cmp_rxn"] = _cur_rxn
        _rxn_row = reactions[reactions["reaction_name"] == _cur_rxn]
        if not _rxn_row.empty:
            _rxn_solvent = str(_rxn_row.iloc[0].get("solvent", ""))
            _rxn_T = _safe_float(_rxn_row.iloc[0].get("T_C"), 25.0)
            _matched, _fluid_warn = _find_fluid_for_reaction(
                _rxn_solvent, _rxn_T, _fluid_names)
            if _matched and _matched != st.session_state.get("cmp_fluid"):
                st.session_state["cmp_fluid"] = _matched
                st.session_state["_sel_cmp_fluid"] = _matched
                if _fluid_warn:
                    st.session_state["_fluid_auto_warn"] = _fluid_warn
                else:
                    st.session_state.pop("_fluid_auto_warn", None)
                st.rerun()

def _sel_idx(lst, key, default=0):
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col1, col2 = st.columns(2)
with col1:
    if not fluids.empty:
        fluid_name = st.selectbox("Fluid system", _fluid_names, index=_sel_idx(_fluid_names, "_sel_cmp_fluid"), key="cmp_fluid")
        st.session_state["_sel_cmp_fluid"] = fluid_name
        fluid = fluids[fluids["fluid_name"] == fluid_name].iloc[0]
        rho = float(fluid["rho_kg_m3"])
        mu = float(fluid["mu_Pa_s"])
        D_mol = float(fluid["D_mol_m2_s"])
        fluid_T_C = _parse_fluid_temp(fluid_name)
    else:
        st.info("No fluids in database – using water defaults.")
        fluid_name = ""
        rho, mu, D_mol = 997.0, 0.00089, 2.3e-9
        fluid_T_C = 25.0
# Show auto-match warning if present
if st.session_state.get("_fluid_auto_warn"):
    st.warning(st.session_state.pop("_fluid_auto_warn"))
with col2:
    if not reactions.empty:
        _rxn_list = reactions["reaction_name"].tolist()
        rxn_name = st.selectbox("Reaction (for Da numbers)", _rxn_list, index=_sel_idx(_rxn_list, "_sel_cmp_rxn"), key="cmp_rxn")
        st.session_state["_sel_cmp_rxn"] = rxn_name
        rxn = reactions[reactions["reaction_name"] == rxn_name].iloc[0]
        t_rxn = float(rxn["t_rxn_s"])
        rxn_k = float(rxn["k_value"])
        rxn_C0 = float(rxn["C0_mol_L"])
        rxn_order = str(rxn.get("order", "1"))
        rxn_T_C = _safe_float(rxn.get("T_C"), 25.0)
        rxn_delta_H = _safe_float(rxn.get("delta_H_kJ_mol"), 0.0)
        # Ensure _prev_cmp_rxn is initialised after selectbox renders,
        # so auto-selection only fires when user actually changes the reaction.
        if "_prev_cmp_rxn" not in st.session_state:
            st.session_state["_prev_cmp_rxn"] = rxn_name
        if t_rxn <= 0:
            k = rxn_k
            C0 = rxn_C0
            order = rxn_order
            if order in ("1", "pseudo-1") and k > 0:
                t_rxn = np.log(2) / k
            elif order in ("2", "pseudo-2") and k * C0 > 0:
                t_rxn = 1.0 / (k * C0)
            else:
                t_rxn = 1.0
    else:
        t_rxn = 1.0
        rxn_k, rxn_C0, rxn_order, rxn_T_C, rxn_delta_H = 0.0, 0.0, "1", 25.0, 0.0

col3, col4 = st.columns(2)
with col3:
    v_s = st.number_input("Superficial gas velocity v_s (m/s)", value=0.0, min_value=0.0,
                          format="%.4f", key="cmp_vs",
                          help="Set > 0 to compute kLa. Typical: 0.001 – 0.05 m/s.")
with col4:
    coal_choice = st.selectbox("Liquid type (for kLa)",
                               ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                               key="cmp_coal")
    is_coalescing = coal_choice.startswith("Coalescing")

st.divider()

# ── Per-reactor correlation mode selector ─────────────────────────────────
st.header("2 · Correlation Source Selection")
corr_modes: dict[str, dict[str, str] | str] = {}  # reactor_name → per-param dict
_any_has_alt = any(has_any_alt_correlations(n) for n in selected_names)
if _any_has_alt:
    with st.expander("Correlation source selection (per reactor)", expanded=False):
        st.caption(
            "Select **Literature**, **ROM**, or **Experimental** per parameter "
            "for each reactor independently."
        )
        for rname in selected_names:
            if has_any_alt_correlations(rname):
                st.subheader(rname)
                corr_modes[rname] = render_correlation_matrix(
                    rname, key_prefix=f"cmp_corr_{rname}",
                )
                st.divider()
            else:
                corr_modes[rname] = "Literature"
else:
    for rname in selected_names:
        corr_modes[rname] = "Literature"

st.divider()

# ── Additional options ─────────────────────────────────────────────────────────
st.header("3 · Additional Options")

# ── Particle options ──────────────────────────────────────────────────────────
include_particles = st.checkbox("Include solid particles", value=False, key="cmp_include_particles")
cmp_rho_p = cmp_d50_um = cmp_phi_p = cmp_X_wt = cmp_S_zw = 0.0
if include_particles and not particles_db.empty:
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        cmp_particle_name = st.selectbox("Particle", particles_db["particle_name"].tolist(), key="cmp_particle")
        cmp_part = particles_db[particles_db["particle_name"] == cmp_particle_name].iloc[0]
        cmp_rho_p = float(cmp_part["rho_p_kg_m3"])
        cmp_d50_um = float(cmp_part["d50_um"])
        cmp_phi_p = float(cmp_part["shape_factor"])
    with pcol2:
        cmp_X_wt = st.number_input("Solids loading X (wt-%)", value=5.0, min_value=0.01,
                                    format="%.2f", key="cmp_Xwt")
    with pcol3:
        cmp_S_zw = st.number_input("Zwietering S", value=5.5, min_value=0.5,
                                    max_value=20.0, format="%.1f", key="cmp_Szw")
elif include_particles and particles_db.empty:
    st.warning("Particle database is empty.")

# ── Heat balance options ──────────────────────────────────────────────────────
_heat_context = f"{fluid_name}|{rxn_name if not reactions.empty else ''}"

# Deactivate heat balance when fluid or reaction selection changes
if st.session_state.get("_heat_context") != _heat_context:
    if st.session_state.get("_heat_active", False):
        st.session_state["_heat_active"] = False
        st.session_state["_heat_stale"] = True
    st.session_state["_heat_context"] = _heat_context

if st.button("🔥 Run Heat Balance"):
    st.session_state["_heat_active"] = True
    st.session_state.pop("_heat_stale", None)
    # Default process temp from reaction temperature; coolant 20 °C below
    _default_T = rxn_T_C if rxn_T_C > 0 else fluid_T_C
    st.session_state["cmp_T_process"] = _default_T
    st.session_state["cmp_T_cool"] = _default_T - 20.0
    st.rerun()

if st.session_state.get("_heat_stale"):
    st.info("Fluid or reaction changed — click **🔥 Run Heat Balance** to update.")

include_heat = st.session_state.get("_heat_active", False)
cmp_T_process = fluid_T_C
cmp_T_coolant = fluid_T_C - 20.0
if include_heat:
    hcol1, hcol2, hcol3 = st.columns(3)
    with hcol1:
        cmp_T_process = st.number_input(
            "Process temperature (°C)",
            format="%.1f", key="cmp_T_process",
            help="Defaults to reaction temperature; adjust as needed",
            step=1.0)
    with hcol2:
        cmp_T_coolant = st.number_input(
            "Coolant temperature (°C)",
            format="%.1f", key="cmp_T_cool",
            help="Jacket coolant inlet temperature",
            step=1.0)
    with hcol3:
        st.markdown(f"**ΔH** = {rxn_delta_H:.0f} kJ/mol  ")
        st.markdown(f"**ΔT** = {cmp_T_process - cmp_T_coolant:.1f} °C")
    if rxn_delta_H == 0:
        st.info("Selected reaction has no enthalpy value – add ΔH in the Reaction Database.")


# ── Compute 4 corners per reactor ────────────────────────────────────────

def _liquid_height(V_litres: float, D_tank: float, H_max: float,
                   dish_type: str = "") -> float:
    """Estimate liquid height for a given fill volume, accounting for bottom dish."""
    return liquid_height_from_volume(V_litres, D_tank, H_max, dish_type)


CORNER_LABELS = [
    "min RPM / max V",
    "max RPM / max V",
    "min RPM / min V",
    "max RPM / min V",
]

envelope_rows: list[dict] = []  # one row per (reactor, corner)
skipped: list[str] = []
reactor_info: dict = {}         # stash geometry for curve precomputation

for rname in selected_names:
    r = reactors[reactors["reactor_name"] == rname].iloc[0]
    D_imp_v  = _safe_float(r.get("D_imp_m"))
    D_tank_v = _safe_float(r.get("D_tank_m"))
    H_max    = _safe_float(r.get("H_m"))
    Np_v     = _safe_float(r.get("Np"), None)
    Nq_v     = _safe_float(r.get("Nq"), None)
    scale    = r.get("scale", "")

    if D_imp_v == 0 or D_tank_v == 0 or H_max == 0:
        skipped.append(rname)
        continue

    # RPM range → rev/s
    rpm_min = _safe_float(r.get("N_rpm_min"))
    rpm_max = _safe_float(r.get("N_rpm_max"))
    n_rps_default = _safe_float(r.get("N_rps"))
    if rpm_max == 0 and n_rps_default > 0:
        rpm_max = n_rps_default * 60
        rpm_min = rpm_max          # single point
    elif rpm_max == 0:
        skipped.append(rname)
        continue
    if rpm_min == 0:
        rpm_min = rpm_max          # single point if min missing

    N_lo = rpm_min / 60.0
    N_hi = rpm_max / 60.0

    # Volume range
    V_max_L = _safe_float(r.get("V_L_max"))
    V_min_L = _safe_float(r.get("V_L_min"))
    V_default = _safe_float(r.get("V_L"))
    # Fallback: compute from geometry
    V_geo = np.pi / 4 * D_tank_v**2 * H_max * 1000  # litres
    if V_max_L == 0:
        V_max_L = V_default if V_default > 0 else V_geo
    if V_min_L == 0:
        V_min_L = V_max_L       # single volume if missing

    reactor_info[rname] = dict(
        D_imp=D_imp_v, D_tank=D_tank_v, H_max=H_max,
        Np=Np_v, Nq=Nq_v,
        N_lo=N_lo, N_hi=N_hi,
        V_max_L=V_max_L, V_min_L=V_min_L,
        rpm_max=rpm_max,
    )

    # Heat-transfer geometry for this reactor
    _r_material = str(r.get("material", ""))
    _r_bottom_dish = str(r.get("bottom_dish", "")) if pd.notna(r.get("bottom_dish")) else ""
    _r_U_override = _safe_float(r.get("U_W_m2K"), 0.0)
    _r_A_override = _safe_float(r.get("A_ht_m2"), 0.0)
    _r_wall_mm = _safe_float(r.get("wall_thickness_mm"), 0.0)
    reactor_info[rname]["material"] = _r_material
    reactor_info[rname]["bottom_dish"] = _r_bottom_dish
    reactor_info[rname]["U_override"] = _r_U_override
    reactor_info[rname]["A_override"] = _r_A_override
    reactor_info[rname]["wall_thickness_mm"] = _r_wall_mm

    corners = [
        ("min RPM / max V", N_lo, V_max_L),
        ("max RPM / max V", N_hi, V_max_L),
        ("min RPM / min V", N_lo, V_min_L),
        ("max RPM / min V", N_hi, V_min_L),
    ]

    for label, N, V_L in corners:
        H_v = _liquid_height(V_L, D_tank_v, H_max, _r_bottom_dish)
        h, _rom_src = compute_reactor_hydro_with_mode(
            corr_modes.get(rname, "Literature"), rname,
            N=N, D_imp=D_imp_v, D_tank=D_tank_v, H=H_v,
            rho=rho, mu=mu, Np=Np_v, Nq=Nq_v,
            v_s=v_s, coalescing=is_coalescing,
            D_mol=D_mol,
        )
        da = compute_damkohler_numbers(
            h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
            kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
        )
        # Particle parameters (if enabled)
        part_vals = {}
        if include_particles and cmp_d50_um > 0:
            _dp = cmp_d50_um * 1e-6
            _nu = mu / rho
            _drho = abs(cmp_rho_p - rho)
            _vt = settling_velocity(_dp, cmp_rho_p, rho, mu, cmp_phi_p)
            _rep = particle_reynolds(_dp, _vt, rho, mu)
            _njs = zwietering_njs(cmp_S_zw, _nu, _dp, _drho, rho, cmp_X_wt, D_imp_v)
            _ksl = solid_liquid_mass_transfer(_dp, _vt, rho, mu, D_mol)
            _njs_rpm = _njs * 60
            _n_over_njs = N / _njs if _njs > 0 else 0.0
            part_vals = {
                "N_js (RPM)": _njs_rpm,
                "N/N_js": _n_over_njs,
                "v_t (m/s)": _vt,
                "Re_p": _rep,
                "k_SL (m/s)": _ksl,
            }
        # Heat balance parameters (if enabled)
        heat_vals = {}
        if include_heat and rxn_delta_H != 0:
            _r_mol_s = reaction_rate_mol_per_s(rxn_order, rxn_k, rxn_C0, V_L)
            _Q_gen = heat_generation_rate(rxn_delta_H, _r_mol_s)
            _A_ht = _r_A_override if _r_A_override > 0 else estimate_jacket_area(D_tank_v, H_v, _r_bottom_dish)
            if _r_U_override > 0:
                _U_ht = _r_U_override
            else:
                _U_ht, _ = estimate_U_detailed(
                    N_rps=N, D_imp=D_imp_v, D_tank=D_tank_v,
                    rho=rho, mu=mu,
                    material=_r_material,
                    wall_thickness_mm=_r_wall_mm,
                    fluid_name=fluid_name,
                )
            _dT = cmp_T_process - cmp_T_coolant
            _Q_cool = heat_removal_capacity(_U_ht, _A_ht, _dT)
            heat_vals = {
                "Q_gen (W)": _Q_gen,
                "Q_cool (W)": _Q_cool,
                "U (W/m²·K)": _U_ht,
                "A_ht (m²)": _A_ht,
                "Q_gen/Q_cool (%)": _Q_gen / _Q_cool * 100 if _Q_cool > 0 else np.inf,
            }
        envelope_rows.append({
            "Reactor": rname,
            "Scale": scale,
            "Corner": label,
            "N (rev/s)": N,
            "RPM": N * 60,
            "RPM_max": rpm_max,
            "V_L": V_L,
            **h,
            **da,
            **part_vals,
            **heat_vals,
        })

if skipped:
    st.warning(f"Skipped reactors with missing geometry/agitation data: {', '.join(skipped)}")

# Show U estimation warnings (once per reactor, at max RPM)
if include_heat and rxn_delta_H != 0:
    _u_warnings_all: list[str] = []
    for rname, info in reactor_info.items():
        if info["U_override"] > 0:
            continue
        _, _u_warns = estimate_U_detailed(
            N_rps=info["N_hi"], D_imp=info["D_imp"], D_tank=info["D_tank"],
            rho=rho, mu=mu,
            material=info["material"],
            wall_thickness_mm=info["wall_thickness_mm"],
            fluid_name=fluid_name,
        )
        if _u_warns:
            _u_warnings_all.append(f"**{rname}**: " + "; ".join(_u_warns))
    if _u_warnings_all:
        with st.expander("ℹ️ U estimation notes", expanded=False):
            for _w in _u_warnings_all:
                st.markdown(f"- {_w}")

if not envelope_rows:
    st.info("No computable reactors in the selection.")
    st.stop()

env_df = pd.DataFrame(envelope_rows)
env_df["RPM_pct"] = env_df["RPM"] / env_df["RPM_max"] * 100

# ── Aggregate min / max per reactor for plotting ─────────────────────────
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
if include_particles and cmp_d50_um > 0:
    PLOT_PARAMS = PLOT_PARAMS + PARTICLE_PARAMS

agg_dict = {p: ["min", "max"] for p in PLOT_PARAMS}
agg_dict["Scale"] = "first"
agg_dict["Volume (L)"] = ["min", "max"]
agg_df = env_df.groupby("Reactor", sort=False).agg(agg_dict)
agg_df.columns = ["_".join(c).strip("_") for c in agg_df.columns]
agg_df = agg_df.reset_index()

# ── Precompute true boundary curves at many RPM points ───────────────────
_N_INTERP = 50

curve_data: dict = {}  # rname → {pct_arr, maxV: {param: arr}, minV: {param: arr}}

# Pre-compute RPM-independent particle quantities (hoisted out of inner loop)
_p7_part_static: dict[str, dict] = {}
if include_particles and cmp_d50_um > 0:
    _dp_p7 = cmp_d50_um * 1e-6
    _nu_p7 = mu / rho
    _drho_p7 = abs(cmp_rho_p - rho)
    _vt_p7 = settling_velocity(_dp_p7, cmp_rho_p, rho, mu, cmp_phi_p)
    _rep_p7 = particle_reynolds(_dp_p7, _vt_p7, rho, mu)
    for rname, info in reactor_info.items():
        _njs_p7 = zwietering_njs(cmp_S_zw, _nu_p7, _dp_p7, _drho_p7, rho, cmp_X_wt, info["D_imp"])
        _ksl_p7 = solid_liquid_mass_transfer(_dp_p7, _vt_p7, rho, mu, D_mol)
        _p7_part_static[rname] = {
            "N_js_rps": _njs_p7, "N_js (RPM)": _njs_p7 * 60,
            "v_t (m/s)": _vt_p7, "Re_p": _rep_p7, "k_SL (m/s)": _ksl_p7,
        }

for rname, info in reactor_info.items():
    N_arr = np.linspace(info["N_lo"], info["N_hi"], _N_INTERP)
    pct_arr = N_arr / info["N_hi"] * 100 if info["N_hi"] > 0 else np.zeros(_N_INTERP)

    curves: dict = {"pct_arr": pct_arr}
    for vol_key, V_L in [("maxV", info["V_max_L"]), ("minV", info["V_min_L"])]:
        H_v = _liquid_height(V_L, info["D_tank"], info["H_max"], info["bottom_dish"])
        param_arrs: dict = {p: np.empty(_N_INTERP) for p in PLOT_PARAMS}

        # Pre-compute RPM-independent heat quantities for this volume
        _Q_gen_p7 = _A_ht_p7 = 0.0
        if include_heat and rxn_delta_H != 0:
            _r_mol_s = reaction_rate_mol_per_s(rxn_order, rxn_k, rxn_C0, V_L)
            _Q_gen_p7 = heat_generation_rate(rxn_delta_H, _r_mol_s)
            _A_ht_p7 = info["A_override"] if info["A_override"] > 0 else estimate_jacket_area(info["D_tank"], H_v, info["bottom_dish"])

        for j, N in enumerate(N_arr):
            h, _ = compute_reactor_hydro_with_mode(
                corr_modes.get(rname, "Literature"), rname,
                N=N, D_imp=info["D_imp"], D_tank=info["D_tank"], H=H_v,
                rho=rho, mu=mu, Np=info["Np"], Nq=info["Nq"],
                v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
            )
            da = compute_damkohler_numbers(
                h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
                kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
            )
            vals = {**h, **da}
            # Particle parameters — use pre-computed statics
            if rname in _p7_part_static:
                _ps = _p7_part_static[rname]
                vals["N_js (RPM)"] = _ps["N_js (RPM)"]
                vals["N/N_js"] = N / _ps["N_js_rps"] if _ps["N_js_rps"] > 0 else 0.0
                vals["v_t (m/s)"] = _ps["v_t (m/s)"]
                vals["Re_p"] = _ps["Re_p"]
                vals["k_SL (m/s)"] = _ps["k_SL (m/s)"]
            # Heat balance — only U depends on RPM
            if include_heat and rxn_delta_H != 0:
                if info["U_override"] > 0:
                    _U_ht = info["U_override"]
                else:
                    _U_ht, _ = estimate_U_detailed(
                        N_rps=N, D_imp=info["D_imp"], D_tank=info["D_tank"],
                        rho=rho, mu=mu,
                        material=info["material"],
                        wall_thickness_mm=info["wall_thickness_mm"],
                        fluid_name=fluid_name,
                    )
                _dT = cmp_T_process - cmp_T_coolant
                _Q_cool = heat_removal_capacity(_U_ht, _A_ht_p7, _dT)
                vals["Q_gen (W)"] = _Q_gen_p7
                vals["Q_cool (W)"] = _Q_cool
                vals["U (W/m²·K)"] = _U_ht
                vals["A_ht (m²)"] = _A_ht_p7
                vals["Q_gen/Q_cool (%)"] = _Q_gen_p7 / _Q_cool * 100 if _Q_cool > 0 else np.inf
            for p in PLOT_PARAMS:
                param_arrs[p][j] = vals.get(p, np.nan)
        curves[vol_key] = param_arrs
    curve_data[rname] = curves

st.divider()

# ── Summary Table ─────────────────────────────────────────────────────────
st.header("4 · Operating Envelope Summary")
st.caption("Each row shows the range across the 4 corner conditions (min/max RPM × min/max volume).")

table_rows = []
for _, a in agg_df.iterrows():
    row: dict = {"Reactor": a["Reactor"], "Scale": a["Scale_first"]}
    row["Volume (L)"] = f"{a['Volume (L)_min']:.1f} – {a['Volume (L)_max']:.1f}"
    for p in PLOT_PARAMS:
        lo, hi = a[f"{p}_min"], a[f"{p}_max"]
        if abs(lo - hi) < 1e-12:
            row[p] = f"{lo:.3g}"
        else:
            row[p] = f"{lo:.3g} – {hi:.3g}"
    table_rows.append(row)

st.dataframe(pd.DataFrame(table_rows), width="content", hide_index=True)

# ── Detail: 4-corner table ───────────────────────────────────────────────
with st.expander("Full 4-corner detail table", expanded=False):
    detail_cols = ["Reactor", "Corner", "N (rev/s)", "V_L", "Re",
                   "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
                   "Micromix time t_E (s)", "Micromix time t_E_local (s)",
                   "Kolmogorov η (µm)",
                   "Avg shear rate (1/s)", "Max shear rate (1/s)",
                   "Avg shear stress (Pa)", "kLa (1/s)", "kLa_surface (1/s)",
                   "Da_macro", "Da_micro", "Da_GL"]
    if include_particles and cmp_d50_um > 0:
        detail_cols += PARTICLE_PARAMS
    if include_heat and rxn_delta_H != 0:
        detail_cols += HEAT_PARAMS
    # Only include columns that exist in the dataframe
    detail_cols = [c for c in detail_cols if c in env_df.columns]
    fmt = {c: "{:.3g}" for c in detail_cols if c not in ("Reactor", "Corner")}
    st.dataframe(env_df[detail_cols].style.format(fmt), width="content", hide_index=True)

st.divider()

# ── Charts: operating envelopes ──────────────────────────────────────────
st.header("5 · Operating Envelope Charts")

with st.expander("Show / hide envelope charts", expanded=True):

    # RPM-to-percentage mapping table
    st.subheader("Stir Speed Reference Table")
    st.caption("Translates percentage of max RPM (chart x-axis) to actual RPM for each reactor.")

    pct_steps = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    rpm_table_rows = []
    for rname in agg_df["Reactor"].tolist():
        rc = env_df[env_df["Reactor"] == rname].iloc[0]
        rpm_max_val = rc["RPM_max"]
        rpm_min_val = rc["RPM"] if rc["Corner"] == "min RPM / max V" else env_df[(env_df["Reactor"] == rname) & (env_df["Corner"] == "min RPM / max V")].iloc[0]["RPM"]
        pct_min = rpm_min_val / rpm_max_val * 100
        row_data: dict = {"Reactor": rname, "RPM min": f"{rpm_min_val:.0f}", "RPM max": f"{rpm_max_val:.0f}",
                           "% min": f"{pct_min:.0f}%"}
        for pct in pct_steps:
            row_data[f"{pct}%"] = f"{rpm_max_val * pct / 100:.0f}"
        rpm_table_rows.append(row_data)

    st.dataframe(pd.DataFrame(rpm_table_rows), width="content", hide_index=True)

    st.caption(
        "Each reactor's operational region is plotted as a filled polygon. "
        "The x-axis is stir speed as a percentage of each vessel's maximum RPM. "
        "Overlapping regions indicate where reactors share comparable conditions."
    )

    # Friendly display names for selected parameters
    _DISPLAY_NAMES: dict[str, str] = {
        "Da_macro": "Macromixing (Da_macro)",
        "Da_micro": "Micromixing (Da_micro)",
        "Da_GL": "Gas-Liquid Mass Transfer (Da_GL)",
        "Q_gen/Q_cool (%)": "Heat Transfer Capacity (Q_gen/Q_cool (%))",
    }
    _display = lambda p: _DISPLAY_NAMES.get(p, p)

    _DEFAULT_PARAMS = ["Da_micro", "Da_macro", "Da_GL", "Q_gen/Q_cool (%)", "Blend time 95% (s)", "P/V (W/L)"]
    _defaults = [p for p in _DEFAULT_PARAMS if p in PLOT_PARAMS]

    params_to_plot = st.multiselect(
        "Parameters to plot",
        PLOT_PARAMS,
        default=_defaults,
        key="env_params",
        format_func=_display,
    )

    # Extended colour palette (one colour per reactor)
    _PALETTE = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf",
    ]

    _N_INTERP = 50  # (also defined earlier for precomputation)

    for param in params_to_plot:
        fig = go.Figure()

        reactor_list = agg_df["Reactor"].tolist()
        for i, rname in enumerate(reactor_list):
            color = _PALETTE[i % len(_PALETTE)]
            corners = env_df[env_df["Reactor"] == rname]

            # 4 corner values (for markers)
            c = corners.set_index("Corner")
            pct_lo = c.loc["min RPM / max V", "RPM_pct"]
            pct_hi = 100.0

            val_minN_maxV = c.loc["min RPM / max V", param]
            val_maxN_maxV = c.loc["max RPM / max V", param]
            val_maxN_minV = c.loc["max RPM / min V", param]
            val_minN_minV = c.loc["min RPM / min V", param]

            # True boundary curves from precomputed data
            curves = curve_data[rname]
            pct_arr = curves["pct_arr"]
            y_maxV = curves["maxV"][param]
            y_minV = curves["minV"][param]

            # Build polygon from actual curve points (follows nonlinear shapes)
            poly_x = np.concatenate([pct_arr, pct_arr[::-1], [pct_arr[0]]])
            poly_y = np.concatenate([y_maxV, y_minV[::-1], [y_maxV[0]]])

            # Filled polygon (shaded region)
            fig.add_trace(go.Scatter(
                x=poly_x,
                y=poly_y,
                fill="toself",
                fillcolor=color,
                opacity=0.20,
                line=dict(color=color, width=1),
                mode="lines",
                name=rname,
                showlegend=True,
                legendgroup=rname,
                hoverinfo="skip",
            ))

            # Corner markers with hover detail
            marker_x = [pct_lo, pct_hi, pct_lo, pct_hi]
            marker_y = [val_minN_maxV, val_maxN_maxV, val_minN_minV, val_maxN_minV]
            marker_labels = ["min RPM / max V", "max RPM / max V",
                             "min RPM / min V", "max RPM / min V"]
            marker_symbols = ["diamond", "circle", "square", "triangle-up"]

            fig.add_trace(go.Scatter(
                x=marker_x,
                y=marker_y,
                mode="markers",
                marker=dict(
                    size=10,
                    color=color,
                    symbol=marker_symbols,
                    line=dict(width=1, color="white"),
                ),
                text=marker_labels,
                customdata=[[rname]] * 4,
                hovertemplate=(
                    "%{customdata[0]}<br>"
                    "%{text}<br>"
                    "RPM %%max = %{x:.0f}%%<br>"
                    + param + " = %{y:.3g}<extra></extra>"
                ),
                showlegend=False,
                legendgroup=rname,
            ))

            # Max-volume boundary (solid line)
            fig.add_trace(go.Scatter(
                x=pct_arr, y=y_maxV,
                mode="lines", line=dict(color=color, width=2),
                name=f"{rname} (max V)",
                showlegend=False, legendgroup=rname,
                hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.3g}<extra>" + rname + " max V</extra>",
            ))
            # Min-volume boundary (dotted line)
            fig.add_trace(go.Scatter(
                x=pct_arr, y=y_minV,
                mode="lines", line=dict(color=color, width=2, dash="dot"),
                name=f"{rname} (min V)",
                showlegend=False, legendgroup=rname,
                hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.3g}<extra>" + rname + " min V</extra>",
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

        # Reference line for N/N_js = 1 (just-suspended)
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

        # Reference line for Q_gen/Q_cool = 100% (insufficient cooling)
        if param == "Q_gen/Q_cool (%)":
            fig.add_shape(
                type="line", x0=0, x1=1, y0=100.0, y1=100.0,
                xref="paper", yref="y",
                line=dict(color="red", width=1.5, dash="dash"),
            )
            fig.add_annotation(
                x=1.0, xref="paper", xanchor="right",
                y=100.0, yref="y", yanchor="bottom", yshift=2,
                text="Q_gen/Q_cool = 100% (cooling limit)",
                font=dict(size=11, color="red"), showarrow=False,
            )

        _param_label = _display(param)
        fig.update_layout(
            title=_param_label,
            xaxis_title="Stir speed (% of vessel max RPM)",
            yaxis_title=_param_label,
            xaxis=dict(range=[0, 105], dtick=10,
                       showspikes=True, spikemode="across",
                       spikethickness=1, spikecolor="grey", spikedash="dot"),
            yaxis=dict(showspikes=True, spikemode="across",
                       spikethickness=1, spikecolor="grey", spikedash="dot"),
            height=500,
            margin=dict(t=50, b=50),
            legend=dict(title="Reactor"),
            hovermode="closest",
        )
        st.plotly_chart(fig, width="content")

    # ── Chart legend ─────────────────────────────────────────────────────
    with st.expander("Chart legend"):
        st.markdown("""
**Regions:** Each reactor's operational envelope is drawn as a filled polygon
spanning its RPM range (as % of max) on the x-axis.

**Boundary lines:**
- **Solid line** = max fill volume edge
- **Dotted line** = min fill volume edge

**Corner markers:**

| Symbol | Corner condition |
|--------|-----------------|
| ◆ Diamond | min RPM / max Volume |
| ● Circle | max RPM / max Volume |
| ■ Square | min RPM / min Volume |
| ▲ Triangle | max RPM / min Volume |

Overlapping shaded regions indicate where two reactors can achieve similar parameter values.
""")

# ── Heat Balance Summary ─────────────────────────────────────────────────
if include_heat and rxn_delta_H != 0:
    st.divider()
    st.header("6 · Heat Balance Summary")
    st.caption(
        f"Reaction: **{rxn_name if not reactions.empty else 'N/A'}** | "
        f"ΔH = {rxn_delta_H:.0f} kJ/mol | "
        f"T_process = {cmp_T_process:.0f} °C | T_coolant = {cmp_T_coolant:.0f} °C | "
        f"ΔT = {cmp_T_process - cmp_T_coolant:.1f} °C"
    )

    # Build summary table: one row per reactor at max-RPM / max-V corner
    heat_summary_rows = []
    for rname in reactor_info:
        corners = env_df[(env_df["Reactor"] == rname) & (env_df["Corner"] == "max RPM / max V")]
        if corners.empty:
            continue
        c = corners.iloc[0]
        Q_gen = c.get("Q_gen (W)", 0)
        Q_cool = c.get("Q_cool (W)", 0)
        U_val = c.get("U (W/m²·K)", 0)
        A_val = c.get("A_ht (m²)", 0)
        ratio_pct = Q_gen / Q_cool * 100 if Q_cool > 0 else np.inf
        heat_summary_rows.append({
            "Reactor": rname,
            "Volume (L)": c["V_L"],
            "U (W/m²·K)": f"{U_val:.0f}",
            "A (m²)": f"{A_val:.3f}",
            "Q_gen (W)": f"{Q_gen:.1f}",
            "Q_cool (W)": f"{Q_cool:.1f}",
            "Q_gen/Q_cool (%)": f"{ratio_pct:.1f}%" if ratio_pct < 10000 else "∞",
            "Assessment": heat_balance_assessment(Q_gen, Q_cool),
        })

    if heat_summary_rows:
        st.dataframe(pd.DataFrame(heat_summary_rows), width="content", hide_index=True)

    # Operating envelope chart for Q_gen/Q_cool (%)
    if heat_summary_rows:
        _heat_param = "Q_gen/Q_cool (%)"
        fig_heat = go.Figure()

        for i, rname in enumerate(reactor_info):
            if rname not in curve_data:
                continue
            color = _PALETTE[i % len(_PALETTE)]
            curves = curve_data[rname]
            pct_arr = curves["pct_arr"]
            y_maxV = curves["maxV"][_heat_param]
            y_minV = curves["minV"][_heat_param]

            # Filled polygon
            poly_x = np.concatenate([pct_arr, pct_arr[::-1], [pct_arr[0]]])
            poly_y = np.concatenate([y_maxV, y_minV[::-1], [y_maxV[0]]])
            fig_heat.add_trace(go.Scatter(
                x=poly_x, y=poly_y,
                fill="toself", fillcolor=color, opacity=0.20,
                line=dict(color=color, width=1), mode="lines",
                name=rname, showlegend=True, legendgroup=rname,
                hoverinfo="skip",
            ))
            # Max-volume boundary (solid)
            fig_heat.add_trace(go.Scatter(
                x=pct_arr, y=y_maxV,
                mode="lines", line=dict(color=color, width=2),
                showlegend=False, legendgroup=rname,
                hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.1f}%<extra>" + rname + " max V</extra>",
            ))
            # Min-volume boundary (dotted)
            fig_heat.add_trace(go.Scatter(
                x=pct_arr, y=y_minV,
                mode="lines", line=dict(color=color, width=2, dash="dot"),
                showlegend=False, legendgroup=rname,
                hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.1f}%<extra>" + rname + " min V</extra>",
            ))

        # 100% cooling limit reference line
        fig_heat.add_shape(
            type="line", x0=0, x1=1, y0=100.0, y1=100.0,
            xref="paper", yref="y",
            line=dict(color="red", width=1.5, dash="dash"),
        )
        fig_heat.add_annotation(
            x=1.0, xref="paper", xanchor="right",
            y=100.0, yref="y", yanchor="bottom", yshift=2,
            text="100% (cooling limit)",
            font=dict(size=11, color="red"), showarrow=False,
        )
        fig_heat.update_layout(
            title="Heat Transfer Capacity – Operating Envelope",
            xaxis_title="Stir speed (% of vessel max RPM)",
            yaxis_title="Q_gen / Q_cool (%)",
            xaxis=dict(range=[0, 105], dtick=10,
                       showspikes=True, spikemode="across",
                       spikethickness=1, spikecolor="grey", spikedash="dot"),
            yaxis=dict(ticksuffix="%",
                       showspikes=True, spikemode="across",
                       spikethickness=1, spikecolor="grey", spikedash="dot"),
            height=500,
            margin=dict(t=50, b=50),
            legend=dict(title="Reactor"),
            hovermode="closest",
        )
        st.plotly_chart(fig_heat, width="content")

st.divider()

# ── Scale-up summary ─────────────────────────────────────────────────────
st.header("7 · Scale-Up Impact Summary")
st.caption("Ratios use midpoint (average of 4 corners) for each parameter, relative to the first selected reactor.")

if len(agg_df) >= 2:
    # Compute midpoints per reactor
    mid_df = env_df.groupby("Reactor", sort=False)[PLOT_PARAMS + ["Volume (L)"]].mean().reset_index()

    _has_heat = include_heat and rxn_delta_H != 0 and "Q_gen/Q_cool (%)" in mid_df.columns

    ref = mid_df.iloc[0]
    summary_rows = []
    for _, row in mid_df.iloc[1:].iterrows():
        entry: dict = {
            "From → To": f"{ref['Reactor']} → {row['Reactor']}",
            "Volume ratio": row["Volume (L)"] / ref["Volume (L)"] if ref["Volume (L)"] > 0 else np.nan,
            "P/V ratio": row["P/V (W/L)"] / ref["P/V (W/L)"] if ref["P/V (W/L)"] > 0 else np.nan,
            "Tip speed ratio": row["Tip speed (m/s)"] / ref["Tip speed (m/s)"] if ref["Tip speed (m/s)"] > 0 else np.nan,
            "Blend time ratio": row["Blend time 95% (s)"] / ref["Blend time 95% (s)"] if ref["Blend time 95% (s)"] > 0 else np.nan,
            "Da_macro ratio": row["Da_macro"] / ref["Da_macro"] if ref["Da_macro"] > 0 else np.nan,
        }
        if _has_heat:
            ref_pct = ref["Q_gen/Q_cool (%)"]
            row_pct = row["Q_gen/Q_cool (%)"]
            delta = row_pct - ref_pct
            if ref_pct > 0 and np.isfinite(ref_pct) and np.isfinite(row_pct):
                if abs(delta) < 1.0:
                    assessment = "≈ Similar"
                elif delta < 0:
                    assessment = f"✅ Improves ({delta:+.1f} pp)"
                else:
                    assessment = f"⚠️ Worse ({delta:+.1f} pp)"
                entry["Q_gen/Q_cool (%)"] = f"{ref_pct:.1f}% → {row_pct:.1f}%"
            else:
                assessment = "N/A"
                entry["Q_gen/Q_cool (%)"] = "N/A"
            entry["Cooling"] = assessment
        summary_rows.append(entry)
    summary_df = pd.DataFrame(summary_rows)
    numeric_cols = summary_df.select_dtypes(include="number").columns
    st.dataframe(summary_df.style.format({c: "{:.2f}" for c in numeric_cols}), width="content")
else:
    st.info("Select at least two reactors to see scale-up ratios.")
