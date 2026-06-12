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
import pathlib
import plotly.graph_objects as go

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, is_known_solvent, resolve_solvent_name,
)
from utils.calculations import (
    compute_reactor_hydro, compute_damkohler_numbers,
    damkohler_sl, solid_liquid_kla,
    settling_velocity, particle_reynolds, zwietering_njs,
    solid_liquid_mass_transfer, particle_suspension_criterion,
    reaction_rate_mol_per_s, heat_generation_rate,
    estimate_jacket_area, estimate_U, estimate_U_detailed,
    heat_removal_capacity, heat_balance_assessment,
    liquid_height_from_volume,
    mesomixing_time, cooling_rate, time_to_cool_or_heat,
    gmb_njs,
)
from utils.rom_registry import (
    compute_reactor_hydro_with_mode,
    get_correlations,
    has_any_alt_correlations,
    PARAM_DISPLAY,
)
from utils.corr_widgets import render_correlation_matrix
from utils.data_helpers import safe_iloc, load_db, safe_float, find_reactor_image, reactor_search_name, build_search_names

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── Colour palette by scale ──────────────────────────────────────────────
SCALE_COLORS = {"Lab": "#1f77b4", "Pilot": "#ff7f0e", "Manufacturing": "#2ca02c"}


def _load(key, fn):
    # Delegate to the shared cached loader for consistency across pages.
    return load_db(key, fn)


def _safe_float(val, default=0.0):
    try:
        v = float(val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


reactors = _load("reactor_db", "reactors.csv")
custom_fluids = _load("fluid_db", "fluids.csv")
reactions = _load("reaction_db", "reactions.csv")
particles_db = _load("particle_db", "particles.csv")

# Build combined fluid list: built-in solvents + custom fluids
_solvent_names = sorted(SOLVENT_DB.keys())
_custom_names = custom_fluids["fluid_name"].tolist() if not custom_fluids.empty else []
_all_fluid_names = _solvent_names + _custom_names

st.title("📈 Reactor Comparison")

if reactors.empty:
    st.warning("Populate the Reactor Database first.")
    st.stop()

# ── Step gating: track confirmed input sections ──────────────────────────
if "_p7_step" not in st.session_state:
    st.session_state["_p7_step"] = 0

def _reset_p7():
    st.session_state["_p7_step"] = 0

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
    on_change=_reset_p7,
    format_func=lambda n: reactor_search_name(reactors, n),
)
st.session_state["_sel_cmp_reactors"] = selected_names

if not selected_names:
    st.info("Select at least one reactor above.")
    st.stop()

# -- Auto-select fluid & temperature when reaction changes -------------------
if not reactions.empty and _all_fluid_names:
    _prev_rxn = st.session_state.get("_prev_cmp_rxn")
    _cur_rxn = st.session_state.get("cmp_rxn", None)
    if _cur_rxn is not None and _cur_rxn != _prev_rxn:
        st.session_state["_prev_cmp_rxn"] = _cur_rxn
        _rxn_row = reactions[reactions["reaction_name"] == _cur_rxn]
        if not _rxn_row.empty:
            _rxn_solvent = str(_rxn_row.iloc[0].get("solvent", ""))
            _rxn_T = _safe_float(_rxn_row.iloc[0].get("T_C"), 25.0)
            _resolved = resolve_solvent_name(_rxn_solvent)
            if _resolved and _resolved in _all_fluid_names:
                if _resolved != st.session_state.get("cmp_fluid"):
                    st.session_state["cmp_fluid"] = _resolved
                    st.session_state["_sel_cmp_fluid"] = _resolved
                # Also set temperature from reaction
                st.session_state["cmp_fluid_T"] = _rxn_T
                st.rerun()

def _sel_idx(lst, key, default=0):
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col1, col2 = st.columns(2)
col3, col4 = st.columns(2)
with col1:
    fluid_name = st.selectbox("Fluid", _all_fluid_names,
                              index=_sel_idx(_all_fluid_names, "_sel_cmp_fluid"),
                              key="cmp_fluid", on_change=_reset_p7)
    st.session_state["_sel_cmp_fluid"] = fluid_name

    _is_solvent = is_known_solvent(fluid_name)
    if _is_solvent:
        fluid_T_C = col3.number_input("Temperature (°C)", value=25.0, step=1.0,
                                     format="%.1f", key="cmp_fluid_T")
        fluid_P_atm = col4.number_input("Pressure (atm)", value=1.0, min_value=0.01,
                                       step=0.1, format="%.2f", key="cmp_fluid_P")
        _fprops = get_properties(fluid_name, fluid_T_C, fluid_P_atm)
        rho = _fprops["rho_kg_m3"]
        mu = _fprops["mu_Pa_s"]
        D_mol = _fprops["D_mol_m2_s"]
        if not _fprops["in_range"]:
            st.warning(f"⚠️ {fluid_T_C:.1f} °C is outside the liquid range "
                       f"({_fprops['mp_C']:.0f} – {_fprops['bp_at_P_C']:.0f} °C) for {fluid_name}.")
    else:
        _cust = safe_iloc(custom_fluids, "fluid_name", fluid_name, "Custom fluid")
        rho = float(_cust["rho_kg_m3"])
        mu = float(_cust["mu_Pa_s"])
        D_mol = float(_cust["D_mol_m2_s"])
        fluid_T_C = 25.0
        st.caption("Custom fluid — fixed properties")

with col2:
    if not reactions.empty:
        _rxn_list = reactions["reaction_name"].tolist()
        rxn_name = st.selectbox("Reaction (for Da numbers)", _rxn_list, index=_sel_idx(_rxn_list, "_sel_cmp_rxn"), key="cmp_rxn", on_change=_reset_p7)
        st.session_state["_sel_cmp_rxn"] = rxn_name
        rxn = safe_iloc(reactions, "reaction_name", rxn_name, "Reaction")
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

col5, col6, col7 = st.columns(3)
with col5:
    v_s = st.number_input("Superficial gas velocity v_s (m/s)", value=0.0, min_value=0.0,
                          format="%.4f", key="cmp_vs",
                          help="Set > 0 to compute kLa. Typical: 0.001 – 0.05 m/s.")
with col6:
    coal_choice = st.selectbox("Liquid type (for kLa)",
                               ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                               key="cmp_coal")
    is_coalescing = coal_choice.startswith("Coalescing")

with col7:
    cmp_T_coolant = st.number_input(
        "Coolant temperature (°C)", value=15.0, step=1.0,
        format="%.1f", key="cmp_T_cool",
        help="Jacket coolant inlet temperature",
    )

# ── Gate 1: confirm system selection ─────────────────────────────────────
if st.session_state["_p7_step"] < 1:
    st.info("👆 Review the selections above, then confirm to continue.")
    if st.button("✅ Confirm selections", key="p7_gate1", type="primary"):
        st.session_state["_p7_step"] = 1
        st.rerun()
    st.stop()

# ── Show iso images of selected reactors ─────────────────────────────────

_iso_imgs = [(rname, find_reactor_image(reactors, rname, "iso")) for rname in selected_names]
_iso_imgs = [(rname, p) for rname, p in _iso_imgs if p is not None]
if _iso_imgs:
    _MAX_PER_ROW = 3
    for _row_start in range(0, len(_iso_imgs), _MAX_PER_ROW):
        _row_slice = _iso_imgs[_row_start : _row_start + _MAX_PER_ROW]
        cols = st.columns(_MAX_PER_ROW)
        for idx, (rname, img_path) in enumerate(_row_slice):
            with cols[idx]:
                st.image(str(img_path), caption=rname, width=300)

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
    st.info("All selected reactors use **Literature** correlations only. "
            "No ROM or Experimental correlations are registered for these reactors.")

# ── Gate 2: confirm correlations ────────────────────────────────────────
if st.session_state["_p7_step"] < 2:
    st.info("👆 Review the correlation sources above, then confirm to continue.")
    if st.button("✅ Confirm correlations", key="p7_gate2", type="primary"):
        st.session_state["_p7_step"] = 2
        st.rerun()
    st.stop()

st.divider()

# ── Additional options ─────────────────────────────────────────────────────────
st.header("3 · Additional Options")

# ── Particle options ──────────────────────────────────────────────────────────
include_particles = st.checkbox("Include solid particles", value=False, key="cmp_include_particles")
cmp_rho_p = cmp_d50_um = cmp_phi_p = cmp_X_wt = 0.0
cmp_X_vol = 0.0
# Per-reactor particle properties & Zwietering / GMB constants
cmp_d50_per: dict[str, float] = {}
cmp_rho_p_per: dict[str, float] = {}
cmp_phi_p_per: dict[str, float] = {}
cmp_S_zw_per: dict[str, float] = {}
cmp_gmb_z_per: dict[str, float] = {}
cmp_C_D_per: dict[str, float] = {}
if include_particles and not particles_db.empty:
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        cmp_particle_name = st.selectbox("Particle (defaults)", particles_db["particle_name"].tolist(), key="cmp_particle")
        cmp_part = safe_iloc(particles_db, "particle_name", cmp_particle_name, "Particle")
        cmp_rho_p = float(cmp_part["rho_p_kg_m3"])
        cmp_d50_um = float(cmp_part["d50_um"])
        cmp_phi_p = float(cmp_part["shape_factor"])
    # Solids loading – single input with auto-conversion
    _cmp_solids_basis = st.radio("Solids loading basis",
                                  ["Mass (wt-%)", "Volume (vol-%)"],
                                  horizontal=True, key="cmp_solids_basis")
    scol1, scol2, _ = st.columns(3)
    if _cmp_solids_basis == "Mass (wt-%)":
        cmp_X_wt = scol1.number_input("Solids loading X (wt-%)", value=5.0,
                                       min_value=0.01, format="%.2f", key="cmp_Xwt",
                                       help="Mass of solids / mass of liquid × 100")
        cmp_X_vol = 100.0 * cmp_X_wt * rho / (cmp_X_wt * rho + 100.0 * cmp_rho_p) if cmp_rho_p > 0 else 0.0
        scol2.metric("Xv (vol-%)", f"{cmp_X_vol:.2f}")
    else:
        cmp_X_vol = scol1.number_input("Solids loading Xv (vol-%)", value=2.0,
                                        min_value=0.01, max_value=99.0,
                                        format="%.2f", key="cmp_Xv",
                                        help="Volume of solids / volume of slurry × 100")
        cmp_X_wt = 100.0 * cmp_rho_p * cmp_X_vol / (rho * (100.0 - cmp_X_vol)) if (100.0 - cmp_X_vol) > 0 and rho > 0 else 0.0
        scol2.metric("X (wt-%)", f"{cmp_X_wt:.2f}")
    st.markdown("**Particle & suspension parameters (per reactor)**")
    st.caption(
        "Particle properties default to the selected database entry above. "
        "Zwietering S, GMB z, and C/D depend on impeller type and geometry. "
        "Adjust any value independently for each reactor."
    )
    for _rname in selected_names:
        with st.expander(f"⚙️ {_rname}", expanded=False):
            _pc1, _pc2, _pc3 = st.columns(3)
            with _pc1:
                cmp_d50_per[_rname] = st.number_input(
                    "d50 (µm)", value=cmp_d50_um, min_value=0.01,
                    format="%.1f", key=f"cmp_d50_{_rname}",
                )
            with _pc2:
                cmp_rho_p_per[_rname] = st.number_input(
                    "ρ_p (kg/m³)", value=cmp_rho_p, min_value=1.0,
                    format="%.0f", key=f"cmp_rhop_{_rname}",
                )
            with _pc3:
                cmp_phi_p_per[_rname] = st.number_input(
                    "Shape factor φ", value=cmp_phi_p, min_value=0.01,
                    max_value=1.0, format="%.2f",
                    key=f"cmp_phi_{_rname}",
                )
            _zc1, _zc2, _zc3 = st.columns(3)
            with _zc1:
                cmp_S_zw_per[_rname] = st.number_input(
                    "Zwietering S", value=5.5, min_value=0.5,
                    max_value=20.0, format="%.1f",
                    key=f"cmp_Szw_{_rname}",
                )
            with _zc2:
                cmp_gmb_z_per[_rname] = st.number_input(
                    "GMB z constant", value=3.0, min_value=0.1,
                    max_value=30.0, format="%.2f",
                    key=f"cmp_gmb_z_{_rname}",
                    help="Geometry constant (impeller-type dependent)",
                )
            with _zc3:
                cmp_C_D_per[_rname] = st.number_input(
                    "C/D (clearance / impeller dia)", value=0.33,
                    min_value=0.01, max_value=2.0, format="%.3f",
                    key=f"cmp_CD_{_rname}",
                    help="Impeller clearance / impeller diameter",
                )
elif include_particles and particles_db.empty:
    st.warning("Particle database is empty.")

# ── Heat balance: auto-compute when ΔH is available ─────────────────────────
cmp_T_process = fluid_T_C
include_heat = rxn_delta_H != 0

# ── Scale-up matching options ────────────────────────────────────────────
st.header("4 · Scale-Up Matching")
include_scaling = st.checkbox(
    "Perform scale-up matching", value=False, key="cmp_include_scaling",
    help="Find equivalent operating conditions across reactors to match a chosen parameter.",
)

# Scaling parameter choices (params that are computable from hydro)
SCALABLE_PARAMS = [
    "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
    "Micromix time t_E (s)", "Micromix time t_E_local (s)",
    "Re", "kLa (1/s)", "kLa_surface (1/s)",
    "Avg shear rate (1/s)", "Max shear rate (1/s)",
    "Kolmogorov η (µm)", "EDCF (W/kg/s)", "Torque/V (N·m/m³)",
    "Froude number",
]

scale_basis_reactor = ""
scale_param = ""
scale_basis_rpm = 0.0
scale_basis_vol = 0.0
scale_solve_for = "RPM (specify volume)"
scale_target_known: dict[str, float] = {}  # reactor_name → known value (V or RPM)

if include_scaling:
    if len(selected_names) < 2:
        st.warning("Select at least 2 reactors for scale-up matching.")
    else:
        sc1, sc2 = st.columns(2)
        with sc1:
            scale_basis_reactor = st.selectbox(
                "Basis reactor", selected_names, key="cmp_basis_reactor",
            )
        with sc2:
            scale_param = st.selectbox(
                "Parameter to hold constant", SCALABLE_PARAMS, key="cmp_scale_param",
            )

        # Basis reactor conditions
        _basis_r = safe_iloc(reactors, "reactor_name", scale_basis_reactor, "Reactor")
        _basis_rpm_max = _safe_float(_basis_r.get("N_rpm_max"))
        _basis_rpm_min = _safe_float(_basis_r.get("N_rpm_min"))
        _basis_V_max = _safe_float(_basis_r.get("V_L_max"))
        _basis_V_min = _safe_float(_basis_r.get("V_L_min"))
        if _basis_rpm_max == 0:
            _basis_rpm_max = _safe_float(_basis_r.get("N_rps")) * 60
        if _basis_rpm_min == 0:
            _basis_rpm_min = _basis_rpm_max
        if _basis_V_max == 0:
            _basis_V_max = _safe_float(_basis_r.get("V_L"))
        if _basis_V_min == 0:
            _basis_V_min = _basis_V_max

        _basis_rpm_mid = (_basis_rpm_min + _basis_rpm_max) / 2.0
        _basis_V_mid = (_basis_V_min + _basis_V_max) / 2.0

        st.markdown(f"**Basis reactor conditions** ({scale_basis_reactor})")
        bc1, bc2 = st.columns(2)
        with bc1:
            scale_basis_rpm = st.number_input(
                "Basis RPM", value=_basis_rpm_mid,
                min_value=0.1, format="%.1f", key="cmp_scale_basis_rpm",
            )
        with bc2:
            scale_basis_vol = st.number_input(
                "Basis volume (L)", value=_basis_V_mid,
                min_value=0.001, format="%.2f", key="cmp_scale_basis_vol",
            )

        scale_solve_for = st.radio(
            "For target reactors, solve for",
            ["RPM (specify volume)", "Volume (specify RPM)"],
            horizontal=True, key="cmp_scale_solve_for",
        )

        st.markdown("**Target reactor conditions**")
        _target_names = [n for n in selected_names if n != scale_basis_reactor]
        for _tname in _target_names:
            _t_r = safe_iloc(reactors, "reactor_name", _tname, "Reactor")
            if scale_solve_for.startswith("RPM"):
                _t_V_max = _safe_float(_t_r.get("V_L_max"))
                _t_V_min = _safe_float(_t_r.get("V_L_min"))
                if _t_V_max == 0:
                    _t_V_max = _safe_float(_t_r.get("V_L"))
                if _t_V_max == 0:
                    _D_t = _safe_float(_t_r.get("D_tank_m"))
                    _H_t = _safe_float(_t_r.get("H_m"))
                    _t_V_max = np.pi / 4 * _D_t**2 * _H_t * 1000
                if _t_V_min == 0:
                    _t_V_min = _t_V_max
                _t_V_mid = (_t_V_min + _t_V_max) / 2.0
                scale_target_known[_tname] = st.number_input(
                    f"{_tname} — fill volume (L)",
                    value=_t_V_mid, min_value=0.001, format="%.2f",
                    key=f"cmp_scale_tv_{_tname}",
                )
            else:
                _t_rpm_max = _safe_float(_t_r.get("N_rpm_max"))
                _t_rpm_min = _safe_float(_t_r.get("N_rpm_min"))
                if _t_rpm_max == 0:
                    _t_rpm_max = _safe_float(_t_r.get("N_rps")) * 60
                if _t_rpm_min == 0:
                    _t_rpm_min = _t_rpm_max
                _t_rpm_mid = (_t_rpm_min + _t_rpm_max) / 2.0
                scale_target_known[_tname] = st.number_input(
                    f"{_tname} — stir speed (RPM)",
                    value=_t_rpm_mid, min_value=0.1, format="%.1f",
                    key=f"cmp_scale_trpm_{_tname}",
                )

# ── Gate 3: confirm options & compute ───────────────────────────────────
if st.session_state["_p7_step"] < 3:
    st.info("👆 Review the additional options above, then confirm to compute results.")
    if st.button("🔬 Confirm & Compute", key="p7_gate3", type="primary"):
        st.session_state["_p7_step"] = 3
        st.rerun()
    st.stop()


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

_spinner_placeholder = st.empty()
_spinner_placeholder.info("⏳ Computing reactor envelopes…")

for rname in selected_names:
    r = safe_iloc(reactors, "reactor_name", rname, "Reactor")
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

    # C/D ratio for GMB Njs (from reactor DB or user per-reactor default)
    _imp1_C = _safe_float(r.get("imp1_clearance_m"), 0.0)
    reactor_info[rname]["C_D_ratio"] = _imp1_C / D_imp_v if D_imp_v > 0 and _imp1_C > 0 else cmp_C_D_per.get(rname, 0.33)

    # Heat-transfer geometry for this reactor
    _r_material = str(r.get("shell_material", ""))
    _r_lining = str(r.get("lining_material", ""))
    _r_bottom_dish = str(r.get("bottom_dish", "")) if pd.notna(r.get("bottom_dish")) else ""
    _r_U_override = _safe_float(r.get("U_W_m2K"), 0.0)
    _r_A_override = _safe_float(r.get("A_ht_m2"), 0.0)
    _r_wall_mm = _safe_float(r.get("wall_thickness_mm"), 0.0)
    reactor_info[rname]["shell_material"] = _r_material
    reactor_info[rname]["lining_material"] = _r_lining
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
        # Particle parameters (if enabled)
        part_vals = {}
        _kLa_SL_corner = 0.0
        if include_particles and cmp_d50_um > 0:
            _dp = cmp_d50_per.get(rname, cmp_d50_um) * 1e-6
            _rho_p_r = cmp_rho_p_per.get(rname, cmp_rho_p)
            _phi_p_r = cmp_phi_p_per.get(rname, cmp_phi_p)
            _nu = mu / rho if rho > 0 else 0.0
            _drho = abs(_rho_p_r - rho)
            _vt = settling_velocity(_dp, _rho_p_r, rho, mu, _phi_p_r)
            _rep = particle_reynolds(_dp, _vt, rho, mu)
            _njs_zw = zwietering_njs(cmp_S_zw_per.get(rname, 5.5), _nu, _dp, _drho, rho, cmp_X_wt, D_imp_v)
            _Np_corner = Np_v if Np_v else 1.27
            _C_D_corner = _safe_float(r.get("imp1_clearance_m")) / D_imp_v if D_imp_v > 0 and _safe_float(r.get("imp1_clearance_m")) > 0 else cmp_C_D_per.get(rname, 0.33)
            _njs_gmb = gmb_njs(cmp_gmb_z_per.get(rname, 3.0), _Np_corner, D_imp_v, _dp,
                               _drho, rho, cmp_X_vol, _C_D_corner)
            _njs = max(_njs_zw, _njs_gmb)
            _eps_kg_corner = h["P/V (W/kg)"] if h["P/V (W/kg)"] > 0 else 0.0
            _v_slip_corner = max(_vt,
                                 (_eps_kg_corner * _dp) ** (1.0 / 3.0) if _eps_kg_corner > 0 else 0.0)
            _ksl = solid_liquid_mass_transfer(_dp, _v_slip_corner, rho, mu, D_mol)
            _phi_s_corner = cmp_X_vol / 100.0
            _kLa_SL_corner = solid_liquid_kla(_ksl, _dp, _phi_s_corner)
            _njs_rpm = _njs * 60
            _n_over_njs = N / _njs if _njs > 0 else 0.0
            part_vals = {
                "N_js Zwietering (RPM)": _njs_zw * 60,
                "N_js GMB (RPM)": _njs_gmb * 60,
                "N_js (RPM)": _njs_rpm,
                "N/N_js": _n_over_njs,
                "v_t (m/s)": _vt,
                "Re_p": _rep,
                "k_SL (m/s)": _ksl,
                "kLa_SL (1/s)": _kLa_SL_corner,
            }
        da = compute_damkohler_numbers(
            h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
            kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
            kLa_SL=_kLa_SL_corner,
        )
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

_spinner_placeholder.empty()

if not envelope_rows:
    st.info("No computable reactors in the selection.")
    st.stop()

env_df = pd.DataFrame(envelope_rows)
env_df["RPM_pct"] = env_df["RPM"] / env_df["RPM_max"] * 100

# ── Aggregate min / max per reactor for plotting ─────────────────────────
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
    _nu_p7 = mu / rho if rho > 0 else 0.0
    for rname, info in reactor_info.items():
        _dp_p7 = cmp_d50_per.get(rname, cmp_d50_um) * 1e-6
        _rho_p_p7 = cmp_rho_p_per.get(rname, cmp_rho_p)
        _phi_p_p7 = cmp_phi_p_per.get(rname, cmp_phi_p)
        _drho_p7 = abs(_rho_p_p7 - rho)
        _vt_p7 = settling_velocity(_dp_p7, _rho_p_p7, rho, mu, _phi_p_p7)
        _rep_p7 = particle_reynolds(_dp_p7, _vt_p7, rho, mu)
        _njs_zw_p7 = zwietering_njs(cmp_S_zw_per.get(rname, 5.5), _nu_p7, _dp_p7, _drho_p7, rho, cmp_X_wt, info["D_imp"])
        _Np_p7 = info["Np"] if info["Np"] else 1.27
        _C_D_p7 = info.get("C_D_ratio", cmp_C_D_per.get(rname, 0.33))
        _njs_gmb_p7 = gmb_njs(cmp_gmb_z_per.get(rname, 3.0), _Np_p7, info["D_imp"], _dp_p7,
                               _drho_p7, rho, cmp_X_vol, _C_D_p7)
        _njs_p7 = max(_njs_zw_p7, _njs_gmb_p7)
        _phi_s_p7 = cmp_X_vol / 100.0
        _p7_part_static[rname] = {
            "N_js_zw_rps": _njs_zw_p7, "N_js_gmb_rps": _njs_gmb_p7,
            "N_js_rps": _njs_p7,
            "N_js Zwietering (RPM)": _njs_zw_p7 * 60,
            "N_js GMB (RPM)": _njs_gmb_p7 * 60,
            "N_js (RPM)": _njs_p7 * 60,
            "v_t (m/s)": _vt_p7, "Re_p": _rep_p7,
            "d_p": _dp_p7, "phi_s": _phi_s_p7,
        }

# ── Cache the envelope sweep on a signature of all its inputs ────────────
# The reactor × volume × RPM sweep below is the most expensive computation on
# the page. Re-running it on every Streamlit rerun (e.g. an unrelated widget
# change) is wasteful, so memoize it in session_state keyed on a hash of every
# input that affects the result.
_L = locals()


def _sig_round(x, ndigits: int = 9):
    try:
        return round(float(x), ndigits)
    except (TypeError, ValueError):
        return x


_p7_env_sig = (
    repr(sorted(reactor_info.items())),
    repr(sorted(corr_modes.items())) if isinstance(corr_modes, dict) else repr(corr_modes),
    repr(sorted(_p7_part_static.items())),
    _sig_round(rho), _sig_round(mu), _sig_round(v_s), bool(is_coalescing),
    _sig_round(D_mol), _sig_round(t_rxn),
    bool(include_heat), _sig_round(rxn_delta_H),
    str(_L.get("rxn_order", "")), _sig_round(_L.get("rxn_k", 0.0)),
    _sig_round(_L.get("rxn_C0", 0.0)),
    _sig_round(_L.get("cmp_T_process", 0.0)), _sig_round(_L.get("cmp_T_coolant", 0.0)),
    str(fluid_name), tuple(PLOT_PARAMS), _N_INTERP,
)
_p7_cached = (st.session_state.get("_p7_env_sig") == _p7_env_sig
              and "_p7_curve_data" in st.session_state)

if _p7_cached:
    curve_data = st.session_state["_p7_curve_data"]
else:
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
                _kLa_SL_env_p7 = 0.0
                _ksl_env_p7 = 0.0
                if rname in _p7_part_static:
                    _ps = _p7_part_static[rname]
                    _eps_kg_p7 = h["P/V (W/kg)"] if h["P/V (W/kg)"] > 0 else 0.0
                    _v_slip_p7 = max(_ps["v_t (m/s)"],
                                     (_eps_kg_p7 * _ps["d_p"]) ** (1.0 / 3.0) if _eps_kg_p7 > 0 else 0.0)
                    _ksl_env_p7 = solid_liquid_mass_transfer(_ps["d_p"], _v_slip_p7, rho, mu, D_mol)
                    _kLa_SL_env_p7 = solid_liquid_kla(_ksl_env_p7, _ps["d_p"], _ps["phi_s"])
                da = compute_damkohler_numbers(
                    h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
                    kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
                    kLa_SL=_kLa_SL_env_p7,
                )
                vals = {**h, **da}
                # Particle parameters — use pre-computed statics + dynamic k_SL/kLa_SL/Da_SL
                if rname in _p7_part_static:
                    _ps = _p7_part_static[rname]
                    vals["N_js Zwietering (RPM)"] = _ps["N_js Zwietering (RPM)"]
                    vals["N_js GMB (RPM)"] = _ps["N_js GMB (RPM)"]
                    vals["N_js (RPM)"] = _ps["N_js (RPM)"]
                    vals["N/N_js"] = N / _ps["N_js_rps"] if _ps["N_js_rps"] > 0 else 0.0
                    vals["v_t (m/s)"] = _ps["v_t (m/s)"]
                    vals["Re_p"] = _ps["Re_p"]
                    vals["k_SL (m/s)"] = _ksl_env_p7
                    vals["kLa_SL (1/s)"] = _kLa_SL_env_p7
                    vals["Da_SL"] = da["Da_SL"]
                # Heat balance — only U depends on RPM
                if include_heat and rxn_delta_H != 0:
                    if info["U_override"] > 0:
                        _U_ht = info["U_override"]
                    else:
                        _U_ht, _ = estimate_U_detailed(
                            N_rps=N, D_imp=info["D_imp"], D_tank=info["D_tank"],
                            rho=rho, mu=mu,
                            material=info["shell_material"],
                            lining_material=info["lining_material"],
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

    st.session_state["_p7_env_sig"] = _p7_env_sig
    st.session_state["_p7_curve_data"] = curve_data

st.divider()

# ── Summary Table ─────────────────────────────────────────────────────────
st.header("5 · Operating Envelope Summary")
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

st.dataframe(pd.DataFrame(table_rows), width='content', hide_index=True)

# ── Detail: 4-corner table ───────────────────────────────────────────────
with st.expander("Full 4-corner detail table", expanded=False):
    detail_cols = ["Reactor", "Corner", "N (rev/s)", "V_L", "Re",
                   "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
                   "Micromix time t_E (s)", "Micromix time t_E_local (s)",
                   "Kolmogorov η (µm)",
                   "Avg shear rate (1/s)", "Max shear rate (1/s)",
                   "Avg shear stress (Pa)", "kLa (1/s)", "kLa_surface (1/s)",
                   "Da_macro", "Da_micro", "Da_GL", "Da_SL"]
    if include_particles and cmp_d50_um > 0:
        detail_cols += PARTICLE_PARAMS
    if include_heat and rxn_delta_H != 0:
        detail_cols += HEAT_PARAMS
    # Deduplicate while preserving order, then keep only columns in the dataframe
    detail_cols = list(dict.fromkeys(c for c in detail_cols if c in env_df.columns))
    fmt = {c: "{:.3g}" for c in detail_cols if c not in ("Reactor", "Corner")}
    st.dataframe(env_df[detail_cols].style.format(fmt), width='content', hide_index=True)

st.divider()

# ── Charts: operating envelopes ──────────────────────────────────────────
st.header("6 · Operating Envelope Charts")

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

    st.dataframe(pd.DataFrame(rpm_table_rows), width='content', hide_index=True)

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
        "Da_SL": "Solid-Liquid Mass Transfer (Da_SL)",
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
        st.plotly_chart(fig, width='content')

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
    st.header("6b · Heat Balance Summary")
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
        st.dataframe(pd.DataFrame(heat_summary_rows), width='content', hide_index=True)

    # Show U estimation warnings (once per reactor, at max RPM)
    _u_warnings_all: list[str] = []
    for rname, info in reactor_info.items():
        if info["U_override"] > 0:
            continue
        _, _u_warns = estimate_U_detailed(
            N_rps=info["N_hi"], D_imp=info["D_imp"], D_tank=info["D_tank"],
            rho=rho, mu=mu,
            material=info["shell_material"],
            lining_material=info["lining_material"],
            wall_thickness_mm=info["wall_thickness_mm"],
            fluid_name=fluid_name,
        )
        if _u_warns:
            _u_warnings_all.append(f"**{rname}**: " + "; ".join(_u_warns))
    if _u_warnings_all:
        with st.expander("ℹ️ U estimation notes", expanded=False):
            for _w in _u_warnings_all:
                st.markdown(f"- {_w}")

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
        st.plotly_chart(fig_heat, width='content')

# ── Scale-up matching computation ────────────────────────────────────────
scaling_results: list[dict] = []
scaling_all_params: list[dict] = []  # full hydro at matched condition

if include_scaling and scale_basis_reactor and scale_param and len(selected_names) >= 2:
    from scipy.optimize import brentq

    st.divider()
    st.header("7 · Scale-Up Matching Results")

    # 1. Compute target value from basis reactor
    _basis_info_r = safe_iloc(reactors, "reactor_name", scale_basis_reactor, "Reactor")
    _b_D_imp = _safe_float(_basis_info_r.get("D_imp_m"))
    _b_D_tank = _safe_float(_basis_info_r.get("D_tank_m"))
    _b_H_max = _safe_float(_basis_info_r.get("H_m"))
    _b_Np = _safe_float(_basis_info_r.get("Np"), None)
    _b_Nq = _safe_float(_basis_info_r.get("Nq"), None)
    _b_bottom_dish = str(_basis_info_r.get("bottom_dish", "")) if pd.notna(_basis_info_r.get("bottom_dish")) else ""
    _b_N = scale_basis_rpm / 60.0
    _b_H = _liquid_height(scale_basis_vol, _b_D_tank, _b_H_max, _b_bottom_dish)

    _basis_hydro, _ = compute_reactor_hydro_with_mode(
        corr_modes.get(scale_basis_reactor, "Literature"), scale_basis_reactor,
        N=_b_N, D_imp=_b_D_imp, D_tank=_b_D_tank, H=_b_H,
        rho=rho, mu=mu, Np=_b_Np, Nq=_b_Nq,
        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
    )
    target_value = _basis_hydro[scale_param]
    st.caption(
        f"**Basis:** {scale_basis_reactor} at {scale_basis_rpm:.0f} RPM, "
        f"{scale_basis_vol:.1f} L  →  **{scale_param} = {target_value:.4g}**"
    )

    # Build basis row
    _basis_row = {
        "Reactor": scale_basis_reactor,
        "Role": "Basis",
        "RPM": scale_basis_rpm,
        "Volume (L)": scale_basis_vol,
        scale_param: target_value,
        "Status": "—",
    }
    scaling_results.append(_basis_row)
    scaling_all_params.append({"Reactor": scale_basis_reactor, "Role": "Basis",
                                "RPM": scale_basis_rpm, "Volume (L)": scale_basis_vol,
                                **_basis_hydro})

    # 2. Solve for each target reactor
    _target_names = [n for n in selected_names if n != scale_basis_reactor]
    for _tname in _target_names:
        _t_r = safe_iloc(reactors, "reactor_name", _tname, "Reactor")
        _t_D_imp = _safe_float(_t_r.get("D_imp_m"))
        _t_D_tank = _safe_float(_t_r.get("D_tank_m"))
        _t_H_max = _safe_float(_t_r.get("H_m"))
        _t_Np = _safe_float(_t_r.get("Np"), None)
        _t_Nq = _safe_float(_t_r.get("Nq"), None)
        _t_bottom_dish = str(_t_r.get("bottom_dish", "")) if pd.notna(_t_r.get("bottom_dish")) else ""

        # RPM operating limits
        _t_rpm_max = _safe_float(_t_r.get("N_rpm_max"))
        _t_rpm_min = _safe_float(_t_r.get("N_rpm_min"))
        if _t_rpm_max == 0:
            _t_rpm_max = _safe_float(_t_r.get("N_rps")) * 60
        if _t_rpm_min == 0:
            _t_rpm_min = 1.0  # absolute floor for solver

        # Volume limits
        _t_V_max = _safe_float(_t_r.get("V_L_max"))
        _t_V_min = _safe_float(_t_r.get("V_L_min"))
        if _t_V_max == 0:
            _t_V_max = _safe_float(_t_r.get("V_L"))
        if _t_V_max == 0:
            _t_V_max = np.pi / 4 * _t_D_tank**2 * _t_H_max * 1000
        if _t_V_min == 0:
            _t_V_min = _t_V_max * 0.1  # 10% fill as floor for solver

        _known = scale_target_known.get(_tname, 0.0)

        if _t_D_imp == 0 or _t_D_tank == 0 or _t_H_max == 0:
            scaling_results.append({
                "Reactor": _tname, "Role": "Target",
                "RPM": np.nan, "Volume (L)": np.nan,
                scale_param: np.nan,
                "Status": "Missing geometry",
            })
            continue

        if scale_solve_for.startswith("RPM"):
            # Known: volume → solve for RPM
            _t_V_L = _known
            _t_H_v = _liquid_height(_t_V_L, _t_D_tank, _t_H_max, _t_bottom_dish)

            def _obj_rpm(rpm, _V=_t_V_L, _H=_t_H_v, _Di=_t_D_imp,
                         _Dt=_t_D_tank, _Np=_t_Np, _Nq=_t_Nq,
                         _name=_tname):
                _N = rpm / 60.0
                h, _ = compute_reactor_hydro_with_mode(
                    corr_modes.get(_name, "Literature"), _name,
                    N=_N, D_imp=_Di, D_tank=_Dt, H=_H,
                    rho=rho, mu=mu, Np=_Np, Nq=_Nq,
                    v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                )
                return h[scale_param] - target_value

            # Try root-finding
            try:
                _lo = max(_t_rpm_min, 0.5)
                _hi = _t_rpm_max * 1.5  # allow slight overshoot to find root
                _f_lo = _obj_rpm(_lo)
                _f_hi = _obj_rpm(_hi)
                if _f_lo * _f_hi > 0:
                    # No sign change → target may not be achievable
                    # Report closest boundary
                    _v_lo_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_lo / 60, D_imp=_t_D_imp, D_tank=_t_D_tank, H=_t_H_v,
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    _v_hi_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_hi / 60, D_imp=_t_D_imp, D_tank=_t_D_tank, H=_t_H_v,
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    # Pick the boundary closer to the target
                    _err_lo = abs(_v_lo_h[scale_param] - target_value)
                    _err_hi = abs(_v_hi_h[scale_param] - target_value)
                    if _err_lo < _err_hi:
                        _best_rpm, _best_val = _lo, _v_lo_h[scale_param]
                        _best_h = _v_lo_h
                    else:
                        _best_rpm, _best_val = _hi, _v_hi_h[scale_param]
                        _best_h = _v_hi_h
                    _in_range = _t_rpm_min <= _best_rpm <= _t_rpm_max
                    scaling_results.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _best_rpm, "Volume (L)": _t_V_L,
                        scale_param: _best_val,
                        "Status": f"Not achievable (closest: {_best_val:.4g})"
                                  + ("" if _in_range else " [outside RPM range]"),
                    })
                    scaling_all_params.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _best_rpm, "Volume (L)": _t_V_L,
                        **_best_h,
                    })
                else:
                    _solved_rpm = brentq(_obj_rpm, _lo, _hi, xtol=0.01, maxiter=200)
                    _in_range = _t_rpm_min <= _solved_rpm <= _t_rpm_max
                    _solved_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_solved_rpm / 60, D_imp=_t_D_imp, D_tank=_t_D_tank, H=_t_H_v,
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    _status = "Matched"
                    if not _in_range:
                        _status = f"Matched (outside operating range {_t_rpm_min:.0f}–{_t_rpm_max:.0f} RPM)"
                    scaling_results.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _solved_rpm, "Volume (L)": _t_V_L,
                        scale_param: _solved_h[scale_param],
                        "Status": _status,
                    })
                    scaling_all_params.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _solved_rpm, "Volume (L)": _t_V_L,
                        **_solved_h,
                    })
            except Exception as _exc:
                scaling_results.append({
                    "Reactor": _tname, "Role": "Target",
                    "RPM": np.nan, "Volume (L)": _t_V_L,
                    scale_param: np.nan,
                    "Status": f"Solver error: {_exc}",
                })

        else:
            # Known: RPM → solve for volume
            _t_rpm = _known
            _t_N = _t_rpm / 60.0

            def _obj_vol(vol_L, _N=_t_N, _Di=_t_D_imp,
                         _Dt=_t_D_tank, _Hmax=_t_H_max, _Np=_t_Np, _Nq=_t_Nq,
                         _name=_tname, _dish=_t_bottom_dish):
                _H = _liquid_height(vol_L, _Dt, _Hmax, _dish)
                h, _ = compute_reactor_hydro_with_mode(
                    corr_modes.get(_name, "Literature"), _name,
                    N=_N, D_imp=_Di, D_tank=_Dt, H=_H,
                    rho=rho, mu=mu, Np=_Np, Nq=_Nq,
                    v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                )
                return h[scale_param] - target_value

            try:
                _lo_v = max(_t_V_min * 0.5, 0.001)
                _hi_v = _t_V_max * 1.2
                _f_lo = _obj_vol(_lo_v)
                _f_hi = _obj_vol(_hi_v)
                if _f_lo * _f_hi > 0:
                    _v_lo_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_t_N, D_imp=_t_D_imp, D_tank=_t_D_tank,
                        H=_liquid_height(_lo_v, _t_D_tank, _t_H_max, _t_bottom_dish),
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    _v_hi_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_t_N, D_imp=_t_D_imp, D_tank=_t_D_tank,
                        H=_liquid_height(_hi_v, _t_D_tank, _t_H_max, _t_bottom_dish),
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    _err_lo = abs(_v_lo_h[scale_param] - target_value)
                    _err_hi = abs(_v_hi_h[scale_param] - target_value)
                    if _err_lo < _err_hi:
                        _best_vol, _best_val = _lo_v, _v_lo_h[scale_param]
                        _best_h = _v_lo_h
                    else:
                        _best_vol, _best_val = _hi_v, _v_hi_h[scale_param]
                        _best_h = _v_hi_h
                    _in_range = _t_V_min <= _best_vol <= _t_V_max
                    scaling_results.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _t_rpm, "Volume (L)": _best_vol,
                        scale_param: _best_val,
                        "Status": f"Not achievable (closest: {_best_val:.4g})"
                                  + ("" if _in_range else " [outside volume range]"),
                    })
                    scaling_all_params.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _t_rpm, "Volume (L)": _best_vol,
                        **_best_h,
                    })
                else:
                    _solved_vol = brentq(_obj_vol, _lo_v, _hi_v, xtol=0.001, maxiter=200)
                    _in_range = _t_V_min <= _solved_vol <= _t_V_max
                    _solved_h, _ = compute_reactor_hydro_with_mode(
                        corr_modes.get(_tname, "Literature"), _tname,
                        N=_t_N, D_imp=_t_D_imp, D_tank=_t_D_tank,
                        H=_liquid_height(_solved_vol, _t_D_tank, _t_H_max, _t_bottom_dish),
                        rho=rho, mu=mu, Np=_t_Np, Nq=_t_Nq,
                        v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
                    )
                    _status = "Matched"
                    if not _in_range:
                        _status = f"Matched (outside operating range {_t_V_min:.1f}–{_t_V_max:.1f} L)"
                    scaling_results.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _t_rpm, "Volume (L)": _solved_vol,
                        scale_param: _solved_h[scale_param],
                        "Status": _status,
                    })
                    scaling_all_params.append({
                        "Reactor": _tname, "Role": "Target",
                        "RPM": _t_rpm, "Volume (L)": _solved_vol,
                        **_solved_h,
                    })
            except Exception as _exc:
                scaling_results.append({
                    "Reactor": _tname, "Role": "Target",
                    "RPM": _t_rpm, "Volume (L)": np.nan,
                    scale_param: np.nan,
                    "Status": f"Solver error: {_exc}",
                })

    # Display scaling results table
    if scaling_results:
        _scale_df = pd.DataFrame(scaling_results)
        st.subheader("Matched Operating Conditions")
        _fmt_cols = {c: "{:.4g}" for c in _scale_df.columns
                     if c not in ("Reactor", "Role", "Status")}
        st.dataframe(_scale_df.style.format(_fmt_cols), width='content', hide_index=True)

    # Full parameter comparison at matched conditions
    if scaling_all_params:
        with st.expander("Full parameter comparison at matched conditions", expanded=False):
            _sp_df = pd.DataFrame(scaling_all_params)
            _show_cols = ["Reactor", "Role", "RPM", "Volume (L)", "Re",
                          "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
                          "Micromix time t_E (s)", "Micromix time t_E_local (s)",
                          "Kolmogorov η (µm)", "Avg shear rate (1/s)",
                          "Max shear rate (1/s)", "Avg shear stress (Pa)",
                          "kLa (1/s)", "kLa_surface (1/s)",
                          "Torque (N·m)", "Torque/V (N·m/m³)",
                          "EDCF (W/kg/s)", "Froude number"]
            _show_cols = [c for c in _show_cols if c in _sp_df.columns]
            _fmt2 = {c: "{:.4g}" for c in _show_cols
                     if c not in ("Reactor", "Role")}
            st.dataframe(_sp_df[_show_cols].style.format(_fmt2),
                         width='content', hide_index=True)

        # Percentage difference relative to basis reactor
        with st.expander("Percentage difference vs. basis reactor", expanded=False):
            _sp_df2 = pd.DataFrame(scaling_all_params)
            _num_cols = [c for c in _show_cols if c not in ("Reactor", "Role")]
            _basis_row = _sp_df2[_sp_df2["Role"] == "Basis"].iloc[0]
            _pct_rows = []
            for _, _row in _sp_df2.iterrows():
                _pct_entry: dict = {"Reactor": _row["Reactor"], "Role": _row["Role"]}
                for c in _num_cols:
                    _b_val = _basis_row.get(c, 0)
                    _t_val = _row.get(c, 0)
                    if _b_val and np.isfinite(_b_val) and _b_val != 0 and np.isfinite(_t_val):
                        _pct_entry[c] = (_t_val - _b_val) / abs(_b_val) * 100
                    else:
                        _pct_entry[c] = np.nan
                _pct_rows.append(_pct_entry)
            _pct_df = pd.DataFrame(_pct_rows)
            _pct_cols = [c for c in _show_cols if c in _pct_df.columns]
            _fmt_pct = {c: "{:+.1f}%" for c in _pct_cols
                        if c not in ("Reactor", "Role")}
            st.dataframe(_pct_df[_pct_cols].style.format(_fmt_pct),
                         width='content', hide_index=True)

st.divider()

# ── Scale-up summary ─────────────────────────────────────────────────────
st.header("8 · Scale-Up Impact Summary")
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
    st.dataframe(summary_df.style.format({c: "{:.2f}" for c in numeric_cols}), width='content')
else:
    st.info("Select at least two reactors to see scale-up ratios.")

st.divider()

# ── Save results per reactor ─────────────────────────────────────────────
st.header("9 · Save Results")
st.caption(
    "Save the max-RPM / max-V corner result for each selected reactor to "
    "Recorded Results (same format as the Mixing Sensitivity page)."
)

if st.button("📌 Save results for all selected reactors", key="cmp_save_all"):
    if "recorded_results" not in st.session_state:
        st.session_state.recorded_results = pd.DataFrame()

    _saved_count = 0
    for rname in selected_names:
        _corners = env_df[
            (env_df["Reactor"] == rname) & (env_df["Corner"] == "max RPM / max V")
        ]
        if _corners.empty:
            continue
        _c = _corners.iloc[0]
        _result_row = {
            "reactor": rname,
            "reaction": rxn_name if not reactions.empty else "",
            "fluid": fluid_name,
            "fluid_T_C": fluid_T_C,
            "RPM": _c["RPM"],
            "Volume (L)": _c["V_L"],
            "Re": _c.get("Re", ""),
            "P/V (W/L)": _c.get("P/V (W/L)", ""),
            "Tip speed (m/s)": _c.get("Tip speed (m/s)", ""),
            "Blend time (s)": _c.get("Blend time 95% (s)", ""),
            "Circulation time (s)": _c.get("Circulation time (s)", ""),
            "Micromix t_E (s)": _c.get("Micromix time t_E (s)", ""),
            "Micromix t_E_local (s)": _c.get("Micromix time t_E_local (s)", ""),
            "Kolmogorov η (µm)": _c.get("Kolmogorov η (µm)", ""),
            "EDCF (W/kg/s)": _c.get("EDCF (W/kg/s)", ""),
            "Torque (N·m)": _c.get("Torque (N·m)", ""),
            "Froude number": _c.get("Froude number", ""),
            "Avg shear rate (1/s)": _c.get("Avg shear rate (1/s)", ""),
            "Max shear rate (1/s)": _c.get("Max shear rate (1/s)", ""),
            "Avg shear stress (Pa)": _c.get("Avg shear stress (Pa)", ""),
            "kLa (1/s)": _c.get("kLa (1/s)", ""),
            "kLa_surface (1/s)": _c.get("kLa_surface (1/s)", ""),
            "t_rxn (s)": t_rxn,
            "Da_macro": _c.get("Da_macro", ""),
            "Da_micro": _c.get("Da_micro", ""),
            "Da_GL": _c.get("Da_GL", ""),
            "Da_SL": _c.get("Da_SL", ""),
            "Assessment": _c.get("Assessment", ""),
        }
        if include_particles and cmp_d50_um > 0:
            _result_row.update({
                "Particle": cmp_particle_name if "cmp_particle_name" in dir() else "",
                "d50 (µm)": cmp_d50_per.get(rname, cmp_d50_um),
                "ρ_p (kg/m³)": cmp_rho_p_per.get(rname, cmp_rho_p),
                "v_t (m/s)": _c.get("v_t (m/s)", ""),
                "Re_p": _c.get("Re_p", ""),
                "N_js (RPM)": _c.get("N_js (RPM)", ""),
                "k_SL (m/s)": _c.get("k_SL (m/s)", ""),
            })
        st.session_state.recorded_results = pd.concat(
            [st.session_state.recorded_results, pd.DataFrame([_result_row])],
            ignore_index=True,
        )
        _saved_count += 1

    _results_csv = DATA_DIR / "recorded_results.csv"
    st.session_state.recorded_results.to_csv(_results_csv, index=False)
    st.success(f"Saved {_saved_count} reactor result(s) to **Recorded Results**.")

# ── Generate PDF Report ──────────────────────────────────────────────────
st.divider()
st.header("10 · Export Report")

if st.button("📥 Export PDF Report", type="primary", key="p7_export_pdf"):
    with st.spinner("Generating PDF…"):
        try:
            from utils.report_builder import build_reactor_comparison_pdf, report_filename
            import pandas as _pd_report
            # Select key params for report charts
            _report_chart_params = [
                "Da_micro", "Da_macro", "Da_GL", "P/V (W/L)",
                "Blend time 95% (s)", "Tip speed (m/s)",
            ]
            if include_heat and rxn_delta_H != 0:
                _report_chart_params.append("Q_gen/Q_cool (%)")
            _report_chart_params = [p for p in _report_chart_params if p in PLOT_PARAMS]

            _p7_snap = {
                "selected_names": selected_names,
                "fluid": fluid_name,
                "fluid_T_C": fluid_T_C,
                "reaction": rxn_name if not reactions.empty else "N/A",
                "t_rxn": t_rxn,
                "env_df": env_df,
                "agg_df": agg_df,
                "reactor_info": reactor_info,
                "include_heat": include_heat and rxn_delta_H != 0,
                "include_particles": include_particles and cmp_d50_um > 0,
                "scaling_results": scaling_results if include_scaling else [],
                "scaling_all_params": scaling_all_params if include_scaling else [],
                "scale_param": scale_param if include_scaling else "",
                "scale_basis_reactor": scale_basis_reactor if include_scaling else "",
                "curve_data": curve_data,
                "report_chart_params": _report_chart_params,
            }
            _pdf_bytes = build_reactor_comparison_pdf(_p7_snap)
            st.session_state["_p7_pdf_bytes"] = _pdf_bytes
            st.session_state["_p7_pdf_name"] = report_filename(
                "Reactor_Comparison", selected_names[0] if selected_names else ""
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_p7_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_p7_pdf_bytes"],
        file_name=st.session_state["_p7_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")
