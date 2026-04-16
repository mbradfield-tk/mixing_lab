"""
Page – Heat Transfer Modeling Tool
===================================
Estimate heating/cooling times, compute U from first principles,
visualize batch temperature profiles, and compare heat-transfer media.
"""

import streamlit as st
import pandas as pd
import numpy as np
import pathlib
import plotly.graph_objects as go

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, list_solvents, is_known_solvent,
)
from utils.calculations import (
    impeller_power,
    reynolds_number,
    liquid_height_from_volume,
    estimate_jacket_area,
    estimate_U_from_resistances,
    heat_removal_capacity,
    time_to_cool_or_heat,
    cooling_rate,
    nusselt_jacket,
    process_side_htc,
    jacket_side_htc,
    batch_temperature_profile,
    batch_temp_profile_variable_jacket,
    NUSSELT_CORRELATIONS,
    HTM_DB,
    WALL_CONDUCTIVITY,
    LINING_CONDUCTIVITY,
    LINING_THICKNESS_DEFAULT,
    JACKET_HTC_DEFAULT,
    FOULING_DEFAULT,
    SOLVENT_THERMAL,
    _lookup_wall_k,
    _lookup_lining_k,
    _lookup_solvent_thermal,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── Load databases ───────────────────────────────────────────────────────

def _load_csv(key, filename, columns):
    if key not in st.session_state:
        p = DATA_DIR / filename
        if p.exists():
            st.session_state[key] = pd.read_csv(p)
        else:
            st.session_state[key] = pd.DataFrame(columns=columns)
    return st.session_state[key]


reactors = _load_csv("reactor_db", "reactors.csv", ["reactor_name"])
custom_fluids = _load_csv("fluid_db", "fluids.csv", ["fluid_name"])

_solvent_names = sorted(SOLVENT_DB.keys())
_custom_names = custom_fluids["fluid_name"].tolist() if not custom_fluids.empty else []
_all_fluid_names = _solvent_names + _custom_names


def _safe(series_val, default):
    try:
        v = float(series_val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


st.title("🔥 Heat Transfer Modeling")

if reactors.empty:
    st.warning("Please populate the Reactor database first.")
    st.stop()

# ── Callbacks to clear dependent widget keys on upstream changes ─────
_HT_REACTOR_KEYS = [
    "ht_Dt", "ht_Di", "ht_N", "ht_Np",
    "ht_wall_mat", "ht_wall_mm", "ht_lining", "ht_lining_mm",
    "ht_A_calc", "ht_A_ov",
]
_HT_FLUID_KEYS = ["ht_rho", "ht_mu", "ht_Cp", "ht_kf"]
_HT_AREA_KEYS = ["ht_A_calc", "ht_A_ov"]

def _on_reactor_change():
    for k in _HT_REACTOR_KEYS:
        st.session_state.pop(k, None)
    st.session_state.pop("_ht_computed", None)
    # Mark that reactor changed so we can force-update volume
    st.session_state["_ht_reactor_changed"] = True

def _on_fluid_change():
    for k in _HT_FLUID_KEYS:
        st.session_state.pop(k, None)
    st.session_state.pop("_ht_computed", None)

def _on_volume_change():
    for k in _HT_AREA_KEYS:
        st.session_state.pop(k, None)

# =====================================================================
# 1 · REACTOR SELECTION
# =====================================================================
st.header("1 · Reactor & Fluid Selection")

col_r, col_f = st.columns(2)
with col_r:
    _reactor_list = reactors["reactor_name"].tolist()
    reactor_name = st.selectbox("Reactor", _reactor_list, key="ht_reactor",
                                 on_change=_on_reactor_change)
with col_f:
    fluid_name = st.selectbox("Process fluid", _all_fluid_names, key="ht_fluid",
                               on_change=_on_fluid_change)

reactor = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
_is_solvent = is_known_solvent(fluid_name)

col_T, col_P = st.columns(2)
with col_T:
    if _is_solvent:
        fluid_T_C = st.number_input("Reference temperature (°C)", value=25.0,
                                     step=1.0, format="%.1f", key="ht_T")
    else:
        st.caption("Custom fluid — fixed properties")
        fluid_T_C = 25.0
with col_P:
    if _is_solvent:
        fluid_P_atm = st.number_input("Pressure (atm)", value=1.0,
                                       min_value=0.01, step=0.1, format="%.2f",
                                       key="ht_P")
    else:
        fluid_P_atm = 1.0

# Compute fluid properties
if _is_solvent:
    _fprops = get_properties(fluid_name, fluid_T_C, fluid_P_atm)
    _rho = _fprops["rho_kg_m3"]
    _mu = _fprops["mu_Pa_s"]
    _Cp = _fprops["Cp_J_per_kgK"]
    _k_fluid = _fprops["k_W_per_mK"]
    _sigma = _fprops.get("surface_tension_N_m", 0.072)
    if not _fprops["in_range"]:
        st.warning(f"⚠️ {fluid_T_C:.1f} °C outside liquid range for {fluid_name}. "
                   "Values are extrapolated.")
else:
    _cust = custom_fluids[custom_fluids["fluid_name"] == fluid_name].iloc[0]
    _rho = float(_cust["rho_kg_m3"])
    _mu = float(_cust["mu_Pa_s"])
    _Cp = float(_cust.get("Cp_J_per_kgK", 4182.0))
    _k_fluid = float(_cust.get("k_W_per_mK", 0.607))
    if np.isnan(_Cp) or _Cp <= 0:
        _Cp = 4182.0
    if np.isnan(_k_fluid) or _k_fluid <= 0:
        _k_fluid = 0.607

# =====================================================================
# 2 · PARAMETER OVERRIDES
# =====================================================================
st.header("2 · Review / Override Parameters")

with st.expander("Reactor geometry & agitation", expanded=False):
    gc1, gc2, gc3 = st.columns(3)
    D_tank = gc1.number_input("D_tank (m)", value=_safe(reactor.get("D_tank_m"), 0.10),
                               format="%.4f", key="ht_Dt")
    D_imp = gc2.number_input("D_imp (m)", value=_safe(reactor.get("D_imp_m"), 0.05),
                              format="%.4f", key="ht_Di")
    _rpm_lo = _safe(reactor.get("N_rpm_min"), 0.0)
    _rpm_hi = _safe(reactor.get("N_rpm_max"), 0.0)
    _rpm_default = (_rpm_lo + _rpm_hi) / 2 if _rpm_lo > 0 and _rpm_hi > 0 else (
        _rpm_hi if _rpm_hi > 0 else _safe(reactor.get("N_rps"), 5.0) * 60)
    N_rpm = gc3.number_input("N (RPM)", value=_rpm_default, step=2.0,
                              format="%.0f", key="ht_N")
    N_rps = N_rpm / 60.0

    gc4, gc5 = st.columns(2)
    Np_in = gc4.number_input("Np (power number)", value=_safe(reactor.get("Np"), 1.27),
                              format="%.2f", key="ht_Np")

# Volume
V_L_min = _safe(reactor.get("V_L_min"), 0.0)
V_L_max = _safe(reactor.get("V_L_max"), 0.0)
V_L_nom = _safe(reactor.get("V_L"), 0.0)
H_max = _safe(reactor.get("H_max_m"), _safe(reactor.get("H_m"), 0.13))
_bottom_dish = str(reactor.get("bottom_dish", "")) if pd.notna(reactor.get("bottom_dish")) else ""

# Compute the correct default volume for this reactor
if V_L_min > 0 and V_L_max > 0:
    V_L_default = round((V_L_min + V_L_max) / 2.0, 2)
elif V_L_nom > 0:
    V_L_default = round(V_L_nom, 2)
else:
    V_L_default = round(np.pi / 4 * D_tank**2 * _safe(reactor.get("H_m"), 0.13) * 1000, 2)

# Force-update volume when reactor changes
if st.session_state.pop("_ht_reactor_changed", False):
    st.session_state["ht_VL"] = V_L_default

# Show volume range
if V_L_min > 0 and V_L_max > 0:
    _vol_range_text = f"Reactor volume range: **{V_L_min:.1f} – {V_L_max:.1f} L** · Default (midpoint): **{V_L_default:.1f} L**"
elif V_L_nom > 0:
    _vol_range_text = f"Nominal volume: **{V_L_nom:.1f} L** (no min/max defined)"
else:
    _vol_range_text = "No volume range defined — computed from geometry."
st.caption(_vol_range_text)

V_L = st.number_input(
    "Liquid volume (L)", value=V_L_default, min_value=0.01,
    format="%.2f", key="ht_VL",
    on_change=_on_volume_change,
)
V_L_m3 = V_L / 1000.0

H = liquid_height_from_volume(V_L, D_tank, H_max, _bottom_dish)

with st.expander("Fluid properties (process side)", expanded=False):
    if _is_solvent:
        st.caption(f"Computed from **{fluid_name}** at **{fluid_T_C:.1f} °C**")
    fc1, fc2, fc3, fc4 = st.columns(4)
    rho = fc1.number_input("ρ (kg/m³)", value=_rho, format="%.1f", key="ht_rho")
    mu = fc2.number_input("μ (Pa·s)", value=_mu, format="%.6f", key="ht_mu")
    Cp = fc3.number_input("Cp (J/(kg·K))", value=_Cp, format="%.1f", key="ht_Cp")
    k_fluid = fc4.number_input("k (W/(m·K))", value=_k_fluid, format="%.4f", key="ht_kf")

# =====================================================================
# 3 · VESSEL WALL & LINING
# =====================================================================
st.header("3 · Vessel Wall & Jacket Configuration")

with st.expander("Wall and lining", expanded=True):
    wc1, wc2, wc3 = st.columns(3)
    _mat_options = list(WALL_CONDUCTIVITY.keys())
    _default_mat = str(reactor.get("material", "")) if pd.notna(reactor.get("material")) else ""
    if not _default_mat:
        _default_mat = str(reactor.get("shell_material", "")) if pd.notna(reactor.get("shell_material")) else ""
    _mat_idx = 0
    for i, m in enumerate(_mat_options):
        if _default_mat.lower() in m.lower() or m.lower() in _default_mat.lower():
            _mat_idx = i
            break
    wall_material = wc1.selectbox("Wall material", _mat_options,
                                   index=_mat_idx, key="ht_wall_mat")
    wall_k = WALL_CONDUCTIVITY.get(wall_material, 16.0)

    _wall_mm_default = _safe(reactor.get("wall_thickness_mm"), 5.0)
    wall_mm = wc2.number_input("Wall thickness (mm)", value=_wall_mm_default,
                                min_value=0.0, format="%.1f", key="ht_wall_mm")
    wall_m = wall_mm / 1000.0

    _lining_options = ["None"] + list(LINING_CONDUCTIVITY.keys())
    _default_lining = str(reactor.get("lining_material", "")) if pd.notna(reactor.get("lining_material")) else ""
    _lining_idx = 0
    for i, ln in enumerate(_lining_options):
        if _default_lining.lower() in ln.lower() or ln.lower() in _default_lining.lower():
            _lining_idx = i
            break
    lining_material = wc3.selectbox("Lining material", _lining_options,
                                     index=_lining_idx, key="ht_lining")

    lining_k = 0.0
    lining_m = 0.0
    if lining_material != "None":
        _linfo = _lookup_lining_k(lining_material)
        if _linfo:
            lining_k, _lt_default = _linfo
            lc1, lc2 = st.columns(2)
            lining_mm = lc1.number_input("Lining thickness (mm)",
                                          value=_lt_default * 1000,
                                          min_value=0.1, format="%.1f",
                                          key="ht_lining_mm")
            lining_m = lining_mm / 1000.0
            lc2.metric("k_lining (W/(m·K))", f"{lining_k:.3f}")

    fouling_R = st.number_input(
        "Fouling resistance R_f (m²·K/W)", value=FOULING_DEFAULT,
        min_value=0.0, format="%.5f", key="ht_fouling",
        help="Typical: 0.0001 (clean) – 0.001 (heavy fouling)",
    )

# =====================================================================
# 4 · HEAT TRANSFER MEDIUM
# =====================================================================
st.header("4 · Heat Transfer Medium (Utility)")

htm_names = list(HTM_DB.keys())
htm_name = st.selectbox("Heat transfer medium", htm_names,
                         index=0, key="ht_htm")
htm = HTM_DB[htm_name]

with st.expander("HTM properties", expanded=False):
    st.caption(f"**{htm_name}** — {htm.get('notes', '')}")
    st.caption(f"Operating range: {htm['T_min_C']:.0f} to {htm['T_max_C']:.0f} °C")

    hm1, hm2, hm3, hm4 = st.columns(4)
    htm_rho = hm1.number_input("ρ_htm (kg/m³)", value=htm["rho_kg_m3"],
                                format="%.1f", key="ht_htm_rho")
    htm_Cp = hm2.number_input("Cp_htm (J/(kg·K))", value=htm["Cp_J_kgK"],
                               format="%.1f", key="ht_htm_Cp")
    htm_mu = hm3.number_input("μ_htm (Pa·s)", value=htm["mu_Pa_s"],
                               format="%.6f", key="ht_htm_mu")
    htm_k = hm4.number_input("k_htm (W/(m·K))", value=htm["k_W_mK"],
                              format="%.4f", key="ht_htm_k")

_has_override = "h_jacket_override" in htm
jc1, jc2, jc3 = st.columns(3)
if _has_override:
    h_o_default = htm["h_jacket_override"]
    jc1.metric("h_o (jacket side)", f"{h_o_default:.0f} W/(m²·K)")
    jc1.caption("Fixed for condensing steam")
    jacket_type = "condensing"
    v_jacket = 0.0
    D_hyd_jacket = 0.05
else:
    v_jacket = jc1.number_input("Jacket velocity (m/s)", value=1.0,
                                 min_value=0.01, format="%.2f", key="ht_vj",
                                 help="Typical: 0.5–2.0 m/s for liquid utilities")
    D_hyd_jacket = jc2.number_input("Jacket hydraulic dia (m)", value=0.05,
                                     min_value=0.005, format="%.3f", key="ht_Dhyd",
                                     help="Annular gap or half-pipe ID")
    jacket_type = jc3.selectbox("Jacket type",
                                 ["Simple jacket", "Half-pipe coil", "Dimple jacket"],
                                 key="ht_jtype")

# Jacket flow rate for variable-jacket model
m_dot_jacket = st.number_input(
    "Jacket mass flow rate (kg/s)", value=1.0, min_value=0.01,
    format="%.2f", key="ht_mdot",
    help="Mass flow rate of utility fluid through the jacket",
)

# =====================================================================
# 5 · NUSSELT CORRELATION SELECTION
# =====================================================================
st.header("5 · Nusselt Correlation (Process Side)")

_corr_names = list(NUSSELT_CORRELATIONS.keys())
nu_corr = st.selectbox("Nusselt correlation", _corr_names,
                        index=0, key="ht_nu_corr")
_corr_info = NUSSELT_CORRELATIONS[nu_corr]
st.caption(f"**{nu_corr}**: {_corr_info['description']}")
st.caption(f"Reference: {_corr_info['ref']}")

mu_wall = st.number_input(
    "μ at wall temperature (Pa·s)", value=0.0, min_value=0.0,
    format="%.6f", key="ht_mu_wall",
    help="0 = assume μ_wall = μ_bulk (no wall viscosity correction)",
)

# =====================================================================
# 6 · HEAT TRANSFER AREA
# =====================================================================
st.header("6 · Heat Transfer Area")

_A_override = _safe(reactor.get("A_ht_m2"), 0.0)
_A_calc = estimate_jacket_area(D_tank, H, _bottom_dish)

ac1, ac2 = st.columns(2)
with ac1:
    A_source = st.radio("Area source",
                         ["Computed from geometry", "Manual / database override"],
                         horizontal=True, key="ht_A_source")
with ac2:
    if A_source.startswith("Computed"):
        A_ht = _A_calc
        st.metric("A_ht (computed)", f"{A_ht:.4f} m²")
        st.caption("Auto-updated from liquid volume and vessel geometry.")
    else:
        A_ht = st.number_input("A_ht (m²)",
                                value=_A_override if _A_override > 0 else _A_calc,
                                min_value=0.001, format="%.4f", key="ht_A_ov")

# =====================================================================
# 7 · OPERATING CONDITIONS
# =====================================================================
st.header("7 · Operating Conditions")

oc1, oc2, oc3 = st.columns(3)
T_start = oc1.number_input("T_start (°C)", value=25.0, step=1.0,
                            format="%.1f", key="ht_Tstart")
T_target = oc2.number_input("T_target (°C)", value=5.0, step=1.0,
                             format="%.1f", key="ht_Ttarget")
T_jacket_in = oc3.number_input("T_jacket inlet (°C)", value=-10.0, step=1.0,
                                format="%.1f", key="ht_Tj")

oc4, oc5 = st.columns(2)
P_agitator_W = 0.0
with oc4:
    include_agitator_heat = st.checkbox("Include agitator power as heat input",
                                         value=True, key="ht_inc_agit")
    if include_agitator_heat:
        P_agitator_W = impeller_power(Np_in, rho, N_rps, D_imp)
        st.caption(f"P_agitator = {P_agitator_W:.2f} W")
with oc5:
    Q_rxn = st.number_input("Additional heat input Q_rxn (W)", value=0.0,
                              step=1.0, format="%.1f", key="ht_Qrxn",
                              help="Exothermic reaction heat (positive = adds heat to batch)")

# Mode
_is_cooling = T_target < T_start
_mode_label = "Cooling" if _is_cooling else "Heating"
st.info(f"**{_mode_label}** from {T_start:.1f} °C → {T_target:.1f} °C "
        f"| Jacket inlet: {T_jacket_in:.1f} °C")

# Feasibility check
if _is_cooling and T_jacket_in >= T_start:
    st.error("Jacket temperature must be below batch start temperature for cooling.")
    st.stop()
if not _is_cooling and T_jacket_in <= T_start:
    st.error("Jacket temperature must be above batch start temperature for heating.")
    st.stop()
if _is_cooling and T_jacket_in >= T_target:
    pass  # OK, can reach target
elif not _is_cooling and T_jacket_in <= T_target:
    pass
else:
    # Check if target is reachable
    if _is_cooling and T_target < T_jacket_in:
        st.error(f"Cannot cool below jacket temperature ({T_jacket_in:.1f} °C). "
                 f"Target {T_target:.1f} °C is unreachable.")
        st.stop()
    if not _is_cooling and T_target > T_jacket_in:
        st.error(f"Cannot heat above jacket temperature ({T_jacket_in:.1f} °C). "
                 f"Target {T_target:.1f} °C is unreachable.")
        st.stop()

# =====================================================================
# 8 · COMPUTE RESULTS
# =====================================================================
if st.button("🔬 Compute Heat Transfer", type="primary", key="ht_compute"):
    st.session_state["_ht_computed"] = True

if st.session_state.get("_ht_computed"):
    st.divider()
    st.header("8 · Results")

    # ── Process-side h_i ──────────────────────────────────────────────
    Re = rho * N_rps * D_imp**2 / mu if mu > 0 else 0.0
    Pr = Cp * mu / k_fluid if k_fluid > 0 else 0.0
    mu_r = mu / mu_wall if mu_wall > 0 else 1.0
    Nu = nusselt_jacket(Re, Pr, mu_r, nu_corr)
    h_i = Nu * k_fluid / D_tank if D_tank > 0 else 0.0

    # ── Jacket-side h_o ──────────────────────────────────────────────
    if _has_override:
        h_o = htm["h_jacket_override"]
    else:
        # Compute using Dittus-Boelter for jacket side
        if htm_mu > 0 and htm_k > 0 and D_hyd_jacket > 0 and v_jacket > 0:
            Re_j = htm_rho * v_jacket * D_hyd_jacket / htm_mu
            Pr_j = htm_Cp * htm_mu / htm_k
            if Re_j < 2300:
                # Laminar: Sieder-Tate approximation
                Nu_j = 3.66 + 0.065 * (D_hyd_jacket / 1.0) * Re_j * Pr_j / (
                    1.0 + 0.04 * ((D_hyd_jacket / 1.0) * Re_j * Pr_j) ** (2.0/3.0))
            else:
                Nu_j = 0.023 * Re_j**0.8 * Pr_j**0.4
            h_o = Nu_j * htm_k / D_hyd_jacket
        else:
            h_o = JACKET_HTC_DEFAULT

    # ── Overall U ─────────────────────────────────────────────────────
    U = estimate_U_from_resistances(
        h_i=h_i, h_o=h_o,
        wall_k=wall_k, wall_thickness_m=wall_m,
        lining_k=lining_k, lining_thickness_m=lining_m,
        fouling=fouling_R,
    )

    # ── Individual resistances ────────────────────────────────────────
    R_i = 1.0 / h_i if h_i > 0 else np.inf
    R_wall = wall_m / wall_k if wall_k > 0 and wall_m > 0 else 0.0
    R_lining = lining_m / lining_k if lining_k > 0 and lining_m > 0 else 0.0
    R_o = 1.0 / h_o if h_o > 0 else np.inf
    R_total = R_i + R_wall + R_lining + R_o + fouling_R

    # ── Metrics ──────────────────────────────────────────────────────
    st.subheader("Heat Transfer Coefficients")

    r1, r2, r3, r4 = st.columns(4)
    r1.metric("h_i (process side)", f"{h_i:.1f} W/(m²·K)")
    r2.metric("h_o (jacket side)", f"{h_o:.1f} W/(m²·K)")
    r3.metric("U (overall)", f"{U:.1f} W/(m²·K)")
    r4.metric("Nu (process side)", f"{Nu:.1f}")

    r5, r6, r7, r8 = st.columns(4)
    r5.metric("Re (impeller)", f"{Re:.0f}")
    r6.metric("Pr (process)", f"{Pr:.1f}")
    r7.metric("A_ht (m²)", f"{A_ht:.4f}")
    r8.metric("P_agitator (W)", f"{P_agitator_W:.2f}")

    # ── Resistance breakdown ──────────────────────────────────────────
    st.subheader("Thermal Resistance Breakdown")

    # Build resistance items with fixed colours and descriptive labels
    _COLOR_MAP = {
        "convection": "#4A90D9",   # blue  – process side
        "wall":       "#6C757D",   # grey  – vessel wall
        "lining":     "#E8A838",   # amber – lining
        "jacket":     "#28A745",   # green – jacket side
        "fouling":    "#DC3545",   # red   – fouling
    }

    _wall_label = (f"Wall – {wall_material}\n"
                   f"({wall_mm:.1f} mm, k = {wall_k:.1f} W/(m·K))")
    _lining_label = ""
    if lining_material != "None" and lining_m > 0:
        _lining_label = (f"Lining – {lining_material}\n"
                         f"({lining_mm:.1f} mm, k = {lining_k:.3f} W/(m·K))")

    # Ordered list: (label, value, colour_key, group)
    _res_items = [
        ("Process film\n(1/h_i)",   R_i,       "convection", "Convective"),
        (_wall_label,               R_wall,    "wall",       "Conductive"),
    ]
    if lining_material != "None" and R_lining > 0:
        _res_items.append((_lining_label, R_lining, "lining", "Conductive"))
    _res_items += [
        ("Jacket film\n(1/h_o)",    R_o,       "jacket",     "Convective"),
        ("Fouling\n(R_f)",          fouling_R, "fouling",    "Other"),
    ]

    # Filter out zero / infinite entries
    _res_items = [(lbl, val, ck, grp) for lbl, val, ck, grp in _res_items
                  if val > 0 and val < np.inf]

    if _res_items:
        _labels = [lbl for lbl, *_ in _res_items]
        _values = [val for _, val, *_ in _res_items]
        _colors = [_COLOR_MAP[ck] for _, _, ck, _ in _res_items]
        _groups = [grp for *_, grp in _res_items]

        fig_res = go.Figure(go.Bar(
            x=_labels,
            y=_values,
            marker_color=_colors,
            text=[f"{v:.5f}" for v in _values],
            textposition="auto",
        ))

        # Add group annotations at the top
        _seen_groups = {}
        for i, grp in enumerate(_groups):
            _seen_groups.setdefault(grp, []).append(i)

        for grp, indices in _seen_groups.items():
            x0, x1 = min(indices) - 0.4, max(indices) + 0.4
            fig_res.add_shape(
                type="rect", xref="x", yref="paper",
                x0=x0, x1=x1, y0=1.0, y1=1.08,
                fillcolor="rgba(200,200,200,0.3)", line_width=0,
            )
            fig_res.add_annotation(
                x=(min(indices) + max(indices)) / 2, y=1.04,
                xref="x", yref="paper",
                text=f"<b>{grp}</b>", showarrow=False,
                font=dict(size=11),
            )

        fig_res.update_layout(
            title="Thermal Resistance Breakdown (m²·K/W)",
            yaxis_title="R (m²·K/W)",
            height=380,
            margin=dict(t=80),
        )
        st.plotly_chart(fig_res, use_container_width=True)

        # Percentage breakdown — show wall and lining separately
        _pct_parts = []
        for lbl, val, ck, _grp in _res_items:
            pct = val / R_total * 100
            # Use short inline label
            short = lbl.split("\n")[0]
            _pct_parts.append(f"{short}: **{pct:.1f}%**")
        st.caption("**Percentage of total resistance:**")
        st.markdown(" · ".join(_pct_parts))

        # Combined wall + lining metric when both present
        if lining_material != "None" and R_lining > 0:
            _R_cond = R_wall + R_lining
            _pct_cond = _R_cond / R_total * 100
            st.caption(f"Wall + Lining combined: **{_R_cond:.5f} m²·K/W** "
                       f"(**{_pct_cond:.1f}%** of total)")

        # Identify controlling resistance
        _max_idx = max(range(len(_values)), key=lambda i: _values[i])
        _max_label = _res_items[_max_idx][0].split("\n")[0]
        _max_pct = _values[_max_idx] / R_total * 100
        st.info(f"**Controlling resistance:** {_max_label} "
                f"({_max_pct:.1f}% of total)")

    # ── Steady-state heat duty ────────────────────────────────────────
    st.subheader("Steady-State Heat Duty")
    _dT_init = abs(T_start - T_jacket_in)
    Q_max = heat_removal_capacity(U, A_ht, _dT_init)
    _dT_dt = cooling_rate(Q_max, P_agitator_W, rho, V_L_m3, Cp)

    ss1, ss2, ss3 = st.columns(3)
    ss1.metric("Q_max (initial)", f"{Q_max:.1f} W")
    ss2.metric("Initial ΔT", f"{_dT_init:.1f} °C")
    ss3.metric("Initial dT/dt", f"{_dT_dt * 60:.3f} °C/min")

    # ── Analytical time estimate (log-mean) ───────────────────────────
    st.subheader("Time Estimates")

    t_log = time_to_cool_or_heat(rho, V_L_m3, Cp, U, A_ht,
                                  T_start, T_target, T_jacket_in)

    te1, te2 = st.columns(2)
    if t_log < np.inf:
        te1.metric("Analytical time (log-mean, no agitator heat)",
                    f"{t_log:.0f} s ({t_log / 60:.1f} min)")
    else:
        te1.metric("Analytical time", "∞ (unreachable)")

    # ── Transient temperature profiles ────────────────────────────────
    st.subheader("Temperature vs. Time Profiles")

    _dt_sim = max(0.5, t_log / 2000) if t_log < np.inf else 1.0
    _t_max_sim = min(t_log * 2.0 if t_log < np.inf else 36000.0, 86400.0)

    # Model 1: Constant jacket temperature (isothermal utility)
    t1, T1 = batch_temperature_profile(
        rho=rho, V_L_m3=V_L_m3, Cp=Cp,
        U=U, A=A_ht,
        T_start=T_start, T_target=T_target, T_jacket=T_jacket_in,
        P_agitator=P_agitator_W, Q_rxn=Q_rxn,
        dt=_dt_sim, t_max=_t_max_sim,
    )

    # Model 2: Variable jacket (finite flow rate)
    t2, T2, Tj2 = batch_temp_profile_variable_jacket(
        rho=rho, V_L_m3=V_L_m3, Cp=Cp,
        U=U, A=A_ht,
        T_start=T_start, T_target=T_target,
        T_jacket_in=T_jacket_in,
        m_dot_jacket=m_dot_jacket,
        Cp_jacket=htm_Cp,
        P_agitator=P_agitator_W, Q_rxn=Q_rxn,
        dt=_dt_sim, t_max=_t_max_sim,
    )

    # Time to target from simulations
    _t_target_1 = t1[-1] if len(t1) > 1 else np.inf
    _t_target_2 = t2[-1] if len(t2) > 1 else np.inf

    te2.metric("Simulated time (const. jacket)",
               f"{_t_target_1:.0f} s ({_t_target_1 / 60:.1f} min)")

    te3, te4 = st.columns(2)
    te3.metric("Simulated time (variable jacket)",
               f"{_t_target_2:.0f} s ({_t_target_2 / 60:.1f} min)")

    # Unit selection
    _time_unit = st.radio("Time axis unit", ["Minutes", "Seconds", "Hours"],
                           horizontal=True, key="ht_tunit")
    _t_factor = {"Seconds": 1.0, "Minutes": 60.0, "Hours": 3600.0}[_time_unit]
    _t_label = _time_unit.lower()

    # Main temperature profile chart
    fig_T = go.Figure()
    fig_T.add_trace(go.Scatter(
        x=t1 / _t_factor, y=T1,
        mode="lines", name="Batch (const. jacket)",
        line=dict(color="#4A90D9", width=2),
    ))
    fig_T.add_trace(go.Scatter(
        x=t2 / _t_factor, y=T2,
        mode="lines", name="Batch (variable jacket)",
        line=dict(color="#E8A838", width=2),
    ))
    fig_T.add_trace(go.Scatter(
        x=t2 / _t_factor, y=Tj2,
        mode="lines", name="Jacket outlet",
        line=dict(color="#28A745", width=1, dash="dash"),
    ))
    # Reference lines
    fig_T.add_hline(y=T_target, line_dash="dot", line_color="red",
                     annotation_text=f"Target {T_target:.1f} °C")
    fig_T.add_hline(y=T_jacket_in, line_dash="dot", line_color="grey",
                     annotation_text=f"Jacket inlet {T_jacket_in:.1f} °C")

    fig_T.update_layout(
        title="Batch Temperature Profile",
        xaxis_title=f"Time ({_t_label})",
        yaxis_title="Temperature (°C)",
        height=500,
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
    )
    st.plotly_chart(fig_T, use_container_width=True)

    # Heat duty vs time
    st.subheader("Heat Duty vs. Time")
    Q_jacket_1 = U * A_ht * (T_jacket_in - T1)
    Q_jacket_2 = np.zeros_like(t2)
    if m_dot_jacket > 0 and htm_Cp > 0:
        NTU = U * A_ht / (m_dot_jacket * htm_Cp)
        eff = 1.0 - np.exp(-NTU)
        Q_jacket_2 = eff * m_dot_jacket * htm_Cp * (T_jacket_in - T2)
    else:
        Q_jacket_2 = U * A_ht * (T_jacket_in - T2)

    fig_Q = go.Figure()
    fig_Q.add_trace(go.Scatter(
        x=t1 / _t_factor, y=np.abs(Q_jacket_1),
        mode="lines", name="|Q_jacket| (const. jacket)",
        line=dict(color="#4A90D9", width=2),
    ))
    fig_Q.add_trace(go.Scatter(
        x=t2 / _t_factor, y=np.abs(Q_jacket_2),
        mode="lines", name="|Q_jacket| (variable jacket)",
        line=dict(color="#E8A838", width=2),
    ))
    fig_Q.update_layout(
        title="Jacket Heat Duty over Time",
        xaxis_title=f"Time ({_t_label})",
        yaxis_title="|Q| (W)",
        height=400,
        hovermode="x unified",
    )
    st.plotly_chart(fig_Q, use_container_width=True)

    # dT/dt vs time
    st.subheader("Cooling / Heating Rate vs. Time")
    _m_Cp = rho * V_L_m3 * Cp
    if _m_Cp > 0:
        dTdt_1 = (Q_jacket_1 + P_agitator_W + Q_rxn) / _m_Cp * 60  # °C/min
        dTdt_2 = (Q_jacket_2 + P_agitator_W + Q_rxn) / _m_Cp * 60
    else:
        dTdt_1 = np.zeros_like(t1)
        dTdt_2 = np.zeros_like(t2)

    fig_rate = go.Figure()
    fig_rate.add_trace(go.Scatter(
        x=t1 / _t_factor, y=dTdt_1,
        mode="lines", name="dT/dt (const. jacket)",
        line=dict(color="#4A90D9", width=2),
    ))
    fig_rate.add_trace(go.Scatter(
        x=t2 / _t_factor, y=dTdt_2,
        mode="lines", name="dT/dt (variable jacket)",
        line=dict(color="#E8A838", width=2),
    ))
    fig_rate.update_layout(
        title="Batch Temperature Rate of Change",
        xaxis_title=f"Time ({_t_label})",
        yaxis_title="dT/dt (°C/min)",
        height=400,
        hovermode="x unified",
    )
    st.plotly_chart(fig_rate, use_container_width=True)

    # ── RPM sensitivity ───────────────────────────────────────────────
    st.subheader("RPM Sensitivity")
    st.caption("Effect of impeller speed on U and time-to-target "
               "(constant jacket temperature model).")

    if _rpm_lo > 0 and _rpm_hi > 0:
        _n_sweep = 30
        rpm_arr = np.linspace(_rpm_lo, _rpm_hi, _n_sweep)
    else:
        rpm_arr = np.linspace(max(N_rpm * 0.3, 10), N_rpm * 1.5, 30)

    U_arr = np.empty_like(rpm_arr)
    t_arr_rpm = np.empty_like(rpm_arr)
    hi_arr = np.empty_like(rpm_arr)

    for i, _rpm in enumerate(rpm_arr):
        _N = _rpm / 60.0
        _Re = rho * _N * D_imp**2 / mu if mu > 0 else 0.0
        _Pr = Cp * mu / k_fluid if k_fluid > 0 else 0.0
        _mu_r = mu / mu_wall if mu_wall > 0 else 1.0
        _Nu = nusselt_jacket(_Re, _Pr, _mu_r, nu_corr)
        _hi = _Nu * k_fluid / D_tank if D_tank > 0 else 0.0
        _U = estimate_U_from_resistances(
            h_i=_hi, h_o=h_o,
            wall_k=wall_k, wall_thickness_m=wall_m,
            lining_k=lining_k, lining_thickness_m=lining_m,
            fouling=fouling_R,
        )
        U_arr[i] = _U
        hi_arr[i] = _hi
        t_arr_rpm[i] = time_to_cool_or_heat(rho, V_L_m3, Cp, _U, A_ht,
                                              T_start, T_target, T_jacket_in)

    fig_rpm = go.Figure()
    fig_rpm.add_trace(go.Scatter(
        x=rpm_arr, y=U_arr,
        mode="lines", name="U (W/(m²·K))",
        line=dict(color="#4A90D9", width=2),
    ))
    fig_rpm.add_trace(go.Scatter(
        x=rpm_arr, y=hi_arr,
        mode="lines", name="h_i (W/(m²·K))",
        line=dict(color="#28A745", width=2, dash="dash"),
    ))
    fig_rpm.add_vline(x=N_rpm, line_dash="dot", line_color="red",
                       annotation_text=f"Current {N_rpm:.0f} RPM")
    fig_rpm.update_layout(
        title="Heat Transfer Coefficients vs. RPM",
        xaxis_title="RPM",
        yaxis_title="Coefficient (W/(m²·K))",
        height=400,
        hovermode="x unified",
    )
    st.plotly_chart(fig_rpm, use_container_width=True)

    fig_trpm = go.Figure()
    fig_trpm.add_trace(go.Scatter(
        x=rpm_arr, y=t_arr_rpm / 60.0,
        mode="lines", name="Time to target (min)",
        line=dict(color="#E8A838", width=2),
    ))
    fig_trpm.add_vline(x=N_rpm, line_dash="dot", line_color="red",
                        annotation_text=f"Current {N_rpm:.0f} RPM")
    fig_trpm.update_layout(
        title="Time to Target vs. RPM",
        xaxis_title="RPM",
        yaxis_title="Time (min)",
        height=400,
        hovermode="x unified",
    )
    st.plotly_chart(fig_trpm, use_container_width=True)

    # ── Nusselt correlation comparison ────────────────────────────────
    st.subheader("Nusselt Correlation Comparison")
    st.caption("Compares all available correlations at the current operating point.")

    _comp_rows = []
    for _cn, _cv in NUSSELT_CORRELATIONS.items():
        _Nu_c = nusselt_jacket(Re, Pr, mu_r, _cn)
        _hi_c = _Nu_c * k_fluid / D_tank if D_tank > 0 else 0.0
        _U_c = estimate_U_from_resistances(
            h_i=_hi_c, h_o=h_o,
            wall_k=wall_k, wall_thickness_m=wall_m,
            lining_k=lining_k, lining_thickness_m=lining_m,
            fouling=fouling_R,
        )
        _t_c = time_to_cool_or_heat(rho, V_L_m3, Cp, _U_c, A_ht,
                                     T_start, T_target, T_jacket_in)
        _comp_rows.append({
            "Correlation": _cn,
            "Nu": f"{_Nu_c:.1f}",
            "h_i (W/(m²·K))": f"{_hi_c:.1f}",
            "U (W/(m²·K))": f"{_U_c:.1f}",
            "Time (min)": f"{_t_c / 60:.1f}" if _t_c < np.inf else "∞",
        })

    st.dataframe(pd.DataFrame(_comp_rows), use_container_width=True, hide_index=True)

    # ── HTM comparison ────────────────────────────────────────────────
    st.subheader("Heat Transfer Media Comparison")
    st.caption("Compares all available media at the current process conditions.")

    _htm_rows = []
    for _hname, _hdata in HTM_DB.items():
        if "h_jacket_override" in _hdata:
            _ho_c = _hdata["h_jacket_override"]
        elif _hdata["mu_Pa_s"] > 0 and _hdata["k_W_mK"] > 0:
            _Re_hc = _hdata["rho_kg_m3"] * v_jacket * D_hyd_jacket / _hdata["mu_Pa_s"]
            _Pr_hc = _hdata["Cp_J_kgK"] * _hdata["mu_Pa_s"] / _hdata["k_W_mK"]
            if _Re_hc < 2300:
                _Nu_hc = 3.66
            else:
                _Nu_hc = 0.023 * _Re_hc**0.8 * _Pr_hc**0.4
            _ho_c = _Nu_hc * _hdata["k_W_mK"] / D_hyd_jacket
        else:
            _ho_c = JACKET_HTC_DEFAULT

        _U_hc = estimate_U_from_resistances(
            h_i=h_i, h_o=_ho_c,
            wall_k=wall_k, wall_thickness_m=wall_m,
            lining_k=lining_k, lining_thickness_m=lining_m,
            fouling=fouling_R,
        )
        _t_hc = time_to_cool_or_heat(rho, V_L_m3, Cp, _U_hc, A_ht,
                                      T_start, T_target, T_jacket_in)
        _in_range = (_hdata["T_min_C"] <= T_jacket_in <= _hdata["T_max_C"])
        _htm_rows.append({
            "Medium": _hname,
            "h_o (W/(m²·K))": f"{_ho_c:.0f}",
            "U (W/(m²·K))": f"{_U_hc:.1f}",
            "Time (min)": f"{_t_hc / 60:.1f}" if _t_hc < np.inf else "∞",
            "T range (°C)": f"{_hdata['T_min_C']:.0f} to {_hdata['T_max_C']:.0f}",
            "In range?": "✅" if _in_range else "⚠️",
        })

    st.dataframe(pd.DataFrame(_htm_rows), use_container_width=True, hide_index=True)

    # ── Summary table ─────────────────────────────────────────────────
    st.subheader("Full Results Summary")
    _summary = {
        "Reactor": reactor_name,
        "Fluid": fluid_name,
        "Volume (L)": V_L,
        "RPM": N_rpm,
        "T_start (°C)": T_start,
        "T_target (°C)": T_target,
        "T_jacket (°C)": T_jacket_in,
        "HTM": htm_name,
        "Nusselt correlation": nu_corr,
        "Re": f"{Re:.0f}",
        "Pr": f"{Pr:.1f}",
        "Nu": f"{Nu:.1f}",
        "h_i (W/(m²·K))": f"{h_i:.1f}",
        "h_o (W/(m²·K))": f"{h_o:.1f}",
        "U (W/(m²·K))": f"{U:.1f}",
        "A_ht (m²)": f"{A_ht:.4f}",
        "Q_max initial (W)": f"{Q_max:.1f}",
        "Initial dT/dt (°C/min)": f"{_dT_dt * 60:.3f}",
        "P_agitator (W)": f"{P_agitator_W:.2f}",
        "Time analytical (min)": f"{t_log / 60:.1f}" if t_log < np.inf else "∞",
        "Time simulated const. jacket (min)": f"{_t_target_1 / 60:.1f}",
        "Time simulated var. jacket (min)": f"{_t_target_2 / 60:.1f}",
        "Wall material": wall_material,
        "Wall thickness (mm)": wall_mm,
        "Lining": lining_material,
        "Fouling R (m²·K/W)": fouling_R,
    }
    st.dataframe(pd.DataFrame([_summary]).T.rename(columns={0: "Value"}),
                 use_container_width=True)
