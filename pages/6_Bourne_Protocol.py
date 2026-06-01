"""
Page 8 – Bourne Protocol
====================================================
Guides the user step-by-step through the Bourne mixing-sensitivity screening
protocol as described by Bourne (2003).

The protocol systematically varies:
  Test 1 – Impeller speed  (does mixing matter at all?)
  Test 2 – Feed rate/time  (mesomixing vs micromixing)
  Test 3 – Feed location   (mesomixing vs macromixing)
  Confirmatory – # feed points / viscosity change

Based on process-response changes at each level, the decision tree identifies
which mixing scale (macro / meso / micro) is controlling.

References
----------
- Bourne, J.R. (2003). Mixing and the Selectivity of Chemical Reactions.
  Org. Process Res. Dev., 7(4), 471-508.
- Sarafinas, A. (2018). Test Process Mixing Sensitivities Using the Bourne
  Protocol. Scientific Update Webinar, 13 Nov 2018.
- Sarafinas, A. & Teich, C.I. (2016). Chapter 13 in Advances in Industrial
  Mixing (Kresta et al., Wiley).
"""

import streamlit as st

import pandas as pd
import numpy as np
import pathlib

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, is_known_solvent,
)
from utils.calculations import (
    compute_reactor_hydro,
    reynolds_number,
    power_number_correlation,
    impeller_power,
    power_per_volume,
    tip_speed,
    micromixing_time_engulfment,
    blend_time_turbulent,
    kolmogorov_length,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

from utils.data_helpers import load_db, safe_float, all_fluid_names, safe_iloc, reactor_search_name

reactors = load_db("reactor_db", "reactors.csv")
custom_fluids = load_db("fluid_db", "fluids.csv")

_all_fluid_names = all_fluid_names(custom_fluids)

# ══════════════════════════════════════════════════════════════════════════
st.title("🅱️ Bourne Protocol – Mixing Sensitivity Screening")
st.markdown("""
Determines **which mixing scale** (macro / meso / micro) controls your
process outcome. Applicable to any semi-batch stirred-tank process with
competitive rate processes.
""")

_BP_IMG = pathlib.Path(__file__).resolve().parent.parent / "images" / "general" / "bourne_protocol_decision_tree.png"
with st.expander("📋 Protocol overview – Decision Tree", expanded=False):
    if _BP_IMG.exists():
        st.image(str(_BP_IMG), caption="Bourne Protocol – Decision Tree")
    else:
        st.info("Decision tree image not found. Generate it from the ⚙️ Admin page.")

# ── Generation counter for reliable widget reset ─────────────────────────
if "_bp_gen" not in st.session_state:
    st.session_state["_bp_gen"] = 0

def _bk(name: str) -> str:
    """Return a widget key namespaced by the current generation counter."""
    return f"bp_{st.session_state['_bp_gen']}_{name}"

def _reset_bourne():
    """Increment generation counter and clear old Bourne Protocol state."""
    old_gen = st.session_state.get("_bp_gen", 0)
    # Remove all keys from the old generation (and non-generation legacy keys)
    for k in list(st.session_state.keys()):
        if k == "_bp_gen" or k == "bp_restart":
            continue
        if k.startswith("bp_") or k.startswith("_sel_bp_"):
            del st.session_state[k]
    st.session_state["_bp_gen"] = old_gen + 1


st.button("🔄 Restart protocol", key="bp_restart", on_click=_reset_bourne)

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 0 – Reactor & Fluid Selection
# ══════════════════════════════════════════════════════════════════════════
st.header("Define Your System")

# Persist selections across page navigations
_reactor_list = reactors["reactor_name"].tolist() if not reactors.empty else []

def _sel_idx(lst, key, default=0):
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col_r, col_f = st.columns(2)

with col_r:
    if not reactors.empty:
        reactor_name = st.selectbox("Reactor", _reactor_list, index=_sel_idx(_reactor_list, "_sel_bp_reactor"), key=_bk("reactor"),
                                    format_func=lambda n: reactor_search_name(reactors, n))
        st.session_state["_sel_bp_reactor"] = reactor_name
        r = safe_iloc(reactors, "reactor_name", reactor_name, "Reactor")
        D_tank = safe_float(r["D_tank_m"], 0.10)
        H = safe_float(r["H_m"], 0.13)
        D_imp = safe_float(r["D_imp_m"], 0.05)
        Np_val = safe_float(r["Np"], 1.27)
        Nq_val = safe_float(r["Nq"], 0.79)
        V_L_min = safe_float(r.get("V_L_min"), 0.0) or safe_float(r["V_L"], 0.0)
        V_L_max = safe_float(r.get("V_L_max"), 0.0) or safe_float(r["V_L"], 0.0)
        V_L_avg = (V_L_min + V_L_max) / 2.0
        # RPM bounds from reactor database
        N_rpm_min = safe_float(r.get("N_rpm_min"), 0.0) or None
        N_rpm_max = safe_float(r.get("N_rpm_max"), 0.0) or None
        # Use average of min/max RPM as the centerpoint speed
        if N_rpm_min is not None and N_rpm_max is not None:
            N_center = (N_rpm_min + N_rpm_max) / 2.0 / 60.0  # rev/s
        else:
            N_center = safe_float(r["N_rps"], 5.0)
        N_rps_min = N_rpm_min / 60.0 if N_rpm_min is not None else None
        N_rps_max = N_rpm_max / 60.0 if N_rpm_max is not None else None
    else:
        st.info("No reactors in database — enter manually below.")
        D_tank = st.number_input("Tank diameter (m)", value=0.10, format="%.4f")
        H = st.number_input("Liquid height (m)", value=0.13, format="%.4f")
        D_imp = st.number_input("Impeller diameter (m)", value=0.05, format="%.4f")
        N_center = st.number_input("Centerpoint speed (rev/s)", value=5.0, format="%.2f")
        Np_val = st.number_input("Power number Np", value=1.27, format="%.2f")
        Nq_val = st.number_input("Pumping number Nq", value=0.79, format="%.2f")
        V_L_avg = np.pi / 4 * D_tank**2 * H * 1000
        V_L_min = V_L_avg
        V_L_max = V_L_avg
        N_rpm_min = None
        N_rpm_max = None
        N_rps_min = None
        N_rps_max = None

col_t, col_p = st.columns(2)

with col_f:
    fluid_name = st.selectbox("Fluid system", _all_fluid_names,
                              index=_sel_idx(_all_fluid_names, "_sel_bp_fluid"),
                              key=_bk("fluid"))
    st.session_state["_sel_bp_fluid"] = fluid_name
    _is_solvent = is_known_solvent(fluid_name)
    if _is_solvent:
        _bp_T_C = col_t.number_input("Temperature (°C)", value=25.0, step=1.0,
                                    format="%.1f", key=_bk("fluid_T"))
        _bp_P_atm = col_p.number_input("Pressure (atm)", value=1.0, min_value=0.01,
                                     step=0.1, format="%.2f", key=_bk("fluid_P"))
        _fprops = get_properties(fluid_name, _bp_T_C, _bp_P_atm)
        rho = _fprops["rho_kg_m3"]
        mu = _fprops["mu_Pa_s"]
        if not _fprops["in_range"]:
            st.warning(f"⚠️ {_bp_T_C:.1f} °C is outside the liquid range "
                       f"({_fprops['mp_C']:.0f} – {_fprops['bp_at_P_C']:.0f} °C) "
                       f"for {fluid_name}. Values are extrapolated.")
    else:
        _cust = safe_iloc(custom_fluids, "fluid_name", fluid_name, "Custom fluid")
        rho = float(_cust["rho_kg_m3"])
        mu = float(_cust["mu_Pa_s"])

nu = mu / rho

# ── Volume selection ──────────────────────────────────────────────────────
st.subheader("Working Volume")
if V_L_min != V_L_max:
    st.caption(f"Reactor volume range: {V_L_min:.2f} – {V_L_max:.2f} L  •  Average: {V_L_avg:.2f} L")

_vol_min = V_L_min if V_L_min > 0 else 0.1
_reactor_tag = reactor_name.replace(" ", "_") if not reactors.empty else "manual"
V_L = st.number_input(
    "Working volume (L)", min_value=_vol_min, value=_vol_min,
    step=1.0, format="%.2f", key=_bk(f"vol_L_{_reactor_tag}"),
    help="Defaults to the minimum volume from the reactor database.",
)
V_m3 = V_L / 1000.0  # m³

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 1 – TEST 1: Impeller Speed  (Does mixing matter?)
# ══════════════════════════════════════════════════════════════════════════
st.header("Test 1 - Impeller Speed")
st.markdown("""
Vary **impeller speed only** (hold feed rate & location constant).\n
Target ≈**100× P/m range** across three agitation speeds.\n
If the response changes → mixing matters; proceed to Test 2.

**Speed selection:** Default centerpoint P/m ≈ **0.2 W/kg**. Set high/low at
**10×** and **0.1×** (speed ratios ≈ 2.15× and 0.46×, since P ∝ N³).
Use plant condition as centerpoint if targeting existing equipment.
""")

st.subheader("Suggested Experimental Conditions")

# Default centerpoint P/m from Sarafinas (2018); allow user to customize or use plant RPM
centerpoint_mode = st.radio(
    "Centerpoint selection method",
    ["Default (0.2 W/kg)", "Custom P/m", "Custom RPM"],
    horizontal=True, key=_bk("ctr_mode"),
)

if centerpoint_mode == "Custom P/m":
    PV_center_wkg = st.number_input("Centerpoint P/m (W/kg)", value=0.2, format="%.4g", key=_bk("pv_ctr"))
    _N_from_pv = (PV_center_wkg * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
    st.write(f"Custom centerpoint: **P/m = {PV_center_wkg:.4g} W/kg** "
             f"at N = {_N_from_pv:.2f} rev/s ({_N_from_pv*60:.0f} RPM)")
elif centerpoint_mode == "Custom RPM":
    _rpm_default = float(N_center * 60) if N_center else 300.0
    _custom_rpm = st.number_input("Centerpoint RPM", value=_rpm_default, min_value=1.0,
                                  format="%.1f", key=_bk("rpm_ctr"))
    _N_custom = _custom_rpm / 60.0
    _P_custom = impeller_power(Np_val, rho, _N_custom, D_imp)
    PV_center_wkg = power_per_volume(_P_custom, V_m3) / rho
    st.write(f"Custom centerpoint: **N = {_custom_rpm:.1f} RPM** → "
             f"**P/m = {PV_center_wkg:.4g} W/kg** "
             f"({power_per_volume(_P_custom, V_m3)/1000:.4g} W/L)")
    if N_rpm_min is not None and N_rpm_max is not None:
        if _custom_rpm < N_rpm_min or _custom_rpm > N_rpm_max:
            st.warning(f"⚠️ {_custom_rpm:.1f} RPM is outside the reactor range "
                       f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM).")
else:
    # Default centerpoint P/m = 0.2 W/kg (Sarafinas, 2018).
    # Check whether the reactor can achieve it within its RPM range.
    _DEFAULT_PV = 0.2  # W/kg
    _N_for_default = (_DEFAULT_PV * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
    _N_for_default_rpm = _N_for_default * 60

    _can_reach_default = True
    if N_rps_min is not None and _N_for_default < N_rps_min:
        _can_reach_default = False
    if N_rps_max is not None and _N_for_default > N_rps_max:
        _can_reach_default = False

    if _can_reach_default:
        PV_center_wkg = _DEFAULT_PV
        st.write(f"Centerpoint default: **P/m = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {_N_for_default_rpm:.0f} RPM)")
        if N_rpm_min is not None and N_rpm_max is not None:
            st.caption(f"Default of 0.2 W/kg is within reactor RPM range "
                       f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM).")
    else:
        # Fall back to average reactor speed
        _P_ctr = impeller_power(Np_val, rho, N_center, D_imp)
        PV_center_wkg = power_per_volume(_P_ctr, V_m3) / rho
        _ctr_rpm = N_center * 60
        st.write(f"Centerpoint from reactor average speed: **P/m = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {N_center:.2f} rev/s ({_ctr_rpm:.0f} RPM)")
        st.caption(f"Default of 0.2 W/kg requires {_N_for_default_rpm:.0f} RPM, "
                   f"which is outside the reactor range "
                   f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM). "
                   f"Using average speed instead.")

# Compute the three test speeds
# P/m ∝ N³  →  N ∝ (P/m)^(1/3)
N_center_calc = (PV_center_wkg * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
N_high_ideal = N_center_calc * 10**(1/3)   # ≈ 2.154×
N_low_ideal = N_center_calc * 0.1**(1/3)   # ≈ 0.464×

# Clamp to reactor min/max RPM bounds
N_low = N_low_ideal
N_high = N_high_ideal
_low_clamped = False
_high_clamped = False

if N_rps_min is not None and N_low < N_rps_min:
    N_low = N_rps_min
    _low_clamped = True
if N_rps_max is not None and N_high > N_rps_max:
    N_high = N_rps_max
    _high_clamped = True

# Build the conditions table
def _hydro_row(label, N):
    Re = reynolds_number(N, D_imp, rho, mu)
    Np = Np_val
    P = impeller_power(Np, rho, N, D_imp)
    eps = power_per_volume(P, V_m3)
    eps_kg = eps / rho                          # W/kg for Kolmogorov & micromixing
    u = tip_speed(N, D_imp)
    t_blend = blend_time_turbulent(Nq_val, V_m3, D_imp, N)
    t_micro = micromixing_time_engulfment(eps_kg, nu)
    eta = kolmogorov_length(nu, eps_kg)
    pv_wkg = eps_kg
    pv_rel = pv_wkg / PV_center_wkg if PV_center_wkg > 0 else 0
    return {
        "Condition": label,
        "Volume (L)": round(V_L, 3),
        # "N (rev/s)": round(N, 3),
        "N (RPM)": round(N * 60, 1),
        "P/m (W/kg)": round(pv_wkg, 4),
        "P/m rel. to center": f"{pv_rel:.2f}×",
        "P/V (W/L)": round(eps / 1000, 4),
        "Re": round(Re, 0),
        "Tip speed (m/s)": round(u, 3),
        "Blend time (s)": round(t_blend, 2),
        "t_E micro (s)": round(t_micro, 5),
        "η Kolmogorov (µm)": round(eta * 1e6, 1),
    }

_low_label = "Low  (min RPM)" if _low_clamped else "Low  (0.1× P/m)"
_high_label = "High (max RPM)" if _high_clamped else "High (10× P/m)"

t1_rows = [
    _hydro_row(_low_label, N_low),
    _hydro_row("Center (1× P/m)", N_center_calc),
    _hydro_row(_high_label, N_high),
]
t1_df = pd.DataFrame(t1_rows)
_t1_basic_cols = ["Condition", "Volume (L)", "N (RPM)", "P/m (W/kg)", "P/m rel. to center", "P/V (W/L)", "Tip speed (m/s)"]
st.dataframe(t1_df[_t1_basic_cols], use_container_width=True, hide_index=True)

with st.expander("Additional hydrodynamic parameters"):
    _t1_detail_cols = ["Condition", "Re", "Blend time (s)", "t_E micro (s)", "η Kolmogorov (µm)"]
    st.dataframe(t1_df[_t1_detail_cols], use_container_width=True, hide_index=True)

pv_ratio = t1_rows[2]["P/m (W/kg)"] / t1_rows[0]["P/m (W/kg)"] if t1_rows[0]["P/m (W/kg)"] > 0 else 0
st.caption(f"P/m ratio (high/low) = **{pv_ratio:.1f}×**  •  Speed ratio (high/low) = **{N_high/N_low:.2f}×**")

if _low_clamped or _high_clamped:
    _msgs = []
    if _low_clamped:
        _actual_low_rel = t1_rows[0]["P/m (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"Low speed clamped to reactor minimum ({N_rpm_min:.1f} RPM) → "
                     f"actual P/m = **{_actual_low_rel:.3f}×** centerpoint "
                     f"(target was 0.1×)")
    if _high_clamped:
        _actual_high_rel = t1_rows[2]["P/m (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"High speed clamped to reactor maximum ({N_rpm_max:.1f} RPM) → "
                     f"actual P/m = **{_actual_high_rel:.3f}×** centerpoint "
                     f"(target was 10×)")
    st.warning("⚠️ **Speed limits applied:**\n\n" + "\n\n".join(_msgs) +
              "\n\nThe achievable P/V range is narrower than the ideal 100× span. "
              "Consider whether this is sufficient to draw conclusions.")

with st.expander("Practical limits to consider"):
    st.markdown("""
    **Min speed:** N_js (solids), N_jd (liquids), flooding (gas)  
    **Max speed:** vortex / surface aeration, mechanical limits, splashing
    """)

# ── Quantitative response entry for Test 1 ──────────────────────────────
st.subheader("Record Test 1 Responses")
#  Sensitivity criterion default to 5%
st.caption("Sensitivity criterion: ≥ 5% relative change from center = sensitive.")

_t1_resp_name = st.text_input("Response metric name", value="Yield (%)",
                              key=_bk("t1_resp_name"))
col_t1a, col_t1b, col_t1c = st.columns(3)
with col_t1a:
    t1_resp_low = st.number_input(f"{_low_label}", value=0.0, format="%.4g",
                                  key=_bk("t1_resp_low"))
with col_t1b:
    t1_resp_ctr = st.number_input("Center (1× P/m)", value=0.0, format="%.4g",
                                  key=_bk("t1_resp_ctr"))
with col_t1c:
    t1_resp_high = st.number_input(f"{_high_label}", value=0.0, format="%.4g",
                                   key=_bk("t1_resp_high"))

# Sensitivity assessment
_SENSITIVITY_THRESHOLD = 5.0  # % relative change (Sarafinas, 2018)

def _assess_sensitivity(resp_values, ref_value, threshold_pct=_SENSITIVITY_THRESHOLD):
    """Return (max_pct_change, is_sensitive, detail_text)."""
    if ref_value == 0:
        # Fall back to absolute range if center is zero
        rng = max(resp_values) - min(resp_values)
        return (rng, rng > 0, f"Absolute range = {rng:.4g} (center = 0; cannot compute relative change)")
    pct_changes = [abs(v - ref_value) / abs(ref_value) * 100 for v in resp_values]
    max_pct = max(pct_changes)
    sensitive = max_pct >= threshold_pct
    return (max_pct, sensitive, pct_changes)

# Button-triggered assessment
if st.button("📊 Assess Test 1 Responses", key=_bk("t1_assess")):
    if t1_resp_ctr == 0.0 and t1_resp_low == 0.0 and t1_resp_high == 0.0:
        st.warning("Enter all three response values before assessing.")
    else:
        t1_max_pct, t1_sensitive, t1_pct_detail = _assess_sensitivity(
            [t1_resp_low, t1_resp_ctr, t1_resp_high], t1_resp_ctr
        )
        st.session_state[_bk("t1_assessed")] = {
            "max_pct": t1_max_pct, "sensitive": t1_sensitive,
            "pct_detail": t1_pct_detail,
            "resp": [t1_resp_low, t1_resp_ctr, t1_resp_high],
            "resp_name": _t1_resp_name,
            "labels": [_low_label, "Center (1× P/m)", _high_label],
        }

# Display results if assessed
if _bk("t1_assessed") not in st.session_state:
    st.info("Enter all three response values, then click **Assess Test 1 Responses**.")
    st.stop()

_t1a = st.session_state[_bk("t1_assessed")]
t1_max_pct = _t1a["max_pct"]
t1_sensitive = _t1a["sensitive"]
t1_pct_detail = _t1a["pct_detail"]

if isinstance(t1_pct_detail, list):
    _t1_res_df = pd.DataFrame({
        "Condition": _t1a["labels"],
        _t1a["resp_name"]: _t1a["resp"],
        "Δ from center (%)": [f"{p:.1f}%" for p in t1_pct_detail],
    })
else:
    _t1_res_df = pd.DataFrame({
        "Condition": _t1a["labels"],
        _t1a["resp_name"]: _t1a["resp"],
    })
st.dataframe(_t1_res_df, use_container_width=True, hide_index=True)
st.caption(f"Maximum relative change from center = **{t1_max_pct:.1f}%**  •  "
           f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**")

if not t1_sensitive:
    st.success(
        f"✅ **Mixing not critical.** Max variation {t1_max_pct:.1f}% "
        f"< {_SENSITIVITY_THRESHOLD:.0f}% threshold. Protocol complete."
    )
else:
    st.warning(
        f"⚠️ **Mixing matters!** {t1_max_pct:.1f}% change "
        f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Proceed to **Test 2**."
    )

# Allow user to override the automatic assessment
t1_override = st.checkbox("Override automatic assessment", key=_bk("t1_override"))
if t1_override:
    t1_result_manual = st.radio(
        "Manual assessment:",
        ["Sensitive – mixing matters", "Not sensitive – mixing not critical"],
        key=_bk("t1_manual"),
        horizontal=True,
    )
    t1_sensitive = t1_result_manual.startswith("Sensitive")

if not t1_sensitive:
    st.stop()

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 2 – TEST 2: Feed Rate / Time  (Meso vs Micro)
# ══════════════════════════════════════════════════════════════════════════
st.header("Test 2 - Feed Rate / Feed Time")
st.markdown("""
Vary **feed time only** (hold speed & location at centerpoint). Test a **9× flow-rate
range** (1/3× and 3× the centerpoint feed time).

Sensitive → mesomixing-controlled. Insensitive → micromixing-controlled.

""")

st.subheader("Suggested Feed Time Conditions")

t_feed_center = st.number_input(
    "Centerpoint feed time (min)",
    min_value=0.1, value=60.0, step=5.0,
    help="Fastest safe & practical plant-scale feed time, or your current lab feed time.",
    key=_bk("tfeed_ctr"),
)

t_feed_fast = t_feed_center / 3.0
t_feed_slow = t_feed_center * 3.0

feed_vol = st.number_input("Total feed volume (mL)", min_value=0.1, value=100.0, key=_bk("feed_vol"))

t2_conditions = pd.DataFrame([
    {
        "Condition": "Fast  (1/3× feed time)",
        "Feed time (min)": round(t_feed_fast, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_fast, 2),
        "Notes": "Highest mesomixing stress — more 'bad stuff' if meso-controlled",
    },
    {
        "Condition": "Center (1× feed time)",
        "Feed time (min)": round(t_feed_center, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_center, 2),
        "Notes": "Baseline condition",
    },
    {
        "Condition": "Slow  (3× feed time)",
        "Feed time (min)": round(t_feed_slow, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_slow, 2),
        "Notes": "Less mesomixing stress — better result if meso-controlled",
    },
])
st.dataframe(t2_conditions, width='content', hide_index=True)
st.caption(f"Flow-rate ratio (fast/slow) = **9×**  •  Hold impeller speed at N = {N_center_calc:.2f} rev/s ({N_center_calc*60:.0f} RPM)")

with st.expander("Understanding mesomixing"):
    st.markdown(r"""
    $$t_{meso} \propto \frac{Q_{feed}}{k \cdot \varepsilon_{loc}}$$

    Higher feed rate → longer plume persistence → more concentration gradients.
    If response is **insensitive** to feed rate → **micromixing-controlled**
    (plume disperses fast; molecular-scale mixing is rate-limiting).
    """)

# ── Quantitative response entry for Test 2 ──────────────────────────────
st.subheader("Record Test 2 Responses")
st.caption("Use the same response metric as Test 1.")

_t2_resp_name = st.text_input("Response metric name", value=_t1_resp_name,
                              key=_bk("t2_resp_name"))
col_t2a, col_t2b, col_t2c = st.columns(3)
with col_t2a:
    t2_resp_fast = st.number_input("Fast (1/3× feed time)", value=0.0,
                                   format="%.4g", key=_bk("t2_resp_fast"))
with col_t2b:
    t2_resp_ctr = st.number_input("Center (1× feed time)", value=0.0,
                                  format="%.4g", key=_bk("t2_resp_ctr"))
with col_t2c:
    t2_resp_slow = st.number_input("Slow (3× feed time)", value=0.0,
                                   format="%.4g", key=_bk("t2_resp_slow"))

# Button-triggered assessment
if st.button("📊 Assess Test 2 Responses", key=_bk("t2_assess")):
    if t2_resp_ctr == 0.0 and t2_resp_fast == 0.0 and t2_resp_slow == 0.0:
        st.warning("Enter all three response values before assessing.")
    else:
        t2_max_pct, t2_sensitive, t2_pct_detail = _assess_sensitivity(
            [t2_resp_fast, t2_resp_ctr, t2_resp_slow], t2_resp_ctr
        )
        st.session_state[_bk("t2_assessed")] = {
            "max_pct": t2_max_pct, "sensitive": t2_sensitive,
            "pct_detail": t2_pct_detail,
            "resp": [t2_resp_fast, t2_resp_ctr, t2_resp_slow],
            "resp_name": _t2_resp_name,
        }

if _bk("t2_assessed") not in st.session_state:
    st.info("Enter all three response values, then click **Assess Test 2 Responses**.")
    st.stop()

_t2a = st.session_state[_bk("t2_assessed")]
t2_max_pct = _t2a["max_pct"]
t2_sensitive = _t2a["sensitive"]
t2_pct_detail = _t2a["pct_detail"]

if isinstance(t2_pct_detail, list):
    _t2_res_df = pd.DataFrame({
        "Condition": ["Fast (1/3×)", "Center (1×)", "Slow (3×)"],
        _t2a["resp_name"]: _t2a["resp"],
        "Δ from center (%)": [f"{p:.1f}%" for p in t2_pct_detail],
    })
else:
    _t2_res_df = pd.DataFrame({
        "Condition": ["Fast (1/3×)", "Center (1×)", "Slow (3×)"],
        _t2a["resp_name"]: _t2a["resp"],
    })
st.dataframe(_t2_res_df, use_container_width=True, hide_index=True)
st.caption(f"Maximum relative change from center = **{t2_max_pct:.1f}%**  •  "
           f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**")

if not t2_sensitive:
    st.info(
        f"🔬 **Micromixing controls.** {t2_max_pct:.1f}% change "
        f"(< {_SENSITIVITY_THRESHOLD:.0f}%). Scale-up: hold **local ε** constant."
    )
    _micro_conclusion = True
else:
    st.warning(
        f"⚠️ **Mesomixing matters!** {t2_max_pct:.1f}% change "
        f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Proceed to **Test 3**."
    )
    _micro_conclusion = False

# Allow user to override
t2_override = st.checkbox("Override automatic assessment", key=_bk("t2_override"))
if t2_override:
    t2_result_manual = st.radio(
        "Manual assessment:",
        ["Sensitive – mesomixing matters", "Not sensitive – micromixing controls"],
        key=_bk("t2_manual"),
        horizontal=True,
    )
    t2_sensitive = t2_result_manual.startswith("Sensitive")
    _micro_conclusion = not t2_sensitive

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 3 – TEST 3: Feed Location  (Meso vs Macro)
# ══════════════════════════════════════════════════════════════════════════

# Defaults when Test 3 is not reached (micromixing controls)
t3_sensitive = False
t3_max_pct = 0.0
_macro_conclusion = False
_meso_conclusion = False

# Compute ε_avg at centerpoint (used in Test 3 and save section)
eps_avg = power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3)
eps_avg_kg = eps_avg / rho  # W/kg for micromixing calculations

if t2_sensitive:
    st.header("Test 3 - Feed Location")
    st.markdown("""
    Vary **feed location only** (hold speed & feed time at centerpoint). Move
    from surface → impeller zone. This changes **local ε** without affecting
    macromixing.

    | Feed location | ε_loc / ε_avg |
    |---|---|
    | Surface | ≈ 0.1 |
    | Sub-surface (mid-tank) | ≈ 1 |
    | Impeller zone | ≈ 2 – 5 |

    Sensitive → **mesomixing**. Insensitive → **macromixing**.
    """)

    t3_locs = pd.DataFrame([
        {
            "Feed Location": "Surface",
            "ε_loc / ε_avg": 0.1,
            "ε_loc (W/m³)": round(0.1 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(0.1 * eps_avg_kg, nu), 5),
        },
        {
            "Feed Location": "Sub-surface (mid-tank)",
            "ε_loc / ε_avg": 1.0,
            "ε_loc (W/m³)": round(1.0 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(1.0 * eps_avg_kg, nu), 5),
        },
        {
            "Feed Location": "Impeller zone",
            "ε_loc / ε_avg": 3.0,
            "ε_loc (W/m³)": round(3.0 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(3.0 * eps_avg_kg, nu), 5),
        },
    ])
    st.dataframe(t3_locs, use_container_width=True, hide_index=True)
    st.caption("Local ε estimated at centerpoint impeller speed. Actual values depend on impeller type and geometry.")

    with st.expander("Practical considerations for feed location"):
        st.markdown("""
        - Subsurface feeds → higher local ε but risk of **back-mixing** and plugging.
        - Feed velocity must exceed local bulk velocity (Jo et al., 1994).
        - *"Subsurface addition — prove it is not needed!"* — Ed Paul
        """)

    # ── Quantitative response entry for Test 3 ──────────────────────────
    st.subheader("Record Test 3 Responses")
    st.caption("Use the same response metric as Tests 1 & 2.")

    _t3_resp_name = st.text_input("Response metric name", value=_t2_resp_name,
                                  key=_bk("t3_resp_name"))
    col_t3a, col_t3b, col_t3c = st.columns(3)
    with col_t3a:
        t3_resp_surf = st.number_input("Surface feed", value=0.0, format="%.4g",
                                       key=_bk("t3_resp_surf"))
    with col_t3b:
        t3_resp_mid = st.number_input("Sub-surface (mid-tank)", value=0.0,
                                      format="%.4g", key=_bk("t3_resp_mid"))
    with col_t3c:
        t3_resp_imp = st.number_input("Impeller zone", value=0.0, format="%.4g",
                                      key=_bk("t3_resp_imp"))

    # Button-triggered assessment
    if st.button("📊 Assess Test 3 Responses", key=_bk("t3_assess")):
        if t3_resp_mid == 0.0 and t3_resp_surf == 0.0 and t3_resp_imp == 0.0:
            st.warning("Enter all three response values before assessing.")
        else:
            t3_max_pct, t3_sensitive, t3_pct_detail = _assess_sensitivity(
                [t3_resp_surf, t3_resp_mid, t3_resp_imp], t3_resp_mid
            )
            st.session_state[_bk("t3_assessed")] = {
                "max_pct": t3_max_pct, "sensitive": t3_sensitive,
                "pct_detail": t3_pct_detail,
                "resp": [t3_resp_surf, t3_resp_mid, t3_resp_imp],
                "resp_name": _t3_resp_name,
            }

    if _bk("t3_assessed") not in st.session_state:
        st.info("Enter all three response values, then click **Assess Test 3 Responses**.")
        st.stop()

    _t3a = st.session_state[_bk("t3_assessed")]
    t3_max_pct = _t3a["max_pct"]
    t3_sensitive = _t3a["sensitive"]
    t3_pct_detail = _t3a["pct_detail"]

    if isinstance(t3_pct_detail, list):
        _t3_res_df = pd.DataFrame({
            "Feed Location": ["Surface", "Sub-surface (mid)", "Impeller zone"],
            _t3a["resp_name"]: _t3a["resp"],
            "Δ from mid-tank (%)": [f"{p:.1f}%" for p in t3_pct_detail],
        })
    else:
        _t3_res_df = pd.DataFrame({
            "Feed Location": ["Surface", "Sub-surface (mid)", "Impeller zone"],
            _t3a["resp_name"]: _t3a["resp"],
        })
    st.dataframe(_t3_res_df, use_container_width=True, hide_index=True)
    st.caption(f"Maximum relative change from mid-tank = **{t3_max_pct:.1f}%**  •  "
               f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**")

    if not t3_sensitive:
        st.info(
            f"🌀 **Macromixing controls.** {t3_max_pct:.1f}% change "
            f"(< {_SENSITIVITY_THRESHOLD:.0f}%). Bulk blending is rate-limiting."
        )
        st.caption("Scale-up: maintain short blend times (hydrofoils, static mixers, multiple impellers).")
        _macro_conclusion = True
        _meso_conclusion = False
    else:
        st.success(
            f"📐 **Mesomixing controls.** {t3_max_pct:.1f}% change "
            f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Feed-plume dispersion is rate-limiting."
        )
        st.caption("Scale-up: constant impeller speed, extended feed time, or multiple feed points.")
        _macro_conclusion = False
        _meso_conclusion = True

    # Allow user to override
    t3_override = st.checkbox("Override automatic assessment", key=_bk("t3_override"))
    if t3_override:
        t3_result_manual = st.radio(
            "Manual assessment:",
            ["Sensitive – mesomixing controls", "Not sensitive – macromixing controls"],
            key=_bk("t3_manual"),
            horizontal=True,
        )
        t3_sensitive = t3_result_manual.startswith("Sensitive")
        _macro_conclusion = not t3_sensitive
        _meso_conclusion = t3_sensitive

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 4 – Confirmatory Experiments
# ══════════════════════════════════════════════════════════════════════════
st.header("4 · Confirmatory Experiments (Optional)")
st.caption("Optional tests to validate the conclusion from Tests 2–3.")

tab_fp, tab_visc = st.tabs(["A – Number of Feed Points", "B – Viscosity Change"])

with tab_fp:
    st.markdown("""
    ### A — Number of Feed Points
    Increase feed points (1 → 2–3) at same speed, feed time, and location.
    More feed points → lower local feed velocity at each point.

    | Outcome | Conclusion |
    |---|---|
    | Insensitive to # feed points | Micromixing controlled |
    | Changes with # feed points | Mesomixing controlled |
    """)

    fp_result = st.radio(
        "Was the process response sensitive to the number of feed points?",
        ["— Select —", "No – insensitive (supports micromixing control)",
         "Yes – changed (supports mesomixing control)"],
        key=_bk("fp_result"),
    )
    if fp_result.startswith("No"):
        st.info("✅ Confirms **micromixing** control.")
    elif fp_result.startswith("Yes"):
        st.info("✅ Confirms **mesomixing** control.")

with tab_visc:
    st.markdown(r"""
    ### B — Viscosity Change
    Change bulk viscosity (e.g. co-solvent ratio) while holding recipe constant.
    $t_E = 17.3 (\nu / \varepsilon)^{1/2}$

    | Outcome | Conclusion |
    |---|---|
    | Response **changes** with viscosity | Micromixing controlled |
    | Response **insensitive** | Not micromixing controlled |

    ⚠️ Large viscosity changes may shift the flow regime — keep Re turbulent.
    """)

    visc_result = st.radio(
        "Was the process response sensitive to viscosity change?",
        ["— Select —", "Yes – changed (supports micromixing control)",
         "No – insensitive (does not support micromixing control)"],
        key=_bk("visc_result"),
    )
    if visc_result.startswith("Yes"):
        st.info("✅ Supports **micromixing** control.")
    elif visc_result.startswith("No"):
        st.info("✅ Suggests micromixing is **not** the controlling mechanism.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 5 – Overall Conclusion & Decision Tree
# ══════════════════════════════════════════════════════════════════════════
st.header("5 · Protocol Summary & Conclusion")

# Build the conclusion from all recorded answers
conclusions = []
scaleup_notes = []

# Test 1 always passed if we reached here
conclusions.append(("Test 1 – Impeller Speed",
                    f"**Sensitive** ({t1_max_pct:.1f}% change) → Mixing matters", "⚠️"))

# Test 2
if not t2_sensitive:
    conclusions.append(("Test 2 – Feed Rate",
                        f"**Insensitive** ({t2_max_pct:.1f}% change) → Micromixing-controlled", "🔬"))
    scaleup_notes.append("Hold **local ε** constant at the feed point on scale-up.")
else:
    conclusions.append(("Test 2 – Feed Rate",
                        f"**Sensitive** ({t2_max_pct:.1f}% change) → Not purely micromixing-controlled", "⚠️"))

# Test 3
if t2_sensitive:  # only reached Test 3 if Test 2 was sensitive
    if t3_sensitive:
        conclusions.append(("Test 3 – Feed Location",
                            f"**Sensitive** ({t3_max_pct:.1f}% change) → Mesomixing-controlled", "📐"))
        scaleup_notes.append("Control feed plume dispersion: constant impeller speed, extended feed time, or multiple feed points.")
    else:
        conclusions.append(("Test 3 – Feed Location",
                            f"**Insensitive** ({t3_max_pct:.1f}% change) → Macromixing-controlled", "🌀"))
        scaleup_notes.append("Maintain short blend times on scale-up (high-efficiency impellers, static mixers).")

# Determine the dominant mixing limitation
if not t2_sensitive:
    dominant = "Micromixing"
    dominant_icon = "🔬"
    dominant_color = "blue"
elif t3_sensitive:
    dominant = "Mesomixing"
    dominant_icon = "📐"
    dominant_color = "orange"
else:
    dominant = "Macromixing"
    dominant_icon = "🌀"
    dominant_color = "red"

# Summary table
st.subheader("Test Results")
for test_name, result, icon in conclusions:
    st.markdown(f"{icon} **{test_name}:** {result}")

st.subheader("Dominant Mixing Limitation")

if dominant == "Micromixing":
    st.success(f"""
    {dominant_icon} **{dominant}** — Molecular-scale engulfment is rate-limiting.

    **Scale-up:** Maintain constant **local ε** at the feed point.
    """)
elif dominant == "Mesomixing":
    st.warning(f"""
    {dominant_icon} **{dominant}** — Feed-plume disintegration is rate-limiting.

    **Scale-up:** Constant impeller speed, extended feed time, or multiple feed points.
    """)
elif dominant == "Macromixing":
    st.error(f"""
    {dominant_icon} **{dominant}** — Bulk blending / circulation is rate-limiting.

    **Scale-up:** Reduce blend time (hydrofoil impellers, multiple impellers, static mixers).
    """)
else:
    st.info("Complete Tests 2 and 3 to determine the controlling mixing scale.")

# Scale-up notes
if scaleup_notes:
    st.subheader("Key Scale-Up Actions")
    for note in scaleup_notes:
        st.markdown(f"- {note}")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 6 – Save Protocol Record
# ══════════════════════════════════════════════════════════════════════════
st.header("6 · Save Protocol Record")

if st.button("📌 Save Bourne Protocol result to Recorded Results", key=_bk("save")):
    result_row = {
        "reactor": reactor_name if not reactors.empty else "Manual entry",
        "fluid": fluid_name,
        "reaction": "Bourne Protocol",
        "Re": round(reynolds_number(N_center_calc, D_imp, rho, mu), 0),
        "P/V (W/L)": round(power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3) / 1000, 4),
        "Tip speed (m/s)": round(tip_speed(N_center_calc, D_imp), 3),
        "Blend time (s)": round(blend_time_turbulent(Nq_val, V_m3, D_imp, N_center_calc), 2),
        "Micromix t_E (s)": round(micromixing_time_engulfment(eps_avg / rho, nu), 5),
        "Kolmogorov η (µm)": round(kolmogorov_length(nu, eps_avg / rho) * 1e6, 1),
        "t_rxn (s)": "N/A (protocol)",
        "Da_macro": "N/A",
        "Da_micro": "N/A",
        "Assessment": f"Bourne Protocol → {dominant}-controlled",
    }
    if "recorded_results" not in st.session_state:
        st.session_state.recorded_results = pd.DataFrame()
    st.session_state.recorded_results = pd.concat(
        [st.session_state.recorded_results, pd.DataFrame([result_row])],
        ignore_index=True,
    )
    results_csv = DATA_DIR / "recorded_results.csv"
    st.session_state.recorded_results.to_csv(results_csv, index=False)
    st.success("Protocol result saved to **Recorded Results**.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 7 – Export Report
# ══════════════════════════════════════════════════════════════════════════
st.header("7 · Export Report")

if st.button("📥 Export PDF Report", type="primary", key=_bk("export_pdf")):
    with st.spinner("Generating PDF…"):
        try:
            from utils.report_builder import build_bourne_protocol_pdf, report_filename

            _bp_snap = {
                "reactor": reactor_name if not reactors.empty else "Manual entry",
                "fluid": fluid_name,
                "V_L": V_L,
                "dominant": dominant,
                "conclusions": conclusions,
                "scaleup_notes": scaleup_notes,
                "t1_conditions": t1_rows,
                "t1_responses": st.session_state.get(_bk("t1_assessed")),
                "t2_conditions": {
                    "N_RPM": round(N_center_calc * 60, 1),
                    "feed_vol_mL": feed_vol,
                    "feed_location": "Held constant (centerpoint)",
                    "rows": [
                        {"Condition": "Fast (1/3x feed time)", "Feed time (min)": round(t_feed_fast, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_fast, 2)},
                        {"Condition": "Center (1x feed time)", "Feed time (min)": round(t_feed_center, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_center, 2)},
                        {"Condition": "Slow (3x feed time)", "Feed time (min)": round(t_feed_slow, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_slow, 2)},
                    ],
                },
                "t2_responses": st.session_state.get(_bk("t2_assessed")),
                "t3_conditions": {
                    "N_RPM": round(N_center_calc * 60, 1),
                    "feed_time_min": round(t_feed_center, 2),
                    "eps_avg_W_m3": round(eps_avg, 2),
                    "rows": [
                        {"Feed Location": "Surface", "eps_loc/eps_avg": 0.1, "eps_loc (W/m3)": round(0.1 * eps_avg, 2)},
                        {"Feed Location": "Sub-surface (mid-tank)", "eps_loc/eps_avg": 1.0, "eps_loc (W/m3)": round(1.0 * eps_avg, 2)},
                        {"Feed Location": "Impeller zone", "eps_loc/eps_avg": 3.0, "eps_loc (W/m3)": round(3.0 * eps_avg, 2)},
                    ],
                },
                "t3_responses": st.session_state.get(_bk("t3_assessed")),
                "centerpoint_metrics": {
                    "N (RPM)": round(N_center_calc * 60, 1),
                    "P/m (W/kg)": round(PV_center_wkg, 4),
                    "Re": round(reynolds_number(N_center_calc, D_imp, rho, mu), 0),
                    "Tip speed (m/s)": round(tip_speed(N_center_calc, D_imp), 3),
                    "Blend time (s)": round(blend_time_turbulent(Nq_val, V_m3, D_imp, N_center_calc), 2),
                    "Micromix t_E (s)": round(micromixing_time_engulfment(eps_avg_kg, nu), 5),
                    "Kolmogorov eta (um)": round(kolmogorov_length(nu, eps_avg_kg) * 1e6, 1),
                },
            }
            _pdf_bytes = build_bourne_protocol_pdf(_bp_snap)
            st.session_state["_p6_pdf_bytes"] = _pdf_bytes
            st.session_state["_p6_pdf_name"] = report_filename(
                "Bourne_Protocol", reactor_name if not reactors.empty else ""
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_p6_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_p6_pdf_bytes"],
        file_name=st.session_state["_p6_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 7 – References
# ══════════════════════════════════════════════════════════════════════════
st.header("References")
st.markdown("""
1. Bourne, J.R. (2003). *Mixing and the Selectivity of Chemical Reactions.*
   Org. Process Res. Dev., 7(4), 471–508.
2. Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). *Handbook of
   Industrial Mixing.* Wiley-Interscience.
3. Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical
   Reactions.* Wiley.
4. Jo, M.C., Penney, W.R. & Fasano, J.B. (1994). Back-mixing into reactor
   feed pipes caused by turbulence in an agitated vessel. AIChE Symp. Ser.,
   90(299).
""")
