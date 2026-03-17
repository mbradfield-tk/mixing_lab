"""
Page 8 – Bourne Protocol (Sarafinas Modification)
====================================================
Guides the user step-by-step through the Bourne mixing-sensitivity screening
protocol as described by Bourne (2003) and modified/extended by
Aaron Sarafinas (2018).

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


# ── Helper: load DB ──────────────────────────────────────────────────────
def _load(key, fn):
    if key not in st.session_state:
        p = DATA_DIR / fn
        st.session_state[key] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    return st.session_state[key]


reactors = _load("reactor_db", "reactors.csv")
custom_fluids = _load("fluid_db", "fluids.csv")

# Build combined fluid list: built-in solvents + custom fluids
_solvent_names = sorted(SOLVENT_DB.keys())
_custom_names = custom_fluids["fluid_name"].tolist() if not custom_fluids.empty else []
_all_fluid_names = _solvent_names + _custom_names

# ══════════════════════════════════════════════════════════════════════════
st.title("🧐 Bourne Protocol – Mixing Sensitivity Screening")
st.markdown("""
The **Bourne Protocol** (Bourne 2003, as extended by Sarafinas 2018) is an
efficient experimental method to determine **which scale of mixing** —
macro, meso, or micro — controls your process outcome.

The protocol is applicable to **any** semi-batch or fed-batch process in a
stirred tank where competitive rate processes (reactions, crystallisation,
precipitation, etc.) can impact the result.

> *"Always assume there is a mixing problem until proven otherwise."*
> — E. L. Paul (2003)
""")

_BP_IMG = pathlib.Path(__file__).resolve().parent.parent / "images" / "general" / "bourne_protocol_decision_tree.png"
with st.expander("📋 Protocol overview – Decision Tree", expanded=False):
    if _BP_IMG.exists():
        st.image(str(_BP_IMG), caption="Bourne Protocol – Decision Tree (Bourne 2003 / Sarafinas 2018)")
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
st.header("0 · Define Your System")

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
        reactor_name = st.selectbox("Reactor", _reactor_list, index=_sel_idx(_reactor_list, "_sel_bp_reactor"), key=_bk("reactor"))
        st.session_state["_sel_bp_reactor"] = reactor_name
        r = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
        D_tank = float(r["D_tank_m"])
        H = float(r["H_m"])
        D_imp = float(r["D_imp_m"])
        Np_val = float(r["Np"])
        Nq_val = float(r["Nq"])
        V_L = float(r["V_L"])
        # RPM bounds from reactor database
        N_rpm_min = float(r["N_rpm_min"]) if pd.notna(r.get("N_rpm_min")) else None
        N_rpm_max = float(r["N_rpm_max"]) if pd.notna(r.get("N_rpm_max")) else None
        # Use average of min/max RPM as the centerpoint speed
        if N_rpm_min is not None and N_rpm_max is not None:
            N_center = (N_rpm_min + N_rpm_max) / 2.0 / 60.0  # rev/s
        else:
            N_center = float(r["N_rps"])
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
        V_L = np.pi / 4 * D_tank**2 * H * 1000
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
        _cust = custom_fluids[custom_fluids["fluid_name"] == fluid_name].iloc[0]
        rho = float(_cust["rho_kg_m3"])
        mu = float(_cust["mu_Pa_s"])

nu = mu / rho
V_m3 = V_L / 1000.0  # m³

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 1 – TEST 1: Impeller Speed  (Does mixing matter?)
# ══════════════════════════════════════════════════════════════════════════
st.header("1 · Test 1 — Impeller Speed (Does Mixing Matter?)")
st.markdown("""
**What to vary:** Impeller speed only (hold feed rate and feed location constant; if any).

**Goal:** Achieve approximately a **100× change in P/V** across three speeds (low /
center / high).  If the process response changes, mixing matters — proceed to
Test 2.  If there is no change, mixing is **not critical** over this range.

**Sarafinas guidelines for choosing speeds:**
- If targeting existing plant equipment → use the **plant condition** as the
  centerpoint.
- Default centerpoint: vessel-average P/V ≈ **0.2 W/kg** (≈ 1 HP / 1000 US gal).
- Set upper and lower P/V at **10×** and **0.1×** the centerpoint.
- For turbulent flow, this corresponds to speed ratios of ≈ **2.15×** and
  **1/2.15×** the centerpoint RPM (since P ∝ N³).
""")

st.subheader("Suggested Experimental Conditions")

use_custom_center = st.checkbox("Set a custom centerpoint P/V", value=False, key=_bk("custom_pv"))

if use_custom_center:
    PV_center_wkg = st.number_input("Centerpoint P/V (W/kg)", value=0.2, format="%.4g", key=_bk("pv_ctr"))
else:
    # Sarafinas recommends 0.2 W/kg as the default centerpoint.
    # Check whether the reactor can achieve it within its RPM range.
    _SARAFINAS_DEFAULT_PV = 0.2  # W/kg
    _N_for_default = (_SARAFINAS_DEFAULT_PV * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
    _N_for_default_rpm = _N_for_default * 60

    _can_reach_default = True
    if N_rps_min is not None and _N_for_default < N_rps_min:
        _can_reach_default = False
    if N_rps_max is not None and _N_for_default > N_rps_max:
        _can_reach_default = False

    if _can_reach_default:
        PV_center_wkg = _SARAFINAS_DEFAULT_PV
        st.write(f"Centerpoint from Sarafinas default: **P/V = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {_N_for_default:.2f} rev/s ({_N_for_default_rpm:.0f} RPM)")
        if N_rpm_min is not None and N_rpm_max is not None:
            st.caption(f"Sarafinas default 0.2 W/kg is within reactor RPM range "
                       f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM).")
    else:
        # Fall back to average reactor speed
        _P_ctr = impeller_power(Np_val, rho, N_center, D_imp)
        PV_center_wkg = power_per_volume(_P_ctr, V_m3) / rho
        _ctr_rpm = N_center * 60
        st.write(f"Centerpoint from reactor average speed: **P/V = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {N_center:.2f} rev/s ({_ctr_rpm:.0f} RPM)")
        st.caption(f"Sarafinas default 0.2 W/kg requires {_N_for_default_rpm:.0f} RPM, "
                   f"which is outside the reactor range "
                   f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM). "
                   f"Using average speed instead.")

# Compute the three test speeds
# P/V ∝ N³  →  N ∝ (P/V)^(1/3)
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
        "N (rev/s)": round(N, 3),
        "N (RPM)": round(N * 60, 1),
        "P/V (W/kg)": round(pv_wkg, 4),
        "P/V rel. to center": f"{pv_rel:.2f}×",
        "P/V (W/L)": round(eps / 1000, 4),
        "Re": round(Re, 0),
        "Tip speed (m/s)": round(u, 3),
        "Blend time (s)": round(t_blend, 2),
        "t_E micro (s)": round(t_micro, 5),
        "η Kolmogorov (µm)": round(eta * 1e6, 1),
    }

_low_label = "Low  (min RPM)" if _low_clamped else "Low  (0.1× P/V)"
_high_label = "High (max RPM)" if _high_clamped else "High (10× P/V)"

t1_rows = [
    _hydro_row(_low_label, N_low),
    _hydro_row("Center (1× P/V)", N_center_calc),
    _hydro_row(_high_label, N_high),
]
t1_df = pd.DataFrame(t1_rows)
st.dataframe(t1_df, use_container_width=True, hide_index=True)

pv_ratio = t1_rows[2]["P/V (W/kg)"] / t1_rows[0]["P/V (W/kg)"] if t1_rows[0]["P/V (W/kg)"] > 0 else 0
st.caption(f"P/V ratio (high/low) = **{pv_ratio:.1f}×**  •  Speed ratio (high/low) = **{N_high/N_low:.2f}×**")

if _low_clamped or _high_clamped:
    _msgs = []
    if _low_clamped:
        _actual_low_rel = t1_rows[0]["P/V (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"Low speed clamped to reactor minimum ({N_rpm_min:.1f} RPM) → "
                     f"actual P/V = **{_actual_low_rel:.3f}×** centerpoint "
                     f"(target was 0.1×)")
    if _high_clamped:
        _actual_high_rel = t1_rows[2]["P/V (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"High speed clamped to reactor maximum ({N_rpm_max:.1f} RPM) → "
                     f"actual P/V = **{_actual_high_rel:.3f}×** centerpoint "
                     f"(target was 10×)")
    st.warning("⚠️ **Speed limits applied:**\n\n" + "\n\n".join(_msgs) +
              "\n\nThe achievable P/V range is narrower than the ideal 100× span. "
              "Consider whether this is sufficient to draw conclusions.")

# Additional practical limits
with st.expander("Practical limits to consider"):
    st.markdown("""
    **Minimum speed constraints:**
    - Solid–liquid systems → just-suspended speed (N_js)
    - Liquid–liquid systems → just-dispersed speed (N_jd)
    - Gas–liquid systems → flooding condition

    **Maximum speed constraints:**
    - Surface aeration / vortex formation
    - Mechanical limits of agitator / seal
    - Splashing at small scale
    """)

# ── Quantitative response entry for Test 1 ──────────────────────────────
st.subheader("Record Test 1 Responses")
st.markdown(
    "Enter the measured process response (e.g. yield %, selectivity %, "
    "impurity level, particle size) at each condition. The app will calculate "
    "the magnitude of change and determine sensitivity using the "
    "**Sarafinas criterion** (≥ 5% relative change = sensitive)."
)

_t1_resp_name = st.text_input("Response metric name", value="Yield (%)",
                              key=_bk("t1_resp_name"))
col_t1a, col_t1b, col_t1c = st.columns(3)
with col_t1a:
    t1_resp_low = st.number_input(f"{_low_label}", value=0.0, format="%.4g",
                                  key=_bk("t1_resp_low"))
with col_t1b:
    t1_resp_ctr = st.number_input("Center (1× P/V)", value=0.0, format="%.4g",
                                  key=_bk("t1_resp_ctr"))
with col_t1c:
    t1_resp_high = st.number_input(f"{_high_label}", value=0.0, format="%.4g",
                                   key=_bk("t1_resp_high"))

# Sensitivity assessment
_SARAFINAS_THRESHOLD = 5.0  # % relative change

def _assess_sensitivity(resp_values, ref_value, threshold_pct=_SARAFINAS_THRESHOLD):
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
            "labels": [_low_label, "Center (1× P/V)", _high_label],
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
           f"Sarafinas threshold = **{_SARAFINAS_THRESHOLD:.0f}%**")

if not t1_sensitive:
    st.success(
        f"✅ **Mixing does not appear to be critical.** Maximum response variation "
        f"({t1_max_pct:.1f}%) is below the {_SARAFINAS_THRESHOLD:.0f}% threshold. "
        f"Protocol complete."
    )
    st.info("Consider re-running the protocol if process conditions change significantly "
            "(concentration, phases, viscosity, reagents).")
else:
    st.warning(
        f"⚠️ **Mixing matters!** Response varied by {t1_max_pct:.1f}% "
        f"(≥ {_SARAFINAS_THRESHOLD:.0f}% threshold). "
        f"Proceed to **Test 2** to identify which mixing scale is responsible."
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
st.header("2 · Test 2 — Feed Rate / Feed Time")
st.markdown("""
**What to vary:** Feed time (or equivalently, feed rate) **only** — keep impeller
speed and feed location constant (use the centerpoint speed from Test 1).

**Goal:** Determine whether the process is in the **mesomixing-controlled**
regime (feed-rate sensitive) or **micromixing-controlled** regime
(feed-rate insensitive).

**Sarafinas guidelines:**
- Use the **fastest safe and practical plant-scale feed rate** as the midpoint.
- Test a **9× range** on volumetric flow rate: feed times of **1/3×** and **3×** the
  centerpoint feed time.
- Safety considerations (exotherms, gas evolution) must bound the range.

> *Changing the feed rate affects ONLY mesomixing.* — Sarafinas (2018)
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
    The **mesomixing time** characterises the disintegration of the fresh feed
    plume in the turbulent environment around the feed point.  Two expressions
    exist depending on whether the impeller-driven flow or the feed jet
    dominates locally:

    $$t_{meso} \propto \frac{Q_{feed}}{k \cdot \varepsilon_{loc}}$$

    Increasing the feed rate **increases** the mesomixing time (the feed plume
    persists longer), pushing the system toward higher local concentration
    gradients and more "bad stuff."

    If the process response is **insensitive** to feed rate, you are in the
    **micromixing-controlled** plateau — the feed plume disperses fast enough
    that the molecular-scale mixing step is rate-limiting.
    """)

# ── Quantitative response entry for Test 2 ──────────────────────────────
st.subheader("Record Test 2 Responses")
st.markdown(
    "Enter the measured process response at each feed-time condition. "
    "Use the **same metric** as Test 1."
)

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
           f"Sarafinas threshold = **{_SARAFINAS_THRESHOLD:.0f}%**")

if not t2_sensitive:
    st.info(
        f"🔬 **Micromixing appears to control the process.** Response varied by only "
        f"{t2_max_pct:.1f}% (< {_SARAFINAS_THRESHOLD:.0f}% threshold). "
        f"The feed plume disintegrates quickly; the molecular-scale engulfment "
        f"step is rate-limiting."
    )
    st.markdown("""
    **Scale-up recommendation (Bourne 2003):** Hold **local energy dissipation**
    (ε_loc) constant where the mixing and reaction occur.

    You may optionally run the **confirmatory experiments** below to strengthen
    this conclusion.
    """)
    _micro_conclusion = True
else:
    st.warning(
        f"⚠️ **Mesomixing (feed-plume dispersion) matters!** Response varied by "
        f"{t2_max_pct:.1f}% (≥ {_SARAFINAS_THRESHOLD:.0f}% threshold). "
        f"Proceed to **Test 3** to distinguish mesomixing from macromixing."
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
    st.header("3 · Test 3 — Feed Location")
    st.markdown("""
    **What to vary:** Feed location **only** — keep impeller speed and feed time at
    the centerpoint values from Tests 1 & 2.

    **Goal:** Move the feed point from a **low-intensity** region (liquid surface) to
    a **high-intensity** region (near the impeller), or vice versa.  This changes
    the **local** ε (and therefore both micro- and mesomixing times) **without
    affecting macromixing**.

    | Feed location | ε_loc / ε_avg (rule of thumb) |
    |---------------|-------------------------------|
    | Surface feed | ≈ 0.1 |
    | Sub-surface (mid-tank) | ≈ 1 |
    | Impeller zone feed | ≈ 2 – 5 |

    > *Changing the feed location affects both micromixing and mesomixing,
    > but not macromixing.* — Sarafinas (2018)
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
        - **Subsurface feeds** deliver the reagent into a higher-intensity mixing
          zone but carry a risk of **feed-pipe back-mixing** and plugging.
        - Design feed velocity to exceed the local bulk fluid velocity to prevent
          back-mixing (Jo et al., 1994; Vicum & Baldyga, 2004).
        - Consider the **plant configuration**: can the plant accommodate subsurface
          feed points?
        - Ed Paul's 3rd rule of mixing: *"Subsurface addition for reaction and
          crystallisation — prove that it is not needed!"*
        """)

    # ── Quantitative response entry for Test 3 ──────────────────────────
    st.subheader("Record Test 3 Responses")
    st.markdown(
        "Enter the measured process response at each feed-location condition. "
        "Use the **same metric** as Tests 1 & 2."
    )

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
               f"Sarafinas threshold = **{_SARAFINAS_THRESHOLD:.0f}%**")

    if not t3_sensitive:
        st.info(
            f"🌀 **Macromixing appears to control the process.** Response varied by only "
            f"{t3_max_pct:.1f}% (< {_SARAFINAS_THRESHOLD:.0f}% threshold). "
            f"The bulk blending / circulation time is the rate-limiting mixing step."
        )
        st.markdown("""
        **Scale-up recommendation (Bourne 2003):**
        - Maintain short blend times on scale-up.
        - Consider **high-efficiency hydrofoil impellers** (A310, A320),
          **static mixers**, or multiple impellers.
        - Note: *"The feed time for a semi-batch reactor is usually so long that
          macromixing is not controlling."* — Bourne (2003).  Consider re-examining
          Test 2 results.
        """)
        _macro_conclusion = True
        _meso_conclusion = False
    else:
        st.success(
            f"📐 **Mesomixing controls the process.** Response varied by "
            f"{t3_max_pct:.1f}% (≥ {_SARAFINAS_THRESHOLD:.0f}% threshold). "
            f"The disintegration of the feed plume in the local mixing environment "
            f"is the rate-limiting mixing step."
        )
        st.markdown("""
        **Scale-up recommendation (Bourne 2003 / Sarafinas 2018):**
        - **Keep impeller speed constant** on scale-up (expensive but effective).
        - Or **extend the feed time** proportionally (may be conservative).
        - Consider **multiple feed points** to reduce local feed intensity.
        - Design the feed system to control local exit velocity and plume
          dispersion.
        """)
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
st.markdown("""
Bourne proposed two optional confirmatory tests to validate the conclusion
from Tests 2 and 3.  These are most useful when the Test 2 result is
ambiguous or when additional confidence is needed before committing to
a scale-up strategy.
""")

tab_fp, tab_visc = st.tabs(["A – Number of Feed Points", "B – Viscosity Change"])

with tab_fp:
    st.markdown("""
    ### Confirmatory Test A — Number of Feed Points

    **Procedure:** At the same impeller speed, feed time, and feed location,
    increase the number of feed points (e.g., from 1 to 2 or 3).

    - **Same total feed time & tube diameter** → each feed point delivers
      1/N_fp of the total flow → lower local feed velocity → lower mesomixing
      time at each feed point.
    - Compare with the equivalent single-feed-point experiment at a longer
      feed time that gives the same local feed velocity.

    **Interpretation:**
    | Outcome | Conclusion |
    |---------|------------|
    | Process response **insensitive** to # feed points | Micromixing controlled |
    | Process response **changes** with # feed points | Mesomixing controlled |
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
    st.markdown("""
    ### Confirmatory Test B — Viscosity Change

    **Procedure:** Change the bulk-fluid viscosity (e.g., co-solvent dilution
    or different solvent ratio) while keeping the recipe otherwise unchanged.

    Viscosity directly affects the micromixing (engulfment) time:

    $$t_E = 17.3 \\left( \\frac{\\nu}{\\varepsilon} \\right)^{1/2}$$

    **Interpretation:**
    | Outcome | Conclusion |
    |---------|------------|
    | Process response **changes** with viscosity | Micromixing controlled |
    | Process response **insensitive** to viscosity | Not micromixing controlled |

    > ⚠️ *Caution:* Large viscosity changes can shift the flow regime from
    > turbulent to transitional → changing ε distribution and bulk blending.
    > Only useful when the regime remains turbulent.
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
    {dominant_icon} **{dominant}** — The molecular-scale engulfment step is rate-limiting.

    The feed plume disperses quickly, but concentration homogenisation at the
    Kolmogorov / Batchelor scale is slow relative to the reaction time.

    **Scale-up strategy:**
    - Maintain constant **local energy dissipation** (ε_loc) at the feed point.
    - ε_loc can be estimated from P/V × (ε_loc / ε_avg) for the feed location.
    - Consider impeller type and feed-point proximity to the impeller.
    """)
elif dominant == "Mesomixing":
    st.warning(f"""
    {dominant_icon} **{dominant}** — Feed-plume disintegration is rate-limiting.

    The fresh feed stream does not break up fast enough before competitive
    rate processes act on the locally high concentrations.

    **Scale-up strategy:**
    - **Keep impeller speed constant** on scale-up (costly in power).
    - Or **extend the feed time** to reduce local feed rate.
    - Use **multiple feed points** to reduce local feed velocity at each point.
    - Design feed nozzle diameter and velocity to optimise plume break-up.
    """)
elif dominant == "Macromixing":
    st.error(f"""
    {dominant_icon} **{dominant}** — Bulk blending / circulation is rate-limiting.

    The vessel contents are not homogenised quickly enough between feed
    additions, leading to large-scale concentration or temperature gradients.

    **Scale-up strategy:**
    - Focus on **blend time reduction**: high-efficiency impellers (hydrofoils),
      multiple impellers, or static mixers.
    - Avoid increasing vessel H/T ratio without compensating with additional
      impellers.
    - Consider continuous-flow alternatives with in-line mixing.
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
# SECTION 7 – References
# ══════════════════════════════════════════════════════════════════════════
st.header("References")
st.markdown("""
1. Bourne, J.R. (2003). *Mixing and the Selectivity of Chemical Reactions.*
   Org. Process Res. Dev., 7(4), 471–508.
2. Sarafinas, A. (2018). *Test Process Mixing Sensitivities Using the Bourne
   Protocol.* Scientific Update Webinar, 13 November 2018.
3. Sarafinas, A. & Teich, C.I. (2016). Chapter 13 in *Advances in Industrial
   Mixing* (Kresta et al., Wiley).
4. Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). *Handbook of
   Industrial Mixing.* Wiley-Interscience.
5. Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical
   Reactions.* Wiley.
6. Jo, M.C., Penney, W.R. & Fasano, J.B. (1994). Back-mixing into reactor
   feed pipes caused by turbulence in an agitated vessel. AIChE Symp. Ser.,
   90(299).
""")
