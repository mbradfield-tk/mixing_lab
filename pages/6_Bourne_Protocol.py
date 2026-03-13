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
fluids = _load("fluid_db", "fluids.csv")

# ══════════════════════════════════════════════════════════════════════════
st.title("🧫 Bourne Protocol – Mixing Sensitivity Screening")
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

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 0 – Reactor & Fluid Selection
# ══════════════════════════════════════════════════════════════════════════
st.header("0 · Define Your System")

# Persist selections across page navigations
_reactor_list = reactors["reactor_name"].tolist() if not reactors.empty else []
_fluid_list = fluids["fluid_name"].tolist() if not fluids.empty else []

def _sel_idx(lst, key, default=0):
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col_r, col_f = st.columns(2)

with col_r:
    if not reactors.empty:
        reactor_name = st.selectbox("Reactor", _reactor_list, index=_sel_idx(_reactor_list, "_sel_bp_reactor"), key="bp_reactor")
        st.session_state["_sel_bp_reactor"] = reactor_name
        r = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
        D_tank = float(r["D_tank_m"])
        H = float(r["H_m"])
        D_imp = float(r["D_imp_m"])
        N_center = float(r["N_rps"])
        Np_val = float(r["Np"])
        Nq_val = float(r["Nq"])
        V_L = float(r["V_L"])
    else:
        st.info("No reactors in database — enter manually below.")
        D_tank = st.number_input("Tank diameter (m)", value=0.10, format="%.4f")
        H = st.number_input("Liquid height (m)", value=0.13, format="%.4f")
        D_imp = st.number_input("Impeller diameter (m)", value=0.05, format="%.4f")
        N_center = st.number_input("Centerpoint speed (rev/s)", value=5.0, format="%.2f")
        Np_val = st.number_input("Power number Np", value=1.27, format="%.2f")
        Nq_val = st.number_input("Pumping number Nq", value=0.79, format="%.2f")
        V_L = np.pi / 4 * D_tank**2 * H * 1000

with col_f:
    if not fluids.empty:
        fluid_name = st.selectbox("Fluid system", _fluid_list, index=_sel_idx(_fluid_list, "_sel_bp_fluid"), key="bp_fluid")
        st.session_state["_sel_bp_fluid"] = fluid_name
        fl = fluids[fluids["fluid_name"] == fluid_name].iloc[0]
        rho = float(fl["rho_kg_m3"])
        mu = float(fl["mu_Pa_s"])
    else:
        st.info("No fluids in database — enter manually.")
        rho = st.number_input("Density ρ (kg/m³)", value=997.0)
        mu = st.number_input("Viscosity μ (Pa·s)", value=0.00089, format="%.6f")

nu = mu / rho
V_m3 = V_L / 1000.0  # m³

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 1 – TEST 1: Impeller Speed  (Does mixing matter?)
# ══════════════════════════════════════════════════════════════════════════
st.header("1 · Test 1 — Impeller Speed (Does Mixing Matter?)")
st.markdown("""
**What to vary:** Impeller speed only (hold feed rate and feed location constant).

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

use_custom_center = st.checkbox("Set a custom centerpoint P/V", value=False, key="bp_custom_pv")

if use_custom_center:
    PV_center_wkg = st.number_input("Centerpoint P/V (W/kg)", value=0.2, format="%.4g", key="bp_pv_ctr")
else:
    # Compute from database reactor at its listed speed
    _P_ctr = impeller_power(Np_val, rho, N_center, D_imp)
    PV_center_wkg = power_per_volume(_P_ctr, V_m3) / rho
    st.write(f"Centerpoint from reactor database: **P/V = {PV_center_wkg:.4g} W/kg** at N = {N_center:.2f} rev/s")

# Compute the three test speeds
# P/V ∝ N³  →  N ∝ (P/V)^(1/3)
PV_high = PV_center_wkg * 10
PV_low = PV_center_wkg * 0.1

N_center_calc = (PV_center_wkg * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
N_high = N_center_calc * 10**(1/3)   # ≈ 2.154×
N_low = N_center_calc * 0.1**(1/3)   # ≈ 0.464×

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
    return {
        "Condition": label,
        "N (rev/s)": round(N, 3),
        "N (RPM)": round(N * 60, 1),
        "P/V (W/kg)": round(eps / rho, 4),
        "P/V (W/L)": round(eps / 1000, 4),
        "Re": round(Re, 0),
        "Tip speed (m/s)": round(u, 3),
        "Blend time (s)": round(t_blend, 2),
        "t_E micro (s)": round(t_micro, 5),
        "η Kolmogorov (µm)": round(eta * 1e6, 1),
    }


t1_rows = [
    _hydro_row("Low  (0.1× P/V)", N_low),
    _hydro_row("Center (1× P/V)", N_center_calc),
    _hydro_row("High (10× P/V)", N_high),
]
t1_df = pd.DataFrame(t1_rows)
st.dataframe(t1_df, width="stretch", hide_index=True)

pv_ratio = t1_rows[2]["P/V (W/kg)"] / t1_rows[0]["P/V (W/kg)"] if t1_rows[0]["P/V (W/kg)"] > 0 else 0
st.caption(f"P/V ratio (high/low) = **{pv_ratio:.0f}×**  •  Speed ratio (high/low) = **{N_high/N_low:.1f}×**")

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

# User records outcome
st.subheader("Record Test 1 Outcome")
t1_result = st.radio(
    "Did the process response change when impeller speed was varied?",
    ["— Select —", "Yes – process response changed", "No – no significant change"],
    key="bp_t1_result",
    horizontal=True,
)

if t1_result.startswith("No"):
    st.success(
        "✅ **Mixing does not appear to be critical** for this process over the "
        "tested range of P/V. Protocol complete."
    )
    st.info("Consider re-running the protocol if process conditions change significantly "
            "(concentration, phases, viscosity, reagents).")
    st.stop()

if t1_result == "— Select —":
    st.info("Run the experiments above, then record the outcome to continue the protocol.")
    st.stop()

st.warning("⚠️ Mixing matters! Proceed to **Test 2** to identify which mixing scale is responsible.")

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
    key="bp_tfeed_ctr",
)

t_feed_fast = t_feed_center / 3.0
t_feed_slow = t_feed_center * 3.0

feed_vol = st.number_input("Total feed volume (mL)", min_value=0.1, value=100.0, key="bp_feed_vol")

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
st.dataframe(t2_conditions, width="stretch", hide_index=True)
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

st.subheader("Record Test 2 Outcome")
t2_result = st.radio(
    "Did the process response change when feed rate was varied?",
    ["— Select —", "Yes – process response changed", "No – no significant change"],
    key="bp_t2_result",
    horizontal=True,
)

if t2_result.startswith("No"):
    st.info(
        "🔬 **Micromixing appears to control the process.**  "
        "The feed plume disintegrates quickly; the molecular-scale engulfment "
        "step is rate-limiting."
    )
    st.markdown("""
    **Scale-up recommendation (Bourne 2003):** Hold **local energy dissipation**
    (ε_loc) constant where the mixing and reaction occur.

    You may optionally run the **confirmatory experiments** below to strengthen
    this conclusion.
    """)
    _micro_conclusion = True
elif t2_result.startswith("Yes"):
    st.warning(
        "⚠️ **Mesomixing (feed-plume dispersion) matters!** Proceed to **Test 3** "
        "to distinguish mesomixing from macromixing."
    )
    _micro_conclusion = False
else:
    st.info("Run the experiments, then record the outcome to continue.")
    st.stop()

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 3 – TEST 3: Feed Location  (Meso vs Macro)
# ══════════════════════════════════════════════════════════════════════════
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

# Compute local ε estimates for the user
eps_avg = power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3)
eps_avg_kg = eps_avg / rho  # W/kg for micromixing calculations

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
st.dataframe(t3_locs, width="stretch", hide_index=True)
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

st.subheader("Record Test 3 Outcome")
t3_result = st.radio(
    "Did the process response change when feed location was varied?",
    ["— Select —", "Yes – process response changed", "No – no significant change"],
    key="bp_t3_result",
    horizontal=True,
)

if t3_result.startswith("No"):
    st.info(
        "🌀 **Macromixing appears to control the process.**  "
        "The bulk blending / circulation time is the rate-limiting mixing step."
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
elif t3_result.startswith("Yes"):
    st.success(
        "📐 **Mesomixing controls the process.**  "
        "The disintegration of the feed plume in the local mixing environment "
        "is the rate-limiting mixing step."
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
else:
    st.info("Run the experiments, then record the outcome to continue.")
    st.stop()

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
        key="bp_fp_result",
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
        key="bp_visc_result",
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
conclusions.append(("Test 1 – Impeller Speed", "**Sensitive** → Mixing matters", "⚠️"))

# Test 2
if t2_result.startswith("Yes"):
    conclusions.append(("Test 2 – Feed Rate", "**Sensitive** → Not purely micromixing-controlled", "⚠️"))
elif t2_result.startswith("No"):
    conclusions.append(("Test 2 – Feed Rate", "**Insensitive** → Micromixing-controlled", "🔬"))
    scaleup_notes.append("Hold **local ε** constant at the feed point on scale-up.")

# Test 3
if t3_result.startswith("Yes"):
    conclusions.append(("Test 3 – Feed Location", "**Sensitive** → Mesomixing-controlled", "📐"))
    scaleup_notes.append("Control feed plume dispersion: constant impeller speed, extended feed time, or multiple feed points.")
elif t3_result.startswith("No"):
    conclusions.append(("Test 3 – Feed Location", "**Insensitive** → Macromixing-controlled", "🌀"))
    scaleup_notes.append("Maintain short blend times on scale-up (high-efficiency impellers, static mixers).")

# Determine the dominant mixing limitation
if t2_result.startswith("No"):
    dominant = "Micromixing"
    dominant_icon = "🔬"
    dominant_color = "blue"
elif t3_result.startswith("Yes"):
    dominant = "Mesomixing"
    dominant_icon = "📐"
    dominant_color = "orange"
elif t3_result.startswith("No"):
    dominant = "Macromixing"
    dominant_icon = "🌀"
    dominant_color = "red"
else:
    dominant = "Undetermined"
    dominant_icon = "❓"
    dominant_color = "gray"

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
# SECTION 6 – Decision Tree Diagram
# ══════════════════════════════════════════════════════════════════════════
st.header("6 · Decision Tree Reference")
st.markdown("""
```
                ┌───────────────────────────────┐
                │  TEST 1: Vary Impeller Speed  │
                │  (~100× change in P/V)        │
                └──────────────┬────────────────┘
                               │
               ┌───────────────┴──────────────-─┐
               │                                │
        No effect                         Effect observed
               │                                │
     ┌─────────▼─────────┐           ┌──────────▼──────────┐
     │ MIXING IS NOT     │           │ TEST 2: Vary Feed   │
     │ CRITICAL — EXIT   │           │ Rate (9× range)     │
     └───────────────────┘           └──────────┬──────────┘
                                                │
                                ┌───────────────┴──────────────-─┐
                                │                                │
                         No effect                         Effect observed
                                │                                │
                   ┌────────────▼────────────┐     ┌────────────-▼───────────┐
                   │ MICROMIXING CONTROLS    │     │ TEST 3: Vary Feed       │
                   │ → Hold local ε constant │     │ Location                │
                   └─────────────────────────┘     └────────────┬───────────-┘
                                                                │
                                                ┌───────────────┴───────────────┐
                                                │                               │
                                         No effect                        Effect observed
                                                │                               │
                                   ┌────────────▼────────────┐     ┌────────────▼────────────┐
                                   │ MACROMIXING CONTROLS    │     │ MESOMIXING CONTROLS     │
                                   │ → Improve blend time    │     │ → Control feed plume    │
                                   └─────────────────────────┘     └─────────────────────────┘
```
""")

# ══════════════════════════════════════════════════════════════════════════
# SECTION 7 – Save Protocol Record
# ══════════════════════════════════════════════════════════════════════════
st.header("7 · Save Protocol Record")

if st.button("📌 Save Bourne Protocol result to Recorded Results", key="bp_save"):
    result_row = {
        "reactor": reactor_name if not reactors.empty else "Manual entry",
        "fluid": fluid_name if not fluids.empty else "Manual entry",
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
# SECTION 8 – References
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
