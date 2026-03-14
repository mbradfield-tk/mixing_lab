"""
Page 10 – Mixing Sensitivity Protocol (Interactive Decision Tree)
==================================================================
Guides the user through a structured assessment to determine which
mixing mechanisms (micro-, meso-, macromixing, mass transfer, heat
transfer) are likely to limit a chemical reaction at scale.

The decision tree proceeds through:
  0. Pre-screening – Bourne Protocol Part 1 quick screen
  1. Kinetics – are reaction kinetics available?
  2. Phases – single-phase vs multi-phase considerations
  3. Competing reactions – mesomixing / feed sensitivity
  4. Heat transfer – exothermicity screening
  5. Mixing-time comparison – Damköhler-based screening
  6. Summary & recommendations
"""

import streamlit as st
import pandas as pd
import numpy as np
import pathlib
import graphviz

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"


# ── helpers ──────────────────────────────────────────────────────────────
def _load(key: str, fn: str) -> pd.DataFrame:
    if key not in st.session_state:
        p = DATA_DIR / fn
        st.session_state[key] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    return st.session_state[key]


def _safe(val, default=0.0):
    try:
        v = float(val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


reactions = _load("reaction_db", "reactions.csv")
fluids = _load("fluid_db", "fluids.csv")

# Session-state key prefix for this page
_PFX = "_msp_"


def _key(name: str) -> str:
    return f"{_PFX}{name}"


def _get(name: str, default=None):
    return st.session_state.get(_key(name), default)


def _set(name: str, value):
    st.session_state[_key(name)] = value


# ── Page header ──────────────────────────────────────────────────────────
st.title("🧭 Mixing Sensitivity Protocol")
st.caption(
    "An interactive decision tree to assess which mixing mechanisms "
    "may limit your reaction at scale."
)

# ── Visual overview of the decision tree ─────────────────────────────────
_DOT_SOURCE = """
    digraph protocol {
        rankdir=TB
        fontname="Arial"
        node [fontname="Arial" fontsize=10 style=filled shape=box
              fillcolor="#F5F5F5" color="#BDBDBD" margin="0.15,0.08"]
        edge [fontname="Arial" fontsize=9 color="#616161"]

        /* start / end nodes */
        START [label="🧭  Start Protocol" shape=Mrecord
               fillcolor="#4CAF50" fontcolor=white fontsize=11]
        SUM   [label="6 · Summary & \\nRecommendations" shape=Mrecord
               fillcolor="#2196F3" fontcolor=white fontsize=11]

        /* Step 0 – Bourne pre-screening */
        BOURNE [label="0 · Bourne Protocol\\nPart 1 suggested pre-screen" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        BOURNE_DO [label="Run Bourne Protocol\\nPart 1 (quick screen)\\n→ vary P/V (i.e stir speed)"
                   fillcolor="#FFF3E0" color="#FFB74D"]
        BOURNE_OK [label="✅ Pre-screen\\ncomplete"
                   fillcolor="#E8F5E9" color="#81C784"]

        /* Step 1 – Kinetics */
        K     [label="1 · Kinetics\\navailable?" shape=diamond
               fillcolor="#E3F2FD" color="#90CAF9"]
        KACT  [label="Measure kinetics\\n& calorimetry\\n→ add to Reaction DB"
               fillcolor="#FFF3E0" color="#FFB74D"]
        KAPPROX [label="Select approximate\\nkinetics from DB\\n(similar reaction)"
                 fillcolor="#FFF8E1" color="#FFD54F"]
        KSEL  [label="Select reaction\\nfrom database"
               fillcolor="#E8F5E9" color="#81C784"]

        /* Step 2 – Phases */
        PH    [label="2 · How many\\nphases?" shape=diamond
               fillcolor="#E3F2FD" color="#90CAF9"]
        PH_OK [label="✅ Mass transfer\\nnot a factor"
               fillcolor="#E8F5E9" color="#81C784"]
        PH_MT [label="⚠️ Interphase mass\\ntransfer may limit\\n→ characterise kLa, k_SL"
               fillcolor="#FFF8E1" color="#FFD54F"]

        /* Step 3 – Competing reactions */
        CR      [label="3 · Competing\\nreactions?" shape=diamond
                 fillcolor="#E3F2FD" color="#90CAF9"]
        MESO    [label="⚠️ Micro/meso-mixing limitation possible\\n→ Bourne Protocol"
                 fillcolor="#FFF8E1" color="#EF9A9A"]
        MESO_OK [label="✅ Micro/meso-mixing\\nnot a factor"
                 fillcolor="#E8F5E9" color="#81C784"]

        /* Step 4 – Heat */
        HT     [label="4 · ΔH data\\navailable?" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        HT_MEASURE [label="Perform calorimetry\\n→ measure ΔH"
                    fillcolor="#FFF3E0" color="#FFB74D"]
        HT_EST [label="Estimate ΔH from\\nsimilar reaction"
                fillcolor="#FFF8E1" color="#FFD54F"]
        HT_CHK [label="|ΔH| ≥ 50\\nkJ/mol?" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        HT_HOT [label="🔴 Heat-sensitive\\n→ run heat balance"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        HT_OK  [label="🟢 Modest\\nexothermicity"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Step 5 – Mixing time */
        MIX    [label="5 · Reaction\\ntime t_rxn" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        MICRO  [label="🔴 Micromixing\\nlikely sensitive"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        MACRO  [label="🟡 Macromixing\\nmay matter at scale"
                fillcolor="#FFF8E1" color="#FFD54F"]
        MIX_OK [label="🟢 Mixing unlikely\\nto be limiting"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Edges — Step 0 */
        START -> BOURNE

        BOURNE -> BOURNE_DO [label="Not done\\nyet"]
        BOURNE -> BOURNE_OK [label="Done /\\nSkip"]
        BOURNE_DO -> BOURNE_OK [label="Complete"]
        BOURNE_OK -> K

        /* Edges — Step 1 */
        K -> KACT    [label="No kinetics"]
        K -> KAPPROX [label="Approximate"]
        K -> KSEL    [label="Yes"]
        KACT -> K    [label="Data obtained\\n→ repeat Step 1"]
        KAPPROX -> PH
        KSEL -> PH

        /* Edges — Step 2 */
        PH -> PH_OK [label="Single\\nliquid"]
        PH -> PH_MT [label="Multi-phase\\n(G / L / S)"]

        PH_OK -> CR
        PH_MT -> CR

        /* Edges — Step 3 */
        CR -> MESO    [label="Yes /\\nNot sure"]
        CR -> MESO_OK [label="No"]

        MESO    -> HT
        MESO_OK -> HT

        /* Edges — Step 4 */
        HT -> HT_MEASURE [label="No ΔH"]
        HT -> HT_EST     [label="Estimate"]
        HT -> HT_CHK     [label="Yes"]
        HT_MEASURE -> HT [label="Data obtained\\n→ repeat Step 4"]
        HT_EST -> HT_CHK
        HT_CHK -> HT_HOT [label="Yes"]
        HT_CHK -> HT_OK  [label="No"]

        /* Edges — Step 5 */
        HT_HOT -> MIX
        HT_OK  -> MIX

        MIX -> MICRO  [label="< 1 s"]
        MIX -> MACRO  [label="1 – 10 s"]
        MIX -> MIX_OK [label="> 10 s"]

        MICRO  -> SUM
        MACRO  -> SUM
        MIX_OK -> SUM
    }
"""

with st.expander("📋 Protocol overview", expanded=True):
    st.graphviz_chart(_DOT_SOURCE, use_container_width=False)

    # Export buttons
    _graph = graphviz.Source(_DOT_SOURCE)
    _col_png, _col_svg, _ = st.columns([1, 1, 4])
    with _col_png:
        st.download_button(
            "⬇️ Download PNG",
            data=_graph.pipe(format="png"),
            file_name="mixing_sensitivity_protocol.png",
            mime="image/png",
        )
    with _col_svg:
        st.download_button(
            "⬇️ Download SVG",
            data=_graph.pipe(format="svg"),
            file_name="mixing_sensitivity_protocol.svg",
            mime="image/svg+xml",
        )

if st.button("🔄 Restart protocol", key=_key("restart")):
    for k in list(st.session_state.keys()):
        if k.startswith(_PFX) and k != _key("restart"):
            del st.session_state[k]
    st.rerun()

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 0 – BOURNE PROTOCOL PRE-SCREENING
# ══════════════════════════════════════════════════════════════════════════
st.header("0 · Bourne Protocol Pre-Screening")

with st.expander("ℹ️ Background — Bourne pre-screening", expanded=False):
    st.markdown(
        "Before diving into the full assessment, a quick **Bourne Protocol "
        "Part 1** experiment can reveal whether the reaction shows *any* "
        "mixing sensitivity.  This involves running the reaction at two or "
        "more impeller speeds and/or feed times and checking for changes "
        "in yield or impurity profile."
    )

bourne_screen = st.radio(
    "Have you performed (or do you want to skip) the Bourne pre-screen?",
    [
        "Already done – results show mixing sensitivity",
        "Already done – no mixing sensitivity observed",
        "Skip for now – proceed with protocol",
    ],
    key=_key("bourne_screen"),
)

if bourne_screen == "Already done – results show mixing sensitivity":
    st.warning(
        "**Mixing sensitivity confirmed experimentally.**  "
        "The remaining steps will help identify *which* mechanisms "
        "(micro-, meso-, macromixing, heat transfer) are responsible."
    )
    _bourne_sensitive = True
elif bourne_screen == "Already done – no mixing sensitivity observed":
    st.success(
        "✅ **No mixing sensitivity observed** in the pre-screen.  "
        "Proceed through the remaining steps to confirm and identify "
        "any latent risks at larger scale."
    )
    _bourne_sensitive = False
else:
    st.info(
        "💡 **Recommendation:** Performing Bourne Part 1 before this "
        "protocol provides a direct experimental answer.  You can run it "
        "from the **🧫 Bourne Protocol** page.\n\n"
        "Proceeding with the theoretical assessment for now."
    )
    _bourne_sensitive = None

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 1 – KINETICS
# ══════════════════════════════════════════════════════════════════════════
st.header("1 · Reaction Kinetics")

with st.expander("ℹ️ Background — why kinetics matter", expanded=False):
    st.markdown(
        "Before assessing mixing sensitivity we need kinetic data "
        "(rate constant *k*, reaction order, concentration *C₀*) and "
        "thermochemical data (enthalpy of reaction *ΔH*).\n\n"
        "The characteristic reaction time $t_{rxn}$ is the benchmark against "
        "which all mixing and transport times are compared via Damköhler "
        "numbers.  Without it, the protocol cannot determine whether mixing "
        "or mass transfer is fast or slow relative to the chemistry."
    )

kinetics_choice = st.radio(
    "Do you have kinetics available in the Reaction Database?",
    [
        "Yes – select from database",
        "Approximate – use a similar reaction as proxy",
        "No – I need to measure them",
    ],
    key=_key("kinetics_choice"),
    index=0 if _get("kinetics_choice") is None else None,
)

if kinetics_choice == "No – I need to measure them":
    st.warning(
        "**Action required before proceeding:**\n\n"
        "1. Conduct **time-course experiments** (e.g. in-situ FTIR / ReactIR, "
        "sampling + HPLC) to determine the rate constant *k* and reaction order.\n"
        "2. Perform **reaction calorimetry** (e.g. RC1, µRC, or Simular) to "
        "measure the heat of reaction *ΔH*.\n"
        "3. Add the results to the **🧪 Reaction Database** page.\n"
        "4. Return here and **select the reaction from the database** to continue."
    )
    st.info(
        "💡 **Tip:** If exact kinetics are unavailable but a structurally "
        "similar reaction is in the database, choose **\"Approximate\"** above "
        "to select it as a proxy and proceed with the assessment."
    )
    st.stop()

_using_approximate = kinetics_choice == "Approximate – use a similar reaction as proxy"

if _using_approximate:
    st.info(
        "🔶 You are using **approximate kinetics** from a similar reaction.  "
        "Results should be treated as indicative — confirm with measured data "
        "when available."
    )

# ── Reaction selection ────────────────────────────────────────────────────
if reactions.empty:
    st.error("The Reaction Database is empty.  Please add at least one reaction first.")
    st.stop()

_rxn_names = reactions["reaction_name"].tolist()
_rxn_idx = _rxn_names.index(st.session_state["_sel_msp_rxn"]) if st.session_state.get("_sel_msp_rxn") in _rxn_names else 0
rxn_name = st.selectbox(
    "Select reaction", _rxn_names, index=_rxn_idx, key=_key("rxn_sel")
)
st.session_state["_sel_msp_rxn"] = rxn_name
rxn = reactions[reactions["reaction_name"] == rxn_name].iloc[0]
rxn_order = str(rxn.get("order", "1"))
rxn_k = _safe(rxn.get("k_value"), 0.0)
rxn_C0 = _safe(rxn.get("C0_mol_L"), 0.0)
rxn_t_rxn = _safe(rxn.get("t_rxn_s"), 0.0)
rxn_delta_H = _safe(rxn.get("delta_H_kJ_mol"), 0.0)
rxn_T = _safe(rxn.get("T_C"), 25.0)
rxn_solvent = str(rxn.get("solvent", ""))

# Compute characteristic reaction time
if rxn_t_rxn > 0:
    t_rxn = rxn_t_rxn
elif rxn_k > 0:
    if rxn_order in ("1", "pseudo-1"):
        t_rxn = np.log(2) / rxn_k
    elif rxn_order in ("2", "pseudo-2") and rxn_C0 > 0:
        t_rxn = 1.0 / (rxn_k * rxn_C0)
    else:
        t_rxn = 1.0
else:
    t_rxn = 0.0

with st.expander("Selected reaction details", expanded=False):
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Order", rxn_order)
    col2.metric("k", f"{rxn_k:.4g}")
    col3.metric("C₀ (mol/L)", f"{rxn_C0:.4g}")
    col4.metric("t_rxn (s)", f"{t_rxn:.4g}")
    col5, col6, col7, _ = st.columns(4)
    col5.metric("T (°C)", f"{rxn_T:.0f}")
    col6.metric("ΔH (kJ/mol)", f"{rxn_delta_H:.1f}" if rxn_delta_H != 0 else "N/A")
    col7.metric("Solvent", rxn_solvent or "—")

if t_rxn <= 0:
    st.error(
        "Cannot determine a characteristic reaction time from the selected "
        "reaction data. Please check *k*, *C₀*, and *t_rxn* in the Reaction Database."
    )
    st.stop()

st.success(f"✅ Kinetics available — characteristic reaction time **t_rxn = {t_rxn:.4g} s**.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 2 – PHASES
# ══════════════════════════════════════════════════════════════════════════
st.header("2 · Phase Assessment")

with st.expander("ℹ️ Background — interphase mass transfer", expanded=False):
    st.markdown(
        "The number of phases present determines whether interphase mass "
        "transfer could be a rate-limiting step.\n\n"
        "In a multi-phase system, the reactant must cross a phase boundary "
        "(e.g. gas → liquid, or dissolution of a solid) before it can react. "
        "Each transport step has a characteristic time scale; if any of these "
        "is longer than the reaction time, the overall rate becomes "
        "transport-limited rather than kinetics-limited."
    )

phases = st.multiselect(
    "Which phases are present in the reaction system?",
    ["Liquid", "Solid", "Gas"],
    default=["Liquid"],
    key=_key("phases"),
)

_multiphase = len(phases) > 1
_has_solid = "Solid" in phases
_has_gas = "Gas" in phases

if _multiphase:
    st.warning("**Multi-phase system detected**")

    _mass_transfer_notes = []
    if _has_gas and not _has_solid:
        _mass_transfer_notes = [
            "**Gas → Liquid mass transfer** can be rate-limiting.  "
            "The volumetric mass-transfer coefficient *kLa* depends on "
            "agitation intensity, gas flow rate, and fluid properties.",
            "At scale, gas hold-up and bubble size distribution change with "
            "impeller design, vessel geometry, and sparger type.",
            "The gas-liquid Damköhler number **Da_GL = 1 / (kLa · t_rxn)** "
            "(= transport time / reaction time) indicates whether mass transfer "
            "or kinetics controls the rate.",
        ]
    elif _has_solid and not _has_gas:
        _mass_transfer_notes = [
            "**Liquid → Solid mass transfer** can be rate-limiting.  "
            "The solid-liquid mass-transfer coefficient *k_SL* depends on "
            "particle size, settling velocity and turbulence.",
            "Particle suspension (Zwietering *N_js*) must be ensured — "
            "unsuspended solids drastically reduce effective surface area.",
            "At scale, maintaining uniform solids distribution may require "
            "higher specific power input or different impeller configuration.",
        ]
    else:
        _mass_transfer_notes = [
            "**Three-phase (gas–liquid–solid) system** — multiple mass-transfer "
            "resistances act in series: gas → liquid → solid surface.",
            "Each transport step has its own rate, and the slowest step controls.",
            "Scale-up is particularly challenging: gas-liquid *kLa*, solid "
            "suspension *N_js*, and liquid-solid *k_SL* may all scale differently.",
        ]

    _mass_transfer_notes.extend([
        "It is **highly recommended** to characterise each vessel's "
        "hydrodynamics (see ⚙️ Mixing Sensitivity and 📊 Reactor Comparison pages) "
        "to understand how each transport mechanism varies with scale.",
        "If the intrinsic reaction is fast relative to any mass-transfer step, "
        "the observed rate will be transport-limited and mixing-sensitive.",
    ])

    for note in _mass_transfer_notes:
        st.markdown(f"- {note}")

    # ── Damköhler numbers for mass transfer ──────────────────────────────
    with st.expander("📐 Damköhler Numbers for Mass Transfer — theory", expanded=False):
        st.markdown(
            "Just as the mixing Damköhler number compares mixing time to "
            "reaction time (Step 5), the **mass-transfer Damköhler number** "
            "compares the interphase transport time to the reaction time.  "
            "When Da_MT ≪ 1 transport is fast (reaction-limited); "
            "when Da_MT ≫ 1 transport is the bottleneck (mass-transfer-limited)."
        )

        if _has_gas:
            st.markdown("##### Gas–Liquid Transport")
            st.markdown(
                "The gas-liquid Damköhler number is defined as:\n\n"
                "$$Da_{GL} = \\frac{1}{k_L a \\; t_{rxn}} = \\frac{t_{\\text{transfer}}}{t_{rxn}}$$\n\n"
                "where:\n"
                "- $k_L a$ = volumetric gas-liquid mass-transfer coefficient (s⁻¹)\n"
                "- $t_{rxn}$ = characteristic reaction time (s)\n"
                "- $1 / k_L a$ = characteristic gas-liquid transport time (s)\n\n"
                "This is the ratio of the **transport time** to the **reaction time**, "
                "consistent with the definition on the 📐 Mixing & Damköhler page."
            )
            st.markdown(
                "| Da_GL range | Regime |\n"
                "|-------------|--------|\n"
                "| Da_GL ≪ 1 | **Reaction-limited** — transport is fast relative to reaction |\n"
                "| Da_GL ≈ 1 | **Transition** — both transport and kinetics matter |\n"
                "| Da_GL ≫ 1 | **Mass-transfer-limited** — gas absorption controls the rate |\n\n"
                "For fast reactions, the **Hatta number** Ha = √($k \\cdot D_A$) / $k_L$ "
                "is also useful: when Ha > 3 the reaction occurs within the liquid film "
                "and enhancement of absorption must be considered."
            )

        if _has_solid:
            st.markdown("##### Liquid–Solid Transport")
            st.markdown(
                "The liquid-solid Damköhler number is defined as:\n\n"
                "$$Da_{LS} = \\frac{1}{k_{SL} \\cdot a_s \\; t_{rxn}} = \\frac{t_{\\text{transfer,LS}}}{t_{rxn}}$$\n\n"
                "where:\n"
                "- $k_{SL}$ = solid-liquid mass-transfer coefficient (m/s)\n"
                "- $a_s$ = specific surface area of the solids (m²/m³), "
                "depends on particle size ($a_s$ ≈ 6 / $d_p$ for spheres)\n"
                "- $1 / (k_{SL} \\cdot a_s)$ = characteristic liquid-solid transport time (s)\n\n"
                "This is the ratio of the **transport time** to the **reaction time**, "
                "consistent with the Da_GL convention above."
            )
            st.markdown(
                "| Da_LS range | Regime |\n"
                "|-------------|--------|\n"
                "| Da_LS ≪ 1 | **Reaction-limited** — dissolution is fast relative to reaction |\n"
                "| Da_LS ≈ 1 | **Transition** — both transport and reaction compete |\n"
                "| Da_LS ≫ 1 | **Mass-transfer-limited** — particle dissolution rate controls |\n\n"
                "Solid-liquid transport also depends on particle suspension quality. "
                "If particles are not fully suspended (below Zwietering $N_{js}$), "
                "the effective surface area drops dramatically."
            )

    if _has_gas:
        _da_gl = 1.0 / (0.05 * t_rxn) if t_rxn > 0 else 0.0  # placeholder: typical kLa ~ 0.05 s⁻¹
        st.caption(
            f"💡 Example: for a typical kLa ≈ 0.05 s⁻¹ and your t_rxn = {t_rxn:.4g} s → "
            f"Da_GL ≈ {_da_gl:.3g}.  {'Transport likely fast enough.' if _da_gl < 1 else 'Mass transfer may limit — measure kLa for your vessel.'}"
        )

    st.info(
        "📌 **Proceed with caution** — mass-transfer limitations are system-dependent. "
        "The remaining steps of this protocol will still assess micro/meso/macromixing "
        "and heat transfer, but keep in mind that interphase transport may dominate.\n\n"
        "👉 Compute Da_GL and/or Da_LS for your specific reactor on the "
        "**⚙️ Mixing Sensitivity** or **📊 Reactor Comparison** pages."
    )
else:
    st.success(
        "✅ **Single liquid phase** — interphase mass transfer is not a factor.  "
        "Micro-, meso-, and macromixing may still affect the reaction."
    )

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 3 – COMPETING REACTIONS / MICRO- & MESOMIXING
# ══════════════════════════════════════════════════════════════════════════
st.header("3 · Competing Reactions — Micro- & Mesomixing")

with st.expander("ℹ️ Background — micromixing & mesomixing theory", expanded=False):
    st.markdown(
        "When multiple reactions compete for the same reagent, incomplete "
        "mixing at two distinct length scales can shift selectivity:\n\n"
        "- **Micromixing** (molecular scale) — the rate at which fluid "
        "elements are homogenised at the Kolmogorov / Batchelor scale.  "
        "Governed by the local turbulent energy dissipation rate *ε*.  "
        "The Engulfment model (E-model, Baldyga & Bourne) describes how "
        "fresh eddies engulf surrounding fluid: *E* ≈ 0.058 √(*ε* / *ν*).  "
        "If the reaction is faster than engulfment, local concentration "
        "gradients at the molecular level control selectivity.\n\n"
        "- **Mesomixing** (feed-point scale) — the rate at which the feed "
        "stream is dispersed into the turbulent bulk by the mean flow and "
        "large-eddy convection.  Feed-point turbulence, feed rate, and "
        "pipe diameter determine the initial plume dilution rate.  "
        "Selectivity is affected when the feed is consumed faster than it "
        "can be dispersed."
    )

competing = st.radio(
    "Are there competing (parallel or consecutive) reactions that could "
    "reduce yield or increase impurity formation?",
    ["Yes", "No", "Not sure"],
    key=_key("competing"),
)

if competing == "Yes":
    st.warning(
        "**Micro- and/or mesomixing sensitivity likely.**\n\n"
        "When competing reactions exist, both molecular-level and "
        "feed-point mixing can influence selectivity."
    )

    # Micromixing sub-section
    with st.expander("🔬 Micromixing effects — theory & details", expanded=False):
        st.markdown(
            "At the smallest scales, the Engulfment rate *E* determines how "
            "fast reactant and surrounding fluid reach molecular homogeneity.  "
            "Key points:\n\n"
            "- The **micromixing time** $t_{micro}$ ≈ 17.3 √(*ν* / *ε*) "
            "(Baldyga & Bourne, 1999) depends on kinematic viscosity *ν* and "
            "local energy dissipation *ε* at the feed point.\n"
            "- When $Da_{micro}$ = $t_{micro}$ / $t_{rxn}$ > 1, the reaction "
            "is faster than micromixing — **selectivity is micromixing-controlled**.\n"
            "- In stirred vessels, *ε* varies by 1–2 orders of magnitude "
            "across the tank; the value at the feed point is what matters.\n"
            "- At larger scale, local *ε* near the impeller typically "
            "**decreases** for constant P/V, making micromixing limitations "
            "worse."
        )

    # Mesomixing sub-section
    with st.expander("🌊 Mesomixing effects — theory & details", expanded=False):
        st.markdown(
            "At the feed-point scale, incomplete turbulent dispersion of the "
            "feed plume creates concentration gradients larger than the "
            "Kolmogorov scale.  Critical factors:\n\n"
            "- **Feed rate / addition time** — slower addition gives the impeller "
            "more time to disperse the feed before local concentrations build up.\n"
            "- **Feed location** — feeding near the impeller (high-turbulence zone) "
            "promotes rapid incorporation; feeding at the surface can create "
            "concentration pockets.\n"
            "- **Number of feed points** — multiple injection points reduce local "
            "super-saturation.\n"
            "- **Feed pipe diameter** — affects the feed jet momentum and initial "
            "mixing at the nozzle.\n"
            "- The **mesomixing time** can be estimated from the turbulent "
            "dispersion model: $t_{meso}$ ≈ 0.15 ($Q_f$ / $k_T$ $ε$)^(1/3), "
            "where $Q_f$ is feed flow rate and $k_T$ is a turbulence constant."
        )

    st.info(
        "👉 **Recommendation:** Use the **🧫 Bourne Protocol** page to "
        "experimentally screen for micro- and mesomixing sensitivity by varying "
        "impeller speed (changes *ε*), feed time, and feed location.  "
        "The Villermaux–Dushman (iodide–iodate) test reaction is a classic "
        "tool to decouple micro- and mesomixing effects (Bourne, 2003)."
    )
    _meso_sensitive = True
elif competing == "Not sure":
    st.info(
        "If you are unsure whether side-reactions occur, consider:\n\n"
        "- Run the reaction at two different impeller speeds **and** two "
        "different addition rates.  If yield or impurity profile changes "
        "with speed, micromixing is involved; if it changes with feed time, "
        "mesomixing is involved.\n"
        "- Check the literature or process-chemistry knowledge for known "
        "by-products or degradation pathways.\n"
        "- If in doubt, treat the system as potentially sensitive and "
        "use the **🧫 Bourne Protocol** to screen."
    )
    _meso_sensitive = True
else:
    st.success(
        "✅ **No competing reactions identified** — micro- and mesomixing "
        "(selectivity-related) are unlikely to be a factor.  "
        "Macromixing (blend time) may still matter at scale."
    )
    _meso_sensitive = False

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 4 – HEAT TRANSFER
# ══════════════════════════════════════════════════════════════════════════
st.header("4 · Heat Transfer Screening")

_has_enthalpy = rxn_delta_H != 0.0

if not _has_enthalpy:
    st.warning(
        "The selected reaction has **no enthalpy (ΔH) data** in the database."
    )

    dh_action = st.radio(
        "How would you like to proceed?",
        [
            "Perform calorimetry – measure ΔH experimentally (recommended)",
            "Estimate ΔH from a similar reaction",
        ],
        key=_key("dh_action"),
    )

    if dh_action == "Perform calorimetry – measure ΔH experimentally (recommended)":
        st.info(
            "**Action required:**\n\n"
            "1. Perform **reaction calorimetry** (e.g. RC1, µRC, or Simular) "
            "to measure the heat of reaction *ΔH*.\n"
            "2. Update the reaction entry in the **🧪 Reaction Database**.\n"
            "3. Return here — the protocol will pick up the new value automatically."
        )
        _heat_sensitive = False
        st.stop()
    else:
        # Estimate: let user pick a similar reaction from the DB
        st.info(
            "Select a reaction from the database whose enthalpy can serve "
            "as a rough estimate for your system."
        )
        _dh_candidates = reactions[
            reactions["delta_H_kJ_mol"].notna()
            & (reactions["delta_H_kJ_mol"] != 0)
        ]
        if _dh_candidates.empty:
            st.error(
                "No reactions with ΔH data found in the database.  "
                "Please add at least one reaction with calorimetry data, "
                "or perform calorimetry on the current reaction."
            )
            _heat_sensitive = False
            st.stop()

        est_rxn_name = st.selectbox(
            "Select reference reaction for ΔH estimate",
            _dh_candidates["reaction_name"].tolist(),
            key=_key("dh_est_rxn"),
        )
        _est_row = _dh_candidates[
            _dh_candidates["reaction_name"] == est_rxn_name
        ].iloc[0]
        rxn_delta_H = _safe(_est_row.get("delta_H_kJ_mol"), 0.0)
        _has_enthalpy = rxn_delta_H != 0.0
        if _has_enthalpy:
            st.warning(
                f"🔶 Using **estimated ΔH = {rxn_delta_H:.1f} kJ/mol** "
                f"from *{est_rxn_name}*.  Treat the heat-transfer assessment "
                f"as indicative — confirm with measured data when available."
            )

if _has_enthalpy:
    abs_dH = abs(rxn_delta_H)
    st.markdown(
        f"The selected reaction has **ΔH = {rxn_delta_H:.1f} kJ/mol** "
        f"({'exothermic' if rxn_delta_H < 0 else 'endothermic'})."
    )

    # Heuristic classification
    if abs_dH >= 100:
        _heat_class = "highly exothermic"
        _heat_color = "error"
    elif abs_dH >= 50:
        _heat_class = "moderately exothermic"
        _heat_color = "warning"
    elif abs_dH >= 20:
        _heat_class = "mildly exothermic"
        _heat_color = "info"
    else:
        _heat_class = "low exothermicity"
        _heat_color = "success"

    # Classification banner
    _heat_msg = f"**{_heat_class.title()}** — |ΔH| = {abs_dH:.1f} kJ/mol"
    if _heat_color == "error":
        st.error(f"🔴 {_heat_msg}")
    elif _heat_color == "warning":
        st.warning(f"🟡 {_heat_msg}")
    elif _heat_color == "info":
        st.info(f"🔵 {_heat_msg}")
    else:
        st.success(f"🟢 {_heat_msg}")

    # Detailed guidance
    if abs_dH >= 50:
        with st.expander("ℹ️ Background — heat-transfer scale-up considerations", expanded=False):
            st.markdown(
                "**Heat-transfer sensitivity is likely**, especially in vessels "
                "with low surface-area-to-volume ratios (i.e. larger scale).\n\n"
                "Key considerations:"
            )
            st.markdown(
                "- **Cooling capacity** scales with jacket area (∝ D²) while "
                "heat generation scales with volume (∝ D³).  "
                "This means large reactors are inherently harder to cool.\n"
                "- **Jacket ΔT** — the available temperature driving force depends "
                "on coolant type and the allowable process temperature range.\n"
                "- **Feed rate control** — for semi-batch reactions, limiting "
                "the addition rate can control the instantaneous heat load.\n"
                "- **Overall U** — wall material, fouling, and agitation intensity "
                "all affect the heat-transfer coefficient."
            )
        st.info(
            "👉 **Recommendation:** Run the heat balance on the "
            "**⚙️ Mixing Sensitivity** or **📊 Reactor Comparison** page "
            "to quantify Q_gen vs Q_cool for your specific reactor(s)."
        )
        _heat_sensitive = True
    else:
        st.markdown(
            "The reaction enthalpy is modest.  Heat transfer is **unlikely "
            "to be limiting** in most reactor configurations, but may become "
            "relevant at very large scale or in poorly cooled vessels.  "
            "Consider running a heat balance if scaling beyond pilot scale."
        )
        _heat_sensitive = False

    # Flag known highly-exothermic reaction types
    _rxn_type = str(rxn.get("type", "")).lower()
    _known_hot = ["grignard", "nitration", "sulfonation", "diazotization",
                   "polymerization", "hydrogenation", "oxidation"]
    _flagged = [t for t in _known_hot if t in _rxn_type]
    if _flagged:
        st.warning(
            f"⚠️ Reaction type **{rxn.get('type', '')}** is commonly "
            f"associated with significant exothermicity.  "
            f"Heat-transfer assessment is strongly recommended regardless "
            f"of the reported ΔH magnitude."
        )
        _heat_sensitive = True

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 5 – MIXING TIME COMPARISON (micro / macro)
# ══════════════════════════════════════════════════════════════════════════
st.header("5 · Mixing Time vs Reaction Time")
with st.expander("ℹ️ Background — Damköhler number interpretation", expanded=False):
    st.markdown(
        "The **Damköhler number** (Da) compares the characteristic mixing time "
        "to the reaction time:\n\n"
        r"$$\mathrm{Da} = \frac{t_{\text{mix}}}{t_{\text{rxn}}}$$" "\n\n"
        "When Da > 1 the reaction is faster than mixing and "
        "the system is **mixing-sensitive** at that scale.  "
        "Different mixing mechanisms have different characteristic times, "
        "so there is a separate Da for micromixing, mesomixing, and macromixing."
    )
    st.markdown(
        "| Da range | Interpretation |\n"
        "|----------|----------------|\n"
        "| Da < 0.01 | Reaction is much slower than mixing — **not sensitive** |\n"
        "| 0.01 – 0.1 | **Likely not sensitive**, but monitor at scale |\n"
        "| 0.1 – 1 | **Potentially sensitive** — mixing and reaction on similar timescales |\n"
        "| 1 – 10 | **Likely sensitive** — mixing limits observed rate |\n"
        "| Da > 10 | **Highly sensitive** — mixing fully controls the process |"
    )

st.markdown(f"**Your reaction time: t_rxn = {t_rxn:.4g} s**")

# Provide heuristic guidance based on t_rxn alone
if t_rxn < 0.1:
    st.error(
        f"🔴 **Very fast reaction** (t_rxn = {t_rxn:.4g} s).  "
        "**Micromixing-sensitive** in most reactor configurations.  "
        "Local turbulent energy dissipation near the impeller determines "
        "the effective mixing rate.  Feed location and impeller tip speed "
        "are critical parameters."
    )
    _micro_likely = True
elif t_rxn < 1.0:
    st.warning(
        f"🟡 **Fast reaction** (t_rxn = {t_rxn:.4g} s).  "
        "Micromixing likely relevant in larger vessels where local ε "
        "at the feed point decreases.  Confirm with Damköhler analysis "
        "on the ⚙️ Mixing Sensitivity page."
    )
    _micro_likely = True
elif t_rxn < 10:
    st.info(
        f"🔵 **Moderate reaction** (t_rxn = {t_rxn:.4g} s).  "
        "Micromixing is less likely to dominate, but macromixing "
        "(blend time) could be relevant in larger vessels.  "
        "Check blend time relative to t_rxn."
    )
    _micro_likely = False
else:
    st.success(
        f"🟢 **Slow reaction** (t_rxn = {t_rxn:.4g} s).  "
        "Mixing is unlikely to limit the reaction in well-agitated vessels.  "
        "Only very large, poorly mixed tanks might approach Da ~ 1."
    )
    _micro_likely = False

with st.expander("ℹ️ Background — macromixing (blend time) at scale", expanded=False):
    st.markdown(
        "### What is macromixing?\n\n"
        "**Macromixing** describes the large-scale convective blending that "
        "distributes feed material throughout the entire vessel.  The characteristic "
        "time is the **95 % blend time** (θ₉₅) — the time required for a tracer "
        "to reach 95 % uniformity after a pulse injection.\n\n"
        "### How blend time scales with vessel size\n\n"
        "For a turbulent, baffled stirred tank the Grenville correlation gives:\n\n"
        r"$$\theta_{95} \;=\; 5.2 \, \frac{V}{N_Q \, N \, D^3}$$" "\n\n"
        "where *V* is the liquid volume, *N_Q* is the impeller pumping number, "
        "*N* is the impeller speed (rev/s), and *D* is the impeller diameter.\n\n"
        "At constant power-per-unit-volume (P/V), blend time scales approximately as:\n\n"
        r"$$\theta_{95} \;\propto\; T^{2/3}$$" "\n\n"
        "where *T* is the tank diameter.  This means a 10× increase in tank diameter "
        "leads to roughly a **4.6× increase** in blend time — a relationship that is "
        "often under-appreciated during scale-up.\n\n"
        "### Typical blend times by scale\n\n"
        "| Scale | Tank diameter | Typical θ₉₅ |\n"
        "|-------|---------------|-------------|\n"
        "| Lab (1–5 L) | 0.1–0.15 m | 1–5 s |\n"
        "| Pilot (50–200 L) | 0.3–0.5 m | 5–20 s |\n"
        "| Production (1–10 m³) | 1.0–2.0 m | 20–60 s |\n"
        "| Large production (10–50 m³) | 2.0–3.5 m | 40–120 s |\n\n"
        "### The macromixing Damköhler number\n\n"
        r"$$\mathrm{Da_{macro}} \;=\; \frac{\theta_{95}}{t_{\text{rxn}}}$$" "\n\n"
        "When Da_macro > 1 the blend time exceeds the reaction time, meaning "
        "reagents react before they are uniformly distributed.  This can cause "
        "local concentration gradients, reduced selectivity, and hot spots.\n\n"
        "### Why it matters\n\n"
        "- **Semi-batch additions:** Feed added at one point reacts locally before "
        "being blended into the bulk — leading to high local stoichiometric excess.\n"
        "- **Feed location:** Top-surface vs. sub-surface addition significantly "
        "affects the effective mixing path length.\n"
        "- **Impeller selection:** High-flow (axial) impellers reduce θ₉₅ more "
        "effectively than high-shear (radial) impellers at the same P/V.\n"
        "- **Scale-up trap:** A reaction insensitive to mixing at lab scale "
        "(θ₉₅ ≈ 2 s, Da_macro ≪ 1) can become sensitive in production "
        "(θ₉₅ ≈ 60 s, Da_macro > 1) — even without any micromixing issues."
    )

st.markdown(
    "**Macromixing** (blend time) typically ranges from ~1 s in lab "
    "reactors to 30–120 s in large production vessels.  If your reaction "
    "time is in that range, macromixing may matter at scale."
)

if t_rxn < 30:
    st.info(
        "👉 **Recommendation:** Compute Damköhler numbers for your specific "
        "reactor(s) on the **⚙️ Mixing Sensitivity** page or compare "
        "multiple vessels on the **📊 Reactor Comparison** page."
    )

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 6 – SUMMARY & RECOMMENDATIONS
# ══════════════════════════════════════════════════════════════════════════
st.header("6 · Summary & Recommendations")

if _using_approximate:
    st.warning(
        "⚠️ This assessment is based on **approximate kinetics** from a similar "
        "reaction. Results are indicative — confirm with measured data."
    )

findings: list[tuple[str, str, str]] = []  # (mechanism, status, detail)

# Bourne pre-screening
if _bourne_sensitive is True:
    findings.append((
        "Bourne pre-screen",
        "🔴 Mixing sensitivity confirmed",
        "Experimental pre-screen showed yield/impurity changes with mixing conditions.",
    ))
elif _bourne_sensitive is False:
    findings.append((
        "Bourne pre-screen",
        "🟢 No sensitivity observed",
        "Experimental pre-screen showed no mixing sensitivity at lab scale.",
    ))
else:
    findings.append((
        "Bourne pre-screen",
        "⚪ Not performed",
        "Consider running Bourne Protocol Part 1 for a direct experimental answer.",
    ))

# Micromixing
if _micro_likely:
    findings.append((
        "Micromixing",
        "🔴 Likely sensitive",
        f"t_rxn = {t_rxn:.4g} s — fast enough that local energy dissipation "
        "controls the mixing rate.",
    ))
else:
    findings.append((
        "Micromixing",
        "🟢 Unlikely",
        f"t_rxn = {t_rxn:.4g} s — reaction is slow relative to typical "
        "micromixing times.",
    ))

# Micro- & Mesomixing (selectivity)
if _meso_sensitive:
    findings.append((
        "Micro/mesomixing (selectivity)",
        "🟡 Potentially sensitive" if competing == "Not sure" else "🔴 Likely sensitive",
        "Competing reactions present — both micromixing (local ε) and "
        "mesomixing (feed dispersion) may affect selectivity.",
    ))
else:
    findings.append((
        "Micro/mesomixing (selectivity)",
        "🟢 Not a factor",
        "No competing reactions identified.",
    ))

# Macromixing
if t_rxn < 60:
    findings.append((
        "Macromixing (blend time)",
        "🟡 Check at scale",
        f"t_rxn = {t_rxn:.4g} s is within the range of blend times in "
        "larger vessels (10–120 s).  Compute Da_macro for your reactor.",
    ))
else:
    findings.append((
        "Macromixing (blend time)",
        "🟢 Unlikely",
        f"t_rxn = {t_rxn:.4g} s is much longer than typical blend times.",
    ))

# Mass transfer
if _multiphase:
    _mt_phases = " + ".join(phases)
    findings.append((
        f"Mass transfer ({_mt_phases})",
        "🟡 System-dependent",
        "Multi-phase system — interphase transport may limit the observed rate.  "
        "Characterise kLa and/or k_SL for each reactor.",
    ))
else:
    findings.append((
        "Mass transfer",
        "🟢 Not applicable",
        "Single liquid phase — no interphase transport.",
    ))

# Heat transfer
if _has_enthalpy and _heat_sensitive:
    findings.append((
        "Heat transfer",
        "🔴 Likely sensitive",
        f"|ΔH| = {abs(rxn_delta_H):.1f} kJ/mol — run a heat balance "
        "to confirm adequate cooling capacity.",
    ))
elif _has_enthalpy and not _heat_sensitive:
    findings.append((
        "Heat transfer",
        "🟢 Manageable",
        f"|ΔH| = {abs(rxn_delta_H):.1f} kJ/mol — modest exothermicity, "
        "unlikely to be limiting in most configurations.",
    ))
else:
    findings.append((
        "Heat transfer",
        "⚪ Unknown",
        "No ΔH data available — consider measuring by calorimetry.",
    ))

# Render summary table
summary_df = pd.DataFrame(findings, columns=["Mechanism", "Status", "Detail"])
for _, row in summary_df.iterrows():
    st.markdown(f"**{row['Mechanism']}** — {row['Status']}")
    st.caption(row["Detail"])

# Next steps
st.subheader("Recommended Next Steps")
_steps = []

if _bourne_sensitive is None:
    _steps.append(
        "Run **Bourne Protocol Part 1** (quick screen) on the "
        "**🧫 Bourne Protocol** page to confirm whether mixing sensitivity "
        "exists experimentally."
    )

if _using_approximate:
    _steps.append(
        "Obtain **measured kinetics and calorimetry** for the actual "
        "reaction to replace the approximate values used in this assessment."
    )

if _micro_likely or (t_rxn < 60):
    _steps.append(
        "Compute **Damköhler numbers** (Da_macro, Da_micro) for your specific "
        "reactor on the **⚙️ Mixing Sensitivity** page."
    )

if _meso_sensitive:
    _steps.append(
        "Run the **🧫 Bourne Protocol** to experimentally screen micro- and "
        "mesomixing sensitivity (vary impeller speed to probe micromixing, "
        "vary feed time/location to probe mesomixing)."
    )

if _multiphase:
    _steps.append(
        "Assess interphase **mass-transfer coefficients** (kLa, k_SL) on "
        "the ⚙️ Mixing Sensitivity or 📊 Reactor Comparison pages."
    )

if _has_enthalpy and _heat_sensitive:
    _steps.append(
        "Run a full **heat balance** (🔥 button on ⚙️ Mixing Sensitivity or "
        "📊 Reactor Comparison) to evaluate Q_gen vs Q_cool."
    )

if not _steps:
    _steps.append(
        "The reaction appears **low risk** for mixing sensitivity.  "
        "Standard scale-up practices should be sufficient, but consider "
        "a quick check on the 📊 Reactor Comparison page for completeness."
    )

for i, step in enumerate(_steps, 1):
    st.markdown(f"{i}. {step}")
