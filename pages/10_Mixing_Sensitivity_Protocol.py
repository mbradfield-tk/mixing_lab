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
  3. Competing reactions
  4. Heat transfer – exothermicity screening
  5. Mixing-time comparison – Damköhler-based screening
  6. Summary & recommendations
"""

import streamlit as st
import pandas as pd
import numpy as np
import pathlib

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

from utils.data_helpers import load_db, safe_float as _safe
from utils.solvent_properties import is_known_solvent, get_properties

reactions = load_db("reaction_db", "reactions.csv")
fluids = load_db("fluid_db", "fluids.csv")

# Session-state key prefix for this page
_PFX = "_msp_"

# ── Generation counter for reliable widget reset ─────────────────────────
if "_msp_gen" not in st.session_state:
    st.session_state["_msp_gen"] = 0


def _key(name: str) -> str:
    return f"{_PFX}{st.session_state['_msp_gen']}_{name}"


def _get(name: str, default=None):
    return st.session_state.get(_key(name), default)


def _set(name: str, value):
    st.session_state[_key(name)] = value


# ── Page header ──────────────────────────────────────────────────────────
st.title("🧭 Reaction Sensitivity Protocol")
st.caption(
    "An interactive decision tree to assess which mixing mechanisms "
    "may limit your reaction at scale."
)

# ── Visual overview of the decision tree ─────────────────────────────────
_IMG_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "general"
_MSP_IMG = _IMG_DIR / "mixing_sensitivity_protocol.png"

with st.expander("📋 Protocol overview", expanded=False):
    if _MSP_IMG.exists():
        st.image(str(_MSP_IMG), caption="Mixing Sensitivity Protocol – Decision Tree")
    else:
        st.info("Protocol overview image not found. Generate it from the 🛠️ Admin page.")

def _reset_msp():
    """Increment generation counter and clear old MSP state."""
    old_gen = st.session_state.get("_msp_gen", 0)
    for k in list(st.session_state.keys()):
        if k == "_msp_gen" or k == _key("restart"):
            continue
        if k.startswith(_PFX):
            del st.session_state[k]
    st.session_state["_msp_gen"] = old_gen + 1


st.button("🔄 Restart protocol", key=_key("restart"), on_click=_reset_msp)

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
    index=None,
)

if bourne_screen is None:
    st.info("👆 Please select an option above to continue.")
    st.stop()

if bourne_screen == "Already done – results show mixing sensitivity":
    st.warning(
        "**Mixing sensitivity confirmed experimentally.**  "
        "The remaining steps will help identify *which* mechanisms "
        "(micro-, meso-, macromixing, heat transfer) are responsible."
    )
    _bourne_sensitive = True

    _bourne_mech_opts = [
        "Micromixing (impeller-speed sensitive)",
        "Mesomixing (feed-rate / feed-time sensitive)",
        "Macromixing (feed-location sensitive)",
    ]
    _bourne_mech_sel = st.multiselect(
        "If the Bourne Protocol identified the controlling scale(s), select them "
        "(leave blank if not yet known):",
        _bourne_mech_opts,
        key=_key("bourne_mech"),
        help=(
            "Bourne Test 1 (impeller speed) probes **micromixing**; "
            "Test 2 (feed rate / time) probes **mesomixing**; "
            "Test 3 (feed location) probes **meso-/macromixing**."
        ),
    )
    _bourne_mechanisms = [m.split("(")[0].strip() for m in _bourne_mech_sel]
    if _bourne_mechanisms:
        st.caption(
            "Identified controlling scale(s): **" + ", ".join(_bourne_mechanisms) +
            "**. These are carried into the summary as experimentally confirmed "
            "sensitivities and drive the conclusions and recommendations below."
        )
    else:
        st.caption(
            "No specific scale selected — run the full Bourne Protocol (Tests 1–3) "
            "to pinpoint micro-, meso-, or macromixing control."
        )
elif bourne_screen == "Already done – no mixing sensitivity observed":
    st.success(
        "✅ **No mixing sensitivity observed** in the pre-screen.  "
        "Proceed through the remaining steps to confirm and identify "
        "any latent risks at larger scale."
    )
    _bourne_sensitive = False
    _bourne_mechanisms = []
else:
    st.info(
        "💡 **Recommendation:** Performing Bourne Part 1 before this "
        "protocol provides a direct experimental answer.  You can run it "
        "from the **🅱️ Bourne Protocol** page.\n\n"
        "Proceeding with the theoretical assessment for now."
    )
    _bourne_sensitive = None
    _bourne_mechanisms = []

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
    index=None,
)

if kinetics_choice is None:
    st.info("👆 Please select an option above to continue.")
    st.stop()

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
        # Characteristic time = reciprocal rate constant (time constant 1/k),
        # the e-folding time used in the Damkohler definition.
        t_rxn = 1.0 / rxn_k
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
st.caption(
    "Convention: t_rxn is the reaction time constant — 1/k (first order) or "
    "1/(k·C₀) (second order) — i.e. the e-folding time used in the Damköhler "
    "number, not the 50 % half-life (which is ~30 % shorter for first order)."
)

# ── Process mode ─────────────────────────────────────────────────────────
is_semi_batch = st.checkbox(
    "Semi-batch (fed-batch) process",
    value=False,
    key=_key("semi_batch"),
    help="Reagent or anti-solvent is dosed during the reaction.",
)

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
    default=[],
    key=_key("phases"),
)

if not phases:
    st.info("👆 Please select at least one phase above to continue.")
    st.stop()

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
        "hydrodynamics (see 🌀 Mixing Assessment and 📈 Reactor Comparison pages) "
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
                "consistent with the Da_GL definition on the 📐 Equations Reference page."
            )
            st.markdown(
                "| Da_GL range | Regime |\n"
                "|-------------|--------|\n"
                "| Da_GL ≪ 1 | **Reaction-limited** — transport is fast relative to reaction |\n"
                "| Da_GL ≈ 1 | **Transition** — both transport and kinetics matter |\n"
                "| Da_GL ≫ 1 | **Mass-transfer-limited** — gas absorption controls the rate |\n\n"
                "For fast reactions, the **Hatta number** Ha = √($k_1 \\cdot D_A$) / $k_L$ "
                "(with $k_1$ the pseudo-first-order rate constant, s⁻¹) is also useful: "
                "when Ha > 3 the reaction occurs within the liquid film "
                "and enhancement of absorption must be considered."
            )

        if _has_solid:
            st.markdown("##### Liquid–Solid Transport")
            st.markdown(
                "The solid-liquid Damköhler number is defined as:\n\n"
                "$$Da_{SL} = \\frac{1}{k_{SL} \\cdot a_s \\; t_{rxn}} = \\frac{t_{\\text{transfer,SL}}}{t_{rxn}}$$\n\n"
                "where:\n"
                "- $k_{SL}$ = solid-liquid mass-transfer coefficient (m/s)\n"
                "- $a_s$ = specific surface area of the solids **per unit liquid "
                "volume** (m²/m³): $a_s = 6\\,\\phi_s / d_p$ for spheres, where "
                "$\\phi_s$ is the volumetric solids fraction\n"
                "- $1 / (k_{SL} \\cdot a_s)$ = characteristic solid-liquid transport time (s)\n\n"
                "This is the ratio of the **transport time** to the **reaction time**, "
                "consistent with the Da_GL convention above and the Da_SL definition "
                "on the 📐 Equations Reference page."
            )
            st.markdown(
                "| Da_SL range | Regime |\n"
                "|-------------|--------|\n"
                "| Da_SL ≪ 1 | **Reaction-limited** — dissolution is fast relative to reaction |\n"
                "| Da_SL ≈ 1 | **Transition** — both transport and reaction compete |\n"
                "| Da_SL ≫ 1 | **Mass-transfer-limited** — particle dissolution rate controls |\n\n"
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
        "**🌀 Mixing Assessment** or **📈 Reactor Comparison** pages."
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
st.header("3 · Competing Reactions")

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
    index=None,
)

if competing is None:
    st.info("👆 Please select an option above to continue.")
    st.stop()

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
            "- Baldyga & Bourne (1999) describe **two** mesomixing time scales:\n"
            "  - *Inertial-convective disintegration* of the feed plume: "
            r"$\tau_S = A\,(\Lambda_C^{2}/\varepsilon)^{1/3}$, "
            "where $\\Lambda_C$ is the feed-plume (or feed-pipe) scale and "
            "$A\\approx 1{-}2$ (this app uses $\\tau_S = 2\\,(d_{feed}^2/\\varepsilon)^{1/3}$).\n"
            "  - *Turbulent dispersion* of the plume by the mean flow: "
            r"$\tau_D = Q_{feed}/(\bar{u}\,D_t)$, "
            "where $Q_{feed}$ is the feed flow rate, $\\bar{u}$ the local mean "
            "velocity and $D_t$ the turbulent diffusivity.\n"
            "  The larger of the two scales controls the mesomixing rate."
        )

    st.info(
        "👉 **Recommendation:** Use the **� Bourne Protocol** page to "
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
        "use the **� Bourne Protocol** to screen."
    )
    _meso_sensitive = True
else:
    st.success(
        "✅ **No competing reactions identified** — micro- and mesomixing "
        "(selectivity-related) are unlikely to be a factor.  "
        "Macromixing (blend time) may still matter at scale."
    )
    _meso_sensitive = False

# Semi-batch overrides: mesomixing is inherently relevant for fed-batch
# processes because the feed plume must be dispersed into the bulk before
# the reagent reacts (Baldyga & Bourne, 1999; Paul et al., 2004).
if is_semi_batch and not _meso_sensitive:
    _meso_sensitive = True
    st.info(
        "\u26a0\ufe0f **Semi-batch process** \u2014 even without competing reactions, "
        "**mesomixing** (feed-plume dispersion) controls local concentration at the "
        "feed point. Selectivity, local heat release, and supersaturation can all be "
        "affected by feed rate, feed location, and local turbulence."
    )

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 4 – HEAT TRANSFER
# ══════════════════════════════════════════════════════════════════════════
st.header("4 · Heat Transfer Screening")

_has_enthalpy = rxn_delta_H != 0.0
_dT_ad = None  # adiabatic temperature rise (K), computed below when data allow

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
        index=None,
    )

    if dh_action is None:
        st.info("👆 Please select an option above to continue.")
        st.stop()

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

    # ── Adiabatic temperature rise (ΔT_ad) — the proper thermal criterion ──
    # ΔH per mole alone does not determine thermal-runaway risk: a large ΔH at
    # low concentration may be benign, while a modest ΔH at high concentration
    # can be hazardous. ΔT_ad = |ΔH|·C0 / (ρ·Cp) is the temperature the batch
    # would reach with no cooling, and is the basis of the Stoessel
    # criticality classification.
    _rho_cp_auto = None  # volumetric heat capacity, kJ/(m³·K)
    if rxn_solvent and is_known_solvent(rxn_solvent):
        try:
            _sp_heat = get_properties(rxn_solvent, rxn_T, 1.0)
            _rho_cp_auto = _sp_heat["rho_kg_m3"] * _sp_heat["Cp_J_per_kgK"] / 1000.0
        except Exception:
            _rho_cp_auto = None

    with st.expander("🌡️ Adiabatic temperature rise (ΔT_ad) — refine the estimate", expanded=True):
        st.caption(
            "ΔT_ad = |ΔH|·C₀ / (ρ·Cp) is the temperature the batch would reach if "
            "all reaction heat were retained (no cooling). It — not ΔH per mole — "
            "is the proper measure of thermal-runaway potential (Stoessel)."
        )
        _rc1, _rc2 = st.columns(2)
        _rho_cp_in = _rc1.number_input(
            "Volumetric heat capacity ρ·Cp (kJ/m³·K)",
            value=float(_rho_cp_auto) if _rho_cp_auto else 1800.0,
            min_value=100.0, step=50.0, key=_key("rho_cp"),
            help="Auto-filled from the reaction solvent when known "
                 "(≈1800 for typical organics, ≈4180 for water).",
        )
        _C0_heat = _rc2.number_input(
            "Limiting-reagent concentration C₀ (mol/L)",
            value=float(rxn_C0) if rxn_C0 > 0 else 1.0,
            min_value=0.0, step=0.1, key=_key("c0_heat"),
            help="Concentration of the reagent that drives the heat release.",
        )
        if _C0_heat > 0 and _rho_cp_in > 0:
            # |ΔH| [kJ/mol] · C0 [mol/L] · 1000 [L/m³] / (ρCp [kJ/m³·K]) = K
            _dT_ad = abs_dH * _C0_heat * 1000.0 / _rho_cp_in
            if _dT_ad >= 200:
                st.error(
                    f"🔴 **ΔT_ad ≈ {_dT_ad:.0f} K** — very high. Loss of cooling could "
                    "drive a runaway; secondary decomposition (MTSR/MTT) must be assessed."
                )
            elif _dT_ad >= 50:
                st.error(
                    f"🔴 **ΔT_ad ≈ {_dT_ad:.0f} K** — high. Strong cooling and feed-rate "
                    "control are required; quantify Q_gen vs Q_cool at scale."
                )
            elif _dT_ad >= 20:
                st.warning(
                    f"🟡 **ΔT_ad ≈ {_dT_ad:.0f} K** — moderate. Manageable with adequate "
                    "cooling, but verify the heat balance at production scale."
                )
            else:
                st.success(
                    f"🟢 **ΔT_ad ≈ {_dT_ad:.0f} K** — low. Thermal runaway is unlikely, "
                    "though local hot spots at the feed point may still occur."
                )
        else:
            st.caption("Enter C₀ and ρ·Cp to estimate ΔT_ad.")

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

    # Detailed guidance — flag heat sensitivity on either a large molar
    # enthalpy or (more rigorously) a significant adiabatic temperature rise.
    _heat_flag = abs_dH >= 50 or (_dT_ad is not None and _dT_ad >= 50)
    if _heat_flag:
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
            "**🌀 Mixing Assessment** or **📈 Reactor Comparison** page "
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
        "on the 🌀 Mixing Assessment page."
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
        "reactor(s) on the **🌀 Mixing Assessment** page or compare "
        "multiple vessels on the **📈 Reactor Comparison** page."
    )

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
    if _bourne_mechanisms:
        findings.append((
            "Bourne pre-screen",
            "🔴 Mixing sensitivity confirmed",
            "Experimental pre-screen showed yield/impurity changes with mixing "
            "conditions. Identified controlling scale(s): "
            f"{', '.join(_bourne_mechanisms)}.",
        ))
    else:
        findings.append((
            "Bourne pre-screen",
            "🔴 Mixing sensitivity confirmed",
            "Experimental pre-screen showed yield/impurity changes with mixing "
            "conditions. Controlling scale not yet identified — run the full "
            "Bourne Protocol (Tests 1–3) to pinpoint micro-, meso-, or macromixing.",
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
    if is_semi_batch and competing == "No":
        findings.append((
            "Mesomixing (feed-plume)",
            "🟡 Semi-batch — check experimentally",
            "No competing reactions, but feed-plume dispersion controls local concentration. "
            "Vary feed rate and feed location to assess (Bourne Tests 2 & 3).",
        ))
    else:
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
        "No competing reactions; batch process (no feed addition).",
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
    _ht_detail = f"|ΔH| = {abs(rxn_delta_H):.1f} kJ/mol"
    if _dT_ad is not None:
        _ht_detail += f", ΔT_ad ≈ {_dT_ad:.0f} K"
    findings.append((
        "Heat transfer",
        "🔴 Likely sensitive",
        f"{_ht_detail} — run a heat balance to confirm adequate cooling capacity.",
    ))
elif _has_enthalpy and not _heat_sensitive:
    _ht_detail = f"|ΔH| = {abs(rxn_delta_H):.1f} kJ/mol"
    if _dT_ad is not None:
        _ht_detail += f", ΔT_ad ≈ {_dT_ad:.0f} K"
    findings.append((
        "Heat transfer",
        "🟢 Manageable",
        f"{_ht_detail} — modest thermal load, "
        "unlikely to be limiting in most configurations.",
    ))
else:
    findings.append((
        "Heat transfer",
        "⚪ Unknown",
        "No ΔH data available — consider measuring by calorimetry.",
    ))

# Semi-batch
if is_semi_batch:
    findings.append((
        "Semi-batch (fed-batch)",
        "🟡 Feed-point sensitive",
        "Mesomixing (feed-plume dispersion) controls local concentration, "
        "heat release, and supersaturation at the feed point. "
        "Sensitive to feed rate, feed location, and local turbulence.",
    ))

# Overall summary verdict
#
# The theoretical assessment (Steps 1–5) estimates *how likely* the system is
# to be mixing sensitive and *which* mechanism is responsible.  The Bourne
# pre-screen, when actually performed, is independent experimental evidence of
# *whether* a sensitivity exists — but Part 1 alone does not reveal which
# mechanism.  It is therefore combined with AND-style logic rather than being
# counted among the mechanism red/yellow flags.
_mech_findings = [f for f in findings if f[0] != "Bourne pre-screen"]
_n_red = sum(1 for _, s, _ in _mech_findings if "🔴" in s)
_n_yellow = sum(1 for _, s, _ in _mech_findings if "🟡" in s)
_red_mechs = [m for m, s, _ in _mech_findings if "🔴" in s]


def _join_mechs(names: list[str]) -> str:
    """Render mechanism names as a readable lower-case phrase."""
    _clean = [n.split("(")[0].strip().lower() for n in names]
    if len(_clean) == 1:
        return _clean[0]
    if len(_clean) == 2:
        return f"{_clean[0]} and {_clean[1]}"
    return ", ".join(_clean[:-1]) + f", and {_clean[-1]}"


if _bourne_sensitive is True:
    # Experimental evidence confirms a sensitivity is present.
    if _bourne_mechanisms:
        # The Bourne Protocol identified the controlling scale(s) directly —
        # this experimental result takes precedence over the theoretical flags.
        _bm_phrase = _join_mechs(_bourne_mechanisms)
        _verdict = (
            f"Mixing sensitivity confirmed -- the Bourne Protocol identified "
            f"{_bm_phrase} as the controlling scale(s). Focus scale-up on this "
            "mechanism (see recommendations below)."
        )
        st.error(
            f"🔴 **Mixing sensitivity confirmed** — the Bourne Protocol identified "
            f"**{_bm_phrase}** as the controlling scale(s). Focus scale-up efforts "
            "on this mechanism (see recommendations below)."
        )
    elif _red_mechs:
        _verdict = (
            f"Mixing sensitivity confirmed -- the reaction may be {_join_mechs(_red_mechs)} "
            "limited, and the Bourne pre-screen confirms a sensitivity is present. "
            "Detailed characterisation is recommended to identify the controlling mechanism."
        )
        st.error(
            f"🔴 **Mixing sensitivity confirmed** — the reaction may be "
            f"**{_join_mechs(_red_mechs)} limited**, and the Bourne pre-screen confirms "
            "a sensitivity is present. Detailed characterisation is recommended to "
            "identify the controlling mechanism."
        )
    elif _n_yellow >= 1:
        _verdict = (
            "Mixing sensitivity confirmed -- the Bourne pre-screen shows a sensitivity is "
            "present. The theoretical assessment did not flag a specific mechanism as likely, "
            "but some items require verification at scale. Further screening is needed to "
            "pinpoint the controlling scale."
        )
        st.error(
            "🔴 **Mixing sensitivity confirmed** — the Bourne pre-screen shows a sensitivity "
            "is present. The theoretical assessment did not flag a specific mechanism as "
            "*likely*, but some items require verification at scale (see below). Further "
            "screening (full Bourne Protocol) is needed to pinpoint the controlling scale."
        )
    else:
        _verdict = (
            "Mixing sensitivity confirmed -- the Bourne pre-screen shows an experimental "
            "sensitivity, although the theoretical assessment flagged no mechanism. Revisit "
            "the inputs (kinetics, phases, feed strategy) and run the full Bourne Protocol "
            "to identify the controlling scale."
        )
        st.error(
            "🔴 **Mixing sensitivity confirmed** — the Bourne pre-screen shows an experimental "
            "sensitivity even though the theoretical assessment flagged no mechanism. Revisit "
            "the inputs (kinetics, phases, feed strategy) and run the full Bourne Protocol "
            "to identify the controlling scale."
        )
elif _bourne_sensitive is False:
    # Experimental evidence showed no sensitivity at lab scale.
    if _red_mechs:
        _verdict = (
            f"Possible scale-dependent sensitivity -- the Bourne pre-screen showed no "
            f"sensitivity at lab scale, but the assessment flags {_join_mechs(_red_mechs)} "
            "as likely to become limiting at larger scale. Confirm with Damkohler analysis "
            "before scale-up."
        )
        st.warning(
            f"🟡 **Possible scale-dependent sensitivity** — the Bourne pre-screen showed "
            f"no sensitivity at lab scale, but the assessment flags **{_join_mechs(_red_mechs)}** "
            "as likely to become limiting at larger scale. Confirm with Damköhler analysis "
            "before scale-up."
        )
    elif _n_yellow >= 1:
        _verdict = (
            "Low mixing sensitivity risk -- the Bourne pre-screen showed no sensitivity and "
            "no mechanism is flagged as likely, though a few items warrant a check at scale."
        )
        st.success(
            "🟢 **Low mixing sensitivity risk** — the Bourne pre-screen showed no sensitivity "
            "and no mechanism is flagged as likely, though a few items warrant a check at scale."
        )
    else:
        _verdict = (
            "Low mixing sensitivity risk -- the Bourne pre-screen showed no sensitivity and "
            "no mixing mechanism is expected to limit this reaction."
        )
        st.success(
            "🟢 **Low mixing sensitivity risk** — the Bourne pre-screen showed no sensitivity "
            "and no mixing mechanism is expected to limit this reaction."
        )
else:
    # Bourne pre-screen not performed — rely on the theoretical assessment alone.
    if _n_red >= 2:
        _verdict = (
            f"High mixing sensitivity risk -- multiple mechanisms ({_join_mechs(_red_mechs)}) "
            "are likely to limit this reaction at scale. Detailed characterisation is strongly "
            "recommended; a Bourne pre-screen would provide direct experimental confirmation."
        )
        st.error(
            f"🔴 **High mixing sensitivity risk** — multiple mechanisms "
            f"(**{_join_mechs(_red_mechs)}**) are likely to limit this reaction at scale. "
            "Detailed characterisation is strongly recommended. A Bourne pre-screen would "
            "provide direct experimental confirmation."
        )
    elif _n_red == 1:
        _verdict = (
            f"Moderate mixing sensitivity risk -- {_join_mechs(_red_mechs)} is likely to be "
            "sensitive. Targeted investigation is recommended; a Bourne pre-screen would "
            "confirm whether a sensitivity is present experimentally."
        )
        st.warning(
            f"🟡 **Moderate mixing sensitivity risk** — **{_join_mechs(_red_mechs)}** is likely "
            "to be sensitive. Targeted investigation is recommended; a Bourne pre-screen would "
            "confirm whether a sensitivity is present experimentally."
        )
    elif _n_yellow >= 1:
        _verdict = (
            "Low-to-moderate mixing sensitivity risk -- no mechanisms are flagged as likely "
            "sensitive, but some require verification at scale. A Bourne pre-screen would "
            "provide a direct experimental answer."
        )
        st.warning(
            "🟡 **Low-to-moderate mixing sensitivity risk** — no mechanisms are flagged as "
            "likely sensitive, but some require verification at scale. A Bourne pre-screen "
            "would provide a direct experimental answer."
        )
    else:
        _verdict = (
            "Low mixing sensitivity risk -- no mixing mechanisms are expected to limit this "
            "reaction under typical operating conditions."
        )
        st.success(
            "🟢 **Low mixing sensitivity risk** — no mixing mechanisms are expected "
            "to limit this reaction under typical operating conditions."
        )

# Render summary table
_summary_rows = []
for mechanism, status, detail in findings:
    _summary_rows.append({
        "Sensitivity Type": mechanism,
        "Finding": status + " — " + detail,
    })

summary_df = pd.DataFrame(_summary_rows)
st.dataframe(summary_df, width='stretch', hide_index=True)

# st.caption(
#     "🔴 Likely sensitive · 🟡 Potentially sensitive / check at scale · "
#     "🟢 Unlikely / not applicable · ⚪ Not assessed"
# )

# Next steps
st.subheader("Recommended Next Steps")
_steps: list[dict[str, str]] = []

if _bourne_sensitive is None:
    _steps.append({
        "Area": "Bourne pre-screen",
        "Action": "Run Bourne Protocol Part 1 (quick screen) on the 🅱️ Bourne Protocol page "
                  "to confirm whether mixing sensitivity exists experimentally.",
    })

# Targeted scale-up actions for the Bourne-identified controlling scale(s)
if _bourne_sensitive is True and _bourne_mechanisms:
    _mech_actions = {
        "Micromixing": (
            "Hold local energy dissipation ε (P/V) constant on scale-up and keep the "
            "feed point in the high-shear zone near the impeller. Confirm with Da_micro "
            "on the 🌀 Mixing Assessment page.",
        ),
        "Mesomixing": (
            "Control feed-plume dispersion: hold local ε constant at the feed point "
            "(match P/V → a *lower* RPM at larger scale, not constant RPM) and cut the "
            "local feed rate — extend feed time, add feed points, and/or use a smaller "
            "feed-pipe diameter.",
        ),
        "Macromixing": (
            "Reduce bulk blend time on scale-up: use high-efficiency (hydrofoil) or "
            "multiple impellers, optimise baffling, or consider static/in-line mixers.",
        ),
    }
    for _m in _bourne_mechanisms:
        if _m in _mech_actions:
            _steps.append({
                "Area": f"{_m} (Bourne-confirmed)",
                "Action": _mech_actions[_m][0],
            })

if _using_approximate:
    _steps.append({
        "Area": "Kinetics",
        "Action": "Obtain measured kinetics and calorimetry for the actual reaction "
                  "to replace the approximate values used in this assessment.",
    })

if _micro_likely or (t_rxn < 60):
    _steps.append({
        "Area": "Damköhler analysis",
        "Action": "Compute Damköhler numbers (Da_macro, Da_micro) for your specific "
                  "reactor on the 🌀 Mixing Assessment page.",
    })

if _meso_sensitive:
    _steps.append({
        "Area": "Micro/mesomixing",
        "Action": "Run the 🧐 Bourne Protocol to experimentally screen micro- and "
                  "mesomixing sensitivity (vary impeller speed to probe micromixing, "
                  "vary feed time/location to probe mesomixing).",
    })

if _multiphase:
    _steps.append({
        "Area": "Mass transfer",
        "Action": "Assess interphase mass-transfer coefficients (kLa, k_SL) on "
                  "the 🌀 Mixing Assessment or 📈 Reactor Comparison pages.",
    })

if _has_enthalpy and _heat_sensitive:
    _steps.append({
        "Area": "Heat transfer",
        "Action": "Run a full heat balance (🔥 button on 🌀 Mixing Assessment or "
                  "📈 Reactor Comparison) to evaluate Q_gen vs Q_cool.",
    })

if is_semi_batch:
    _steps.append({
        "Area": "Semi-batch / mesomixing",
        "Action": "Run the full 🧐 Bourne Protocol: vary impeller speed (micromixing), "
                  "feed rate / addition time (mesomixing), and feed location "
                  "(meso/macromixing) to identify the dominant mixing scale.",
    })

if not _steps:
    _steps.append({
        "Area": "General",
        "Action": "The reaction appears low risk for mixing sensitivity. "
                  "Standard scale-up practices should be sufficient, but consider "
                  "a quick check on the 🌀 Mixing Assessment page for completeness.",
    })

steps_df = pd.DataFrame(_steps)
st.dataframe(steps_df, width='stretch', hide_index=True)

# ── Generate PDF Report ──────────────────────────────────────────────────
st.divider()
st.header("7 · Export Report")

if st.button("📥 Export PDF Report", type="primary", key="p10_export_pdf"):
    with st.spinner("Generating PDF…"):
        try:
            from utils.report_builder import build_protocol_pdf, report_filename

            # Overall verdict text (_verdict) is computed in the summary section
            # above and already reflects the Bourne pre-screen AND-style logic.

            # Bourne result text
            if _bourne_sensitive is True:
                _bourne_txt = "Mixing sensitivity confirmed experimentally"
                if _bourne_mechanisms:
                    _bourne_txt += (
                        " — controlling scale(s): " + ", ".join(_bourne_mechanisms)
                    )
            elif _bourne_sensitive is False:
                _bourne_txt = "No sensitivity observed"
            else:
                _bourne_txt = "Not performed"

            _p10_snap = {
                "reaction": rxn_name,
                "t_rxn": t_rxn,
                "rxn_delta_H": rxn_delta_H,
                "dT_ad": _dT_ad,
                "phases": phases,
                "findings": findings,
                "next_steps": _steps,
                "bourne_result": _bourne_txt,
                "competing": competing if competing is not None else "Not assessed",
                "overall_verdict": _verdict,
                "using_approximate": _using_approximate,
                "is_semi_batch": is_semi_batch,
            }
            _pdf_bytes = build_protocol_pdf(_p10_snap)
            st.session_state["_p10_pdf_bytes"] = _pdf_bytes
            st.session_state["_p10_pdf_name"] = report_filename(
                "Sensitivity_Protocol", rxn_name
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_p10_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_p10_pdf_bytes"],
        file_name=st.session_state["_p10_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")
