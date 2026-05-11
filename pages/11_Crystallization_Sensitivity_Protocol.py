"""
Page 11 – Crystallization Mixing-Sensitivity Protocol
======================================================
An interactive decision-tree that assesses which mixing mechanisms
(micro-, meso-, macromixing, mass transfer, heat transfer) are likely
to limit a pharmaceutical crystallization at scale.

The protocol mirrors the structure of Page 10 (Mixing Sensitivity
Protocol for chemical reactions) but replaces reaction-kinetic
benchmarks with crystallization-specific characteristic times:

  0. Pre-screening  – Bourne Protocol quick screen
  1. Crystallization parameters – select or enter system data
  2. Crystallization type & supersaturation assessment
  3. Nucleation / growth competition – micro- & mesomixing
  4. Heat transfer screening
  5. Damköhler comparison – mixing times vs crystallization times
  6. Summary & recommendations

References
----------
- Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). Handbook of
  Industrial Crystallization, 3rd ed. Cambridge Univ. Press.
  (esp. Chapters 8, 10, 12, 13)
- Baldyga, J. & Bourne, J.R. (1999). Turbulent Mixing and
  Chemical Reactions. Wiley.
- Green, D.A. (2019). Ch. 10 "Crystallizer Mixing" in Myerson.
- Black, S.N. (2019). Ch. 13 "Crystallization in the Pharmaceutical
  Industry" in Myerson.
"""

import streamlit as st
import pandas as pd
import numpy as np
import pathlib

from utils.calculations import (
    reynolds_number,
    impeller_power,
    power_per_volume,
    tip_speed,
    blend_time_turbulent,
    micromixing_time_engulfment,
    kolmogorov_length,
)
from utils.solvent_properties import SOLVENT_DB, get_properties, is_known_solvent

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

from utils.data_helpers import load_db, safe_float as _safe, all_fluid_names, safe_iloc

crystallizations = load_db("crystallization_db", "crystallizations.csv")
reactors = load_db("reactor_db", "reactors.csv")
custom_fluids = load_db("fluid_db", "fluids.csv")

_all_fluid_names = all_fluid_names(custom_fluids)

# ── session-state key helpers ────────────────────────────────────────────
_PFX = "_csp_"
if "_csp_gen" not in st.session_state:
    st.session_state["_csp_gen"] = 0


def _key(name: str) -> str:
    return f"{_PFX}{st.session_state['_csp_gen']}_{name}"


def _get(name: str, default=None):
    return st.session_state.get(_key(name), default)


def _reset_csp():
    old_gen = st.session_state.get("_csp_gen", 0)
    for k in list(st.session_state.keys()):
        if k == "_csp_gen" or k == _key("restart"):
            continue
        if k.startswith(_PFX):
            del st.session_state[k]
    st.session_state["_csp_gen"] = old_gen + 1


# ═══════════════════════════════════════════════════════════════════════
# PAGE HEADER
# ═══════════════════════════════════════════════════════════════════════
st.title("💎 Crystallization Sensitivity Protocol")
st.caption("Interactive decision tree for crystallization mixing-sensitivity assessment.")

# ── Work-in-progress gate ────────────────────────────────────────────
st.warning(
    "🚧 **Work in Progress** — This page is under active development "
    "and not yet ready for use. Check back soon!"
)
st.stop()

# ── Visual overview of the decision tree ─────────────────────────────────
_IMG_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "general"
_CSP_IMG = _IMG_DIR / "crystallization_sensitivity_protocol.png"

with st.expander("📋 Protocol overview", expanded=False):
    if _CSP_IMG.exists():
        st.image(str(_CSP_IMG), caption="Crystallization Sensitivity Protocol – Decision Tree")
    else:
        st.info("Protocol overview image not found. Generate it from the 🛠️ Admin page.")

st.button("🔄 Restart protocol", key=_key("restart"), on_click=_reset_csp)

with st.expander("ℹ️ About this protocol", expanded=False):
    st.markdown(
        "Adapts the Mixing Sensitivity Protocol (Page 10) for crystallization. "
        "The single $t_{rxn}$ is replaced by **multiple crystallization time scales**:\n\n"
        "- **Induction time** $t_{ind}$ — delay to first nuclei\n"
        "- **Nucleation time** $t_N$ — time to consume supersaturation via nucleation\n"
        "- **Growth time** $t_G$ — time for crystal growth to deplete supersaturation\n\n"
        "These are compared against mixing time scales "
        "($t_E$, $t_{meso}$, $\\theta_{95}$) via Damköhler numbers.\n\n"
        "**Refs:** Myerson *et al.* (2019); Baldyga & Bourne (1999); "
        "Green (2019) Ch. 10; Black (2019) Ch. 13."
    )

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 0 – BOURNE PROTOCOL PRE-SCREENING
# ═══════════════════════════════════════════════════════════════════════
st.header("0 · Bourne Protocol Pre-Screening")

with st.expander("ℹ️ Background — Bourne pre-screening for crystallization", expanded=False):
    st.markdown(
        "Relevant **response metrics** for crystallization Bourne screening:\n\n"
        "- Mean particle size / chord length (FBRM)\n"
        "- Particle count · Polymorph form (XRPD, Raman)\n"
        "- CSD span · Filterability\n\n"
        "Run at 2+ impeller speeds / feed rates and check for changes. "
        "See **🅱️  Bourne Protocol** page for the full procedure."
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
    st.info("Select an option above to continue.")
    st.stop()

if bourne_screen == "Already done – results show mixing sensitivity":
    st.warning(
        "**Mixing sensitivity confirmed.** Remaining steps identify which mechanisms are responsible."
    )
    _bourne_sensitive = True
elif bourne_screen == "Already done – no mixing sensitivity observed":
    st.success(
        "✅ **No sensitivity observed.** Proceed to confirm and identify latent risks at larger scale."
    )
    _bourne_sensitive = False
else:
    st.info(
        "💡 Consider running **Bourne Part 1** first (🅱️  Bourne Protocol page, "
        "response = particle size or polymorph). Proceeding with theoretical assessment."
    )
    _bourne_sensitive = None

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 1 – CRYSTALLIZATION PARAMETERS
# ═══════════════════════════════════════════════════════════════════════
st.header("1 · Crystallization Parameters")

with st.expander("ℹ️ Background — why these parameters matter", expanded=False):
    st.markdown(
        "Key parameters and their roles:\n\n"
        "- **σ** — supersaturation driving force\n"
        "- **$t_{ind}$** — induction time (delay to nucleation)\n"
        "- **$k_g$, $g$** — growth rate constant and order\n"
        "- **MSZW** — max sustainable σ before spontaneous nucleation\n\n"
        "Measure experimentally (turbidity, FBRM/PVM) or estimate from literature."
    )

# ── select from database or enter manually ────────────────────────────
cryst_choice = st.radio(
    "How would you like to provide crystallization data?",
    [
        "Select from Crystallization Database",
        "Enter parameters manually",
    ],
    key=_key("cryst_choice"),
    index=None,
)

if cryst_choice is None:
    st.info("Select an option above to continue.")
    st.stop()

if cryst_choice == "Select from Crystallization Database":
    if crystallizations.empty:
        st.error(
            "Crystallization Database is empty. Add entries to "
            "`data/crystallizations.csv` or enter manually."
        )
        st.stop()

    _cx_names = crystallizations["crystallization_name"].tolist()
    _cx_idx = (
        _cx_names.index(st.session_state["_sel_csp_cx"])
        if st.session_state.get("_sel_csp_cx") in _cx_names
        else 0
    )
    cx_name = st.selectbox(
        "Select crystallization system",
        _cx_names,
        index=_cx_idx,
        key=_key("cx_sel"),
    )
    st.session_state["_sel_csp_cx"] = cx_name
    cx = safe_iloc(crystallizations, "crystallization_name", cx_name, "Crystallization")

    cryst_type = str(cx.get("type", "Cooling"))
    solute_name = str(cx.get("solute", ""))
    solvent_name = str(cx.get("solvent", ""))
    anti_solvent_name = str(cx.get("anti_solvent", ""))
    c_sat_ref = _safe(cx.get("c_sat_ref_g_L"), 0.0)
    T_ref = _safe(cx.get("T_ref_C"), 25.0)
    dc_dT = _safe(cx.get("dc_dT_g_L_K"), 0.0)
    sigma_target = _safe(cx.get("sigma_target"), 0.3)
    MSZW = _safe(cx.get("MSZW_C"), 10.0)
    t_ind = _safe(cx.get("t_ind_s"), 0.0)
    k_g = _safe(cx.get("k_g_m_s"), 0.0)
    g_order = _safe(cx.get("g_order"), 1.5)
    delta_H_cryst = _safe(cx.get("delta_H_cryst_kJ_mol"), 0.0)
    has_polymorphs = str(cx.get("polymorphs", "No")).lower() not in ("no", "nan", "")
    polymorph_info = str(cx.get("polymorphs", "No"))
    is_seeded = str(cx.get("seeded", "No")).lower().startswith("y")
    target_L = _safe(cx.get("target_L_um"), 100.0)
    cx_notes = str(cx.get("notes", ""))

    with st.expander("Selected crystallization details", expanded=False):
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Type", cryst_type)
        c2.metric("Solute", solute_name)
        c3.metric("Solvent", solvent_name)
        c4.metric("Anti-solvent", anti_solvent_name if anti_solvent_name != "nan" else "—")
        c5, c6, c7, c8 = st.columns(4)
        c5.metric("c_sat (g/L)", f"{c_sat_ref:.1f}")
        c6.metric("σ target", f"{sigma_target:.2f}")
        c7.metric("MSZW (°C)", f"{MSZW:.1f}")
        c8.metric("t_ind (s)", f"{t_ind:.4g}")
        c9, c10, c11, c12 = st.columns(4)
        c9.metric("k_g (m/s)", f"{k_g:.2e}")
        c10.metric("g (order)", f"{g_order:.1f}")
        c11.metric("ΔH_cryst (kJ/mol)", f"{delta_H_cryst:.1f}")
        c12.metric("Polymorphs", polymorph_info)
        if cx_notes and cx_notes != "nan":
            st.caption(cx_notes)

else:
    # Manual entry
    st.subheader("Manual Parameter Entry")
    cryst_type = st.selectbox(
        "Crystallization type",
        ["Cooling", "Anti-solvent", "Reactive", "Evaporative", "pH-shift"],
        key=_key("m_type"),
    )
    solute_name = st.text_input("Solute name", value="API", key=_key("m_solute"))

    col_m1, col_m2 = st.columns(2)
    with col_m1:
        solvent_name = st.text_input("Solvent", value="Ethanol", key=_key("m_solvent"))
    with col_m2:
        anti_solvent_name = st.text_input(
            "Anti-solvent (if applicable)", value="", key=_key("m_antisolv")
        )

    col_m3, col_m4, col_m5 = st.columns(3)
    with col_m3:
        c_sat_ref = st.number_input(
            "c_sat at T_ref (g/L)", value=100.0, format="%.1f", key=_key("m_csat")
        )
    with col_m4:
        T_ref = st.number_input("T_ref (°C)", value=25.0, format="%.1f", key=_key("m_Tref"))
    with col_m5:
        dc_dT = st.number_input(
            "dc/dT (g/L/K)", value=2.0, format="%.2f", key=_key("m_dcdT"),
            help="Slope of solubility curve near T_ref.",
        )

    col_m6, col_m7, col_m8 = st.columns(3)
    with col_m6:
        sigma_target = st.number_input(
            "Target σ (dimensionless)", value=0.3, format="%.3f", key=_key("m_sigma"),
            help="σ = (c − c_sat) / c_sat",
        )
    with col_m7:
        MSZW = st.number_input(
            "Metastable zone width (°C)", value=10.0, format="%.1f",
            key=_key("m_mszw"),
        )
    with col_m8:
        t_ind = st.number_input(
            "Induction time t_ind (s)", value=30.0, format="%.4g",
            key=_key("m_tind"),
            help="Time to first detected nucleation at target σ.",
        )

    col_m9, col_m10, col_m11 = st.columns(3)
    with col_m9:
        k_g = st.number_input(
            "Growth rate constant k_g (m/s)",
            value=1.0e-8, format="%.2e", key=_key("m_kg"),
        )
    with col_m10:
        g_order = st.number_input(
            "Growth order g", value=1.5, format="%.1f", key=_key("m_g"),
        )
    with col_m11:
        delta_H_cryst = st.number_input(
            "ΔH_cryst (kJ/mol)", value=-25.0, format="%.1f", key=_key("m_dH"),
        )

    col_m12, col_m13, col_m14 = st.columns(3)
    with col_m12:
        _poly_choice = st.selectbox(
            "Known polymorphs?", ["No", "Yes"], key=_key("m_poly")
        )
        has_polymorphs = _poly_choice == "Yes"
        polymorph_info = _poly_choice
    with col_m13:
        _seed_choice = st.selectbox("Seeded?", ["Yes", "No"], key=_key("m_seed"))
        is_seeded = _seed_choice == "Yes"
    with col_m14:
        target_L = st.number_input(
            "Target crystal size (µm)", value=100.0, format="%.0f", key=_key("m_L"),
        )

# ── validate minimum data ─────────────────────────────────────────────
if t_ind <= 0 and k_g <= 0:
    st.error(
        "At least one of **t_ind** or **k_g** must be > 0 to proceed.  "
        "Please check the crystallization data."
    )
    st.stop()

# ── derived characteristic times ──────────────────────────────────────
# Growth time estimate: t_G ≈ target_L / (k_g · σ^g)
# This is the time for a crystal to grow to the target size.
if k_g > 0 and sigma_target > 0:
    G_linear = k_g * sigma_target**g_order  # m/s linear growth rate
    t_G = (target_L * 1e-6) / G_linear if G_linear > 0 else np.inf
else:
    G_linear = 0.0
    t_G = np.inf

st.success(
    f"✅ Crystallization parameters loaded.\n\n"
    f"- **Induction time:** $t_{{ind}}$ = {t_ind:.4g} s\n"
    f"- **Growth time to {target_L:.0f} µm:** $t_G$ = {t_G:.4g} s  "
    f"(G = {G_linear:.3e} m/s at σ = {sigma_target:.3f})"
)

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 2 – CRYSTALLIZATION TYPE ASSESSMENT
# ═══════════════════════════════════════════════════════════════════════
st.header("2 · Crystallization Type & Supersaturation Profile")

with st.expander("ℹ️ Background — how crystallization type affects mixing sensitivity", expanded=False):
    st.markdown(
        "| Type | σ generation | Dominant mixing concern |\n"
        "|------|-------------|------------------------|\n"
        "| **Cooling** | ΔT | Macromixing — T-field uniformity |\n"
        "| **Anti-solvent** | Poor-solvent addition | Meso/micromixing — feed plume σ spike |\n"
        "| **Reactive** | Reaction → insoluble product | Micromixing — local σ |\n"
        "| **Evaporative** | Solvent removal | Macromixing — boiling zone σ gradient |\n"
        "| **pH-shift** | Acid/base addition | Meso/micromixing — local pH plume |\n\n"
        "*Anti-solvent and reactive are most mixing-sensitive — σ generated locally at the "
        "feed point before mixing disperses it.* — Myerson (2019), Ch. 10"
    )

_TYPE_RISK = {
    "Cooling": ("low", "🟢"),
    "Anti-solvent": ("high", "🔴"),
    "Reactive": ("high", "🔴"),
    "Evaporative": ("moderate", "🟡"),
    "pH-shift": ("high", "🔴"),
}
_risk_level, _risk_icon = _TYPE_RISK.get(cryst_type, ("moderate", "🟡"))

if _risk_level == "high":
    st.warning(
        f"{_risk_icon} **{cryst_type} crystallization** — **high mixing sensitivity**. "
        f"σ generated locally at feed point; micro/mesomixing controls the σ profile."
    )
elif _risk_level == "moderate":
    st.info(
        f"{_risk_icon} **{cryst_type} crystallization** — **moderate sensitivity**. "
        f"σ gradients can develop between generation zone and bulk."
    )
else:
    st.success(
        f"{_risk_icon} **{cryst_type} crystallization** — **lower sensitivity**. "
        f"Macromixing (T uniformity) is the primary concern."
    )

# ── Supersaturation & MSZW assessment ─────────────────────────────────
st.subheader("Supersaturation & Metastable Zone")

if sigma_target > 0 and MSZW > 0:
    # Ratio of target σ to approximate max σ from MSZW
    sigma_max_approx = dc_dT * MSZW / c_sat_ref if c_sat_ref > 0 else sigma_target * 3
    sigma_ratio = sigma_target / sigma_max_approx if sigma_max_approx > 0 else 0

    st.markdown(
        f"- Target supersaturation: **σ = {sigma_target:.3f}**\n"
        f"- MSZW ≈ {MSZW:.1f} °C → estimated max σ ≈ {sigma_max_approx:.3f}\n"
        f"- **σ / σ_max ≈ {sigma_ratio:.2f}**"
    )

    if sigma_ratio > 0.7:
        st.warning(
            "⚠️ Near metastable limit — local σ spikes could trigger "
            "**uncontrolled nucleation** or polymorph change."
        )
        _sigma_risk = True
    elif sigma_ratio > 0.4:
        st.info(
            "Mid-range of metastable zone. Local σ spikes from feed could exceed limit."
        )
        _sigma_risk = True
    else:
        st.success("Well within metastable zone — σ fluctuations unlikely to cause problems.")
        _sigma_risk = False
else:
    _sigma_risk = cryst_type in ("Anti-solvent", "Reactive", "pH-shift")
    st.info("Insufficient data to assess σ / σ_max ratio.")

_feed_sensitive = cryst_type in ("Anti-solvent", "Reactive", "pH-shift")

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 3 – NUCLEATION / GROWTH COMPETITION (micro- & mesomixing)
# ═══════════════════════════════════════════════════════════════════════
st.header("3 · Nucleation / Growth Competition — Micro- & Mesomixing")

with st.expander("ℹ️ Background — the crystallization analogy to competing reactions", expanded=False):
    st.markdown(
        "Crystallization analogy to competing reactions: **nucleation** ($B \\propto \\sigma^n$, $n$ = 2–10) "
        "vs **growth** ($G \\propto \\sigma^g$, $g$ = 1–2). Since $n \\gg g$, "
        "poor mixing → local σ spikes → excess nucleation → fines, broad CSD, wrong polymorph.\n\n"
        "### Micromixing\n"
        "$$t_E = 17.3 \\sqrt{\\nu / \\varepsilon}$$\n\n"
        "When $t_E > t_{ind}$, nucleation occurs in the unmixed feed plume → "
        "**micromixing-controlled** CSD.\n\n"
        "### Mesomixing\n"
        "Feed-plume dispersion time sets local σ. Higher feed rate → "
        "longer mesomixing → higher local σ → more fines."
    )

# ── Feed-sensitive or not ─────────────────────────────────────────────
if _feed_sensitive:
    st.warning(
        f"**{cryst_type}** — feed addition; micro/mesomixing controls local σ and nucleation/growth balance."
    )

    st.markdown(
        "**Risk factors:** $t_{ind}$ < 10 s · σ > 0.5 · polymorphism · oiling-out tendency"
    )

    _micro_meso_sensitive = True
else:
    # Cooling or evaporative — still evaluate
    if t_ind < 10:
        st.warning(
            f"Short $t_{{ind}}$ = {t_ind:.4g} s — even cooling crystallization can be "
            "sensitive to local T/σ non-uniformity from imperfect macromixing."
        )
        _micro_meso_sensitive = True
    elif sigma_target > 0.5:
        st.info("High σ — nucleation kinetics may be influenced by local mixing conditions.")
        _micro_meso_sensitive = True
    else:
        st.success(
            "✅ **Low risk** — moderate σ and reasonable $t_{ind}$; typically not micromixing-controlled."
        )
        _micro_meso_sensitive = False

# ── Polymorphism ──────────────────────────────────────────────────────
if has_polymorphs:
    st.warning(
        f"⚠️ **Polymorphism** ({polymorph_info}) — local σ spikes can nucleate a "
        "**metastable form** (Ostwald's rule). May or may not convert to stable form."
    )
    _polymorph_risk = True
else:
    _polymorph_risk = False

# ── Seeding ───────────────────────────────────────────────────────────
if is_seeded:
    st.info(
        "✅ **Seeded** — seed surface area consumes σ via growth, reducing "
        "(not eliminating) mixing sensitivity. *Rule of thumb:* 2–5× less sensitive than unseeded."
    )
else:
    st.info("⚠️ **Unseeded** — all nucleation is primary/secondary; **more sensitive** to local σ.")

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 4 – HEAT TRANSFER
# ═══════════════════════════════════════════════════════════════════════
st.header("4 · Heat Transfer Screening")

_has_enthalpy = delta_H_cryst != 0.0
_heat_sensitive = False

with st.expander("ℹ️ Background — heat effects in crystallization", expanded=False):
    st.markdown(
        "Crystallization is typically **exothermic**.\n\n"
        "- **Cooling:** Poor jacket HT → cold spots → local σ spikes → encrustation\n"
        "- **Anti-solvent:** Different stream temperatures → coupled mass/heat transfer\n"
        "- **Reactive:** Reaction + crystallization enthalpy can be very significant\n\n"
        "At scale, lower SA/V makes heat removal harder → slower addition or lower T."
    )

if _has_enthalpy:
    abs_dH = abs(delta_H_cryst)
    st.markdown(
        f"**ΔH_cryst = {delta_H_cryst:.1f} kJ/mol** "
        f"({'exothermic' if delta_H_cryst < 0 else 'endothermic'})"
    )
    if abs_dH >= 40:
        st.warning(
            f"🔴 **High |ΔH|** = {abs_dH:.1f} kJ/mol — heat removal may limit rates at scale."
        )
        _heat_sensitive = True
    elif abs_dH >= 20:
        st.info(
            f"🟡 **Moderate |ΔH|** = {abs_dH:.1f} kJ/mol — unlikely to limit, but check at scale."
        )
    else:
        st.success(f"🟢 **Low |ΔH|** = {abs_dH:.1f} kJ/mol.")

    if cryst_type == "Cooling":
        st.markdown(
            "Also consider: jacket ΔT limits · wall encrustation · required cooling rate vs UA."
        )
else:
    st.info("No ΔH_cryst data. Measure by calorimetry (RC1, µRC) or estimate via van 't Hoff.")

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 5 – DAMKÖHLER COMPARISON
# ═══════════════════════════════════════════════════════════════════════
st.header("5 · Mixing Time vs Crystallization Time (Damköhler Analysis)")

with st.expander("ℹ️ Background — Damköhler numbers for crystallization", expanded=False):
    st.markdown(
        "$$Da = \\frac{t_{mix}}{t_{cryst}}$$\n\n"
        "| Da | Interpretation |\n"
        "|----------|----------------|\n"
        "| < 0.01 | **Not sensitive** |\n"
        "| 0.01–0.1 | **Likely not sensitive** — monitor at scale |\n"
        "| 0.1–1 | **Potentially sensitive** |\n"
        "| 1–10 | **Likely sensitive** — mixing limits process |\n"
        "| > 10 | **Highly sensitive** — mixing controls outcome |\n\n"
        "Comparisons: $Da_{micro} = t_E / t_{ind}$ · "
        "$Da_{macro} = \\theta_{95} / t_G$ · "
        "$Da_{meso}$ (qualitative, feed-sensitive processes)"
    )

st.subheader("Heuristic Assessment (no reactor selected)")

st.markdown(
    f"**Your times:** $t_{{ind}}$ = {t_ind:.4g} s · $t_G$ = {t_G:.4g} s\n\n"
    "Heuristic Da based on typical mixing times at each scale:"
)

# Typical mixing time ranges by scale
_SCALE_DATA = [
    ("Lab (1–5 L)", 0.001, 0.05, 2, 5),
    ("Pilot (50–200 L)", 0.01, 0.2, 8, 20),
    ("Production (1–10 m³)", 0.05, 0.5, 25, 60),
    ("Large (10–50 m³)", 0.1, 1.0, 50, 120),
]

_da_rows = []
for scale, te_lo, te_hi, tb_lo, tb_hi in _SCALE_DATA:
    da_micro_lo = te_lo / t_ind if t_ind > 0 else 0
    da_micro_hi = te_hi / t_ind if t_ind > 0 else 0
    da_macro_lo = tb_lo / t_G if t_G > 0 and t_G != np.inf else 0
    da_macro_hi = tb_hi / t_G if t_G > 0 and t_G != np.inf else 0

    def _flag(da):
        if da > 1:
            return "🔴"
        elif da > 0.1:
            return "🟡"
        elif da > 0.01:
            return "🔵"
        else:
            return "🟢"

    _da_rows.append({
        "Scale": scale,
        "t_E range (s)": f"{te_lo:.3f} – {te_hi:.3f}",
        "Da_micro (t_E/t_ind)": f"{da_micro_lo:.3f} – {da_micro_hi:.3f}  {_flag(da_micro_hi)}",
        "θ₉₅ range (s)": f"{tb_lo} – {tb_hi}",
        "Da_macro (θ₉₅/t_G)": f"{da_macro_lo:.3f} – {da_macro_hi:.3f}  {_flag(da_macro_hi)}",
    })

da_df = pd.DataFrame(_da_rows)
st.dataframe(da_df, use_container_width=True, hide_index=True)
st.caption(
    "🟢 Da < 0.01 (not sensitive)  •  🔵 0.01–0.1 (monitor)  •  "
    "🟡 0.1–1 (potentially sensitive)  •  🔴 Da > 1 (likely sensitive)"
)

# Overall micro assessment
_micro_likely = False
if t_ind > 0 and t_ind < 1.0:
    st.error(
        f"🔴 **Very fast nucleation** ($t_{{ind}}$ = {t_ind:.4g} s) — "
        "**micromixing-sensitive**; local ε at feed point dominates CSD."
    )
    _micro_likely = True
elif t_ind > 0 and t_ind < 10.0:
    st.warning(
        f"🟡 **Fast nucleation** ($t_{{ind}}$ = {t_ind:.4g} s) — micromixing relevant at pilot/production scale."
    )
    _micro_likely = True
elif t_ind > 0 and t_ind < 60.0:
    st.info(
        f"🔵 **Moderate nucleation** ($t_{{ind}}$ = {t_ind:.4g} s) — "
        "micromixing less likely; macromixing may matter at production scale."
    )
else:
    st.success(
        f"🟢 **Slow nucleation** ($t_{{ind}}$ = {t_ind:.4g} s) — mixing unlikely to limit in well-agitated vessels."
    )

# Macro assessment
_macro_likely = False
if t_G < np.inf and t_G < 120:
    st.info(
        f"$t_G$ = {t_G:.4g} s — within blend-time range of larger vessels; "
        f"macromixing may affect σ uniformity and growth rate distribution."
    )
    if t_G < 30:
        _macro_likely = True

# ── Optional: compute Da for a specific reactor ──────────────────────
st.subheader("(Optional) Compute Da for a Specific Reactor")

if reactors.empty:
    st.info("No reactors in database. Add via 🧪 Reactor Database page.")
else:
    _reactor_list = reactors["reactor_name"].tolist()
    reactor_name = st.selectbox(
        "Select reactor (optional)",
        ["— Skip —"] + _reactor_list,
        key=_key("reactor_sel"),
    )

    if reactor_name != "— Skip —":
        r = safe_iloc(reactors, "reactor_name", reactor_name, "Reactor")
        D_tank = _safe(r["D_tank_m"], 0.10)
        D_imp = _safe(r["D_imp_m"], 0.05)
        Np_val = _safe(r.get("Np"), 1.0)
        Nq_val = _safe(r.get("Nq"), 0.79)
        N_rpm_min = _safe(r.get("N_rpm_min"), 0.0) or None
        N_rpm_max = _safe(r.get("N_rpm_max"), 0.0) or None

        V_L_min = _safe(r.get("V_L_min"), 0.0) or _safe(r["V_L"], 0.0)
        V_L_max = _safe(r.get("V_L_max"), 0.0) or _safe(r["V_L"], 0.0)
        V_L_avg = (V_L_min + V_L_max) / 2.0

        col_rv, col_rn = st.columns(2)
        with col_rv:
            V_L = st.number_input(
                "Working volume (L)", min_value=V_L_min, value=V_L_avg,
                step=1.0, format="%.1f", key=_key("r_vol"),
            )
        with col_rn:
            if N_rpm_min is not None and N_rpm_max is not None:
                N_rpm = st.number_input(
                    "Impeller speed (RPM)",
                    min_value=float(N_rpm_min),
                    max_value=float(N_rpm_max),
                    value=(N_rpm_min + N_rpm_max) / 2.0,
                    step=1.0, format="%.0f", key=_key("r_rpm"),
                )
            else:
                N_rpm = st.number_input(
                    "Impeller speed (RPM)",
                    min_value=1.0, value=float(r.get("N_rps", 3.0)) * 60,
                    step=1.0, format="%.0f", key=_key("r_rpm"),
                )

        # Get fluid properties for kinematic viscosity
        fluid_sel = st.selectbox(
            "Fluid for viscosity", _all_fluid_names, key=_key("r_fluid")
        )
        _is_solv = is_known_solvent(fluid_sel)
        if _is_solv:
            _fp = get_properties(fluid_sel, T_ref, 1.0)
            rho = _fp["rho_kg_m3"]
            mu = _fp["mu_Pa_s"]
        else:
            _cf = safe_iloc(custom_fluids, "fluid_name", fluid_sel, "Custom fluid")
            rho = float(_cf["rho_kg_m3"])
            mu = float(_cf["mu_Pa_s"])

        nu = mu / rho
        N_rps = N_rpm / 60.0
        V_m3 = V_L / 1000.0

        Re = reynolds_number(N_rps, D_imp, rho, mu)
        P = impeller_power(Np_val, rho, N_rps, D_imp)
        eps = power_per_volume(P, V_m3)
        eps_kg = eps / rho
        t_E = micromixing_time_engulfment(eps_kg, nu)
        theta_95 = blend_time_turbulent(Nq_val, V_m3, D_imp, N_rps)
        eta_k = kolmogorov_length(nu, eps_kg)

        Da_micro = t_E / t_ind if t_ind > 0 else 0
        Da_macro = theta_95 / t_G if t_G > 0 and t_G != np.inf else 0

        st.markdown(f"**Reactor: {reactor_name}** at {N_rpm:.0f} RPM, {V_L:.0f} L")

        rc1, rc2, rc3, rc4 = st.columns(4)
        rc1.metric("Re", f"{Re:,.0f}")
        rc2.metric("P/V (W/kg)", f"{eps_kg:.4g}")
        rc3.metric("t_E (s)", f"{t_E:.4g}")
        rc4.metric("θ₉₅ (s)", f"{theta_95:.2f}")

        rc5, rc6, rc7, rc8 = st.columns(4)
        rc5.metric("η Kolmogorov (µm)", f"{eta_k * 1e6:.1f}")
        rc6.metric("Da_micro", f"{Da_micro:.4g}")
        rc7.metric("Da_macro", f"{Da_macro:.4g}")
        rc8.metric("Tip speed (m/s)", f"{tip_speed(N_rps, D_imp):.2f}")

        def _da_status(da, name):
            if da > 1:
                return f"🔴 **{name} = {da:.3g}** — mixing-sensitive"
            elif da > 0.1:
                return f"🟡 **{name} = {da:.3g}** — potentially sensitive"
            elif da > 0.01:
                return f"🔵 **{name} = {da:.3g}** — monitor at scale"
            else:
                return f"🟢 **{name} = {da:.3g}** — not sensitive"

        st.markdown(_da_status(Da_micro, "Da_micro (t_E / t_ind)"))
        st.markdown(_da_status(Da_macro, "Da_macro (θ₉₅ / t_G)"))

        if _feed_sensitive:
            st.info(
                "📌 **Feed-sensitive:** also consider feed location (near impeller vs surface), "
                "feed rate, and number of feed points."
            )

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 6 – SUMMARY & RECOMMENDATIONS
# ═══════════════════════════════════════════════════════════════════════
st.header("6 · Summary & Recommendations")

findings: list[tuple[str, str, str]] = []

# Bourne pre-screening
if _bourne_sensitive is True:
    findings.append((
        "Bourne pre-screen",
        "🔴 Mixing sensitivity confirmed",
        "CSD / polymorph changes observed with mixing conditions.",
    ))
elif _bourne_sensitive is False:
    findings.append((
        "Bourne pre-screen",
        "🟢 No sensitivity observed",
        "No sensitivity at lab scale.",
    ))
else:
    findings.append((
        "Bourne pre-screen",
        "⚪ Not performed",
        "Consider Bourne Part 1 (response = particle size or polymorph).",
    ))

# Crystallization type
findings.append((
    f"Crystallization type ({cryst_type})",
    f"{_risk_icon} {_risk_level.title()} inherent risk",
    {
        "high": "Feed-based supersaturation generation → micro/mesomixing directly controls local σ.",
        "moderate": "σ gradients can develop between generation zone and bulk.",
        "low": "Macromixing (temperature uniformity) is the primary concern.",
    }.get(_risk_level, ""),
))

# Supersaturation
if _sigma_risk:
    findings.append((
        "Supersaturation profile",
        "🟡 Elevated risk",
        f"σ = {sigma_target:.3f}; {'near' if sigma_target / (dc_dT * MSZW / c_sat_ref if c_sat_ref > 0 and MSZW > 0 else 1) > 0.5 else 'mid-range of'} metastable limit.",
    ))
else:
    findings.append((
        "Supersaturation profile",
        "🟢 Low risk",
        f"σ = {sigma_target:.3f}; within metastable zone.",
    ))

# Micromixing
if _micro_likely:
    findings.append((
        "Micromixing (nucleation control)",
        "🔴 Likely sensitive",
        f"t_ind = {t_ind:.4g} s — fast relative to typical micromixing times.",
    ))
else:
    findings.append((
        "Micromixing (nucleation control)",
        "🟢 Unlikely",
        f"t_ind = {t_ind:.4g} s — slow; micromixing not rate-limiting.",
    ))

# Micro/mesomixing (feed-related)
if _micro_meso_sensitive and _feed_sensitive:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🔴 Likely sensitive",
        "Feed generates local σ spike; feed-point mixing controls CSD.",
    ))
elif _micro_meso_sensitive:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🟡 Potentially sensitive",
        "Short t_ind or high σ suggests sensitivity even without direct feed.",
    ))
else:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🟢 Not a primary factor",
        "No direct feed or slow nucleation.",
    ))

# Macromixing
if _macro_likely:
    findings.append((
        "Macromixing (blend time)",
        "🟡 Check at scale",
        f"t_G = {t_G:.4g} s — within blend-time range at larger scale.",
    ))
elif t_G < 300 and t_G != np.inf:
    findings.append((
        "Macromixing (blend time)",
        "🔵 Monitor",
        f"t_G = {t_G:.4g} s — could matter in very large vessels.",
    ))
else:
    findings.append((
        "Macromixing (blend time)",
        "🟢 Unlikely",
        f"t_G = {t_G:.4g} s — much longer than typical blend times.",
    ))

# Polymorphism
if _polymorph_risk:
    findings.append((
        "Polymorphism",
        "⚠️ Risk present",
        f"Known polymorphs ({polymorph_info}); σ spikes → metastable form risk.",
    ))
else:
    findings.append((
        "Polymorphism",
        "🟢 Low risk",
        "No known polymorphs.",
    ))

# Heat transfer
if _heat_sensitive:
    findings.append((
        "Heat transfer",
        "🔴 Likely limiting",
        f"|ΔH| = {abs(delta_H_cryst):.1f} kJ/mol — check cooling capacity at scale.",
    ))
elif _has_enthalpy:
    findings.append((
        "Heat transfer",
        "🟢 Manageable",
        f"|ΔH| = {abs(delta_H_cryst):.1f} kJ/mol — modest.",
    ))
else:
    findings.append((
        "Heat transfer",
        "⚪ Unknown",
        "No ΔH data — consider measuring.",
    ))

# Seeding
if is_seeded:
    findings.append((
        "Seeding",
        "🟢 Seeded process",
        "Seed surface area reduces nucleation reliance → lower mixing sensitivity.",
    ))
else:
    findings.append((
        "Seeding",
        "🟡 Unseeded",
        "Primary/secondary nucleation only → higher sensitivity.",
    ))

# ── Render findings ───────────────────────────────────────────────────
st.subheader("Assessment Results")
for mech, status, detail in findings:
    st.markdown(f"**{mech}** — {status}")
    st.caption(detail)

# ── Dominant risk ─────────────────────────────────────────────────────
st.subheader("Dominant Mixing Limitation")

_any_red = any("🔴" in s for _, s, _ in findings)
_any_yellow = any("🟡" in s for _, s, _ in findings)

if _micro_likely and _feed_sensitive:
    st.error(
        f"🔬 **Micro/Mesomixing** — fast nucleation ($t_{{ind}}$ = {t_ind:.4g} s) + "
        f"feed-based σ ({cryst_type}) → **highly mixing-sensitive**.\n\n"
        "**Scale-up:** maintain local ε at feed point · feed near impeller · "
        "reduce feed rate or add feed points · consider continuous crystallization."
    )
elif _micro_likely:
    st.warning(
        "🔬 **Micromixing** — local ε controls outcome.\n\n"
        "**Scale-up:** hold local ε constant · consider seeding · ensure adequate RPM at scale."
    )
elif _feed_sensitive and _micro_meso_sensitive:
    st.warning(
        "🌊 **Mesomixing** — feed-plume dispersion controls local σ.\n\n"
        "**Scale-up:** extend feed time · sub-surface feed near impeller · "
        "multiple feed points · maintain pumping rate."
    )
elif _macro_likely:
    st.info(
        "🌀 **Macromixing** — blend time may approach growth time at scale → σ non-uniformity.\n\n"
        "**Scale-up:** high-efficiency axial impellers · multiple impellers for tall vessels · "
        "monitor θ₉₅ vs t_G."
    )
elif not _any_red and not _any_yellow:
    st.success("✅ **Low mixing sensitivity.** Standard scale-up practices sufficient; confirm at pilot scale.")
else:
    st.info("⚠️ **Moderate risk** — some mechanisms may become relevant at scale. Run Da analysis for your reactor.")

# ── Recommended next steps ────────────────────────────────────────────
st.subheader("Recommended Next Steps")
_steps: list[str] = []

if _bourne_sensitive is None:
    _steps.append(
        "Run **Bourne Part 1** (🅱️  Bourne Protocol page) — response = particle size or polymorph."
    )

if _micro_likely or _micro_meso_sensitive:
    _steps.append(
        "Compute **Da numbers** for your reactor(s) (selector above or 🌀 Mixing Assessment page)."
    )

if _feed_sensitive:
    _steps.append(
        "Test **feed location** and **feed rate** at 2–3 conditions (Bourne Tests 2 & 3)."
    )

if _polymorph_risk:
    _steps.append(
        "Confirm **polymorph** under scale-representative mixing (in-situ Raman or XRPD)."
    )

if not is_seeded and _micro_meso_sensitive:
    _steps.append("Consider **seeding** to reduce nucleation reliance and mixing sensitivity.")

if _heat_sensitive:
    _steps.append(
        "Run **heat balance** on target vessel (🌀 Mixing Assessment or 📈 Reactor Comparison)."
    )

if not _steps:
    _steps.append("**Low risk** — standard scale-up practices sufficient.")

for i, step in enumerate(_steps, 1):
    st.markdown(f"{i}. {step}")

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# REFERENCES
# ═══════════════════════════════════════════════════════════════════════
st.header("References")
st.markdown("""
1. Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). *Handbook of Industrial
   Crystallization*, 3rd ed. Cambridge University Press.
2. Green, D.A. (2019). Ch. 10: Understanding and Modeling Crystallizer Mixing
   and Suspension Flow. In Myerson *et al.* (2019).
3. Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical
   Reactions.* Wiley.
4. Karpinski, P.H. & Baldyga, J. (2019). Ch. 8: Precipitation Processes.
   In Myerson *et al.* (2019).
5. Black, S.N. (2019). Ch. 13: Crystallization in the Pharmaceutical
   Industry. In Myerson *et al.* (2019).
6. Baldyga, J., Bourne, J.R. & Hearn, S.J. (1997). Interaction between
   Chemical Reactions and Mixing on Various Scales. *Chem. Eng. Sci.*,
   52(4), 457–466.
7. Mersmann, A. (2001). *Crystallization Technology Handbook*, 2nd ed.
   Marcel Dekker.
8. am Ende, D.J. (2019). *Chemical Engineering in the Pharmaceutical
   Industry: R&D to Manufacturing*, 2nd ed. Wiley.
""")
