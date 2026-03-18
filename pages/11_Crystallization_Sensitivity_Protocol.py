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


crystallizations = _load("crystallization_db", "crystallizations.csv")
reactors = _load("reactor_db", "reactors.csv")
custom_fluids = _load("fluid_db", "fluids.csv")

_solvent_names = sorted(SOLVENT_DB.keys())
_custom_names = (
    custom_fluids["fluid_name"].tolist() if not custom_fluids.empty else []
)
_all_fluid_names = _solvent_names + _custom_names

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
st.caption(
    "An interactive decision tree to assess which mixing mechanisms "
    "may limit a pharmaceutical crystallization at scale."
)

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
        "This protocol adapts the **Mixing Sensitivity Protocol** "
        "(Page 10) for crystallization processes.  The key difference "
        "is that the single characteristic reaction time $t_{rxn}$ is "
        "replaced by **multiple crystallization time scales**:\n\n"
        "- **Induction time** $t_{ind}$ — the delay between achieving "
        "supersaturation and detecting the first nuclei.\n"
        "- **Nucleation time** $t_N$ — the time to generate enough nuclei "
        "to consume the supersaturation.\n"
        "- **Growth time** $t_G$ — the time for crystal growth to deplete "
        "the supersaturation.\n\n"
        "These are compared against the same mixing time scales "
        "(micromixing $t_E$, mesomixing $t_{meso}$, macromixing $\\theta_{95}$) "
        "via Damköhler numbers to determine whether mixing limits the "
        "process.\n\n"
        "**References:** Myerson *et al.* (2019) *Handbook of Industrial "
        "Crystallization*, 3rd ed.; Baldyga & Bourne (1999); "
        "Green, D.A. (2019) Ch. 10; Black, S.N. (2019) Ch. 13."
    )

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 0 – BOURNE PROTOCOL PRE-SCREENING
# ═══════════════════════════════════════════════════════════════════════
st.header("0 · Bourne Protocol Pre-Screening")

with st.expander("ℹ️ Background — Bourne pre-screening for crystallization", expanded=False):
    st.markdown(
        "The Bourne Protocol pre-screen applies to crystallization just "
        "as it does to chemical reactions.  The relevant **response "
        "metrics** for crystallization are:\n\n"
        "- Mean particle size / chord length (FBRM)\n"
        "- Particle count or nucleation count\n"
        "- Polymorph form (XRPD, Raman)\n"
        "- Particle size distribution width (span)\n"
        "- Filterability / cake resistance\n\n"
        "Run the crystallization at two or more impeller speeds and/or "
        "feed rates and check whether **any** of these metrics change.  "
        "Use the **🧐 Bourne Protocol** page for the full experimental "
        "protocol."
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
        "from the **🧐 Bourne Protocol** page (use mean particle size "
        "or polymorph purity as the response metric).\n\n"
        "Proceeding with the theoretical assessment for now."
    )
    _bourne_sensitive = None

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 1 – CRYSTALLIZATION PARAMETERS
# ═══════════════════════════════════════════════════════════════════════
st.header("1 · Crystallization Parameters")

with st.expander("ℹ️ Background — why these parameters matter", expanded=False):
    st.markdown(
        "The mixing-sensitivity assessment requires characteristic time "
        "scales for the crystallization process.  These depend on:\n\n"
        "- **Supersaturation** ($\\sigma$) — the driving force for both "
        "nucleation and growth.\n"
        "- **Induction time** ($t_{ind}$) — the delay before nucleation "
        "is detected; inversely related to $J$ (nucleation rate).\n"
        "- **Growth rate constant** ($k_g$) and order ($g$) — govern how "
        "fast crystals grow once nucleated.\n"
        "- **Metastable zone width (MSZW)** — the maximum sustainable "
        "supersaturation before spontaneous nucleation.\n\n"
        "These can be measured experimentally (e.g., turbidity induction "
        "time, FBRM/PVM growth studies) or estimated from literature "
        "for similar compounds."
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
    st.info("👆 Please select an option above to continue.")
    st.stop()

if cryst_choice == "Select from Crystallization Database":
    if crystallizations.empty:
        st.error(
            "The Crystallization Database is empty.  Add entries to "
            "`data/crystallizations.csv` or choose **Enter parameters manually**."
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
    cx = crystallizations[crystallizations["crystallization_name"] == cx_name].iloc[0]

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
        "The method by which supersaturation is generated determines "
        "which mixing scales are most critical:\n\n"
        "| Type | How σ is generated | Dominant mixing concern |\n"
        "|------|-------------------|------------------------|\n"
        "| **Cooling** | Temperature reduction | Macromixing — uniformity of T field |\n"
        "| **Anti-solvent** | Addition of a poor solvent | Meso- & micromixing — feed plume σ spike |\n"
        "| **Reactive** | Chemical reaction generates insoluble product | Micromixing — instantaneous local σ |\n"
        "| **Evaporative** | Solvent removal at boiling surface | Macromixing — boiling zone σ gradient |\n"
        "| **pH-shift** | Acid/base addition changes solubility | Meso- & micromixing — local pH plume |\n\n"
        "*Anti-solvent and reactive crystallizations are the most "
        "mixing-sensitive because supersaturation is generated at the "
        "feed point, creating locally very high σ in the feed plume "
        "before mixing disperses it.* — Myerson (2019), Ch. 10"
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
        f"{_risk_icon} **{cryst_type} crystallization** — inherently "
        f"**high mixing sensitivity**.  Supersaturation is generated "
        f"locally at the feed point; micro- and mesomixing directly "
        f"control the supersaturation profile experienced by nucleating "
        f"and growing crystals."
    )
elif _risk_level == "moderate":
    st.info(
        f"{_risk_icon} **{cryst_type} crystallization** — **moderate "
        f"mixing sensitivity**.  Supersaturation gradients can develop "
        f"between the generation zone and the bulk."
    )
else:
    st.success(
        f"{_risk_icon} **{cryst_type} crystallization** — **lower "
        f"mixing sensitivity** (compared to anti-solvent or reactive).  "
        f"Macromixing (temperature uniformity) is the primary concern."
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
            "⚠️ Operating close to the metastable limit.  "
            "Any local σ spike from poor mixing could trigger "
            "**uncontrolled primary nucleation** or a polymorph change."
        )
        _sigma_risk = True
    elif sigma_ratio > 0.4:
        st.info(
            "Operating in the mid-range of the metastable zone.  "
            "Local σ spikes from feed addition could push regions "
            "above the metastable limit."
        )
        _sigma_risk = True
    else:
        st.success(
            "Operating well within the metastable zone.  "
            "Small local σ fluctuations are unlikely to cause problems."
        )
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
        "In a chemical reaction, micro- and mesomixing affect **selectivity** "
        "when competing reactions exist.  In crystallization, the analogous "
        "competition is between **nucleation** and **crystal growth**:\n\n"
        "- Both consume supersaturation.\n"
        "- **Nucleation rate** $B \\propto \\sigma^n$ ($n$ typically 2–10 for "
        "primary nucleation).\n"
        "- **Growth rate** $G \\propto \\sigma^g$ ($g$ typically 1–2).\n\n"
        "Because $n \\gg g$, nucleation is much more sensitive to local σ "
        "than growth is.  **Poor mixing creates local σ spikes that "
        "disproportionately accelerate nucleation** → excess fines, "
        "broader CSD, and potentially wrong polymorph.\n\n"
        "This is precisely analogous to a fast side-reaction consuming "
        "reagent before the desired slow reaction can proceed.\n\n"
        "### Micromixing\n"
        "At the Kolmogorov/Batchelor scale, the engulfment rate "
        "$E = 0.058 \\sqrt{\\varepsilon/\\nu}$ determines how fast the "
        "feed solution is homogenised.  The micromixing time:\n\n"
        "$$t_E = 17.3 \\sqrt{\\nu / \\varepsilon}$$\n\n"
        "When $t_E > t_{ind}$, nucleation occurs in the unmixed feed "
        "plume at locally high σ → **micromixing-controlled** CSD.\n\n"
        "### Mesomixing\n"
        "At the feed-plume scale, the turbulent dispersion time "
        "determines how quickly the fresh feed is diluted into the bulk.  "
        "For anti-solvent or reactive crystallization, mesomixing directly "
        "sets the local σ in the feed zone.  Increasing the feed rate "
        "increases the mesomixing time → higher local σ → more fines."
    )

# ── Feed-sensitive or not ─────────────────────────────────────────────
if _feed_sensitive:
    st.warning(
        f"**{cryst_type} crystallization** involves continuous feed addition.  "
        "Micro- and mesomixing at the feed point directly control the local "
        "supersaturation and therefore the nucleation / growth balance."
    )

    st.markdown(
        "**Key risk factors for micro/mesomixing sensitivity:**\n"
        "- Short induction time ($t_{ind}$ < 10 s)\n"
        "- High target supersaturation ($\\sigma > 0.5$)\n"
        "- Known polymorphism\n"
        "- Oiling-out tendency at high σ\n"
    )

    _micro_meso_sensitive = True
else:
    # Cooling or evaporative — still evaluate
    if t_ind < 10:
        st.warning(
            f"Short induction time ($t_{{ind}}$ = {t_ind:.4g} s).  "
            "Even in cooling crystallization, rapid nucleation kinetics can "
            "make the process sensitive to local temperature (and hence σ) "
            "non-uniformity caused by imperfect macromixing."
        )
        _micro_meso_sensitive = True
    elif sigma_target > 0.5:
        st.info(
            "High target supersaturation — nucleation kinetics may be fast "
            "enough to be influenced by local mixing conditions."
        )
        _micro_meso_sensitive = True
    else:
        st.success(
            "✅ **Low mixing-sensitivity risk** from nucleation/growth "
            "competition.  Cooling crystallization with moderate σ and "
            "a reasonable induction time is typically not micromixing-controlled."
        )
        _micro_meso_sensitive = False

# ── Polymorphism ──────────────────────────────────────────────────────
if has_polymorphs:
    st.warning(
        f"⚠️ **Polymorphism detected** ({polymorph_info}).  "
        "Local supersaturation spikes from poor mixing can nucleate a "
        "**metastable polymorph** (Ostwald's rule of stages).  "
        "The metastable form nucleates faster at high σ but may convert "
        "to the stable form over time — or may not, depending on kinetics."
    )
    _polymorph_risk = True
else:
    _polymorph_risk = False

# ── Seeding ───────────────────────────────────────────────────────────
if is_seeded:
    st.info(
        "✅ **Seeded process** — seed crystals provide surface area "
        "that consumes supersaturation via growth rather than nucleation.  "
        "This reduces (but does not eliminate) mixing sensitivity, "
        "especially if the seed loading is adequate.\n\n"
        "*Rule of thumb:* A well-seeded process is 2–5× less sensitive "
        "to mixing than an unseeded one at the same σ."
    )
else:
    st.info(
        "⚠️ **Unseeded process** — all nucleation is primary or "
        "secondary.  The process is **more sensitive** to local σ "
        "conditions set by mixing."
    )

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 4 – HEAT TRANSFER
# ═══════════════════════════════════════════════════════════════════════
st.header("4 · Heat Transfer Screening")

_has_enthalpy = delta_H_cryst != 0.0
_heat_sensitive = False

with st.expander("ℹ️ Background — heat effects in crystallization", expanded=False):
    st.markdown(
        "Crystallization is **exothermic** ($\\Delta H_{cryst}$ < 0 for most "
        "systems).  Additionally:\n\n"
        "- **Cooling crystallization:** The cooling rate must be controlled "
        "to maintain supersaturation within the metastable zone.  Poor "
        "jacket heat transfer at scale can create cold spots → local σ "
        "spikes → wall encrustation and excessive nucleation.\n"
        "- **Anti-solvent crystallization:** The anti-solvent is often at a "
        "different temperature; mixing of streams involves both mass and "
        "heat transfer.\n"
        "- **Reactive crystallization:** Reaction enthalpy adds to "
        "crystallization enthalpy; this can be very significant.\n\n"
        "At scale, the **surface-area-to-volume ratio** decreases, making "
        "it harder to remove heat.  This can force slower addition rates "
        "or lower batch temperatures, altering the σ profile."
    )

if _has_enthalpy:
    abs_dH = abs(delta_H_cryst)
    st.markdown(
        f"**ΔH_cryst = {delta_H_cryst:.1f} kJ/mol** "
        f"({'exothermic' if delta_H_cryst < 0 else 'endothermic'})"
    )
    if abs_dH >= 40:
        st.warning(
            f"🔴 **High crystallization enthalpy** — |ΔH| = {abs_dH:.1f} kJ/mol.  "
            "Heat removal may limit cooling rates and feed rates at scale."
        )
        _heat_sensitive = True
    elif abs_dH >= 20:
        st.info(
            f"🟡 **Moderate crystallization enthalpy** — |ΔH| = {abs_dH:.1f} kJ/mol.  "
            "Unlikely to be limiting in most configurations, but check at scale."
        )
    else:
        st.success(
            f"🟢 **Low crystallization enthalpy** — |ΔH| = {abs_dH:.1f} kJ/mol."
        )

    if cryst_type == "Cooling":
        st.markdown(
            "For cooling crystallization, also consider:\n"
            "- Jacket ΔT limitations at scale\n"
            "- Wall encrustation from local cold spots\n"
            "- Required cooling rate vs. available UA"
        )
else:
    st.info(
        "No ΔH_cryst data provided.  Consider measuring by calorimetry "
        "(e.g., RC1, µRC) or estimating from the van 't Hoff equation."
    )

st.divider()

# ═══════════════════════════════════════════════════════════════════════
# STEP 5 – DAMKÖHLER COMPARISON
# ═══════════════════════════════════════════════════════════════════════
st.header("5 · Mixing Time vs Crystallization Time (Damköhler Analysis)")

with st.expander("ℹ️ Background — Damköhler numbers for crystallization", expanded=False):
    st.markdown(
        "Following Baldyga & Bourne (1999), the controlling mechanism is "
        "identified by comparing characteristic mixing times to "
        "crystallization times:\n\n"
        "$$Da = \\frac{t_{mix}}{t_{cryst}}$$\n\n"
        "| Da range | Interpretation |\n"
        "|----------|----------------|\n"
        "| Da < 0.01 | Crystallization much slower than mixing — **not sensitive** |\n"
        "| 0.01 – 0.1 | **Likely not sensitive**, but monitor at scale |\n"
        "| 0.1 – 1 | **Potentially sensitive** — on similar timescales |\n"
        "| 1 – 10 | **Likely sensitive** — mixing limits process |\n"
        "| Da > 10 | **Highly sensitive** — mixing fully controls outcome |\n\n"
        "Three comparisons are made:\n"
        "- $Da_{micro} = t_E / t_{ind}$ — micromixing vs induction time\n"
        "- $Da_{macro} = \\theta_{95} / t_G$ — macromixing vs growth time\n"
        "- For feed-sensitive processes: $Da_{meso}$ is evaluated qualitatively"
    )

st.subheader("Heuristic Assessment (no reactor selected)")

st.markdown(
    "Before computing Da for a specific reactor, here is a heuristic "
    "assessment based on typical mixing times at different scales.\n\n"
    f"**Your crystallization times:**\n"
    f"- Induction time: $t_{{ind}}$ = {t_ind:.4g} s\n"
    f"- Growth time: $t_G$ = {t_G:.4g} s"
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
        f"🔴 **Very fast nucleation** ($t_{{ind}}$ = {t_ind:.4g} s).  "
        "**Micromixing-sensitive** in most configurations.  "
        "Local energy dissipation at the feed point dominates the "
        "supersaturation profile and CSD."
    )
    _micro_likely = True
elif t_ind > 0 and t_ind < 10.0:
    st.warning(
        f"🟡 **Fast nucleation** ($t_{{ind}}$ = {t_ind:.4g} s).  "
        "Micromixing likely relevant in pilot and production vessels."
    )
    _micro_likely = True
elif t_ind > 0 and t_ind < 60.0:
    st.info(
        f"🔵 **Moderate nucleation** ($t_{{ind}}$ = {t_ind:.4g} s).  "
        "Micromixing less likely to dominate, but macromixing "
        "(blend time) may matter at production scale."
    )
else:
    st.success(
        f"🟢 **Slow nucleation** ($t_{{ind}}$ = {t_ind:.4g} s).  "
        "Mixing is unlikely to limit crystallization in well-agitated vessels."
    )

# Macro assessment
_macro_likely = False
if t_G < np.inf and t_G < 120:
    st.info(
        f"Growth time $t_G$ = {t_G:.4g} s is within the range of blend "
        f"times in larger vessels.  Macromixing may affect uniformity of "
        f"supersaturation and therefore growth rate distribution."
    )
    if t_G < 30:
        _macro_likely = True

# ── Optional: compute Da for a specific reactor ──────────────────────
st.subheader("(Optional) Compute Da for a Specific Reactor")

if reactors.empty:
    st.info("No reactors in the database. Add reactors on the 🧪 Reactor Database page.")
else:
    _reactor_list = reactors["reactor_name"].tolist()
    reactor_name = st.selectbox(
        "Select reactor (optional)",
        ["— Skip —"] + _reactor_list,
        key=_key("reactor_sel"),
    )

    if reactor_name != "— Skip —":
        r = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
        D_tank = float(r["D_tank_m"])
        D_imp = float(r["D_imp_m"])
        Np_val = float(r["Np"]) if pd.notna(r.get("Np")) else 1.0
        Nq_val = float(r["Nq"]) if pd.notna(r.get("Nq")) else 0.79
        N_rpm_min = float(r["N_rpm_min"]) if pd.notna(r.get("N_rpm_min")) else None
        N_rpm_max = float(r["N_rpm_max"]) if pd.notna(r.get("N_rpm_max")) else None

        V_L_min = float(r["V_L_min"]) if pd.notna(r.get("V_L_min")) else float(r["V_L"])
        V_L_max = float(r["V_L_max"]) if pd.notna(r.get("V_L_max")) else float(r["V_L"])
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
            _cf = custom_fluids[custom_fluids["fluid_name"] == fluid_sel].iloc[0]
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
                "📌 For **feed-sensitive** crystallizations, also consider:\n"
                "- **Feed location** — feed near the impeller (high ε) vs. surface\n"
                "- **Feed rate** — slower addition → lower local σ\n"
                "- **Number of feed points** — multiple points reduce local σ"
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
        "Experimental pre-screen showed CSD / polymorph changes with mixing conditions.",
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
        "Consider running Bourne Protocol Part 1 (response = mean particle size or polymorph).",
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
        f"σ_target = {sigma_target:.3f}; operating {'close to' if sigma_target / (dc_dT * MSZW / c_sat_ref if c_sat_ref > 0 and MSZW > 0 else 1) > 0.5 else 'in mid-range of'} metastable limit.",
    ))
else:
    findings.append((
        "Supersaturation profile",
        "🟢 Low risk",
        f"σ_target = {sigma_target:.3f}; well within metastable zone.",
    ))

# Micromixing
if _micro_likely:
    findings.append((
        "Micromixing (nucleation control)",
        "🔴 Likely sensitive",
        f"t_ind = {t_ind:.4g} s — fast nucleation relative to typical micromixing times.",
    ))
else:
    findings.append((
        "Micromixing (nucleation control)",
        "🟢 Unlikely",
        f"t_ind = {t_ind:.4g} s — slow enough that micromixing is not rate-limiting.",
    ))

# Micro/mesomixing (feed-related)
if _micro_meso_sensitive and _feed_sensitive:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🔴 Likely sensitive",
        "Feed addition generates local σ spike; mixing at the feed point controls CSD.",
    ))
elif _micro_meso_sensitive:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🟡 Potentially sensitive",
        "Short induction time or high σ suggests sensitivity, even without direct feed addition.",
    ))
else:
    findings.append((
        "Micro/mesomixing (feed plume)",
        "🟢 Not a primary factor",
        "No direct feed addition or slow nucleation kinetics.",
    ))

# Macromixing
if _macro_likely:
    findings.append((
        "Macromixing (blend time)",
        "🟡 Check at scale",
        f"t_G = {t_G:.4g} s — within blend time range of larger vessels.",
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
        f"Known polymorphs ({polymorph_info}); local σ spikes can nucleate metastable form.",
    ))
else:
    findings.append((
        "Polymorphism",
        "🟢 Low risk",
        "No known polymorphs or single confirmed form.",
    ))

# Heat transfer
if _heat_sensitive:
    findings.append((
        "Heat transfer",
        "🔴 Likely limiting",
        f"|ΔH_cryst| = {abs(delta_H_cryst):.1f} kJ/mol — check cooling capacity at scale.",
    ))
elif _has_enthalpy:
    findings.append((
        "Heat transfer",
        "🟢 Manageable",
        f"|ΔH_cryst| = {abs(delta_H_cryst):.1f} kJ/mol — modest enthalpy.",
    ))
else:
    findings.append((
        "Heat transfer",
        "⚪ Unknown",
        "No ΔH_cryst data — consider measuring.",
    ))

# Seeding
if is_seeded:
    findings.append((
        "Seeding",
        "🟢 Seeded process",
        "Seed surface area reduces reliance on primary nucleation → lower mixing sensitivity.",
    ))
else:
    findings.append((
        "Seeding",
        "🟡 Unseeded",
        "All nucleation is primary/secondary → higher mixing sensitivity.",
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
        "🔬 **Micromixing / Mesomixing** — The combination of fast nucleation "
        f"($t_{{ind}}$ = {t_ind:.4g} s) and feed-based supersaturation generation "
        f"({cryst_type}) makes this crystallization **highly mixing-sensitive**.\n\n"
        "**Scale-up strategy:**\n"
        "- Maintain or increase **local ε** at the feed point.\n"
        "- Feed near the impeller tip (high-turbulence zone).\n"
        "- Reduce feed rate or use multiple feed points.\n"
        "- Consider in-line or continuous crystallization (e.g., impinging jets, "
        "static mixers, MSMPR cascade)."
    )
elif _micro_likely:
    st.warning(
        "🔬 **Micromixing** — Fast nucleation kinetics mean local energy "
        "dissipation controls the outcome.\n\n"
        "**Scale-up strategy:**\n"
        "- Hold **local ε** at the feed/nucleation zone constant.\n"
        "- Consider seeding to reduce nucleation dependence.\n"
        "- Ensure adequate impeller speed at scale."
    )
elif _feed_sensitive and _micro_meso_sensitive:
    st.warning(
        "🌊 **Mesomixing** — Feed-plume dispersion controls local σ.\n\n"
        "**Scale-up strategy:**\n"
        "- Extend feed time or reduce feed rate.\n"
        "- Use sub-surface feed near the impeller.\n"
        "- Consider multiple feed points or a distributor ring.\n"
        "- Maintain impeller pumping rate to ensure rapid plume dilution."
    )
elif _macro_likely:
    st.info(
        "🌀 **Macromixing** — Blend time may become comparable to growth time "
        "at scale, causing σ non-uniformity.\n\n"
        "**Scale-up strategy:**\n"
        "- Use high-efficiency axial impellers (hydrofoils).\n"
        "- Consider multiple impellers for tall vessels.\n"
        "- Monitor blend time vs. growth time at scale."
    )
elif not _any_red and not _any_yellow:
    st.success(
        "✅ **Low mixing sensitivity overall.**  Standard scale-up practices "
        "should be sufficient.  Monitor at pilot scale for confirmation."
    )
else:
    st.info(
        "⚠️ **Moderate risk** — some mixing mechanisms may become relevant "
        "at scale.  Targeted Damköhler analysis for your specific reactor "
        "is recommended."
    )

# ── Recommended next steps ────────────────────────────────────────────
st.subheader("Recommended Next Steps")
_steps: list[str] = []

if _bourne_sensitive is None:
    _steps.append(
        "Run **Bourne Protocol Part 1** (🧐 Bourne Protocol page) using "
        "**mean particle size** or **polymorph form** as the response metric."
    )

if _micro_likely or _micro_meso_sensitive:
    _steps.append(
        "Compute **Damköhler numbers** for your specific reactor(s) "
        "(use the reactor selector above, or the 🌀 Mixing Assessment page)."
    )

if _feed_sensitive:
    _steps.append(
        "Evaluate **feed location** and **feed rate** experimentally — "
        "run the reaction at 2–3 feed positions and feed times "
        "(Bourne Tests 2 & 3)."
    )

if _polymorph_risk:
    _steps.append(
        "Confirm **polymorphic form** under mixing conditions representative "
        "of the target scale (in-situ Raman or XRPD)."
    )

if not is_seeded and _micro_meso_sensitive:
    _steps.append(
        "Consider **seeding** to reduce reliance on primary nucleation "
        "and lower mixing sensitivity."
    )

if _heat_sensitive:
    _steps.append(
        "Run a **heat balance** on the target vessel to confirm adequate "
        "cooling capacity (🔥 on 🌀 Mixing Assessment or 📈 Reactor Comparison)."
    )

if not _steps:
    _steps.append(
        "The crystallization appears **low risk** for mixing sensitivity.  "
        "Standard scale-up practices should be sufficient."
    )

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
