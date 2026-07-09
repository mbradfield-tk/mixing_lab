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
import re

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


# ── Consistent section-outcome helpers ───────────────────────────────────
# Every step closes with the same colour-coded, three-part layout so the user
# reads each section the same way:
#   1. Assessment      — the section's conclusion (colour encodes severity)
#   2. Warning         — a specific risk to act on (amber), only when needed
#   3. Recommendation  — a concrete next action (blue), only when needed
_ASSESS_ICON = {"critical": "🔴", "warning": "🟡", "caution": "🟡", "ok": "🟢"}


def _assessment(kind: str, text: str) -> None:
    """Render a section's assessment conclusion using a traffic-light colour.

    🔴 critical (red) · 🟡 warning / caution (amber) · 🟢 ok (green).
    """
    msg = f"{_ASSESS_ICON.get(kind, '🟡')} **Assessment:** {text}"
    if kind == "critical":
        st.error(msg)
    elif kind == "ok":
        st.success(msg)
    else:  # "warning" and "caution" → traffic-light amber
        st.warning(msg)


def _section_warning(text: str) -> None:
    """Render a specific risk callout (consistent amber styling)."""
    st.warning(f"⚠️ **Warning:** {text}")


def _recommendation(text: str) -> None:
    """Render a recommended action (consistent blue styling)."""
    st.info(f"👉 **Recommendation:** {text}")


# Scale-up actions keyed by the Bourne-identified controlling mixing scale.
# Shared by the full summary (Step 6) and the Bourne-only final output (Step 1).
_BOURNE_MECH_ACTIONS = {
    "Micromixing": (
        "Hold local energy dissipation ε (P/V) constant on scale-up and keep the "
        "feed point in the high-shear zone near the impeller. Confirm with Da_micro "
        "on the 🌀 Vessel Assessment page."
    ),
    "Mesomixing": (
        "Control feed-plume dispersion: hold local ε constant at the feed point "
        "(match P/V → a lower RPM at larger scale, not constant RPM) and cut the "
        "local feed rate — extend feed time, add feed points, and/or use a smaller "
        "feed-pipe diameter."
    ),
    "Macromixing": (
        "Reduce bulk blend time on scale-up: use high-efficiency (hydrofoil) or "
        "multiple impellers, optimise baffling, or consider static/in-line mixers."
    ),
}


def _bourne_sensitive_kpi_phrase(parsed) -> str:
    """Return a comma-joined list of the KPIs the Bourne import flagged sensitive."""
    _names: list[str] = []
    _seen = set()
    if parsed:
        for _n in (1, 2, 3):
            _sk = (parsed.get(f"test{_n}_sensitive_kpis") or "").strip()
            for _entry in (x.strip() for x in _sk.split(";")):
                if not _entry:
                    continue
                _name = re.sub(r"\s*\((?:[\d.]+%|qualitative)\)\s*$", "", _entry).strip()
                if _name and _name not in _seen:
                    _seen.add(_name)
                    _names.append(_name)
    return ", ".join(_names) if _names else "the tracked response(s)"


# Purpose of each Bourne test, used to phrase recommendations about the tests
# that still need to be run to resolve the controlling mixing scale.
_TEST_PURPOSE = {1: "impeller speed", 2: "feed rate/time", 3: "feed location"}


def _bourne_remaining_tests(done_tests, needed=(2, 3)) -> list[int]:
    """Tests still needed (among ``needed``) to resolve the controlling scale."""
    return [t for t in needed if t not in (done_tests or [])]


def _fmt_tests(nums) -> str:
    """Render test numbers as e.g. 'Test 3' or 'Tests 2 and 3'."""
    if not nums:
        return ""
    if len(nums) == 1:
        return f"Test {nums[0]}"
    return "Tests " + " and ".join(str(n) for n in nums)


def _fmt_test_purposes(nums) -> str:
    """Render the purposes of the given tests, e.g. 'feed location'."""
    return " and ".join(_TEST_PURPOSE[n] for n in nums if n in _TEST_PURPOSE)


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
        "The **Bourne Protocol** experimentally reveals whether — and how — a "
        "reaction is mixing-sensitive by varying impeller speed (Test 1), feed "
        "rate / time (Test 2), and feed location (Test 3).\n\n"
        "Run it on the **🅱️ Bourne Protocol** page, use its "
        "**📤 Export results for the Sensitivity Protocol** button, and import the "
        "resulting CSV below. The findings from each completed test are parsed "
        "automatically to drive this assessment."
    )

# Import the machine-readable results CSV exported from the Bourne Protocol page.
_bourne_up = st.file_uploader(
    "Import Bourne Protocol results (CSV)",
    type=["csv"],
    key=_key("bourne_import"),
    help="Generate this file on the 🅱️ Bourne Protocol page via "
         "'📤 Export results for the Sensitivity Protocol'.",
)
if _bourne_up is not None:
    try:
        _bdf = pd.read_csv(_bourne_up, dtype=str, keep_default_na=False)
        _bdict = (
            dict(zip(_bdf["field"], _bdf["value"]))
            if {"field", "value"} <= set(_bdf.columns) else {}
        )
        if _bdict.get("record_type") == "bourne_results":
            st.session_state[_key("bourne_parsed")] = _bdict
        else:
            st.error(
                "That CSV is not a Bourne Protocol results export. Use the "
                "**📤 Export results for the Sensitivity Protocol** button on the "
                "🅱️ Bourne Protocol page."
            )
    except Exception as exc:
        st.error(f"Could not read the file: {exc}")

_bourne_parsed = st.session_state.get(_key("bourne_parsed"))

# Persisted views of the imported Bourne results, surfaced later in the Summary
# section and the PDF report (empty when nothing was imported).
_bourne_test_details: list[dict] = []
_bourne_meta: dict = {}
_bourne_done_tests: list[int] = []

_skip_bourne = st.checkbox(
    "I haven't run the Bourne Protocol — skip this step",
    key=_key("bourne_skip"),
)

if _bourne_parsed is None and not _skip_bourne:
    st.info("👆 Import a Bourne Protocol results CSV above, or tick the box to skip.")
    st.stop()

if _bourne_parsed is not None and not _skip_bourne:
    # ── Provenance of the imported results ───────────────────────────────
    _meta_bits = []
    for _lbl, _fld in [
        ("Project", "project_name"), ("Step", "step_number"),
        ("Unit op", "unit_operation"), ("Reactor", "reactor"), ("Fluid", "fluid"),
    ]:
        _val = (_bourne_parsed.get(_fld) or "").strip()
        if _val:
            _meta_bits.append(f"**{_lbl}:** {_val}")
            _bourne_meta[_lbl] = _val
    if _meta_bits:
        st.caption("Imported Bourne results — " + "  •  ".join(_meta_bits))

    # ── Per-test completion & findings ───────────────────────────────────
    # Only tests that were actually assessed are present in the CSV, so only
    # those appear here (e.g. an unrun Test 3 is not listed).
    _probe_names = {
        1: "Test 1 — Impeller speed (micromixing screen)",
        2: "Test 2 — Feed rate / time (meso vs micro)",
        3: "Test 3 — Feed location (meso vs macro)",
    }
    _test_rows = []
    _done_tests = []
    for _n in (1, 2, 3):
        _assessed = (
            _bourne_parsed.get(f"test{_n}_assessed")
            or _bourne_parsed.get(f"test{_n}_completed")  # back-compat
            or "no"
        ) == "yes"
        if not _assessed:
            continue
        _done_tests.append(_n)
        _finding = (_bourne_parsed.get(f"test{_n}_finding") or "").strip()
        _sens_kpis = (_bourne_parsed.get(f"test{_n}_sensitive_kpis") or "").strip()
        _test_rows.append({
            "Test": _probe_names[_n],
            "Finding": _finding or "—",
            "Sensitive KPI(s)": _sens_kpis if _sens_kpis else "None (no KPI over threshold)",
        })
    if _test_rows:
        st.dataframe(pd.DataFrame(_test_rows), width='stretch', hide_index=True)
    _bourne_test_details = list(_test_rows)
    _bourne_done_tests = list(_done_tests)
    _not_done = [_n for _n in (1, 2, 3) if _n not in _done_tests]
    if _not_done:
        st.caption(
            "Not performed: " + ", ".join(f"Test {_n}" for _n in _not_done) + "."
        )

    # ── Derive sensitivity + controlling mechanism from the import ───────
    _overall = (_bourne_parsed.get("overall_sensitive") or "unknown").strip()
    _dominant = (_bourne_parsed.get("dominant_mechanism") or "").strip()

    if _overall == "yes":
        _bourne_sensitive = True
    elif _overall == "no":
        _bourne_sensitive = False
    else:
        _bourne_sensitive = None

    _bourne_mechanisms = (
        [_dominant] if _dominant in ("Micromixing", "Mesomixing", "Macromixing")
        else []
    )

    # ── Section outcome ──────────────────────────────────────────────────
    if _bourne_sensitive is True:
        if _bourne_mechanisms:
            _assessment(
                "critical",
                f"Bourne Protocol confirmed **mixing sensitivity** — the controlling "
                f"scale is **{_bourne_mechanisms[0].lower()}**. This is carried into "
                "the summary as an experimentally confirmed result.",
            )
        else:
            _assessment(
                "critical",
                "Bourne Protocol confirmed **mixing sensitivity**, but the controlling "
                "scale is not yet fully resolved.",
            )
            _rem = _bourne_remaining_tests(_done_tests)
            if _rem:
                _recommendation(
                    f"Complete **{_fmt_tests(_rem)}** of the Bourne Protocol "
                    f"({_fmt_test_purposes(_rem)}) to pinpoint whether micro-, meso-, "
                    "or macromixing controls, then re-export and re-import."
                )
            else:
                _recommendation(
                    "Re-run the Bourne Protocol decision tree to resolve the "
                    "controlling scale, then re-export and re-import."
                )
    elif _bourne_sensitive is False:
        _assessment(
            "ok",
            "Bourne Protocol showed **no mixing sensitivity** at lab scale (the Test 1 "
            "response was insensitive to impeller speed). The remaining steps confirm "
            "this and check for latent risks at larger scale.",
        )
    else:
        _assessment(
            "caution",
            "Bourne results imported, but **Test 1 was not completed** — experimental "
            "mixing sensitivity is undetermined.",
        )
        _recommendation(
            "Complete at least **Test 1** on the 🅱️ Bourne Protocol page, then re-export "
            "and re-import the results here."
        )
else:
    # User chose to skip the Bourne import.
    _assessment(
        "caution",
        "Bourne pre-screen skipped — proceeding with the theoretical assessment for now.",
    )
    _recommendation(
        "Running the Bourne Protocol and importing its results here provides a direct "
        "experimental answer. You can run it from the **🅱️ Bourne Protocol** page."
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
        "Approximate – use a similar reaction as proxy (from database)",
        "No – I need to measure them",
    ],
    key=_key("kinetics_choice"),
    index=None,
)

if kinetics_choice is None:
    st.info("👆 Please select an option above to continue.")
    st.stop()

if kinetics_choice == "No – I need to measure them":
    # If the Bourne Protocol was executed, offer to finish here and use those
    # experimental results as the final output (the kinetics-based timescale /
    # Damköhler analysis cannot be done without kinetics).
    _bourne_available = bool(_bourne_test_details) or _bourne_sensitive is not None
    _finish_bourne = False
    if _bourne_available:
        _finish_bourne = st.checkbox(
            "✅ Finish now and use the Bourne Protocol results as the final output "
            "(the kinetics-based timescale / Damköhler analysis will be skipped)",
            key=_key("finish_bourne_only"),
        )

    if _finish_bourne:
        st.header("Bourne Protocol Result — Final Output")

        _bo_caveat = (
            "This result is based on the **Bourne Protocol experiments only**. The "
            "kinetics-based **timescale (Damköhler) analysis has not been performed**, "
            "so the conclusions are valid **only for the operating ranges tested** in "
            "the Bourne Protocol. Mixing sensitivities may still arise under different "
            "conditions (larger scale, different feed rate, viscosity, concentration, "
            "or temperature)."
        )
        st.warning("⚠️ **Caveat:** " + _bo_caveat)

        _bo_kpi_phrase = _bourne_sensitive_kpi_phrase(_bourne_parsed)
        _bo_rem = _bourne_remaining_tests(_bourne_done_tests)
        _bo_rem_phrase = _fmt_tests(_bo_rem) or "the remaining Bourne Protocol test(s)"

        # Concise verdict + status line (mirrors the Bourne Protocol conclusion).
        if _bourne_sensitive is True:
            if _bourne_mechanisms:
                _bo_status = f"Mixing sensitivity confirmed — controlling scale: {_bourne_mechanisms[0]}"
                st.error(
                    f"🔴 **Mixing sensitivity confirmed** by the Bourne Protocol — the "
                    f"controlling scale is **{_bourne_mechanisms[0].lower()}**. Affected "
                    f"response(s): {_bo_kpi_phrase}."
                )
            else:
                _bo_status = "Mixing sensitivity confirmed — controlling scale not yet resolved"
                st.error(
                    f"🔴 **Mixing sensitivity confirmed** by the Bourne Protocol (affected "
                    f"response(s): {_bo_kpi_phrase}), but the controlling scale is not yet "
                    f"fully resolved — complete {_bo_rem_phrase} to pinpoint micro-, meso-, "
                    "or macromixing."
                )
        elif _bourne_sensitive is False:
            _bo_status = "No mixing sensitivity observed over the tested ranges"
            st.success(
                "🟢 **No mixing sensitivity observed** by the Bourne Protocol over the "
                "ranges tested."
            )
        else:
            _bo_status = "Undetermined — Test 1 not completed"
            st.warning(
                "🟡 **Undetermined** — the Bourne results were imported but Test 1 was not "
                "completed, so mixing sensitivity cannot be concluded."
            )

        # Per-test experimental findings
        if _bourne_test_details:
            st.subheader("Bourne Protocol Experimental Findings")
            if _bourne_meta:
                st.caption("  •  ".join(f"**{_k}:** {_v}" for _k, _v in _bourne_meta.items()))
            st.dataframe(
                pd.DataFrame(_bourne_test_details), width='stretch', hide_index=True
            )

        # Recommendations
        st.subheader("Recommended Next Steps")
        _bo_steps: list[dict[str, str]] = []
        for _m in _bourne_mechanisms:
            if _m in _BOURNE_MECH_ACTIONS:
                _bo_steps.append({
                    "Area": f"{_m} (Bourne-confirmed)",
                    "Action": _BOURNE_MECH_ACTIONS[_m],
                })
        if _bourne_sensitive is True and not _bourne_mechanisms and _bo_rem:
            _bo_steps.append({
                "Area": "Resolve controlling scale",
                "Action": f"Complete Bourne {_fmt_tests(_bo_rem)} "
                          f"({_fmt_test_purposes(_bo_rem)}) to identify whether micro-, "
                          "meso-, or macromixing controls.",
            })
        _bo_steps.append({
            "Area": "Kinetics / timescale analysis",
            "Action": "Measure reaction kinetics (k, order, C₀) and calorimetry (ΔH), "
                      "then re-run this protocol for the full Damköhler timescale "
                      "analysis and scale-dependent risk assessment.",
        })
        st.dataframe(pd.DataFrame(_bo_steps), width='stretch', hide_index=True)

        # ── PDF export (Bourne-only) ─────────────────────────────────────
        st.divider()
        st.subheader("Export Report")
        if st.button("📥 Export Bourne Result PDF", type="primary",
                     key="p10_bourne_only_pdf"):
            with st.spinner("Generating PDF…"):
                try:
                    from utils.report_builder import build_protocol_pdf, report_filename

                    _bo_label = (
                        _bourne_meta.get("Project") or _bourne_meta.get("Reactor")
                        or "Bourne Protocol result"
                    )
                    _bo_findings = []
                    if _bourne_sensitive is True:
                        _detail = f"Affected response(s): {_bo_kpi_phrase}."
                        _detail += (
                            f" Controlling scale: {_bourne_mechanisms[0]}."
                            if _bourne_mechanisms else
                            f" Controlling scale not yet resolved (complete {_bo_rem_phrase})."
                        )
                        _bo_findings.append(
                            ("Bourne pre-screen", "🔴 Mixing sensitivity confirmed", _detail)
                        )
                    elif _bourne_sensitive is False:
                        _bo_findings.append((
                            "Bourne pre-screen", "🟢 No sensitivity observed",
                            "No mixing sensitivity observed over the ranges tested.",
                        ))
                    _bo_snap = {
                        "reaction": _bo_label,
                        "t_rxn": None,
                        "rxn_delta_H": 0.0,
                        "phases": [],
                        "findings": _bo_findings,
                        "next_steps": _bo_steps,
                        "bourne_result": _bo_status,
                        "bourne_meta": _bourne_meta,
                        "bourne_tests": _bourne_test_details,
                        "bourne_mechanism": _bourne_mechanisms[0] if _bourne_mechanisms else "",
                        "competing": "Not assessed (no kinetics)",
                        "overall_verdict": "",
                        "using_approximate": False,
                        "is_semi_batch": False,
                        "bourne_only": True,
                        "caveat": _bo_caveat,
                    }
                    st.session_state["_p10_pdf_bytes"] = build_protocol_pdf(_bo_snap)
                    st.session_state["_p10_pdf_name"] = report_filename(
                        "Bourne_Result", _bo_label
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

        st.stop()

    _assessment(
        "warning",
        "Kinetics are not yet available — they must be measured before this "
        "assessment can be completed.",
    )
    _recommendation(
        "1. Conduct **time-course experiments** (e.g. in-situ FTIR / ReactIR, "
        "sampling + HPLC) to determine the rate constant *k* and reaction order.\n\n"
        "2. Perform **reaction calorimetry** (e.g. RC1, µRC, or Simular) to "
        "measure the heat of reaction *ΔH*.\n\n"
        "3. Add the results to the **🧪 Reaction Database** page.\n\n"
        "4. Return here and **select the reaction from the database** to continue.\n\n"
        "If exact kinetics are unavailable but a structurally similar reaction is "
        "in the database, choose **\"Approximate\"** above to use it as a proxy."
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

# ── Rendered kinetics equation ───────────────────────────────────────────
# Show the rate law first in general (symbolic) form, then with the actual
# rate constant substituted, and finally the characteristic reaction time —
# so the user can see exactly how t_rxn was obtained. The concentration
# dependency reflects the actual reaction order rather than a generic f(C).
st.markdown("**Kinetic model**")

# Rate constant symbol (primed for pseudo/lumped constants).
_k_sym = "k'" if rxn_order.startswith("pseudo") else "k"

# Parse the numeric reaction order (handles "1", "2", "0", "1.5", "pseudo-2").
try:
    _ord_val = float(rxn_order.replace("pseudo-", ""))
except ValueError:
    _ord_val = None

# Concentration-dependency term C^n for the rate law.
if _ord_val == 0:
    _conc_term = ""            # zero order: rate = k
elif _ord_val == 1:
    _conc_term = "C"
elif _ord_val is not None and _ord_val == int(_ord_val):
    _conc_term = rf"C^{{{int(_ord_val)}}}"
elif _ord_val is not None:
    _conc_term = rf"C^{{{_ord_val:g}}}"
else:
    _conc_term = "f(C)"        # unknown order

# Units of the rate constant: (mol/L)^(1-n)·s^-1.
if _ord_val == 1:
    _k_units = r"\ \mathrm{s^{-1}}"
elif _ord_val == 2:
    _k_units = r"\ \mathrm{L\,mol^{-1}\,s^{-1}}"
elif _ord_val == 0:
    _k_units = r"\ \mathrm{mol\,L^{-1}\,s^{-1}}"
else:
    _k_units = ""

# General rate law and the same with the rate constant substituted.
if _conc_term:
    st.latex(rf"-\frac{{dC}}{{dt}} = {_k_sym}\,{_conc_term}")
    st.latex(rf"-\frac{{dC}}{{dt}} = {rxn_k:.4g}{_k_units} \cdot {_conc_term}")
else:
    # Zero order: rate is independent of concentration.
    st.latex(rf"-\frac{{dC}}{{dt}} = {_k_sym}")
    st.latex(rf"-\frac{{dC}}{{dt}} = {rxn_k:.4g}{_k_units}")

# Characteristic reaction time.
if rxn_t_rxn > 0:
    st.latex(rf"t_{{rxn}} = {t_rxn:.4g}\ \text{{s}} \quad\text{{(specified directly)}}")
elif rxn_order in ("1", "pseudo-1"):
    st.latex(
        rf"t_{{rxn}} = \frac{{1}}{{{_k_sym}}} "
        rf"= \frac{{1}}{{{rxn_k:.4g}\ \text{{s}}^{{-1}}}} = {t_rxn:.4g}\ \text{{s}}"
    )
elif rxn_order in ("2", "pseudo-2"):
    st.latex(
        rf"t_{{rxn}} = \frac{{1}}{{{_k_sym}\,C_0}} "
        rf"= \frac{{1}}{{{rxn_k:.4g}\times{rxn_C0:.4g}}} = {t_rxn:.4g}\ \text{{s}}"
    )
else:
    st.latex(rf"t_{{rxn}} = {t_rxn:.4g}\ \text{{s}} \quad\text{{(order = {rxn_order})}}")

if _using_approximate:
    _assessment(
        "warning",
        "Approximate kinetics — the characteristic reaction time shown above "
        "is based on a proxy reaction, not measured data. All downstream conclusions "
        "(including whether the reaction is slow or fast) are only valid if the proxy "
        "kinetics match the true reaction.",
    )
else:
    _assessment(
        "ok",
        "Kinetics available — characteristic reaction time shown above.",
    )
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
        "hydrodynamics (see 🌀 Vessel Assessment and 📈 Vessel Comparison pages) "
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

    # ── Section outcome: assessment → recommendation ──
    _assessment(
        "warning",
        "Multi-phase system detected — interphase mass transfer may limit the "
        "observed rate.",
    )
    _recommendation(
        "Mass-transfer limitations are system-dependent — the remaining steps still "
        "assess micro/meso/macromixing and heat transfer, but keep in mind that "
        "interphase transport may dominate. Compute Da_GL and/or Da_LS for your "
        "specific reactor on the **🌀 Vessel Assessment** or **📈 Vessel Comparison** "
        "pages."
    )
else:
    _assessment(
        "ok",
        "Single liquid phase — interphase mass transfer is not a factor. Micro-, "
        "meso-, and macromixing may still affect the reaction.",
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
    # Micromixing sub-section
    with st.expander("🔬 Micromixing effects — theory & details", expanded=False):
        st.markdown(
            "At the smallest scales, the Engulfment rate *E* determines how "
            "fast reactant and surrounding fluid reach molecular homogeneity.  "
            "Key points:\n\n"
            "- Baldyga & Bourne (1999) distinguish **two** micromixing time "
            "scales (handbook Table 8.1):\n"
            "  - *Viscous-convective* engulfment: "
            r"$t_{E} \approx 17.3\,\sqrt{\nu/\varepsilon}$"
            " — the engulfment time that usually controls micromixing.\n"
            "  - *Viscous-diffusive* deformation/diffusion: "
            r"$t_{DS} \approx 0.5\,\sqrt{\nu/\varepsilon}\,\ln(Sc)$, "
            r"with Schmidt number $Sc = \nu/D$ — this becomes relevant for "
            "viscous or high-$Sc$ liquids, where molecular diffusion is slow.\n"
            "  The **larger** of the two scales controls the micromixing rate.\n"
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
            "with $A \\approx 1.2$ (handbook Table 8.1).  Here $\\Lambda_C$ is the "
            "**integral scale of turbulence** at the feed point, which in practice "
            "is approximated by the feed-pipe diameter $d_{feed}$, so this app "
            r"evaluates $\tau_S = 1.2\,(d_{feed}^2/\varepsilon)^{1/3}$ on the "
            "🌀 Vessel Assessment and 📈 Vessel Comparison pages.\n"
            "  - *Turbulent dispersion* of the plume by the mean flow: "
            r"$\tau_D = Q_{feed}/(\bar{u}\,D_t)$, "
            "where $Q_{feed}$ is the feed flow rate, $\\bar{u}$ the local mean "
            "velocity and $D_t$ the turbulent diffusivity.\n"
            "  The larger of the two scales controls the mesomixing rate."
        )

    # ── Section outcome: assessment → recommendation ──
    _assessment(
        "warning",
        "Competing reactions present — micro- and/or mesomixing sensitivity is "
        "likely, as both molecular-level and feed-point mixing can influence "
        "selectivity.",
    )
    _recommendation(
        "Use the **🅱️ Bourne Protocol** page to experimentally screen for micro- "
        "and mesomixing sensitivity by varying impeller speed (changes *ε*), feed "
        "time, and feed location. The Villermaux–Dushman (iodide–iodate) test "
        "reaction is a classic tool to decouple micro- and mesomixing effects "
        "(Bourne, 2003)."
    )
    _meso_sensitive = True
elif competing == "Not sure":
    _assessment(
        "caution",
        "Competing reactions uncertain — treat the system as potentially "
        "sensitive until confirmed.",
    )
    _recommendation(
        "Resolve the uncertainty by:\n\n"
        "- Running the reaction at two different impeller speeds **and** two "
        "different addition rates. A change in yield or impurity profile with "
        "speed implies micromixing; a change with feed time implies mesomixing.\n\n"
        "- Checking the literature or process-chemistry knowledge for known "
        "by-products or degradation pathways.\n\n"
        "- When in doubt, use the **🅱️ Bourne Protocol** to screen."
    )
    _meso_sensitive = True
else:
    _assessment(
        "ok",
        "No competing reactions identified — micro- and mesomixing "
        "(selectivity-related) are unlikely to be a factor. Macromixing (blend "
        "time) may still matter at scale.",
    )
    _meso_sensitive = False

# Semi-batch overrides: mesomixing is inherently relevant for fed-batch
# processes because the feed plume must be dispersed into the bulk before
# the reagent reacts (Baldyga & Bourne, 1999; Paul et al., 2004).
if is_semi_batch and not _meso_sensitive:
    _meso_sensitive = True
    _section_warning(
        "Semi-batch process — even without competing reactions, **mesomixing** "
        "(feed-plume dispersion) controls local concentration at the feed point. "
        "Selectivity, local heat release, and supersaturation can all be affected "
        "by feed rate, feed location, and local turbulence."
    )

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# STEP 4 – HEAT TRANSFER
# ══════════════════════════════════════════════════════════════════════════
st.header("4 · Heat Transfer Screening")

_has_enthalpy = rxn_delta_H != 0.0
_dT_ad = None  # adiabatic temperature rise (K), computed below when data allow

if not _has_enthalpy:
    _section_warning(
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
        _assessment(
            "caution",
            "Heat-transfer risk cannot be evaluated without ΔH — calorimetry is "
            "required to complete this section.",
        )
        _recommendation(
            "1. Perform **reaction calorimetry** (e.g. RC1, µRC, or Simular) to "
            "measure the heat of reaction *ΔH*.\n\n"
            "2. Update the reaction entry in the **🧪 Reaction Database**.\n\n"
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
        f"The selected reaction is "
        f"**{'exothermic' if rxn_delta_H < 0 else 'endothermic'}**:"
    )
    st.latex(rf"\Delta H = {rxn_delta_H:.1f}\ \mathrm{{kJ\,mol^{{-1}}}}")

    # ── Adiabatic temperature rise (ΔT_ad) — the proper thermal criterion ──
    # ΔH per mole alone does not determine thermal-runaway risk: a large ΔH at
    # low concentration may be benign, while a modest ΔH at high concentration
    # can be hazardous. ΔT_ad = |ΔH|·C0 / (ρ·Cp) is the temperature *rise* the
    # batch would undergo with no cooling (final temperature ≈ T0 + ΔT_ad), and
    # is a key input to the Stoessel criticality classification.
    _rho_cp_auto = None  # volumetric heat capacity, kJ/(m³·K)
    if rxn_solvent and is_known_solvent(rxn_solvent):
        try:
            _sp_heat = get_properties(rxn_solvent, rxn_T, 1.0)
            _rho_cp_auto = _sp_heat["rho_kg_m3"] * _sp_heat["Cp_J_per_kgK"] / 1000.0
        except Exception:
            _rho_cp_auto = None

    with st.expander("Adiabatic temperature rise (ΔT_ad) — refine the estimate", expanded=True):
        st.caption(
            "ΔT_ad = |ΔH|·C₀ / (ρ·Cp) is the temperature *rise* the batch would "
            "undergo if all reaction heat were retained (no cooling); the final "
            "temperature would be roughly T₀ + ΔT_ad. This rise — not ΔH per mole — "
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
            st.latex(
                rf"\Delta T_{{ad}} = \frac{{|\Delta H|\,C_0}}{{\rho C_p}} "
                rf"\approx {_dT_ad:.0f}\ \mathrm{{K}}"
            )
            if _dT_ad >= 200:
                st.error(
                    "🔴 **Very high** — loss of cooling could "
                    "drive a runaway; the risk of triggering secondary (decomposition) "
                    "reactions must be assessed (compare the MTSR against the "
                    "decomposition onset T_D24, per the Stoessel criticality framework)."
                )
            elif _dT_ad >= 50:
                st.error(
                    "🔴 **High** — strong cooling and feed-rate "
                    "control are required; quantify Q_gen vs Q_cool at scale."
                )
            elif _dT_ad >= 20:
                st.warning(
                    "🟡 **Moderate** — manageable with adequate "
                    "cooling, but verify the heat balance at production scale."
                )
            else:
                st.success(
                    "🟢 **Low** — thermal runaway is unlikely, "
                    "though local hot spots at the feed point may still occur."
                )
        else:
            st.caption("Enter C₀ and ρ·Cp to estimate ΔT_ad.")

    # Heuristic classification (traffic-light severity)
    if abs_dH >= 100:
        _heat_class = "highly exothermic"
        _heat_kind = "critical"
    elif abs_dH >= 50:
        _heat_class = "moderately exothermic"
        _heat_kind = "warning"
    elif abs_dH >= 20:
        _heat_class = "mildly exothermic"
        _heat_kind = "caution"
    else:
        _heat_class = "low exothermicity"
        _heat_kind = "ok"

    # Heat sensitivity is flagged on either a large molar enthalpy or (more
    # rigorously) a significant adiabatic temperature rise.
    _heat_flag = abs_dH >= 50 or (_dT_ad is not None and _dT_ad >= 50)

    # Background detail (shown before the section outcome)
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
    else:
        st.markdown(
            "The reaction enthalpy is modest.  Heat transfer is **unlikely "
            "to be limiting** in most reactor configurations, but may become "
            "relevant at very large scale or in poorly cooled vessels.  "
            "Consider running a heat balance if scaling beyond pilot scale."
        )

    # Flag known highly-exothermic reaction types
    _rxn_type = str(rxn.get("type", "")).lower()
    _known_hot = ["grignard", "nitration", "sulfonation", "diazotization",
                   "polymerization", "hydrogenation", "oxidation"]
    _flagged = [t for t in _known_hot if t in _rxn_type]

    # ── Section outcome: assessment → recommendation → warning ──
    st.latex(rf"|\Delta H| = {abs_dH:.1f}\ \mathrm{{kJ\,mol^{{-1}}}}")
    _assessment(_heat_kind, f"**{_heat_class.title()}** reaction.")
    _heat_sensitive = _heat_flag
    if _heat_flag:
        _recommendation(
            "Run the heat balance on the **🌀 Vessel Assessment** or "
            "**📈 Vessel Comparison** page to quantify Q_gen vs Q_cool for your "
            "specific reactor(s)."
        )
    if _flagged:
        _section_warning(
            f"Reaction type **{rxn.get('type', '')}** is commonly associated with "
            "significant exothermicity. Heat-transfer assessment is strongly "
            "recommended regardless of the reported ΔH magnitude."
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

st.markdown("**Your reaction time:**")
st.latex(rf"t_{{rxn}} = {t_rxn:.4g}\ \mathrm{{s}}")

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

# ── Section outcome: assessment → recommendation → warning ──
if t_rxn < 0.1:
    _assessment(
        "critical",
        "**Very fast reaction** — micromixing-sensitive in "
        "most reactor configurations. Local turbulent energy dissipation near the "
        "impeller determines the effective mixing rate; feed location and impeller "
        "tip speed are critical parameters.",
    )
    _micro_likely = True
elif t_rxn < 1.0:
    _assessment(
        "warning",
        "**Fast reaction** — micromixing likely relevant in "
        "larger vessels where local ε at the feed point decreases. Confirm with "
        "Damköhler analysis on the 🌀 Vessel Assessment page.",
    )
    _micro_likely = True
elif t_rxn < 10:
    _assessment(
        "caution",
        "**Moderate reaction** — micromixing is less likely "
        "to dominate, but macromixing (blend time) could be relevant in larger "
        "vessels. Check blend time relative to t_rxn.",
    )
    _micro_likely = False
else:
    _assessment(
        "ok",
        "**Slow reaction** — mixing is unlikely to limit the "
        "reaction in well-agitated vessels. Only very large, poorly mixed tanks "
        "might approach Da ~ 1.",
    )
    _micro_likely = False

if t_rxn < 30:
    _recommendation(
        "Compute Damköhler numbers for your specific reactor(s) on the "
        "**🌀 Vessel Assessment** page or compare multiple vessels on the "
        "**📈 Vessel Comparison** page."
    )

if is_semi_batch:
    _section_warning(
        "Semi-batch process — three mixing scales act in series at the feed point:\n\n"
        "1. **Macromixing** (θ₉₅) — bulk blending of added reagent into the vessel\n\n"
        "2. **Mesomixing** ($t_{meso}$) — turbulent dispersion of the feed plume; "
        "controls local concentration before molecular homogenisation\n\n"
        "3. **Micromixing** ($t_E$) — engulfment at the Kolmogorov scale\n\n"
        "To distinguish these experimentally, run the full **🅱️ Bourne Protocol**: "
        "vary impeller speed (micromixing), feed rate / addition time (mesomixing), "
        "and feed location (meso- and macromixing)."
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

# Actual KPIs the imported Bourne Protocol flagged as sensitive, so findings
# cite the responses that were really impacted instead of a generic
# "yield/impurity" statement.
_bourne_kpi_phrase = _bourne_sensitive_kpi_phrase(_bourne_parsed)
# Tests still needed to resolve the controlling scale (so recommendations don't
# ask to re-run tests that were already completed).
_bourne_rem_tests = _bourne_remaining_tests(_bourne_done_tests)
_bourne_rem_phrase = _fmt_tests(_bourne_rem_tests)
_bourne_rem_action = (
    f"complete {_bourne_rem_phrase} of the Bourne Protocol" if _bourne_rem_phrase
    else "re-run the Bourne Protocol decision tree"
)

# Bourne pre-screening
if _bourne_sensitive is True:
    if _bourne_mechanisms:
        findings.append((
            "Bourne pre-screen",
            "🔴 Mixing sensitivity confirmed",
            f"Experimental pre-screen showed {_bourne_kpi_phrase} changed with "
            "mixing conditions. Identified controlling scale(s): "
            f"{', '.join(_bourne_mechanisms)}.",
        ))
    else:
        findings.append((
            "Bourne pre-screen",
            "🔴 Mixing sensitivity confirmed",
            f"Experimental pre-screen showed {_bourne_kpi_phrase} changed with "
            f"mixing conditions. Controlling scale not yet identified — {_bourne_rem_action} "
            "to pinpoint micro-, meso-, or macromixing.",
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
            "but some items require verification at scale. To pinpoint the controlling scale, "
            f"{_bourne_rem_action}."
        )
        st.error(
            "🔴 **Mixing sensitivity confirmed** — the Bourne pre-screen shows a sensitivity "
            "is present. The theoretical assessment did not flag a specific mechanism as "
            "*likely*, but some items require verification at scale (see below). To pinpoint "
            f"the controlling scale, {_bourne_rem_action}."
        )
    else:
        _verdict = (
            "Mixing sensitivity confirmed -- the Bourne pre-screen shows an experimental "
            "sensitivity, although the theoretical assessment flagged no mechanism. Revisit "
            f"the inputs (kinetics, phases, feed strategy) and {_bourne_rem_action} "
            "to identify the controlling scale."
        )
        st.error(
            "🔴 **Mixing sensitivity confirmed** — the Bourne pre-screen shows an experimental "
            "sensitivity even though the theoretical assessment flagged no mechanism. Revisit "
            f"the inputs (kinetics, phases, feed strategy) and {_bourne_rem_action} "
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

# Surface the imported Bourne Protocol experimental findings alongside the
# theoretical assessment so the summary reflects the actual test results.
if _bourne_test_details:
    st.subheader("Bourne Protocol Experimental Findings")
    if _bourne_meta:
        st.caption("  •  ".join(f"**{_k}:** {_v}" for _k, _v in _bourne_meta.items()))
    st.dataframe(
        pd.DataFrame(_bourne_test_details), width='stretch', hide_index=True
    )
    if _bourne_mechanisms:
        st.caption(
            f"Controlling scale identified by the Bourne Protocol: "
            f"**{_bourne_mechanisms[0]}**."
        )

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
            "on the 🌀 Vessel Assessment page.",
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

# Bourne confirmed a sensitivity but did not resolve the controlling scale:
# recommend only the test(s) that still need to be run.
if _bourne_sensitive is True and not _bourne_mechanisms and _bourne_rem_tests:
    _steps.append({
        "Area": "Resolve controlling scale",
        "Action": f"Complete Bourne {_fmt_tests(_bourne_rem_tests)} "
                  f"({_fmt_test_purposes(_bourne_rem_tests)}) to identify whether "
                  "micro-, meso-, or macromixing controls, then re-export and re-import.",
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
                  "reactor on the 🌀 Vessel Assessment page.",
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
                  "the 🌀 Vessel Assessment or 📈 Vessel Comparison pages.",
    })

if _has_enthalpy and _heat_sensitive:
    _steps.append({
        "Area": "Heat transfer",
        "Action": "Run a full heat balance (🔥 button on 🌀 Vessel Assessment or "
                  "📈 Vessel Comparison) to evaluate Q_gen vs Q_cool.",
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
                  "a quick check on the 🌀 Vessel Assessment page for completeness.",
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
                "bourne_meta": _bourne_meta,
                "bourne_tests": _bourne_test_details,
                "bourne_mechanism": _bourne_mechanisms[0] if _bourne_mechanisms else "",
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
