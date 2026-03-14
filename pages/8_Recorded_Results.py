"""
Page 7 – Recorded Results
==========================
View, filter, and export saved mixing-sensitivity results.
"""

import streamlit as st

import pandas as pd
import pathlib

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
RESULTS_CSV = DATA_DIR / "recorded_results.csv"


def _load_results() -> pd.DataFrame:
    if RESULTS_CSV.exists():
        return pd.read_csv(RESULTS_CSV)
    return pd.DataFrame()


if "recorded_results" not in st.session_state:
    st.session_state.recorded_results = _load_results()

st.title("📋 Recorded Results")

df = st.session_state.recorded_results

if df.empty:
    st.info("No results recorded yet.  Use the **Mixing Sensitivity Workflow** page to compute and save results.")
    st.stop()

# ── Filters ──────────────────────────────────────────────────────────────
st.header("Filter Results")
col1, col2, col3 = st.columns(3)

with col1:
    if "reactor" in df.columns:
        reactors_filt = st.multiselect("Reactor", df["reactor"].dropna().unique().tolist())
    else:
        reactors_filt = []
with col2:
    if "reaction" in df.columns:
        reactions_filt = st.multiselect("Reaction", df["reaction"].dropna().unique().tolist())
    else:
        reactions_filt = []
with col3:
    if "fluid" in df.columns:
        fluids_filt = st.multiselect("Fluid", df["fluid"].dropna().unique().tolist())
    else:
        fluids_filt = []

filtered = df.copy()
if reactors_filt:
    filtered = filtered[filtered["reactor"].isin(reactors_filt)]
if reactions_filt:
    filtered = filtered[filtered["reaction"].isin(reactions_filt)]
if fluids_filt:
    filtered = filtered[filtered["fluid"].isin(fluids_filt)]

# ── Display ──────────────────────────────────────────────────────────────
st.header("Results Table")

format_dict = {}
for col in filtered.columns:
    if filtered[col].dtype in ("float64", "float32"):
        format_dict[col] = "{:.4g}"

st.dataframe(filtered.style.format(format_dict), use_container_width=False)

st.caption(f"Showing {len(filtered)} of {len(df)} records.")

# ── Assessment summary ───────────────────────────────────────────────────
if "Assessment" in filtered.columns and not filtered.empty:
    st.header("Assessment Summary")

    # Count occurrences of key phrases
    assessments = filtered["Assessment"].astype(str)
    n_limited = assessments.str.contains("mixing-limited|Mixing-sensitive", case=False).sum()
    n_potential = assessments.str.contains("Potentially sensitive", case=False).sum()
    n_safe = len(assessments) - n_limited - n_potential

    c1, c2, c3 = st.columns(3)
    c1.metric("🟢 Reaction-limited", n_safe)
    c2.metric("🟡 Potentially sensitive", n_potential)
    c3.metric("🔴 Mixing-sensitive / limited", n_limited)

# ── Export ────────────────────────────────────────────────────────────────
st.header("Export")
st.download_button(
    "⬇️ Download results (CSV)",
    data=filtered.to_csv(index=False).encode("utf-8"),
    file_name="mixing_lab_results.csv",
    mime="text/csv",
)

# ── Delete ────────────────────────────────────────────────────────────────
st.header("Manage Records")

if st.button("🗑️ Clear all recorded results", type="secondary"):
    st.session_state.recorded_results = pd.DataFrame()
    if RESULTS_CSV.exists():
        RESULTS_CSV.unlink()
    st.success("All recorded results cleared.")
    st.rerun()
