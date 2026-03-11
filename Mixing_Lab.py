"""
Mixing Lab – Hydrodynamics & Mixing-Sensitivity Assessment Tool
================================================================
Main entry point.  Run with:  streamlit run mixing_lab.py
"""

import streamlit as st
from utils.sidebar import render_sidebar

st.set_page_config(
    page_title="Mixing Lab",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

render_sidebar()

# ── Home page content (shown when user clicks "Mixing Lab" in nav) ────────
st.title("🔬 Mixing Lab")
st.subheader("Hydrodynamics & Mixing-Sensitivity Assessment Tool")

st.markdown("""
Welcome to **Mixing Lab** – an interactive tool for assessing hydrodynamics and
mixing sensitivities across different reactor scales.

### What you can do

| Page | Purpose |
|------|---------|
| **Reactor Database** | Browse, add, and edit reactor geometries from lab to manufacturing scale |
| **Reaction Database** | Manage reactions with kinetic parameters (rate constants, orders, temps) |
| **Fluid Database** | Define or import fluid properties (density, viscosity, diffusivity) |
| **Particle Database** | Define solid-particle properties (density, PSD, shape factor) for suspension calculations |
| **Mixing Sensitivity** | Calculate hydrodynamic parameters and Damköhler numbers to estimate mixing sensitivity |
| **Bourne Protocol** | Step-by-step Bourne mixing-sensitivity screening protocol (Sarafinas modification) |
| **Reactor Comparison** | Side-by-side comparison of hydrodynamics across selected reactors |
| **Recorded Results** | Save, review, and export results for specific reactor/reaction combinations |
| **Equations** | Reference for all correlations and equations used |

---

### Quick-start guide

1. **Verify your databases** – Check that reactor, reaction, and fluid databases are populated.
2. **Select a system** – Choose a reactor, reaction, and fluid on the *Mixing Sensitivity* page.
3. **Review results** – Examine Damköhler numbers and sensitivity flags.
4. **Compare scales** – Use the *Reactor Comparison* page to evaluate scale-up impacts.
5. **Record & export** – Save important results for documentation.
""")
