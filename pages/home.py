"""Home page content for Mixing Lab."""

import streamlit as st

st.title("🔬 Mixing Lab")
st.subheader("Hydrodynamics & Mixing-Sensitivity Assessment Tool")

st.markdown("""
Welcome to **Mixing Lab** – an interactive tool for assessing hydrodynamics and
mixing sensitivities across different reactor scales.

### What you can do

| Page | Purpose |
|------|---------|
| **⚗️ Reactor Database** | Browse, add, and edit reactor geometries from lab to manufacturing scale |
| **🧪 Reaction Database** | Manage reactions with kinetic parameters (rate constants, orders, temps) |
| **💧 Fluid Database** | Define or import fluid properties (density, viscosity, diffusivity) |
| **❇️ Particle Database** | Define solid-particle properties (density, PSD, shape factor) |
| **🧭 Sensitivity Protocol** | Interactive decision tree to assess which mixing mechanisms may limit your reaction at scale |
| **⚙️ Mixing Sensitivity** | Calculate hydrodynamic parameters and Damköhler numbers to estimate mixing sensitivity |
| **🧫 Bourne Protocol** | Step-by-step Bourne mixing-sensitivity screening protocol (Sarafinas modification) |
| **📊 Reactor Comparison** | Side-by-side comparison of hydrodynamics across selected reactors |
| **📋 Recorded Results** | Save, review, and export results for specific reactor/reaction combinations |
| **📐 Equations** | Reference pages for all correlations and equations used (hydrodynamics, mixing, mass transfer, heat balance, particles, scale-up, ROM) |

---

### Quick-start guide

1. **Populate your databases** – Add reactors, reactions, fluids, and particles via the Database pages (or import reactors from the Admin page).
2. **Screen for mixing sensitivity** – Use the *Sensitivity Protocol* to walk through a structured decision tree and identify which mixing mechanisms may limit your process.
3. **Run the Bourne Protocol** – If mixing sensitivity is suspected, use the *Bourne Protocol* page to experimentally determine whether micro-, meso-, or macromixing controls your process.
4. **Analyse a single system** – On the *Mixing Sensitivity* page, select a reactor, fluid, reaction, and operating conditions to compute hydrodynamic parameters, Damköhler numbers, and mixing-sensitivity flags.
5. **Compare across scales** – Use the *Reactor Comparison* page to evaluate operating envelopes and scale-up impacts across multiple reactors side-by-side.
6. **Record & export** – Save results from any analysis page to *Recorded Results* for documentation and comparison.
""")
