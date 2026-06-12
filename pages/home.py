"""Home page content for Mixing Lab."""

import streamlit as st

st.title("🔬 Mixing Lab")
st.badge("v1.0.0", color="blue")
st.subheader("Hydrodynamics & Mixing-Sensitivity Assessment Tool")

st.markdown("""
Welcome to **Mixing Lab** – an interactive tool for assessing reaction mixing sensitivities,
            hydrodynamics and scale dependencies.

### What you can do

| Group | Page | Purpose |
|-------|------|---------|
| **Databases** | **⚗️ Reactor Database** | Browse, add, and edit reactor geometries with 2D/3D schematics; use existing reactors as templates |
| | **🧪 Reaction Database** | Manage reactions with kinetic parameters (rate constants, orders, concentrations, ΔH) |
| | **💧 Fluid Database** | Browse built-in solvents, explore temperature-dependent properties, and manage custom fluids |
| | **❇️ Particle Database** | Define solid-particle properties (density, PSD, shape factor) |
| **Analysis** | **🧭 Reaction Sensitivity Protocol** | Interactive decision tree to assess which mixing mechanisms may limit your reaction at scale |
| | **💎 Crystallization Sensitivity Protocol** | Decision tree for crystallization mixing sensitivity (nucleation, growth, anti-solvent) |
| | **🌀 Mixing Assessment** | Reactor–reaction–fluid analysis: compute hydrodynamics, Damköhler numbers, heat balance, solids suspension & operating envelopes |
| | **🅱️ Bourne Protocol** | Step-by-step Bourne/Sarafinas mixing-sensitivity screening – vary speed, feed rate, and feed location to identify the controlling mixing scale |
| | **📈 Reactor Comparison** | Side-by-side comparison of hydrodynamics and operating envelopes across multiple reactors |
| **Tools** | **🔄 Unit Converter** | Quick conversion between common engineering units |
| | **🔥 Heat Transfer** | Estimate heating/cooling times, compute overall U from first principles, visualize batch temperature profiles, and compare heat-transfer media |
| **Results** | **📋 Recorded Results** | Save, review, filter, and export results for specific reactor/reaction combinations |
| **Reference** | **📐 Equations Reference** | Reference pages for all correlations (hydrodynamics, mixing, mass transfer, heat balance, particles, scale-up, ROM, crystallisation) |
| **Admin** | **🛠️ Admin Tools** | Import reactors from CSV/Excel, convert transposed master files to row-per-reactor format |
| | **🔧 ROM / Correlation Fitting** | Upload measured or CFD data, fit reduced-order models, and register correlations for use throughout Mixing Lab |

---

### Quick-start guide

1. **Populate your databases** – Add reactors, reactions, fluids, and particles via the Database pages (or import reactors from the Admin page). Use the template feature to quickly duplicate reactor geometries.
2. **Screen for mixing sensitivity** – Use the *Reaction Sensitivity Protocol* or *Crystallization Sensitivity Protocol* to walk through a structured decision tree and identify which mixing mechanisms may limit your process.
3. **Run the Bourne Protocol** – If mixing sensitivity is suspected, use the *Bourne Protocol* page to experimentally determine whether micro-, meso-, or macromixing controls your process.
4. **Analyse a single system** – On the *Mixing Assessment* page, select a reactor, fluid, reaction, and operating conditions. Confirm each step to compute hydrodynamic parameters, Damköhler numbers, heat balance, solid-particle suspension, and full operating envelopes.
5. **Compare across scales** – Use the *Reactor Comparison* page to evaluate operating envelopes and scale-up impacts across multiple reactors side-by-side.
6. **Generate a report** – Most analysis pages (*Mixing Assessment*, *Reactor Comparison*, the sensitivity protocols, and *Bourne Protocol*) can export a branded PDF with metrics, sensitivity assessment, and envelope charts directly from the page.
7. **Record & export** – Save results from any analysis page to *Recorded Results* for documentation and comparison.
""")
