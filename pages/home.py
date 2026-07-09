"""Home page content for Mixing Lab."""

import streamlit as st

st.title("Mixing Lab")
st.caption("Hydrodynamics & Mixing-Sensitivity Assessment Tool")
st.caption("""
Welcome to **Mixing Lab** - an interactive tool for assessing reaction mixing sensitivities,
            hydrodynamics and scale dependencies.""")

st.markdown("""
### App Contents:

| Group | Page | Purpose |
|-------|------|---------|
| **Databases** | **⚗️ Vessels** | Browse, add, and edit reactor geometries with 2D/3D schematics; use existing reactors as templates |
| | **🧪 Reactions** | Manage reactions with kinetic parameters (rate constants, orders, concentrations, ΔH) |
| | **💧 Fluids** | Browse built-in solvents, explore temperature-dependent properties, and manage custom fluids |
| | **❇️ Particles** | Define solid-particle properties (density, PSD, shape factor) |
| **Mixing Analysis** | **🅱️ Bourne Protocol** | Step-by-step Bourne/Sarafinas mixing-sensitivity screening – vary agitation speed, feed rate, and feed location to identify the controlling mixing scale |
| | **🧭 Reaction Sensitivity Protocol** | Interactive workflow to identify potential mixing sensitivities of a reaction process and determine which mixing mechanisms may limit the reaction at scale |
| | **🚧 Crystallization Sensitivity Protocol (WIP)** | Interactive workflow to identify potential mixing sensitivities of a crystallization process and determine which mixing mechanisms may limit the reaction at scale |
| | **🌀 Vessel Assessment** | Reactor–reaction–fluid analysis: compute hydrodynamics, Damköhler numbers, heat balance, solids suspension & operating envelopes |
| | **⚖️ Vessel Comparison** | Side-by-side comparison of hydrodynamics and operating envelopes across multiple reactors |
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
4. **Analyse a single system** – On the *Vessel Assessment* page, select a reactor, fluid, reaction, and operating conditions. Confirm each step to compute hydrodynamic parameters, Damköhler numbers, heat balance, solid-particle suspension, and full operating envelopes.
5. **Compare across scales** – Use the *Vessel Comparison* page to evaluate operating envelopes and scale-up impacts across multiple reactors side-by-side.
6. **Generate a report** – Most analysis pages (*Vessel Assessment*, *Vessel Comparison*, the sensitivity protocols, and *Bourne Protocol*) can export a branded PDF with metrics, sensitivity assessment, and envelope charts directly from the page.
7. **Record & export** – Save results from any analysis page to *Recorded Results* for documentation and comparison.
""")

st.divider()
st.badge("v1.0.0", color="blue")
st.badge("Updated: July 9, 2026", color="green")