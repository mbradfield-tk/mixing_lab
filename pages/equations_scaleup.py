"""Equations – Scale-Up Rules & Solvent Properties."""

import streamlit as st

st.title("📐 Scale-Up & Solvent Properties")

# ── Scale-up rules ───────────────────────────────────────────────────────
st.header("Common Scale-Up Rules")

st.markdown("Given geometric similarity ($D/T$ = const), the large-scale speed $N_L$ is related to the small-scale $N_S$:")

st.subheader("Constant tip speed")
st.latex(r"N_L = N_S \frac{D_S}{D_L}")

st.subheader("Constant P/V")
st.latex(r"N_L = N_S \left( \frac{D_S}{D_L} \right)^{2/3}")

st.subheader("Constant Re")
st.latex(r"N_L = N_S \left( \frac{D_S}{D_L} \right)^{2}")

st.markdown("""
---
**General references:**
- Paul, E.L., Atiemo-Obeng, V.A. and Kresta, S.M. (2004). *Handbook of Industrial Mixing.* Wiley-Interscience.
- Baldyga, J. and Bourne, J.R. (1999). *Turbulent Mixing and Chemical Reactions.* Wiley.
- Grenville, R.K. and Nienow, A.W. (2004). Blending of miscible liquids. In *Handbook of Industrial Mixing*, Ch. 9.
""")

# ── Solvent Properties (T-dependent) ────────────────────────────────────
st.header("Temperature-Dependent Solvent Properties")
st.markdown("""
The Fluid Database's *Solvent Properties (T)* tab provides physical properties
at any liquid-phase temperature for common pharmaceutical solvents, using the
correlations below.
""")

st.subheader("Density")
st.latex(r"\rho(T) = \rho_{25} + \frac{d\rho}{dT} \cdot (T - 25) \quad [\mathrm{kg/m^3}]")
st.markdown("Linear fit anchored at the known 25 °C value.  $T$ in °C, slope $d\\rho/dT$ (typically negative) from published data.")

st.subheader("Dynamic Viscosity (Arrhenius)")
st.latex(r"\mu(T) = \mu_{25} \, \exp\!\left[ \frac{E_a}{R} \left( \frac{1}{T_K} - \frac{1}{298.15} \right) \right]")
st.markdown("""
$T_K$ = temperature in Kelvin, $E_a$ = activation energy for viscous flow
(J/mol), $R$ = 8.314 J/(mol·K).  Anchored at the known $\\mu_{25}$ so the
reference value is recovered exactly.  $E_a$ values from published data
(Perry's, CRC Handbook, DIPPR).
""")

st.subheader("Surface Tension")
st.latex(r"\sigma(T) = \sigma_{25} + \frac{d\sigma}{dT} \cdot (T - 25) \quad [\mathrm{N/m}]")
st.markdown("Linear fit anchored at the known 25 °C value.  Clamped to $\\sigma \\ge 0$.")

st.subheader("Molecular Diffusivity (Stokes-Einstein Scaling)")
st.latex(r"D(T) = D_{\mathrm{ref}} \; \frac{T_K}{298.15} \; \frac{\mu_{\mathrm{ref}}}{\mu(T)}")
st.markdown("""
Scaled from a reference value $D_{\\mathrm{ref}}$ measured at 25 °C using the
Stokes-Einstein relation ($D \\propto T / \\mu$).  This captures the dominant
temperature dependence through the viscosity term.

**Data sources:** Yaws' Handbook, Perry's Chemical Engineers' Handbook (9th ed.),
CRC Handbook of Chemistry and Physics, DIPPR correlations.
""")
