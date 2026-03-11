"""Equations – Mass Transfer & kLa."""

import streamlit as st

st.title("📐 Mass Transfer Equations")

# ── kLa – sparged ────────────────────────────────────────────────────────
st.header("13. kLa – Volumetric Mass-Transfer Coefficient (Gas-Liquid)")
st.latex(r"k_L a = C_1 \left( \frac{P}{V} \right)^{C_2} v_s^{C_3}")
st.markdown("""
Van 't Riet (1979) correlation for aerated stirred tanks.

| System | $C_1$ | $C_2$ | $C_3$ |
|--------|-------|-------|-------|
| Coalescing (pure liquids) | 0.026 | 0.4 | 0.5 |
| Non-coalescing (electrolytes) | 0.002 | 0.7 | 0.2 |

| Symbol | Description | Units |
|--------|-------------|-------|
| $P/V$ | Gassed specific power input | W/m³ |
| $v_s$ | Superficial gas velocity | m/s |
| $k_L a$ | Volumetric mass-transfer coefficient | 1/s |

**Note:** $v_s = 0$ (no sparging) → $k_L a = 0$.

**Reference:** van 't Riet, K. (1979). Review of measuring methods and results in non-viscous
gas-liquid mass transfer in stirred vessels. *Ind. Eng. Chem. Process Des. Dev.*, 18(3), 357–364.
""")

# ── kLa – headspace / free-surface ───────────────────────────────────────
st.header("13b. kLa – Free-Surface (Headspace) Mass Transfer")
st.latex(r"k_L = 0.4 \; D_{mol}^{\,1/2} \left( \frac{\varepsilon}{\nu} \right)^{1/4}")
st.latex(r"a_{\text{surface}} = \frac{\pi / 4 \; D_T^2}{V}")
st.latex(r"k_L a_{\text{surface}} = k_L \cdot a_{\text{surface}}")
st.markdown("""
Lamont & Scott (1970) small-eddy surface-renewal model.  Applicable when there
is **no gas sparging** and mass transfer occurs only through the flat liquid
surface exposed to the headspace gas.

| Symbol | Description | Units |
|--------|-------------|-------|
| $D_{mol}$ | Molecular diffusivity of dissolved gas | m²/s |
| $\\varepsilon$ | Mean specific energy dissipation rate | W/kg (= m²/s³) |
| $\\nu$ | Kinematic viscosity | m²/s |
| $D_T$ | Tank inside diameter | m |
| $V$ | Liquid volume | m³ |
| $a_{\\text{surface}}$ | Specific interfacial area (free surface) | 1/m |

**Typical D_mol values for common gases in water at 25 °C:**

| Gas | $D_{mol}$ (m²/s) |
|-----|------------------|
| O₂  | 2.1 × 10⁻⁹ |
| CO₂ | 1.9 × 10⁻⁹ |
| N₂  | 1.9 × 10⁻⁹ |
| H₂  | 4.5 × 10⁻⁹ |

**Reference:** Lamont, J. C. and Scott, D. S. (1970). An eddy cell model of mass
transfer into the surface of a turbulent liquid. *AIChE J.*, 16(4), 513–519.
""")
