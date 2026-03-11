"""Equations – Mixing & Damköhler Numbers."""

import streamlit as st

st.title("📐 Mixing & Damköhler Equations")

# ── Blend Time ───────────────────────────────────────────────────────────
st.header("6. Macro-Blend Time (95 %)")
st.latex(r"\theta_{95} = \frac{5.2 \, V}{N_Q \, N \, D^3}")
st.markdown("""
Grenville correlation for turbulent blending in a baffled stirred tank.

This represents the time to achieve 95 % uniformity on the bulk (macro) scale.

**Reference:** Grenville, R. K. (1992). *Blending of Viscous Newtonian and Non-Newtonian Fluids.* Ph.D. thesis, Cranfield University.
""")

# ── Micromixing Time ─────────────────────────────────────────────────────
st.header("7. Engulfment Micromixing Time")
st.latex(r"t_E = 17.3 \left( \frac{\nu}{\varepsilon} \right)^{1/2}")
st.markdown("""
Baldyga & Bourne engulfment model – characterises the time for mixing at the
smallest (molecular) scales.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\\nu$ | Kinematic viscosity ($\\mu / \\rho$) | m²/s |
| $\\varepsilon$ | Mean energy dissipation rate (P/V) | W/m³ or W/kg |

**Reference:** Baldyga, J. and Bourne, J. R. (1999). *Turbulent Mixing and Chemical Reactions.* Wiley.
""")

# ── Damköhler Numbers ───────────────────────────────────────────────────
st.header("11. Damköhler Numbers")

st.subheader("Macro Damköhler Number")
st.latex(r"Da_{macro} = \frac{\theta_{blend}}{t_{rxn}}")

st.subheader("Micro Damköhler Number")
st.latex(r"Da_{micro} = \frac{t_E}{t_{rxn}}")

st.markdown("""
| $Da$ range | Interpretation |
|-----------|----------------|
| $Da \\ll 0.01$ | Reaction-limited – insensitive to mixing |
| $0.01 < Da < 0.1$ | Likely insensitive |
| $0.1 < Da < 1$ | Potentially sensitive |
| $1 < Da < 10$ | Mixing-sensitive |
| $Da \\gg 10$ | Strongly mixing-limited |

**Reference:** Paul, E.L., Atiemo-Obeng, V.A. and Kresta, S.M. (2004). *Handbook of Industrial Mixing.* Wiley.
""")

# ── Gas-Liquid Damköhler ─────────────────────────────────────────────────
st.header("11b. Gas-Liquid Damköhler Number")

st.latex(r"Da_{GL} = \frac{1}{k_L a \; t_{rxn}} = \frac{t_{\text{transfer}}}{t_{rxn}}")

st.markdown(r"""
Compares the characteristic **gas-liquid mass-transfer time** $\;(1/k_L a)\;$ to
the reaction time.  The effective $k_L a$ is the larger of the sparged or
free-surface value.

| $Da_{GL}$ range | Interpretation |
|----------------|----------------|
| $Da_{GL} \ll 0.01$ | Mass transfer is fast — reaction-limited |
| $0.01 < Da_{GL} < 0.1$ | Likely insensitive to G-L transfer |
| $0.1 < Da_{GL} < 1$ | Potentially transfer-limited |
| $1 < Da_{GL} < 10$ | Transfer-limited |
| $Da_{GL} \gg 10$ | Strongly transfer-limited |

When $k_L a = 0$ (no gas phase present), $Da_{GL}$ is set to zero (not
applicable).

The interpretation follows the same convention as $Da_{macro}$ and
$Da_{micro}$: higher values indicate the system is **more limited** by
the transport process relative to reaction.
""")

# ── Reaction time ────────────────────────────────────────────────────────
st.header("12. Characteristic Reaction Time")

st.subheader("First-order or pseudo-first-order")
st.latex(r"t_{rxn} = \frac{\ln 2}{k}")

st.subheader("Second-order")
st.latex(r"t_{rxn} = \frac{1}{k \, C_0}")

st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $k$ | Rate constant | 1/s (1st order) or L/(mol·s) (2nd order) |
| $C_0$ | Initial concentration | mol/L |
""")

# ── Local micromixing time ───────────────────────────────────────────────
st.header("14. Local Micromixing Time (Impeller Zone)")
st.latex(r"t_{E,\text{local}} = 17.3 \left( \frac{\nu}{\varepsilon_{max}} \right)^{1/2}")
st.markdown("""
This is the same Baldyga & Bourne engulfment model as Equation 7, but evaluated
using the **local maximum** energy dissipation rate $\\varepsilon_{max}$ near the impeller
instead of the volume-averaged $\\varepsilon$.

Because $\\varepsilon_{max}$ can be 1–2 orders of magnitude larger than $\\varepsilon_{avg}$,
$t_{E,\\text{local}}$ is correspondingly shorter.  This estimate is relevant when
reagents are **fed directly into the impeller discharge zone**.

**Reference:** Baldyga, J. and Bourne, J. R. (1999). *Turbulent Mixing and Chemical Reactions.* Wiley.
""")
