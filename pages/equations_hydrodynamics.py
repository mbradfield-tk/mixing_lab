"""Equations – Hydrodynamics & Shear."""

import streamlit as st

st.title("📐 Hydrodynamics & Shear Equations")

# ── Impeller Reynolds Number ─────────────────────────────────────────────
st.header("1. Impeller Reynolds Number")
st.latex(r"Re = \frac{\rho \, N \, D^2}{\mu}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\rho$ | Fluid density | kg/m³ |
| $N$ | Impeller rotational speed | rev/s |
| $D$ | Impeller diameter | m |
| $\\mu$ | Dynamic viscosity | Pa·s |

**Regime classification:** Re < 10 laminar; 10 < Re < 10 000 transitional; Re > 10 000 turbulent.
""")

# ── Power Number & Power ─────────────────────────────────────────────────
st.header("2. Power Number & Impeller Power")
st.latex(r"P = N_p \, \rho \, N^3 \, D^5")
st.markdown("""
$N_p$ is the **power number**, a dimensionless constant that depends on impeller geometry and Re.

| Impeller Type | $N_p$ (turbulent) |
|---------------|-------------------|
| Rushton turbine (6-blade) | 5.0 |
| Pitched-blade turbine (45°, 4-blade) | 1.27 |
| Retreat-curve impeller | 0.3 – 0.5 |
| A310 / A320 hydrofoil | 0.3 |
| Magnetic stir bar | 0.3 – 0.5 (approx.) |

For laminar flow: $N_p \\approx 70 / Re$.
""")

# ── Power per Unit Volume ────────────────────────────────────────────────
st.header("3. Power per Unit Volume (Specific Energy Dissipation)")
st.latex(r"\varepsilon = \frac{P}{V}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Impeller power draw | W |
| $V$ | Liquid volume | m³ |
| $\\varepsilon$ | Mean specific energy dissipation rate | W/m³ |
""")

# ── Tip Speed ────────────────────────────────────────────────────────────
st.header("4. Impeller Tip Speed")
st.latex(r"u_{tip} = \pi \, N \, D")
st.markdown("Tip speed is a common scale-up criterion and relates to shear at the impeller.")

# ── Pumping Rate ─────────────────────────────────────────────────────────
st.header("5. Pumping Rate")
st.latex(r"Q = N_Q \, N \, D^3")
st.markdown("""
$N_Q$ is the **pumping number** (flow number).  Typical values:

| Impeller | $N_Q$ |
|----------|-------|
| Rushton turbine | 0.72 |
| Pitched-blade turbine (down-pumping) | 0.79 |
| A310 hydrofoil | 0.56 |
""")

# ── Kolmogorov Scale ─────────────────────────────────────────────────────
st.header("8. Kolmogorov Length Scale")
st.latex(r"\eta = \left( \frac{\nu^3}{\varepsilon} \right)^{1/4}")
st.markdown("""
The smallest eddy size in turbulent flow.  Below this scale, viscous dissipation dominates.
Typical values: 10 – 100 µm in stirred tanks.
""")

# ── Batchelor Scale ──────────────────────────────────────────────────────
st.header("9. Batchelor Length Scale")
st.latex(r"\lambda_B = \frac{\eta}{\sqrt{Sc}} = \eta \left( \frac{D_{mol}}{\nu} \right)^{1/2}")
st.markdown("""
The scale below which molecular diffusion homogenises concentration.

| Symbol | Description |
|--------|-------------|
| $Sc = \\nu / D_{mol}$ | Schmidt number |
| $D_{mol}$ | Molecular diffusivity (m²/s) |
""")

# ── Local max dissipation ────────────────────────────────────────────────
st.header("10. Maximum Local Energy Dissipation Rate")
st.latex(r"\varepsilon_{max} \approx C \cdot \frac{P}{\rho \, D^3}, \quad C \approx 1 - 5")
st.markdown("""
The energy dissipation rate near the impeller can be 1–2 orders of magnitude
higher than the mean.  We use $C = 3$ as a representative estimate
(Kresta & Wood, 1993).
""")

# ── Shear rate and stress ────────────────────────────────────────────────
st.header("15. Average Shear Rate (Camp-Stein)")
st.latex(r"\dot{\gamma}_{avg} = \sqrt{\frac{P}{\mu \, V}}")
st.markdown("""
Derived from equating the volume-averaged viscous dissipation to the power input:

$$\\mu \\, \\dot{\\gamma}^2 = \\frac{P}{V}$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Power input | W |
| $\\mu$ | Dynamic viscosity | Pa·s |
| $V$ | Liquid volume | m³ |
| $\\dot{\\gamma}_{avg}$ | RMS average shear rate | 1/s |

**Reference:** Camp, T. R. and Stein, P. C. (1943). Velocity gradients and internal work in
fluid motion. *J. Boston Soc. Civil Eng.*, 30, 219–237.
""")

st.header("16. Maximum Shear Rate (Impeller Zone)")
st.latex(r"\dot{\gamma}_{max} = \sqrt{\frac{\varepsilon_{max}}{\nu}}")
st.markdown("""
Evaluated using $\\varepsilon_{max}$ from Equation 10.  This gives the peak
shear rate in the impeller discharge zone.
""")

st.subheader("Shear Stress")
st.latex(r"\tau = \mu \, \dot{\gamma}")
st.markdown("""
For Newtonian fluids, shear stress (Pa) is simply viscosity multiplied by
the shear rate.  Both $\\tau_{avg}$ and $\\tau_{max}$ are reported.
""")
