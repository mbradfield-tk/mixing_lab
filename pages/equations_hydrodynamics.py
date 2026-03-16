"""Equations – Hydrodynamics & Shear."""

import streamlit as st

st.title("📐 Hydrodynamics & Shear Equations")

# ── Impeller Reynolds Number ─────────────────────────────────────────────
st.header("Impeller Reynolds Number")
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
st.header("Power Number & Impeller Power")
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
st.header("Power per Unit Volume (Specific Energy Dissipation)")
st.latex(r"\varepsilon = \frac{P}{V}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Impeller power draw | W |
| $V$ | Liquid volume | m³ |
| $\\varepsilon$ | Mean specific energy dissipation rate | W/m³ |
""")

# ── Tip Speed ────────────────────────────────────────────────────────────
st.header("Impeller Tip Speed")
st.latex(r"u_{tip} = \pi \, N \, D")
st.markdown("Tip speed is a common scale-up criterion and relates to shear at the impeller.")

# ── Pumping Rate ─────────────────────────────────────────────────────────
st.header("Pumping Rate")
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
st.header("Kolmogorov Length Scale")
st.latex(r"\eta = \left( \frac{\nu^3}{\varepsilon} \right)^{1/4}")
st.markdown("""
The smallest eddy size in turbulent flow.  Below this scale, viscous dissipation dominates.
Typical values: 10 – 100 µm in stirred tanks.

**Important:** $\\varepsilon$ here is the **mass-specific** energy dissipation rate in **W/kg** (= m²/s³),
i.e. $\\varepsilon = P / (\\rho\\,V)$, not $P/V$ in W/m³.
""")

# ── Batchelor Scale ──────────────────────────────────────────────────────
st.header("Batchelor Length Scale")
st.latex(r"\lambda_B = \frac{\eta}{\sqrt{Sc}} = \eta \left( \frac{D_{mol}}{\nu} \right)^{1/2}")
st.markdown("""
The scale below which molecular diffusion homogenises concentration.

| Symbol | Description |
|--------|-------------|
| $Sc = \\nu / D_{mol}$ | Schmidt number |
| $D_{mol}$ | Molecular diffusivity (m²/s) |

**Note:** $\\varepsilon$ must be in **W/kg** (= m²/s³) for this formula.
""")

# ── Local max dissipation ────────────────────────────────────────────────
st.header("Maximum Local Energy Dissipation Rate")
st.latex(r"\varepsilon_{max} \approx C \cdot \frac{P}{\rho \, D^3}, \quad C \approx 1 - 5")
st.markdown("""
The energy dissipation rate near the impeller can be 1–2 orders of magnitude
higher than the mean.  We use $C = 3$ as a representative estimate
(Kresta & Wood, 1993).
""")

# ── Shear rate and stress ────────────────────────────────────────────────
st.header("Average Shear Rate (Camp-Stein)")
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

st.header("Maximum Shear Rate (Impeller Zone)")
st.latex(r"\dot{\gamma}_{max} = \sqrt{\frac{\varepsilon_{max}}{\nu}}")
st.markdown("""
Evaluated using $\\varepsilon_{max}$ from the Maximum Local Energy Dissipation Rate equation above.  This gives the peak
shear rate in the impeller discharge zone.
""")

st.subheader("Shear Stress")
st.latex(r"\tau = \mu \, \dot{\gamma}")
st.markdown("""
For Newtonian fluids, shear stress (Pa) is simply viscosity multiplied by
the shear rate.  Both $\\tau_{avg}$ and $\\tau_{max}$ are reported.
""")

# ── Circulation Time ─────────────────────────────────────────────────────
st.header("Circulation Time")
st.latex(r"t_c = \frac{V}{N_q \, N \, D^3}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $V$ | Liquid volume | m³ |
| $N_q$ | Pumping (flow) number | — |
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |

The mean time for a fluid element to complete one loop through the
impeller zone.  Related to blend time via $\\theta_{95} \\approx 5.2 \\, t_c$.

**Reference:** Eq 2.14 (PDF reference).
""")

# ── Torque ────────────────────────────────────────────────────────────────
st.header("Impeller Torque")
st.latex(r"\Lambda = \frac{P}{2 \pi N}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Impeller power draw | W |
| $N$ | Impeller speed | rev/s |
| $\\Lambda$ | Torque | N·m |

**Torque per unit volume** $\\Lambda/V$ provides a scale-independent measure
of mechanical stress on the fluid, useful for comparing across scales.

**Reference:** Eq 2.20–2.21 (PDF reference).
""")

# ── EDCF ──────────────────────────────────────────────────────────────────
st.header("Energy Dissipation Circulation Function (EDCF)")
st.latex(r"\text{EDCF} = \frac{\varepsilon_{max}}{t_c}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\varepsilon_{max}$ | Maximum local energy dissipation rate | W/kg |
| $t_c$ | Circulation time | s |
| EDCF | Mixing intensity × exposure frequency | W/kg/s |

EDCF combines the **intensity** of the impeller discharge zone with the
**frequency** at which fluid passes through it.  A higher EDCF indicates
more vigorous overall mixing — particularly relevant for shear-sensitive
and fast-reacting systems.

**Reference:** Eq 2.23 (PDF reference); Nienow (1997).
""")

# ── Froude Number ─────────────────────────────────────────────────────────
st.header("Froude Number")
st.latex(r"Fr = \frac{N^2 \, D}{g}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |
| $g$ | Gravitational acceleration | m/s² |

Ratio of inertial to gravitational forces.  Important for predicting
free-surface vortex formation in unbaffled or partially baffled vessels.
At high Fr the liquid surface can draw down to the impeller.

**Reference:** Eq 2.30 (PDF reference).
""")
