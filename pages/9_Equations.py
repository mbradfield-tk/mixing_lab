"""
Page 6 – Equations Summary
===========================
Reference for every correlation and equation used in Mixing Lab.
"""

import streamlit as st
from utils.sidebar import render_sidebar
render_sidebar()

st.title("📐 Equations Summary")
st.markdown("All equations used in the Mixing Lab calculations are listed below with references.")

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

# ── Scale-up rules ───────────────────────────────────────────────────────
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

# ── Scale-up rules ───────────────────────────────────────────────────────
st.header("17. Common Scale-Up Rules")

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
st.header("18. Temperature-Dependent Solvent Properties")
st.markdown("""
The Fluid Database's *Solvent Properties (T)* tab provides physical properties
at any liquid-phase temperature for common pharmaceutical solvents, using the
correlations below.
""")

st.subheader("18a. Density")
st.latex(r"\rho(T) = \rho_{25} + \frac{d\rho}{dT} \cdot (T - 25) \quad [\mathrm{kg/m^3}]")
st.markdown("Linear fit anchored at the known 25 °C value.  $T$ in °C, slope $d\\rho/dT$ (typically negative) from published data.")

st.subheader("18b. Dynamic Viscosity (Arrhenius)")
st.latex(r"\mu(T) = \mu_{25} \, \exp\!\left[ \frac{E_a}{R} \left( \frac{1}{T_K} - \frac{1}{298.15} \right) \right]")
st.markdown("""
$T_K$ = temperature in Kelvin, $E_a$ = activation energy for viscous flow
(J/mol), $R$ = 8.314 J/(mol·K).  Anchored at the known $\\mu_{25}$ so the
reference value is recovered exactly.  $E_a$ values from published data
(Perry's, CRC Handbook, DIPPR).
""")

st.subheader("18c. Surface Tension")
st.latex(r"\sigma(T) = \sigma_{25} + \frac{d\sigma}{dT} \cdot (T - 25) \quad [\mathrm{N/m}]")
st.markdown("Linear fit anchored at the known 25 °C value.  Clamped to $\\sigma \\ge 0$.")

st.subheader("18d. Molecular Diffusivity (Stokes-Einstein Scaling)")
st.latex(r"D(T) = D_{\mathrm{ref}} \; \frac{T_K}{298.15} \; \frac{\mu_{\mathrm{ref}}}{\mu(T)}")
st.markdown("""
Scaled from a reference value $D_{\\mathrm{ref}}$ measured at 25 °C using the
Stokes-Einstein relation ($D \\propto T / \\mu$).  This captures the dominant
temperature dependence through the viscosity term.

**Data sources:** Yaws' Handbook, Perry's Chemical Engineers' Handbook (9th ed.),
CRC Handbook of Chemistry and Physics, DIPPR correlations.
""")

# =====================================================================
# Solid-Particle Hydrodynamics
# =====================================================================

st.header("19. Terminal Settling Velocity")
st.latex(r"v_{t,\text{Stokes}} = \frac{d_p^2 \, g \, \Delta\rho}{18 \, \mu} \cdot \phi")
st.markdown(r"""
**Stokes regime** ($Re_p < 0.1$): valid for small or slow-settling particles.
$\phi$ is the sphericity (shape factor, 0–1).

For intermediate Reynolds numbers ($0.1 < Re_p < 1000$) the
**Schiller-Naumann** drag correction is applied iteratively:

$$C_D = \frac{24}{Re_p}\left(1 + 0.15\,Re_p^{0.687}\right)$$

$$v_t = \sqrt{\frac{4\,g\,d_p\,\Delta\rho}{3\,C_D\,\rho_L}} \cdot \phi$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $d_p$ | Particle diameter (typically D50) | m |
| $g$ | Gravitational acceleration | m/s² |
| $\Delta\rho$ | $\lvert\rho_p - \rho_L\rvert$ | kg/m³ |
| $\mu$ | Dynamic viscosity | Pa·s |
| $\phi$ | Sphericity / shape factor | – |

**Reference:** Clift, R., Grace, J.R. and Weber, M.E. (1978).
*Bubbles, Drops, and Particles.* Academic Press.
""")

st.header("20. Particle Reynolds Number")
st.latex(r"Re_p = \frac{\rho_L \, v_t \, d_p}{\mu}")
st.markdown("Classifies the flow regime around a settling particle.")

st.header("21. Zwietering Just-Suspended Speed")
st.latex(r"N_{js} = S \, \nu^{0.1} \, d_p^{0.2} \left(\frac{g\,\Delta\rho}{\rho_L}\right)^{0.45} X^{0.13} \, D^{-0.85}")
st.markdown(r"""
The Zwietering (1958) correlation gives the minimum impeller speed to achieve
**complete off-bottom suspension** of solid particles.

| Symbol | Description | Units |
|--------|-------------|-------|
| $S$ | Zwietering constant (geometry-dependent) | – |
| $\nu$ | Kinematic viscosity | m²/s |
| $d_p$ | Mass-mean particle diameter | m |
| $g$ | Gravitational acceleration | m/s² |
| $\Delta\rho$ | $\lvert\rho_p - \rho_L\rvert$ | kg/m³ |
| $\rho_L$ | Liquid density | kg/m³ |
| $X$ | Solids loading | wt-% |
| $D$ | Impeller diameter | m |

**Typical S values:**

| Impeller / geometry | S |
|-------------------|---|
| PBT (45° down-pumping), D/T ≈ 0.4 | 4.5 – 6.5 |
| Rushton turbine, D/T ≈ 0.33 | 7 – 9 |
| Hydrofoil (A310), D/T ≈ 0.4 | 3 – 5 |

**Reference:** Zwietering, Th.N. (1958). Suspending of solid particles in
liquid by agitators. *Chem. Eng. Sci.*, 8(3–4), 244–253.
""")

st.header("22. Solid-Liquid Mass Transfer (Ranz-Marshall)")
st.latex(r"Sh = 2 + 0.6 \, Re_p^{1/2} \, Sc^{1/3}")
st.latex(r"k_{SL} = \frac{Sh \cdot D_{mol}}{d_p}")
st.markdown(r"""
The Ranz-Marshall correlation estimates the convective mass-transfer
coefficient from the liquid to a suspended particle surface.

| Symbol | Description | Units |
|--------|-------------|-------|
| $Sh$ | Sherwood number | – |
| $Re_p$ | Particle Reynolds number | – |
| $Sc = \nu / D_{mol}$ | Schmidt number | – |
| $D_{mol}$ | Molecular diffusivity | m²/s |
| $d_p$ | Particle diameter | m |
| $k_{SL}$ | Solid-liquid mass transfer coefficient | m/s |

**Reference:** Ranz, W.E. and Marshall, W.R. (1952). Evaporation from drops.
*Chem. Eng. Prog.*, 48(3), 141–146; 48(4), 173–180.
""")

st.header("23. Suspension Quality Assessment")
st.markdown(r"""
The ratio $N / N_{js}$ indicates the quality of particle suspension:

| $N / N_{js}$ | Assessment |
|-------------|------------|
| $< 0.7$ | Poorly suspended |
| $0.7 – 1.0$ | Partially suspended |
| $1.0 – 1.3$ | Just suspended |
| $> 1.3$ | Fully / homogeneously suspended |

Operating at $N \ge N_{js}$ ensures complete off-bottom suspension.
Higher speeds improve homogeneity but energy cost increases as $N^3$.
""")
