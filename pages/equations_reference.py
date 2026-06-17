"""Equations – Combined Reference Page.

All equation reference sections consolidated into a single page,
each within an expandable section.
"""

import streamlit as st
from utils.rom_registry import (
    get_all_correlations,
    get_all_registered_reactors,
    PARAM_DISPLAY,
    SUPPORTED_PARAMS,
)

st.title("📐 Equations Reference")

st.markdown("""
Reference correlations and equations used throughout Mixing Lab.
Expand each section below for details.

---
""")

# ══════════════════════════════════════════════════════════════════════════
# 1. HYDRODYNAMICS & SHEAR
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Hydrodynamics & Shear**", expanded=False):

    # ── Impeller Reynolds Number
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

**Reference:** Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Wiley, Ch. 6.
""")

    # ── Power Number & Power
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

**Reference:** Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Wiley, Ch. 6.
""")

    # ── Power per Unit Volume
    st.header("Power per Unit Volume (Specific Energy Dissipation)")
    st.latex(r"\varepsilon = \frac{P}{V}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Impeller power draw | W |
| $V$ | Liquid volume | m³ |
| $\\varepsilon$ | Mean specific energy dissipation rate | W/m³ |

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")

    # ── Tip Speed
    st.header("Impeller Tip Speed")
    st.latex(r"u_{tip} = \pi \, N \, D")
    st.markdown("""
Tip speed is a common scale-up criterion and relates to shear at the impeller.

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")

    # ── Pumping Rate
    st.header("Pumping Rate")
    st.latex(r"Q = N_Q \, N \, D^3")
    st.markdown("""
$N_Q$ is the **pumping number** (flow number).  Typical values:

| Impeller | $N_Q$ |
|----------|-------|
| Rushton turbine | 0.72 |
| Pitched-blade turbine (down-pumping) | 0.79 |
| A310 hydrofoil | 0.56 |

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")

    # ── Kolmogorov Scale
    st.header("Kolmogorov Length Scale")
    st.latex(r"\eta = \left( \frac{\nu^3}{\varepsilon} \right)^{1/4}")
    st.markdown("""
The smallest eddy size in turbulent flow.  Below this scale, viscous dissipation dominates.
Typical values: 10 – 100 µm in stirred tanks.

**Important:** $\\varepsilon$ here is the **mass-specific** energy dissipation rate in **W/kg** (= m²/s³),
i.e. $\\varepsilon = P / (\\rho\\,V)$, not $P/V$ in W/m³.

**Reference:** Kolmogorov, A.N. (1991). The local structure of turbulence in incompressible viscous fluid for very large Reynolds numbers. [*Proc. R. Soc. Lond. A*](https://doi.org/10.1098/rspa.1991.0075), 434, 9–13. (English translation of the 1941 *Dokl. Akad. Nauk SSSR*, 30, 301–305.)
""")

    # ── Batchelor Scale
    st.header("Batchelor Length Scale")
    st.latex(r"\lambda_B = \frac{\eta}{\sqrt{Sc}} = \eta \left( \frac{D_{mol}}{\nu} \right)^{1/2}")
    st.markdown("""
The scale below which molecular diffusion homogenises concentration.

| Symbol | Description |
|--------|-------------|
| $Sc = \\nu / D_{mol}$ | Schmidt number |
| $D_{mol}$ | Molecular diffusivity (m²/s) |

**Note:** $\\varepsilon$ must be in **W/kg** (= m²/s³) for this formula.

**Reference:** Batchelor, G.K. (1959). Small-scale variation of convected quantities like temperature in turbulent fluid. [*J. Fluid Mech.*](https://doi.org/10.1017/S0022112059000362), 5(1), 113–133.
""")

    # ── Local max dissipation
    st.header("Maximum Local Energy Dissipation Rate")
    st.latex(r"\varepsilon_{max} \approx C \cdot \frac{P}{\rho \, D^3}, \quad C \approx 1 - 5")
    st.markdown("""
The energy dissipation rate near the impeller can be 1–2 orders of magnitude
higher than the mean.  We use $C = 3$ as a representative estimate.

**Reference:** Kresta, S.M. & Wood, P.E. (1993). The flow field produced by a pitched blade turbine: Characterization of the turbulence and estimation of the dissipation rate. [*Chem. Eng. Sci.*](https://doi.org/10.1016/0009-2509(93)80346-R), 48(10), 1761–1774.
""")

    # ── Shear rate and stress
    st.header("Average Shear Rate (Camp–Stein)")
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

**Reference:** Camp, T.R. & Stein, P.C. (1943). Velocity gradients and internal work in fluid motion. [*J. Boston Soc. Civil Eng.*](https://en.wikipedia.org/wiki/Camp%E2%80%93Stein_equation), 30, 219–237.
""")

    st.header("Maximum Shear Rate (Impeller Zone)")
    st.latex(r"\dot{\gamma}_{max} = \sqrt{\frac{\varepsilon_{max}}{\nu}}")
    st.markdown("""
Evaluated using $\\varepsilon_{max}$ from the Maximum Local Energy Dissipation Rate equation above.  This gives the peak
shear rate in the impeller discharge zone.

**Reference:** Kresta & Wood (1993), as above.
""")

    st.subheader("Shear Stress")
    st.latex(r"\tau = \mu \, \dot{\gamma}")
    st.markdown("""
For Newtonian fluids, shear stress (Pa) is simply viscosity multiplied by
the shear rate.  Both $\\tau_{avg}$ and $\\tau_{max}$ are reported.

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")

    # ── Circulation Time
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

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 9.
""")

    # ── Torque
    st.header("Impeller Torque")
    st.latex(r"\Lambda = \frac{P}{2 \pi N}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Impeller power draw | W |
| $N$ | Impeller speed | rev/s |
| $\\Lambda$ | Torque | N·m |

**Torque per unit volume** $\\Lambda/V$ provides a scale-independent measure
of mechanical stress on the fluid.

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")

    # ── EDCF
    st.header("Energy Dissipation Circulation Function (EDCF)")
    st.latex(r"\text{EDCF} = \frac{\varepsilon_{max}}{t_c}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\varepsilon_{max}$ | Maximum local energy dissipation rate | W/kg |
| $t_c$ | Circulation time | s |
| EDCF | Mixing intensity × exposure frequency | W/kg/s |

EDCF combines the **intensity** of the impeller discharge zone with the
**frequency** at which fluid passes through it.

**Reference:** Nienow, A.W. (1997). On impeller circulation and mixing effectiveness in the turbulent flow regime. [*Chem. Eng. Sci.*](https://doi.org/10.1016/S0009-2509(97)00072-9), 52(15), 2557–2565.
""")

    # ── Froude Number
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

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 6.
""")


# ══════════════════════════════════════════════════════════════════════════
# 2. MIXING & DAMKÖHLER
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Mixing & Damköhler Numbers**", expanded=False):

    # ── Blend Time
    st.header("Macro-Blend Time (95 %)")
    st.latex(r"\theta_{95} = \frac{5.2 \, V}{N_Q \, N \, D^3}")
    st.markdown("""
Grenville correlation for turbulent blending in a baffled stirred tank.
This represents the time to achieve 95 % uniformity on the bulk (macro) scale.

**Reference:** Grenville, R.K. (1992). *Blending of Viscous Newtonian and Non-Newtonian Fluids.* Ph.D. thesis, Cranfield University. Also: Grenville, R.K. & Nienow, A.W. (2004). Blending of miscible liquids. In [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 9.
""")

    # ── Micromixing Time
    st.header("Engulfment Micromixing Time")
    st.latex(r"t_E = 17.3 \left( \frac{\nu}{\varepsilon} \right)^{1/2}")
    st.markdown("""
Baldyga & Bourne engulfment model – characterises the time for mixing at the
smallest (molecular) scales.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\\nu$ | Kinematic viscosity ($\\mu / \\rho$) | m²/s |
| $\\varepsilon$ | Mass-specific energy dissipation rate $P/(\\rho V)$ | W/kg (= m²/s³) |

**Reference:** Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
""")

    # ── Damköhler Numbers
    st.header("Damköhler Numbers")

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

**Reference:** Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Wiley, Ch. 13.
""")

    # ── Gas-Liquid Damköhler
    st.header("Gas-Liquid Damköhler Number")
    st.latex(r"Da_{GL} = \frac{1}{k_L a \; t_{rxn}} = \frac{t_{\text{transfer}}}{t_{rxn}}")
    st.markdown(r"""
Compares the characteristic **gas-liquid mass-transfer time** $\;(1/k_L a)\;$ to
the reaction time.

| $Da_{GL}$ range | Interpretation |
|----------------|----------------|
| $Da_{GL} \ll 0.01$ | Mass transfer is fast — reaction-limited |
| $0.01 < Da_{GL} < 0.1$ | Likely insensitive to G-L transfer |
| $0.1 < Da_{GL} < 1$ | Potentially transfer-limited |
| $1 < Da_{GL} < 10$ | Transfer-limited |
| $Da_{GL} \gg 10$ | Strongly transfer-limited |

When $k_L a = 0$ (no gas phase present), $Da_{GL}$ is set to zero (not applicable).

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 13.
""")

    # ── Solid-Liquid Damköhler
    st.header("Solid-Liquid Damköhler Number")
    st.latex(r"Da_{SL} = \frac{1}{k_L a_{SL} \; t_{rxn}} = \frac{t_{\text{SL transfer}}}{t_{rxn}}")
    st.markdown(r"""
Compares the characteristic **solid-liquid mass-transfer time** $\;(1/k_L a_{SL})\;$
to the reaction time, where $k_L a_{SL}$ is the volumetric solid-liquid
mass-transfer coefficient:

$$k_L a_{SL} = k_{SL} \cdot a_s, \quad a_s = \frac{6\,\phi_s}{d_p}$$

The slip velocity used in the Ranz-Marshall correlation is the larger of the
terminal settling velocity and the Kolmogorov-scale turbulent slip velocity:

$$v_{\text{slip}} = \max\!\bigl(v_t,\; (\varepsilon \cdot d_p)^{1/3}\bigr)$$

| $Da_{SL}$ range | Interpretation |
|----------------|----------------|
| $Da_{SL} \ll 0.01$ | Mass transfer is fast — reaction-limited |
| $0.01 < Da_{SL} < 0.1$ | Likely insensitive to S-L transfer |
| $0.1 < Da_{SL} < 1$ | Potentially transfer-limited |
| $1 < Da_{SL} < 10$ | Transfer-limited |
| $Da_{SL} \gg 10$ | Strongly transfer-limited |

When $k_L a_{SL} = 0$ (no solid phase present), $Da_{SL}$ is set to zero.

| Symbol | Description | Units |
|--------|-------------|-------|
| $k_{SL}$ | Solid-liquid mass-transfer coefficient | m/s |
| $a_s$ | Specific particle surface area | 1/m |
| $\phi_s$ | Volumetric solids fraction | — |
| $d_p$ | Particle diameter | m |
| $v_t$ | Terminal settling velocity | m/s |
| $\varepsilon$ | Specific energy dissipation rate | W/kg |
| $t_{rxn}$ | Characteristic reaction time | s |

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 13; Ranz, W.E. & Marshall, W.R. (1952). [*Chem. Eng. Prog.*](https://en.wikipedia.org/wiki/Ranz%E2%80%93Marshall_correlation), 48, 141–146.
""")

    # ── Reaction time
    st.header("Characteristic Reaction Time")

    st.subheader("First-order or pseudo-first-order")
    st.latex(r"t_{rxn} = \frac{1}{k}")

    st.subheader("Second-order")
    st.latex(r"t_{rxn} = \frac{1}{k \, C_0}")

    st.markdown("""
The characteristic reaction time is taken as the **reaction time constant**
(the e-folding time, i.e. the time for the limiting reagent to fall to $1/e$ of
its initial value), consistent with the Damköhler-number definition
$Da = t_{mix} / t_{rxn}$.  Note this is **not** the 50 % half-life
$t_{1/2} = \\ln 2 / k$, which is ~30 % shorter for a first-order reaction.

| Symbol | Description | Units |
|--------|-------------|-------|
| $k$ | Rate constant | 1/s (1st order) or L/(mol·s) (2nd order) |
| $C_0$ | Initial concentration | mol/L |

**Reference:** Levenspiel, O. (1999). [*Chemical Reaction Engineering*](https://openlibrary.org/isbn/047125424X), 3rd ed. Wiley, Ch. 3.
""")

    # ── Local micromixing time
    st.header("Local Micromixing Time (Impeller Zone)")
    st.latex(r"t_{E,\text{local}} = 17.3 \left( \frac{\nu}{\varepsilon_{max}} \right)^{1/2}")
    st.markdown("""
Same Baldyga & Bourne engulfment model, but evaluated using the **local maximum**
energy dissipation rate $\\varepsilon_{max}$ near the impeller.  Relevant when
reagents are **fed directly into the impeller discharge zone**.

**Reference:** Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
""")

    # ── Mesomixing Time Scales
    st.header("Mesomixing Time Scales")
    st.markdown(r"""
Mesomixing acts at the scale of the feed plume / feed-pipe diameter and has
**two** characteristic time scales (Baldyga & Bourne, 1999; Myerson 2019, Table 8.1).
The slower of the two governs feed-plume dispersion.
""")

    st.subheader("Inertial-Convective Disintegration")
    st.latex(r"\tau_S = A \left( \frac{\Lambda_C^{2}}{\varepsilon} \right)^{1/3}, \quad A \approx 1.2")
    st.markdown(r"""
Break-up of the feed plume by inertial-convective disintegration of the large
eddies. $\Lambda_C$ is the integral scale of the velocity fluctuations, set by
the feed-pipe / plume size. **Mixing Lab implements this scale** using the
feed-pipe diameter $d_{\text{feed}}$ as the length scale with $A = 2$:

$$t_{\text{meso}} = 2 \left( \frac{d_{\text{feed}}^2}{\varepsilon} \right)^{1/3}$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $\Lambda_C$ | Integral scale of velocity fluctuations | m |
| $d_{\text{feed}}$ | Feed-pipe internal diameter | m |
| $\varepsilon$ | Local energy dissipation at the feed point | W/kg (= m²/s³) |
""")

    st.subheader("Turbulent Dispersion")
    st.latex(r"\tau_D = \frac{Q_{\text{feed}}}{\bar{u}\, D_t}")
    st.markdown(r"""
Turbulent diffusion of the feed plume into the bulk. This scale **grows with
feed rate** $Q_{\text{feed}}$ — the lever probed by Test 2 of the Bourne Protocol.

| Symbol | Description | Units |
|--------|-------------|-------|
| $Q_{\text{feed}}$ | Volumetric feed rate | m³/s |
| $\bar{u}$ | Local mean velocity at the feed point | m/s |
| $D_t$ | Turbulent diffusivity | m²/s |

Ordering of the mixing time scales: $\;t_{\text{blend}} > \tau_S,\;\tau_D > t_E$.

**Reference:** Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley, Ch. 3; Myerson, A.S. *et al.* (Eds.) (2019). [*Handbook of Industrial Crystallization*](https://doi.org/10.1017/9781139026949), 3rd ed., Cambridge Univ. Press, Ch. 8 (Table 8.1).
""")


# ══════════════════════════════════════════════════════════════════════════
# 3. MASS TRANSFER & kLa
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Mass Transfer & kLa**", expanded=False):

    # ── kLa – sparged
    st.header("kLa – Volumetric Mass-Transfer Coefficient (Gas-Liquid)")
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

**Reference:** van 't Riet, K. (1979). Review of measuring methods and results in non-viscous gas-liquid mass transfer in stirred vessels. [*Ind. Eng. Chem. Process Des. Dev.*](https://doi.org/10.1021/i260071a001), 18(3), 357–364.
""")

    # ── kLa – headspace / free-surface
    st.header("kLa – Free-Surface (Headspace) Mass Transfer")
    st.latex(r"k_L = 0.4 \; D_{mol}^{\,1/2} \left( \frac{\varepsilon}{\nu} \right)^{1/4}")
    st.latex(r"a_{\text{surface}} = \frac{\pi / 4 \; D_T^2}{V}")
    st.latex(r"k_L a_{\text{surface}} = k_L \cdot a_{\text{surface}}")
    st.markdown("""
Lamont & Scott (1970) small-eddy surface-renewal model.

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

**Reference:** Lamont, J.C. & Scott, D.S. (1970). An eddy cell model of mass transfer into the surface of a turbulent liquid. [*AIChE J.*](https://doi.org/10.1002/aic.690160410), 16(4), 513–519.
""")

    # ── Gas Holdup
    st.header("Gas Holdup (Hughmark)")
    st.latex(r"\varepsilon_G = 0.505 \; v_s^{0.47} \left( \frac{P}{V} \right)^{0.4} \left( \frac{\mu}{\sigma} \right)^{0.08}")
    st.markdown(r"""
Simplified Hughmark (1967) / Yawalkar correlation for gas void fraction.

| Symbol | Description | Units |
|--------|-------------|-------|
| $v_s$ | Superficial gas velocity | m/s |
| $P/V$ | Gassed power per unit volume | W/m³ |
| $\mu$ | Liquid dynamic viscosity | Pa·s |
| $\sigma$ | Liquid surface tension | N/m |
| $\varepsilon_G$ | Gas holdup (volume fraction) | — |

**Reference:** Hughmark, G.A. (1967). Holdup and mass transfer in bubble columns. [*Ind. Eng. Chem. Process Des. Dev.*](https://doi.org/10.1021/i260022a011), 6(2), 218–220.
""")

    # ── Sauter Mean Bubble Diameter
    st.header("Sauter Mean Bubble Diameter (Calderbank)")
    st.latex(r"d_{32} = 4.15 \; \frac{\sigma^{0.6}}{(P/V)^{0.4} \; \rho_L^{0.2}} + 0.0009")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\sigma$ | Surface tension | N/m |
| $P/V$ | Power per unit volume | W/m³ |
| $\rho_L$ | Liquid density | kg/m³ |
| $d_{32}$ | Sauter mean bubble diameter | m |

The first term captures breakup by turbulence; the constant 0.0009 m
represents the minimum stable bubble size.

**Reference:** Calderbank, P.H. (1958). Physical rate processes in industrial fermentation. [*Trans. Inst. Chem. Engrs.*](https://scholar.google.com/scholar?q=Calderbank+1958+physical+rate+processes+industrial+fermentation), 36, 443–463.
""")

    # ── Gas Flooding Speed
    st.header("Gas Flooding Speed")
    st.latex(r"Fl_G = \frac{Q_{gas}}{N \, D^3}")
    st.latex(r"N_{flood} = \frac{Q_{gas}}{Fl_{crit} \, D^3}, \quad Fl_{crit} \approx 0.035")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $Q_{gas}$ | Volumetric gas flow rate | m³/s |
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |
| $Fl_{crit}$ | Critical flooding number | — |
| $N_{flood}$ | Minimum speed for gas dispersion | rev/s |

Operating at $N / N_{flood} > 1$ ensures the gas is dispersed.
Recommended: $N / N_{flood} > 1.3$ for reliable dispersion.

**Reference:** Nienow, A.W., Wisdom, D.J. & Middleton, J.C. (1977). The effect of scale and geometry on flooding, recirculation and power in gassed stirred vessels. *Proceedings 2nd European Conference on Mixing*, F1-1.
""")

    # ── Liquid-Liquid section
    st.divider()
    st.subheader("Liquid-Liquid Equations")

    # ── Miscibility Screening
    st.header("Miscibility Screening")
    st.markdown(r"""
**Two-tier approach:**

**1. Known-pair lookup** — Experimental miscibility data from standard references
(Perry's, CRC Handbook, Merck Index).

**2. Hansen Solubility Parameters (fallback):**
""")
    st.latex(r"R_a = \sqrt{4\,(\delta_{d1}-\delta_{d2})^2 + (\delta_{p1}-\delta_{p2})^2 + (\delta_{h1}-\delta_{h2})^2}")
    st.markdown(r"""
| $R_a$ (MPa$^{0.5}$) | Interpretation |
|----------------------|----------------|
| < 15 | Likely **miscible** |
| 15 – 25 | Partially miscible / borderline |
| > 25 | Likely **immiscible** |

**Reference:** Hansen, C.M. (2007). [*Hansen Solubility Parameters: A User's Handbook*](https://doi.org/10.1201/9781420006834), 2nd ed. CRC Press.
""")

    # ── Weber Number
    st.header("Impeller Weber Number (Liquid-Liquid)")
    st.latex(r"We = \frac{\rho_c \, N^2 \, D^3}{\sigma_{LL}}")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\rho_c$ | Continuous-phase density | kg/m³ |
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |
| $\sigma_{LL}$ | Interfacial tension between two liquids | N/m |

**Reference:** Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 12.
""")

    # ── Sauter Mean Drop Diameter
    st.header("Sauter Mean Drop Diameter (Hinze-Kolmogorov)")
    st.latex(r"\frac{d_{32}}{D} = C_1 \, We^{-0.6} (1 + C_2 \, \phi_d)")
    st.markdown(r"""
| Symbol | Description | Typical value |
|--------|-------------|---------------|
| $C_1$ | Constant | 0.053 |
| $C_2$ | Holdup correction | 3.0 |
| $\phi_d$ | Dispersed-phase volume fraction | 0 – 0.5 |
| $D$ | Impeller diameter | m |

**Reference:** Hinze, J.O. (1955). Fundamentals of the hydrodynamic mechanism of splitting in dispersion processes. [*AIChE J.*](https://doi.org/10.1002/aic.690010303), 1(3), 289–295; Chen, H.T. & Middleman, S. (1967). [*AIChE J.*](https://doi.org/10.1002/aic.690130523), 13(5), 989–995.
""")

    # ── Minimum Dispersion Speed
    st.header("Minimum Dispersion Speed (Skelland & Seksaria)")
    st.latex(r"N_{min} = C \left( \frac{\sigma_{LL}}{\rho_c \, D^3} \right)^{1/2} (1 + 2.5 \, \phi_d)")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $C$ | Geometry constant (≈ 1.03 for Rushton) | — |
| $\sigma_{LL}$ | Interfacial tension | N/m |
| $\rho_c$ | Continuous-phase density | kg/m³ |
| $D$ | Impeller diameter | m |
| $\phi_d$ | Dispersed-phase volume fraction | — |

**Reference:** Skelland, A.H.P. & Seksaria, R. (1978). Minimum impeller speeds for liquid-liquid dispersion. [*Ind. Eng. Chem. Process Des. Dev.*](https://doi.org/10.1021/i260065a011), 17(1), 56–61.
""")

    # ── LL Mass Transfer
    st.header("Liquid-Liquid Mass-Transfer Coefficient")
    st.latex(r"v_{slip} = (\varepsilon \, d_{32})^{1/3}")
    st.latex(r"Re_d = \frac{\rho_c \, v_{slip} \, d_{32}}{\mu_c}")
    st.latex(r"Sh = 2 + 0.6 \, Re_d^{1/2} \, Sc^{1/3}")
    st.latex(r"k_{LL} = \frac{Sh \cdot D_{mol}}{d_{32}}")
    st.markdown(r"""
Ranz-Marshall type correlation applied to liquid-liquid systems.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\varepsilon$ | Mass-specific energy dissipation rate | W/kg |
| $d_{32}$ | Sauter mean drop diameter | m |
| $v_{slip}$ | Turbulent slip velocity | m/s |
| $Re_d$ | Drop Reynolds number | — |
| $Sc = \mu_c / (\rho_c \, D_{mol})$ | Schmidt number | — |
| $D_{mol}$ | Molecular diffusivity of solute | m²/s |
| $k_{LL}$ | Liquid-liquid mass-transfer coefficient | m/s |

$$a_{LL} = \frac{6 \, \phi_d}{d_{32}}, \quad k_L a_{LL} = k_{LL} \cdot a_{LL}$$

**Reference:** Calderbank, P.H. & Moo-Young, M.B. (1961). The continuous phase heat and mass-transfer properties of dispersions. [*Chem. Eng. Sci.*](https://doi.org/10.1016/0009-2509(61)80028-8), 16(1–2), 39–54.
""")

    # ── Phase Separation
    st.header("Phase Separation Time Estimate")
    st.latex(r"v_{drop} = \frac{\Delta\rho \, g \, d_{32}^2}{18 \, \mu_c}")
    st.latex(r"t_{sep} \approx \frac{H}{v_{drop}}")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\Delta\rho$ | $\lvert\rho_d - \rho_c\rvert$ | kg/m³ |
| $g$ | Gravitational acceleration | m/s² |
| $d_{32}$ | Sauter mean drop diameter | m |
| $\mu_c$ | Continuous-phase viscosity | Pa·s |
| $H$ | Liquid height | m |

| Separation time | Assessment |
|-----------------|------------|
| < 1 min | Rapid — unstable dispersion |
| 1 – 10 min | Moderate |
| 10 – 60 min | Slow — reasonably stable |
| > 1 h | Very stable dispersion |

**Reference:** Clift, R., Grace, J.R. & Weber, M.E. (1978). [*Bubbles, Drops, and Particles*](https://openlibrary.org/isbn/012176950X), Academic Press.
""")


# ══════════════════════════════════════════════════════════════════════════
# 4. HEAT BALANCE
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Heat Balance**", expanded=False):

    # ── Heat generation
    st.header("Reaction Rate (Molar)")

    st.subheader("First-order or pseudo-first-order")
    st.latex(r"r = k \, C_0 \quad [\mathrm{mol/(L \cdot s)}]")

    st.subheader("Second-order or pseudo-second-order")
    st.latex(r"r = k \, C_0^2 \quad [\mathrm{mol/(L \cdot s)}]")

    st.latex(r"r_{\text{total}} = r \times V_L \quad [\mathrm{mol/s}]")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $k$ | Rate constant | 1/s (1st) or L/(mol·s) (2nd) |
| $C_0$ | Initial concentration | mol/L |
| $V_L$ | Liquid volume | L |

**Reference:** Levenspiel, O. (1999). [*Chemical Reaction Engineering*](https://openlibrary.org/isbn/047125424X), 3rd ed. Wiley.
""")

    st.header("Heat Generation Rate")
    st.latex(r"Q_{\text{gen}} = \lvert \Delta H_{rxn} \rvert \times r_{\text{total}}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\Delta H_{rxn}$ | Enthalpy of reaction (negative = exothermic) | kJ/mol |
| $r_{\\text{total}}$ | Total molar reaction rate | mol/s |
| $Q_{\\text{gen}}$ | Rate of heat release | W |

**Reference:** Levenspiel (1999), Ch. 9.
""")

    # ── Adiabatic temperature rise
    st.header("Adiabatic Temperature Rise")
    st.latex(r"\Delta T_{ad} = \frac{\lvert \Delta H_{rxn} \rvert \, C_0}{\rho \, C_p}")
    st.markdown("""
The temperature the batch would reach if **all** reaction heat were retained
(no cooling).  Unlike the molar enthalpy alone, $\\Delta T_{ad}$ accounts for
reagent concentration and the thermal mass of the mixture, and is the basis of
the Stoessel criticality classification used to rank thermal-runaway risk.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\\Delta H_{rxn}$ | Enthalpy of reaction | J/mol (or kJ/mol × 1000) |
| $C_0$ | Limiting-reagent concentration | mol/m³ (= mol/L × 1000) |
| $\\rho \\, C_p$ | Volumetric heat capacity of the mixture | J/(m³·K) |
| $\\Delta T_{ad}$ | Adiabatic temperature rise | K |

**Indicative severity bands:** ΔT_ad < 20 K low · 20–50 K moderate ·
50–200 K high · > 200 K very high (assess MTSR and secondary decomposition).

**Reference:** Stoessel, F. (2020). [*Thermal Safety of Chemical Processes: Risk Assessment and Process Design*](https://doi.org/10.1002/9783527697854), 2nd ed. Wiley-VCH, Ch. 2–3.
""")

    # ── Heat removal
    st.header("Heat Removal Capacity (Jacket)")
    st.latex(r"Q_{\text{cool}} = U \, A \, \Delta T")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $U$ | Overall heat-transfer coefficient | W/(m²·K) |
| $A$ | Heat-transfer area (jacket + bottom) | m² |
| $\\Delta T$ | $T_{\\text{process}} - T_{\\text{coolant}}$ | °C or K |
| $Q_{\\text{cool}}$ | Maximum steady-state heat removal | W |

**Reference:** Incropera, F.P. *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), 6th ed. Wiley.
""")

    # ── Jacket area estimation
    st.header("Jacket Heat-Transfer Area Estimation")
    st.latex(r"A_{\text{cyl}} = \pi \, D_T \, H")
    st.latex(r"A_{\text{total}} = A_{\text{cyl}} + f_{dish} \cdot \frac{\pi}{4} D_T^2")
    st.markdown("""
| Bottom Dish Type | $f_{dish}$ |
|-----------------|------------|
| Flat | 1.00 |
| 2:1 Elliptical | 1.09 |
| Torispherical / DIN | 1.06 |
| Conical (~60°) | 1.20 |

**Reference:** DIN 28131 (Jacketed agitated vessels).
""")

    # ── U estimation
    st.header("Overall Heat-Transfer Coefficient Estimation")

    st.subheader("Simple Material-Based Estimate (Fallback)")
    st.markdown("""
| Wall Material | $U$ range (W/m²·K) |
|--------------|-------------------|
| Glass-lined | 100 – 250 |
| Stainless steel | 200 – 500 |
| Hastelloy | 200 – 450 |
| Carbon steel | 150 – 350 |

**Reference:** Perry, R.H. & Green, D.W. (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), 9th ed. McGraw-Hill, Section 11.
""")

    st.subheader("Individual Resistances Approach")
    st.latex(r"\frac{1}{U} = \frac{1}{h_i} + \frac{x_w}{k_w} + \frac{1}{h_o} + R_f")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $h_i$ | Process-side (internal) heat-transfer coefficient | W/(m²·K) |
| $x_w$ | Wall thickness | m |
| $k_w$ | Wall thermal conductivity | W/(m·K) |
| $h_o$ | Utility-side (jacket) heat-transfer coefficient | W/(m²·K) |
| $R_f$ | Fouling resistance | m²·K/W |

**Reference:** Incropera *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 3.
""")

    st.subheader("Process-Side Coefficient (Nusselt Correlation)")
    st.latex(r"\mathrm{Nu} = C \cdot \mathrm{Re}^{2/3} \cdot \mathrm{Pr}^{1/3} \cdot \left(\frac{\mu}{\mu_w}\right)^{0.14}")
    st.latex(r"h_i = \frac{\mathrm{Nu} \cdot k_f}{D_T}")
    st.latex(r"\mathrm{Re} = \frac{\rho \, N \, D_i^2}{\mu} \qquad \mathrm{Pr} = \frac{C_p \, \mu}{k_f}")
    st.markdown("""
| Symbol | Description | Typical value |
|--------|-------------|---------------|
| $C$ | Correlation constant | 0.36 (radial turbine, baffled) |
| $N$ | Impeller speed | rev/s |
| $D_i$ | Impeller diameter | m |
| $D_T$ | Tank diameter | m |
| $k_f$ | Fluid thermal conductivity | W/(m·K) |
| $C_p$ | Fluid specific heat capacity | J/(kg·K) |
| $\\mu / \\mu_w$ | Viscosity ratio (bulk / wall) | ≈ 1 for moderate ΔT |

**Reference:** DIN 28131; Seider, E.N. & Tate, G.E. (1936). Heat transfer and pressure drop of liquids in tubes. [*Ind. Eng. Chem.*](https://doi.org/10.1021/ie50324a027), 28(12), 1429–1435.
""")

    st.subheader("Wall Thermal Conductivity")
    st.markdown("""
| Material | $k_w$ (W/(m·K)) |
|----------|-----------------|
| Stainless steel (316/304) | 16 |
| Hastelloy C-276 | 12 |
| Carbon steel | 50 |
| Glass lining | 1.0 |
| Titanium | 22 |
| Inconel | 15 |

For **glass-lined** vessels:

$$R_{wall} = \\frac{x_{glass}}{k_{glass}} + \\frac{x_{steel}}{k_{steel}}$$

**Reference:** Perry & Green (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), Table 11-7.
""")

    st.subheader("Jacket-Side Coefficient")
    st.markdown("""
| Jacket Type | Typical $h_o$ (W/(m²·K)) |
|-------------|--------------------------|
| Simple jacket (water/glycol) | 1500 |
| Half-pipe coil | 2500 |
| Dimple jacket | 1200 |

**Reference:** Perry & Green (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), Section 11.
""")

    st.subheader("Fouling Resistance")
    st.markdown("""
| Condition | $R_f$ (m²·K/W) |
|-----------|----------------|
| Clean (fresh equipment) | 0.0001 |
| Moderate (pharmaceutical) | 0.0002 |
| Heavy fouling | 0.001 |

Default value: **0.0002 m²·K/W** (moderate pharmaceutical process).

**Reference:** Perry & Green (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), Table 11-9.
""")

    # ── Heat balance assessment
    st.header("Heat Balance Assessment")
    st.latex(r"\text{Ratio} = \frac{Q_{\text{gen}}}{Q_{\text{cool}}}")
    st.markdown("""
| $Q_{gen} / Q_{cool}$ | Assessment |
|---------------------|------------|
| < 0.25 | Easily manageable |
| 0.25 – 0.50 | Comfortable margin |
| 0.50 – 0.75 | Moderate – monitor closely |
| 0.75 – 1.00 | ⚠️ Tight – limited safety margin |
| > 1.00 | 🔴 Insufficient cooling capacity |
""")

    # ── Cooling Rate
    st.header("Cooling Rate")
    st.latex(r"\frac{dT}{dt} = \frac{Q_{\text{cool}} - P_{\text{agitator}}}{\rho \, V \, C_p}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $Q_{\\text{cool}}$ | Jacket cooling capacity | W |
| $P_{\\text{agitator}}$ | Impeller power draw | W |
| $\\rho$ | Fluid density | kg/m³ |
| $V$ | Liquid volume | m³ |
| $C_p$ | Specific heat capacity | J/(kg·K) |

The instantaneous rate of temperature change of the batch, accounting
for the heat input from agitator dissipation.

**Reference:** Incropera *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 5.
""")

    # ── Time to Cool / Heat
    st.header("Time to Cool or Heat (Batch)")
    st.latex(r"t = \frac{\rho \, V \, C_p}{U \, A} \, \ln \frac{T_{\text{start}} - T_{\text{jacket}}}{T_{\text{end}} - T_{\text{jacket}}}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\rho$ | Fluid density | kg/m³ |
| $V$ | Liquid volume | m³ |
| $C_p$ | Specific heat capacity | J/(kg·K) |
| $U$ | Overall heat-transfer coefficient | W/(m²·K) |
| $A$ | Heat-transfer area | m² |
| $T_{\\text{start}}$ | Initial batch temperature | °C |
| $T_{\\text{end}}$ | Target batch temperature | °C |
| $T_{\\text{jacket}}$ | Jacket / coolant temperature | °C |

Standard logarithmic solution for a well-mixed batch with constant $U$, $A$,
and jacket temperature.

**Reference:** Incropera *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 5.
""")


# ══════════════════════════════════════════════════════════════════════════
# 5. SOLID PARTICLES
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Solid-Particle Hydrodynamics**", expanded=False):

    # ── Terminal Settling Velocity
    st.header("Terminal Settling Velocity")
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

**Reference:** Clift, R., Grace, J.R. & Weber, M.E. (1978). [*Bubbles, Drops, and Particles*](https://openlibrary.org/isbn/012176950X), Academic Press; Schiller, L. & Naumann, A. (1935). *Z. Ver. Dtsch. Ing.*, 77, 318–320.
""")

    st.header("Particle Reynolds Number")
    st.latex(r"Re_p = \frac{\rho_L \, v_t \, d_p}{\mu}")
    st.markdown("""
Classifies the flow regime around a settling particle.

**Reference:** Clift *et al.* (1978).
""")

    st.header("Zwietering Just-Suspended Speed")
    st.latex(r"N_{js} = S \, \nu^{0.1} \, d_p^{0.2} \left(\frac{g\,\Delta\rho}{\rho_L}\right)^{0.45} X^{0.13} \, D^{-0.85}")
    st.markdown(r"""
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

**Reference:** Zwietering, Th.N. (1958). Suspending of solid particles in liquid by agitators. [*Chem. Eng. Sci.*](https://doi.org/10.1016/0009-2509(58)85031-9), 8(3–4), 244–253.
""")

    # ── GMB Njs
    st.header("Grenville, Mak & Brown (GMB) Just-Suspended Speed")
    st.latex(r"N_{js} = z \, Po^{-1/3} \, D^{-2/3} \left(\frac{g\,\Delta\rho}{\rho_L}\right)^{0.45} X_v^{0.154} \, d_p^{0.167} \left(\frac{C}{D}\right)^{0.1}")
    st.markdown(r"""
An alternative to Zwietering that uses the **power number** $Po$ and
**impeller clearance ratio** $C/D$ explicitly.

| Symbol | Description | Units |
|--------|-------------|-------|
| $z$ | Geometry constant (impeller-type dependent) | – |
| $Po$ | Power number ($N_p$) | – |
| $D$ | Impeller diameter | m |
| $g$ | Gravitational acceleration | m/s² |
| $\Delta\rho$ | $\lvert\rho_p - \rho_L\rvert$ | kg/m³ |
| $\rho_L$ | Liquid density | kg/m³ |
| $X_v$ | Volume percent of solids | vol-% |
| $d_p$ | Particle diameter | m |
| $C/D$ | Impeller clearance / impeller diameter | – |

Note: Zwietering uses mass-based solids loading $X$ (wt-%), while GMB uses volume-based $X_v$ (vol-%).

**Reference:** Grenville, R.K., Mak, A.T.C. & Brown, D.A.R. (2015). Suspension of solid particles in vessels agitated by axial flow impellers. [*Chem. Eng. Res. Des.*](https://doi.org/10.1016/j.cherd.2015.07.009), 100, 282–291.
""")

    st.header("Solid-Liquid Mass Transfer (Ranz-Marshall)")
    st.latex(r"Sh = 2 + 0.6 \, Re_p^{1/2} \, Sc^{1/3}")
    st.latex(r"k_{SL} = \frac{Sh \cdot D_{mol}}{d_p}")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $Sh$ | Sherwood number | – |
| $Re_p$ | Particle Reynolds number | – |
| $Sc = \nu / D_{mol}$ | Schmidt number | – |
| $D_{mol}$ | Molecular diffusivity | m²/s |
| $d_p$ | Particle diameter | m |
| $k_{SL}$ | Solid-liquid mass transfer coefficient | m/s |

**Reference:** Ranz, W.E. & Marshall, W.R. (1952). Evaporation from drops. [*Chem. Eng. Prog.*](https://en.wikipedia.org/wiki/Ranz%E2%80%93Marshall_correlation), 48(3), 141–146; 48(4), 173–180.
""")

    st.header("Suspension Quality Assessment")
    st.markdown(r"""
| $N / N_{js}$ | Assessment |
|-------------|------------|
| $< 0.7$ | Poorly suspended |
| $0.7 – 1.0$ | Partially suspended |
| $1.0 – 1.3$ | Just suspended |
| $> 1.3$ | Fully / homogeneously suspended |

Operating at $N \ge N_{js}$ ensures complete off-bottom suspension.

**Reference:** Zwietering (1958); Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 10.
""")

    # ── Archimedes Number
    st.header("Archimedes Number")
    st.latex(r"Ar = \frac{g \, d_p^3 \, \rho_L \, \Delta\rho}{\mu^2}")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $g$ | Gravitational acceleration | m/s² |
| $d_p$ | Particle diameter | m |
| $\rho_L$ | Liquid density | kg/m³ |
| $\Delta\rho$ | $\lvert\rho_p - \rho_L\rvert$ | kg/m³ |
| $\mu$ | Dynamic viscosity | Pa·s |

Characterises the balance between gravitational and viscous forces for
a particle in a fluid.

**Reference:** Clift *et al.* (1978). [*Bubbles, Drops, and Particles*](https://openlibrary.org/isbn/012176950X), Academic Press.
""")

    # ── Solid-Liquid kLa
    st.header("Solid-Liquid Volumetric Mass Transfer (kLa)")
    st.latex(r"a_s = \frac{6 \, \phi_s}{d_p}")
    st.latex(r"k_L a_{SL} = k_{SL} \cdot a_s")
    st.markdown(r"""
| Symbol | Description | Units |
|--------|-------------|-------|
| $a_s$ | Specific surface area of particles per unit liquid volume | 1/m |
| $\phi_s$ | Volumetric solids fraction | – |
| $d_p$ | Particle diameter | m |
| $k_{SL}$ | Solid-liquid mass transfer coefficient (from Ranz-Marshall) | m/s |
| $k_L a_{SL}$ | Volumetric mass transfer coefficient | 1/s |

**Reference:** Ranz & Marshall (1952); Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 10.
""")


# ══════════════════════════════════════════════════════════════════════════
# 6. SCALE-UP & SOLVENT PROPERTIES
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Scale-Up & Solvent Properties**", expanded=False):

    # ── Scale-up rules
    st.header("Common Scale-Up Rules")
    st.markdown("Given geometric similarity ($D/T$ = const), the large-scale speed $N_L$ is related to the small-scale $N_S$:")

    st.subheader("Constant tip speed")
    st.latex(r"N_L = N_S \frac{D_S}{D_L}")

    st.subheader("Constant P/V")
    st.latex(r"N_L = N_S \left( \frac{D_S}{D_L} \right)^{2/3}")

    st.subheader("Constant Re")
    st.latex(r"N_L = N_S \left( \frac{D_S}{D_L} \right)^{2}")

    st.markdown("""
**References:**
- Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Wiley.
- Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
- Grenville, R.K. & Nienow, A.W. (2004). Blending of miscible liquids. In *Handbook of Industrial Mixing*, Ch. 9.
""")

    # ── Solvent Properties (T-dependent)
    st.header("Temperature-Dependent Solvent Properties")
    st.markdown("""
The Fluid Database's *Solvent Properties (T)* tab provides physical properties
at any liquid-phase temperature for common pharmaceutical solvents.
""")

    st.subheader("Density")
    st.latex(r"\rho(T) = \rho_{25} + \frac{d\rho}{dT} \cdot (T - 25) \quad [\mathrm{kg/m^3}]")
    st.markdown("Linear fit anchored at the known 25 °C value.  $T$ in °C, slope $d\\rho/dT$ (typically negative) from published data.")

    st.subheader("Dynamic Viscosity (Arrhenius)")
    st.latex(r"\mu(T) = \mu_{25} \, \exp\!\left[ \frac{E_a}{R} \left( \frac{1}{T_K} - \frac{1}{298.15} \right) \right]")
    st.markdown("""
$T_K$ = temperature in Kelvin, $E_a$ = activation energy for viscous flow
(J/mol), $R$ = 8.314 J/(mol·K).

**Reference:** Perry, R.H. & Green, D.W. (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), 9th ed. McGraw-Hill; CRC Handbook of Chemistry and Physics; DIPPR correlations.
""")

    st.subheader("Surface Tension")
    st.latex(r"\sigma(T) = \sigma_{25} + \frac{d\sigma}{dT} \cdot (T - 25) \quad [\mathrm{N/m}]")
    st.markdown("Linear fit anchored at the known 25 °C value.  Clamped to $\\sigma \\ge 0$.")

    st.subheader("Molecular Diffusivity (Stokes-Einstein Scaling)")
    st.latex(r"D(T) = D_{\mathrm{ref}} \; \frac{T_K}{298.15} \; \frac{\mu_{\mathrm{ref}}}{\mu(T)}")
    st.markdown("""
Scaled from a reference value $D_{\\mathrm{ref}}$ at 25 °C using the
Stokes-Einstein relation ($D \\propto T / \\mu$).

**Reference:** Yaws' Handbook; Perry's (2019); CRC Handbook; DIPPR correlations.
""")

    # ── Antoine Equation
    st.header("Pressure-Dependent Boiling Point (Antoine Equation)")

    st.subheader("Vapor Pressure")
    st.latex(r"\log_{10}\!\left(P_{\mathrm{mmHg}}\right) = A - \frac{B}{C + T\;(°C)}")
    st.markdown("""
where $A$, $B$, $C$ are substance-specific constants.
""")

    st.subheader("Boiling Point at Pressure P")
    st.latex(r"T_{\mathrm{bp}}\;(°C) = \frac{B}{A - \log_{10}(P_{\mathrm{mmHg}})} - C")
    st.markdown("""
At 1 atm ($P$ = 760 mmHg), this recovers the normal boiling point.

**Conversion:** 1 atm = 760 mmHg = 101.325 kPa.

**Reference:** [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/); Yaws' Handbook of Vapor Pressure (Antoine Constants).
""")


# ══════════════════════════════════════════════════════════════════════════
# 7. ROM & EXPERIMENTAL CORRELATIONS
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**ROM & Experimental Correlations**", expanded=False):

    st.markdown("""
In addition to the standard **literature correlations** used throughout
Mixing Lab, reactor-specific correlations can be registered from:

- **ROM (Reduced Order Models)** – fitted from CFD parametric studies or
  surrogate models.
- **Experimental** – fitted directly to measured data (e.g. kLa from
  gassing-out experiments, blend time from PLIF/decolourisation).

These override the literature correlation for a specific parameter when
the corresponding mode is selected on the **Mixing Sensitivity** or
**Reactor Comparison** pages.
""")

    # ── Supported parameters
    st.subheader("Supported Parameters")
    st.markdown(
        "The following parameters can be overridden by ROM or Experimental "
        "correlations:\n\n"
        + "\n".join(f"- {PARAM_DISPLAY[p]}" for p in SUPPORTED_PARAMS)
    )

    # ── Registered correlations
    all_corrs = get_all_correlations()
    registered_reactors = get_all_registered_reactors()

    if not registered_reactors:
        st.info("No reactor-specific correlations have been registered yet.")
    else:
        st.subheader("Registered Correlations")

        for reactor_name in registered_reactors:
            corrs = all_corrs[reactor_name]
            st.markdown(f"**🔧 {reactor_name}**")

            for corr in corrs:
                label = f"{corr.corr_type} · {PARAM_DISPLAY.get(corr.param, corr.param)}"
                with st.container():
                    st.markdown(f"***{label}***")
                    st.markdown(f"**{corr.name}**")
                    st.latex(corr.latex)

                    if corr.description:
                        st.markdown(corr.description)

                    if corr.input_params:
                        rows = "\n".join(
                            f"| `{k}` | {v} |" for k, v in corr.input_params.items()
                        )
                        st.markdown(
                            "| Input | Description |\n"
                            "|-------|-------------|\n"
                            + rows
                        )

                    st.caption(f"Source: {corr.source}")

    st.markdown(
        "💡 **Adding new correlations:** Register them in "
        "`utils/rom_registry.py` using the `register()` function.  "
        "They will automatically appear here and become selectable on "
        "the calculation pages."
    )


# ══════════════════════════════════════════════════════════════════════════
# 8. CRYSTALLIZATION
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Crystallization**", expanded=False):

    # ── Supersaturation
    st.header("Supersaturation")

    st.subheader("Relative Supersaturation")
    st.latex(r"\sigma = \frac{c - c_{sat}}{c_{sat}}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $c$ | Actual solute concentration | g/L (or mol/L) |
| $c_{sat}$ | Equilibrium saturation concentration at temperature $T$ | g/L (or mol/L) |
| $\\sigma$ | Relative supersaturation (dimensionless) | — |

**Reference:** Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). [*Handbook of Industrial Crystallization*](https://doi.org/10.1017/9781139026949), 3rd ed. Cambridge University Press, Ch. 2.
""")

    st.subheader("Supersaturation Ratio")
    st.latex(r"S = \frac{c}{c_{sat}} = 1 + \sigma")
    st.markdown("""
An alternative definition commonly used in nucleation theory.
$S = 1$ at equilibrium; $S > 1$ indicates a supersaturated solution.

**Reference:** Myerson *et al.* (2019), Ch. 2.
""")

    st.subheader("Solubility Temperature Dependence (Linear Approximation)")
    st.latex(r"c_{sat}(T) \approx c_{sat}(T_{ref}) + \frac{dc}{dT} \cdot (T - T_{ref})")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $dc/dT$ | Slope of the solubility curve near $T_{ref}$ | g/(L·K) |
| $T_{ref}$ | Reference temperature for solubility data | °C |

For cooling crystallization, MSZW can be converted to a maximum supersaturation:

$$\\sigma_{max} \\approx \\frac{(dc/dT) \\cdot \\text{MSZW}}{c_{sat}(T_{ref})}$$

**Reference:** Myerson *et al.* (2019), Ch. 2.
""")

    # ── MSZW
    st.header("Metastable Zone Width (MSZW)")
    st.latex(r"\text{MSZW} = T_{sat} - T_{nuc}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $T_{sat}$ | Saturation (equilibrium) temperature | °C |
| $T_{nuc}$ | Temperature at which nucleation is first detected | °C |
| MSZW | Metastable zone width | °C |

**Reference:** Myerson *et al.* (2019), Ch. 5.
""")

    # ── Primary Nucleation
    st.header("Primary Nucleation Rate")

    st.subheader("Classical Nucleation Theory (CNT)")
    st.latex(r"J = A \exp\!\left( -\frac{16 \pi \gamma^3 v_m^2}{3 k_B^3 T^3 (\ln S)^2} \right)")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $J$ | Nucleation rate | #/(m³·s) |
| $A$ | Pre-exponential kinetic factor | #/(m³·s) |
| $\\gamma$ | Crystal–solution interfacial energy | J/m² |
| $v_m$ | Molecular volume of solute in the crystal | m³ |
| $k_B$ | Boltzmann constant ($1.381 \\times 10^{-23}$) | J/K |
| $T$ | Absolute temperature | K |
| $S$ | Supersaturation ratio $c/c_{sat}$ | — |

**Reference:** Myerson *et al.* (2019), Ch. 3; Mersmann, A. (2001). [*Crystallization Technology Handbook*](https://doi.org/10.1201/9780203908280), 2nd ed. Marcel Dekker, Ch. 2.
""")

    st.subheader("Power-Law Approximation")
    st.latex(r"B = k_b \, \sigma^{\,n}")
    st.markdown("""
| Symbol | Description | Typical range |
|--------|-------------|---------------|
| $B$ | Nucleation rate | #/(m³·s) |
| $k_b$ | Nucleation rate constant | system-specific |
| $n$ | Nucleation order | 2–10 (primary), 1–3 (secondary) |

Because $n \\gg g$ (the growth order, typically 1–2), **nucleation is
far more sensitive to local $\\sigma$ than growth**.

**Reference:** Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
""")

    # ── Induction Time
    st.header("Induction Time")
    st.latex(r"t_{ind} \propto \frac{1}{J} \propto \exp\!\left( +\frac{16 \pi \gamma^3 v_m^2}{3 k_B^3 T^3 (\ln S)^2} \right)")
    st.markdown("""
The induction time $t_{ind}$ is the experimentally observed delay between
achieving supersaturation and detecting the first nuclei.

**Reference:** Myerson *et al.* (2019), Ch. 5; Kashchiev, D. & van Rosmalen, G.M. (2003). Review: Nucleation in solutions revisited. [*Cryst. Res. Technol.*](https://doi.org/10.1002/crat.200310070), 38(7–8), 555–574.
""")

    # ── Crystal Growth
    st.header("Crystal Growth Rate")

    st.subheader("Power-Law Growth Model")
    st.latex(r"G = k_g \, \sigma^{\,g}")
    st.markdown("""
| Symbol | Description | Typical range |
|--------|-------------|---------------|
| $G$ | Linear growth rate (face velocity) | m/s |
| $k_g$ | Growth rate constant | m/s (system-specific) |
| $g$ | Growth order | 1–2 |

**Reference:** Myerson *et al.* (2019), Ch. 6.
""")

    st.subheader("Growth Time (Characteristic)")
    st.latex(r"t_G = \frac{L_{target}}{G} = \frac{L_{target}}{k_g \, \sigma^{\,g}}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $L_{target}$ | Target crystal size | m (or µm) |
| $t_G$ | Time to grow a crystal to the target size | s |

**Reference:** Green, D.A. (2019). Ch. 10 in Myerson *et al.* (2019).
""")

    # ── Damköhler Numbers for Crystallization
    st.header("Damköhler Numbers for Crystallization")
    st.markdown("""
The Damköhler number framework from Baldyga & Bourne (1999) adapted for
crystallization by replacing $t_{rxn}$ with crystallization time scales.
""")

    st.subheader("Micro Damköhler Number (Nucleation)")
    st.latex(r"Da_{micro} = \frac{t_E}{t_{ind}}")
    st.markdown("""
| Symbol | Description |
|--------|-------------|
| $t_E$ | Engulfment micromixing time $= 17.3 \\sqrt{\\nu / \\varepsilon}$ |
| $t_{ind}$ | Induction time at target $\\sigma$ |

**Reference:** Baldyga & Bourne (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
""")

    st.subheader("Macro Damköhler Number (Growth)")
    st.latex(r"Da_{macro} = \frac{\theta_{95}}{t_G}")
    st.markdown("""
| Symbol | Description |
|--------|-------------|
| $\\theta_{95}$ | 95% macro-blend time $= 5.2\\,V / (N_Q\\,N\\,D^3)$ |
| $t_G$ | Growth time to target crystal size |

**Reference:** Baldyga & Bourne (1999); Paul *et al.* (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Ch. 13.
""")

    st.subheader("Interpretation")
    st.markdown("""
| $Da$ range | Interpretation |
|-----------|----------------|
| $Da < 0.01$ | Crystallization much slower than mixing — **not sensitive** |
| $0.01 < Da < 0.1$ | **Likely not sensitive**, but monitor at scale |
| $0.1 < Da < 1$ | **Potentially sensitive** — similar timescales |
| $1 < Da < 10$ | **Likely sensitive** — mixing limits process |
| $Da > 10$ | **Highly sensitive** — mixing fully controls outcome |
""")

    # ── Nucleation vs Growth Competition
    st.header("Nucleation–Growth Competition")
    st.latex(r"B \propto \sigma^{\,n} \qquad (n \sim 2{-}10)")
    st.latex(r"G \propto \sigma^{\,g} \qquad (g \sim 1{-}2)")
    st.markdown("""
Because $n \\gg g$, nucleation accelerates disproportionately with local
supersaturation.  Poor mixing → excess fines, smaller mean crystal size,
broader CSD, and potentially incorrect polymorph.

**Reference:** Baldyga & Bourne (1999); Green, D.A. (2019), Ch. 10 in Myerson *et al.* (2019).
""")

    # ── Ostwald's Rule
    st.header("Ostwald's Rule of Stages")
    st.markdown("""
The **least stable (metastable) polymorph** nucleates first due to lower
interfacial energy $\\gamma$ and hence lower nucleation barrier.

**Mixing implication:** Poor mixing → local $\\sigma$ spikes → metastable
polymorph nucleation → potential quality failure.

**Reference:** Myerson *et al.* (2019), Ch. 3.
""")

    # ── Crystallization Enthalpy
    st.header("Crystallization Enthalpy")
    st.latex(r"Q_{cryst} = |\Delta H_{cryst}| \times \dot{n}_{cryst}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\Delta H_{cryst}$ | Enthalpy of crystallization (usually exothermic) | kJ/mol |
| $\\dot{n}_{cryst}$ | Molar crystallization rate | mol/s |
| $Q_{cryst}$ | Heat generation from crystallization | W |

| $|\\Delta H_{cryst}|$ | Risk level |
|-----------------------|-----------|
| < 20 kJ/mol | 🟢 Low |
| 20–40 kJ/mol | 🟡 Moderate |
| > 40 kJ/mol | 🔴 High — heat removal may limit rates |

**Reference:** Myerson *et al.* (2019), Ch. 5; Mersmann (2001). [*Crystallization Technology Handbook*](https://doi.org/10.1201/9780203908280).
""")

    st.subheader("Van 't Hoff Estimation of ΔH_cryst")
    st.latex(r"\frac{d \ln c_{sat}}{d(1/T)} = -\frac{\Delta H_{cryst}}{R}")
    st.markdown("""
| Symbol | Description |
|--------|-------------|
| $R$ | Gas constant (8.314 J/(mol·K)) |
| $c_{sat}$ | Saturation concentration at temperature $T$ |

**Reference:** Myerson *et al.* (2019), Ch. 2.
""")

    # ── Seeding
    st.header("Seeding — Effect on Mixing Sensitivity")
    st.markdown("""
Seed crystals provide pre-existing surface area $A_{seed}$ that consumes
supersaturation via secondary nucleation and growth.

$$\\frac{d\\sigma}{dt} \\propto -k_g \\, \\sigma^g \\, A_{seed}$$

**Rule of thumb:** A well-seeded process is 2–5× less sensitive to mixing.

**Reference:** Black, S.N. (2019). Ch. 13 in Myerson *et al.* (2019). [*Handbook of Industrial Crystallization*](https://doi.org/10.1017/9781139026949).
""")

    # ── Mixing-Sensitive Crystallization Types
    st.header("Mixing Sensitivity by Crystallization Type")
    st.markdown("""
| Type | How $\\sigma$ is generated | Dominant mixing concern | Inherent risk |
|------|--------------------------|------------------------|---------------|
| **Cooling** | Temperature reduction | Macromixing — uniformity of $T$ field | 🟢 Low |
| **Anti-solvent** | Addition of a poor solvent | Meso- & micromixing — feed plume $\\sigma$ spike | 🔴 High |
| **Reactive** | Reaction produces insoluble product | Micromixing — instantaneous local $\\sigma$ | 🔴 High |
| **Evaporative** | Solvent removal at boiling surface | Macromixing — boiling zone $\\sigma$ gradient | 🟡 Moderate |
| **pH-shift** | Acid/base addition changes solubility | Meso- & micromixing — local pH plume | 🔴 High |

**Reference:** Green, D.A. (2019), Ch. 10 in Myerson *et al.* (2019). [*Handbook of Industrial Crystallization*](https://doi.org/10.1017/9781139026949).
""")

    # ── References
    st.header("Key References")
    st.markdown("""
1. Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). [*Handbook of Industrial Crystallization*](https://doi.org/10.1017/9781139026949), 3rd ed. Cambridge University Press.
2. Baldyga, J. & Bourne, J.R. (1999). [*Turbulent Mixing and Chemical Reactions*](https://openlibrary.org/isbn/0471981710), Wiley.
3. Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). [*Handbook of Industrial Mixing*](https://doi.org/10.1002/0471451452), Wiley.
4. Mersmann, A. (2001). [*Crystallization Technology Handbook*](https://doi.org/10.1201/9780203908280), 2nd ed. Marcel Dekker.
5. Kashchiev, D. & van Rosmalen, G.M. (2003). [*Cryst. Res. Technol.*](https://doi.org/10.1002/crat.200310070), 38(7–8), 555–574.
""")


# ══════════════════════════════════════════════════════════════════════════
# HEAT TRANSFER TOOL
# ══════════════════════════════════════════════════════════════════════════
with st.expander("**Heat Transfer Tool**", expanded=False):

    st.markdown("""
Equations used in the **Heat Transfer Modeling** page (Tools → Heat Transfer).
This page builds on the basic heat-balance equations above with multiple
Nusselt correlations, jacket-side modeling, transient simulations, and
NTU-effectiveness analysis.
""")

    # ── Process-side Nusselt correlations
    st.header("Process-Side Nusselt Correlations")
    st.latex(r"\mathrm{Nu} = C \cdot \mathrm{Re}^{a} \cdot \mathrm{Pr}^{b} \cdot \left(\frac{\mu}{\mu_w}\right)^{c}")
    st.latex(r"h_i = \frac{\mathrm{Nu} \cdot k_f}{D_T}")
    st.markdown("""
| Correlation | $C$ | $a$ | $b$ | $c$ | Application | Reference |
|-------------|-----|-----|-----|-----|-------------|-----------|
| DIN 28131 (standard) | 0.36 | 2/3 | 1/3 | 0.14 | General-purpose, baffled vessels with turbine impellers | DIN 28131:1979 |
| Chilton–Drew–Jebens | 0.36 | 2/3 | 1/3 | 0.14 | Classic correlation for jacketed stirred vessels | Chilton, Drew, Jebens (1944) *Ind. Eng. Chem.* 36(6):510 |
| Lehrer (anchor/helical) | 0.54 | 2/3 | 1/3 | 0.14 | Anchor and helical ribbon impellers | Lehrer (1970) *Chem. Eng. Sci.* 25:1397 |
| Stein–Schmidt (high Re) | 0.50 | 2/3 | 1/3 | 0.14 | High-Re turbulent regimes (Re > 40 000) | Stein & Schmidt (1993) *Chem. Eng. Process.* 32:305 |
| Brooks–Su (Retreat Blade) | 0.33 | 2/3 | 1/3 | 0.14 | Retreat-blade impellers in glass-lined vessels | Brooks & Su (1959) |
| Nagata (paddle) | 0.36 | 2/3 | 1/3 | 0.18 | Paddle impellers with stronger wall viscosity correction | Nagata (1975) *Mixing: Principles and Applications* |

Where:

$$\\mathrm{Re} = \\frac{\\rho \\, N \\, D_i^2}{\\mu} \\qquad \\mathrm{Pr} = \\frac{C_p \\, \\mu}{k_f}$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $C, a, b, c$ | Correlation constants | — |
| $N$ | Impeller speed | rev/s |
| $D_i$ | Impeller diameter | m |
| $D_T$ | Tank diameter | m |
| $k_f$ | Fluid thermal conductivity | W/(m·K) |
| $C_p$ | Specific heat capacity | J/(kg·K) |
| $\\mu$ | Bulk viscosity | Pa·s |
| $\\mu_w$ | Viscosity at wall temperature | Pa·s |
""")

    # ── Overall U from individual resistances
    st.header("Overall Heat-Transfer Coefficient (Resistance Model)")
    st.latex(r"\frac{1}{U} = \frac{1}{h_i} + \frac{x_w}{k_w} + \frac{x_l}{k_l} + \frac{1}{h_o} + R_f")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $h_i$ | Process-side heat-transfer coefficient | W/(m²·K) |
| $x_w$ | Wall thickness | m |
| $k_w$ | Wall thermal conductivity | W/(m·K) |
| $x_l$ | Lining thickness | m |
| $k_l$ | Lining thermal conductivity | W/(m·K) |
| $h_o$ | Jacket-side (utility) heat-transfer coefficient | W/(m²·K) |
| $R_f$ | Fouling resistance | m²·K/W |

When a lining is present (glass, PTFE, etc.), both wall and lining resistances
are summed. Typical lining conductivities:

| Lining | $k_l$ (W/(m·K)) | Default thickness (mm) |
|--------|-----------------|----------------------|
| Glass | 1.0 | 1.5 |
| PTFE / Teflon / PFA | 0.25 | 2.0 |
| PVDF | 0.19 | 3.0 |
| Rubber | 0.16 | 6.0 |
| Epoxy | 0.20 | 3.0 |

**Reference:** Incropera, F.P. *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 3.
""")

    # ── Jacket-side h_o
    st.header("Jacket-Side Heat-Transfer Coefficient")

    st.subheader("Turbulent Flow — Dittus-Boelter Equation")
    st.latex(r"\mathrm{Nu}_j = 0.023 \cdot \mathrm{Re}_j^{0.8} \cdot \mathrm{Pr}_j^{0.4} \qquad (\mathrm{Re}_j > 2300)")
    st.latex(r"h_o = \frac{\mathrm{Nu}_j \cdot k_j}{D_h}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\mathrm{Re}_j$ | Jacket Reynolds number $= \\rho_j \\, v_j \\, D_h / \\mu_j$ | — |
| $\\mathrm{Pr}_j$ | Jacket Prandtl number $= C_{p,j} \\, \\mu_j / k_j$ | — |
| $v_j$ | Jacket fluid velocity | m/s |
| $D_h$ | Hydraulic diameter of jacket annulus | m |
| $k_j$ | Jacket fluid thermal conductivity | W/(m·K) |

**Reference:** Dittus, F.W. & Boelter, L.M.K. (1930). *Univ. Calif. Publ. Eng.* 2:443; Incropera *et al.* (2007), Ch. 8.
""")

    st.subheader("Laminar Flow")
    st.latex(r"\mathrm{Nu}_j = 3.66 + \frac{0.065 \,(D_h / L)\, \mathrm{Re}_j \, \mathrm{Pr}_j}{1 + 0.04 \,\left[(D_h / L)\, \mathrm{Re}_j \, \mathrm{Pr}_j\right]^{2/3}} \qquad (\mathrm{Re}_j < 2300)")
    st.markdown("""
Combined developing / fully-developed laminar flow correlation.

**Reference:** Sieder, E.N. & Tate, G.E. (1936). *Ind. Eng. Chem.* 28(12):1429.
""")

    st.subheader("Condensing Steam")
    st.markdown("""
For condensing steam, the jacket-side coefficient is very high
($h_o \\approx$ 6 000 – 12 000 W/(m²·K)) and is treated as a fixed value
rather than computed from a Nusselt correlation.

**Default:** $h_o$ = 8 000 W/(m²·K) for condensing steam.
""")

    # ── Impeller power
    st.header("Impeller Power (Agitator Heat Input)")
    st.latex(r"P = N_p \, \rho \, N^3 \, D_i^5")
    st.markdown("""
All mechanical energy dissipated by the impeller is ultimately converted to
heat in the fluid. This is added as a heat source in the transient energy balance.

| Symbol | Description | Units |
|--------|-------------|-------|
| $N_p$ | Power number | — |
| $\\rho$ | Fluid density | kg/m³ |
| $N$ | Impeller speed | rev/s |
| $D_i$ | Impeller diameter | m |
| $P$ | Impeller power draw | W |
""")

    # ── Steady-state heat duty
    st.header("Steady-State Heat Duty")
    st.latex(r"Q = U \, A \, \Delta T")
    st.latex(r"\Delta T = |T_{\text{batch}} - T_{\text{jacket}}|")
    st.markdown("""
The maximum instantaneous heat removal/addition rate at any point in time.
As the batch temperature approaches the jacket temperature, $\\Delta T$ and
therefore $Q$ decrease exponentially.

**Reference:** Incropera *et al.* (2007), Ch. 3.
""")

    # ── Instantaneous cooling/heating rate
    st.header("Instantaneous Cooling / Heating Rate")
    st.latex(r"\frac{dT}{dt} = \frac{Q_{\text{jacket}} - P_{\text{agitator}} - Q_{\text{rxn}}}{\rho \, V \, C_p}")
    st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $Q_{\\text{jacket}}$ | Heat transfer through jacket (positive = cooling) | W |
| $P_{\\text{agitator}}$ | Impeller power dissipation (always heats) | W |
| $Q_{\\text{rxn}}$ | Reaction heat (positive = exothermic, heats) | W |
| $\\rho \\, V \\, C_p$ | Thermal mass of batch | J/K |
""")

    # ── Analytical time (log-mean)
    st.header("Analytical Batch Time (Log-Mean)")
    st.latex(r"t = \frac{\rho \, V \, C_p}{U \, A} \, \ln \frac{T_{\text{start}} - T_{\text{jacket}}}{T_{\text{end}} - T_{\text{jacket}}}")
    st.markdown("""
Assumes constant $U$, $A$, and jacket temperature, with no agitator heat or
reaction heat. This is the standard analytical solution for an
exponentially-decaying driving force.

**Reference:** Incropera *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 5.
""")

    # ── Transient model: constant jacket
    st.header("Transient Model — Constant Jacket Temperature (Euler Integration)")
    st.latex(r"\frac{dT}{dt} = \frac{U \, A \,(T_{\text{jacket}} - T) + P_{\text{agitator}} + Q_{\text{rxn}}}{\rho \, V \, C_p}")
    st.latex(r"T(t + \Delta t) = T(t) + \frac{dT}{dt} \cdot \Delta t")
    st.markdown("""
A forward Euler integration of the energy balance at each time step.

**Sign convention:**
- $U \\cdot A \\cdot (T_{\\text{jacket}} - T)$: positive when jacket is hotter (heating), negative when cooler (cooling)
- $P_{\\text{agitator}}$: always positive (adds heat to batch)
- $Q_{\\text{rxn}}$: positive = exothermic (adds heat)

Integration stops when $T$ reaches $T_{\\text{target}}$ or $t$ reaches $t_{\\text{max}}$.
""")

    # ── Transient model: variable jacket (NTU-effectiveness)
    st.header("Transient Model — Variable Jacket Temperature (NTU-Effectiveness)")
    st.latex(r"\mathrm{NTU} = \frac{U \, A}{\dot{m}_j \, C_{p,j}}")
    st.latex(r"\varepsilon = 1 - e^{-\mathrm{NTU}}")
    st.latex(r"Q_{\text{jacket}} = \varepsilon \, \dot{m}_j \, C_{p,j} \, (T_{j,\text{in}} - T_{\text{batch}})")
    st.latex(r"T_{j,\text{out}} = T_{j,\text{in}} + \frac{Q_{\text{jacket}}}{\dot{m}_j \, C_{p,j}}")
    st.markdown("""
When the jacket fluid has a finite flow rate, the outlet temperature differs
from the inlet. The NTU-effectiveness method models this:

| Symbol | Description | Units |
|--------|-------------|-------|
| $\\mathrm{NTU}$ | Number of Transfer Units | — |
| $\\varepsilon$ | Heat exchanger effectiveness | — |
| $\\dot{m}_j$ | Jacket mass flow rate | kg/s |
| $C_{p,j}$ | Jacket fluid specific heat capacity | J/(kg·K) |
| $T_{j,\\text{in}}$ | Jacket inlet temperature (constant) | °C |
| $T_{j,\\text{out}}$ | Jacket outlet temperature (varies) | °C |

At each time step, $Q_{\\text{jacket}}$ is computed using the current batch
temperature and the effectiveness, then the batch temperature is updated
via the same Euler integration as the constant-jacket model.

**Key insight:** Lower jacket flow rates → higher NTU → higher effectiveness
but also higher jacket outlet temperature, which reduces the effective $\\Delta T$.

**Reference:** Incropera *et al.* (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), Ch. 11 (Heat Exchangers).
""")

    # ── Heat transfer media
    st.header("Heat Transfer Media Database")
    st.markdown("""
The tool includes a built-in database of common heat transfer fluids:

| Medium | T range (°C) | ρ (kg/m³) | Cp (J/(kg·K)) | μ (Pa·s) | k (W/(m·K)) |
|--------|-------------|-----------|---------------|----------|-------------|
| Water | 5 – 95 | 997 | 4182 | 8.9×10⁻⁴ | 0.607 |
| Water-Glycol (50/50) | −30 – 105 | 1075 | 3350 | 3.5×10⁻³ | 0.400 |
| Syltherm 800 | −40 – 400 | 935 | 1630 | 9.6×10⁻⁴ | 0.135 |
| Syltherm HF | −73 – 260 | 873 | 1590 | 5.0×10⁻⁴ | 0.100 |
| Dowtherm A | 15 – 400 | 1056 | 1590 | 2.5×10⁻³ | 0.138 |
| Dowtherm Q | −35 – 330 | 993 | 1660 | 2.2×10⁻³ | 0.114 |
| Steam (condensing) | 100 – 250 | — | — | — | — |
| Therminol 66 | −3 – 345 | 1005 | 1680 | 3.0×10⁻³ | 0.118 |
| Marlotherm SH | −5 – 350 | 1040 | 1630 | 2.8×10⁻³ | 0.120 |
| Brine (CaCl₂ 25%) | −40 – 60 | 1230 | 2810 | 4.5×10⁻³ | 0.540 |

Properties are at approximately 25 °C (nominal). All values are user-overridable.

**Sources:** Dow Chemical product data sheets; Perry & Green (2019), Section 11.
""")

    # ── Jacket area estimation
    st.header("Wetted Heat-Transfer Area")
    st.latex(r"A_{\text{total}} = A_{\text{dish}} + \pi \, D_T \, (H - h_{\text{dish}})")
    st.markdown("""
The wetted heat-transfer area is computed from the liquid height and vessel
geometry. When liquid partially fills the bottom dish, only the wetted
fraction of the dish area is counted.

| Bottom Dish Type | $h_{\\text{dish}}$ | $A_{\\text{dish,full}}$ |
|-----------------|-------------------|----------------------|
| Flat | 0 | $\\pi/4 \\cdot D_T^2$ |
| 2:1 Elliptical | $D_T / 4$ | $1.09 \\cdot \\pi/4 \\cdot D_T^2$ |
| Torispherical / DIN | $0.194 \\cdot D_T$ | $1.06 \\cdot \\pi/4 \\cdot D_T^2$ |
| Conical (45°) | $D_T / 2$ | $1.20 \\cdot \\pi/4 \\cdot D_T^2$ |

**Reference:** DIN 28131 (Jacketed agitated vessels).
""")

    # ── References
    st.header("Key References")
    st.markdown("""
1. Incropera, F.P., DeWitt, D.P., Bergman, T.L. & Lavine, A.S. (2007). [*Fundamentals of Heat and Mass Transfer*](https://openlibrary.org/isbn/0471457280), 6th ed. Wiley.
2. Perry, R.H. & Green, D.W. (2019). [*Perry's Chemical Engineers' Handbook*](https://openlibrary.org/isbn/0071834087), 9th ed. McGraw-Hill.
3. DIN 28131:1979 — Agitated vessels; heat transfer at the wall of agitated vessels.
4. Chilton, T.H., Drew, T.B. & Jebens, R.H. (1944). Heat transfer coefficients in agitated vessels. *Ind. Eng. Chem.* 36(6):510–516.
5. Dittus, F.W. & Boelter, L.M.K. (1930). Heat transfer in automobile radiators of the tubular type. *Univ. Calif. Publ. Eng.* 2:443–461.
6. Lehrer, I.H. (1970). Jacket-side Nusselt number. *Chem. Eng. Sci.* 25:1397.
7. Stein, W.A. & Schmidt, M. (1993). Heat transfer for agitated vessels with large impellers. *Chem. Eng. Process.* 32:305.
8. Nagata, S. (1975). *Mixing: Principles and Applications.* Wiley.
""")
