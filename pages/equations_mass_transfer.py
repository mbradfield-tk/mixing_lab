"""Equations – Mass Transfer & kLa."""

import streamlit as st

st.title("📐 Mass Transfer Equations")

# ── kLa – sparged ────────────────────────────────────────────────────────
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

**Reference:** van 't Riet, K. (1979). Review of measuring methods and results in non-viscous
gas-liquid mass transfer in stirred vessels. *Ind. Eng. Chem. Process Des. Dev.*, 18(3), 357–364.
""")

# ── kLa – headspace / free-surface ───────────────────────────────────────
st.header("kLa – Free-Surface (Headspace) Mass Transfer")
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

# ── Gas Holdup ───────────────────────────────────────────────────────────
st.header("Gas Holdup (Hughmark)")
st.latex(r"\varepsilon_G = 0.505 \; v_s^{0.47} \left( \frac{P}{V} \right)^{0.4} \left( \frac{\mu}{\sigma} \right)^{0.08}")
st.markdown(r"""
Simplified Hughmark (1967) / Yawalkar correlation for gas void fraction
in a mechanically agitated vessel.

| Symbol | Description | Units |
|--------|-------------|-------|
| $v_s$ | Superficial gas velocity | m/s |
| $P/V$ | Gassed power per unit volume | W/m³ |
| $\mu$ | Liquid dynamic viscosity | Pa·s |
| $\sigma$ | Liquid surface tension | N/m |
| $\varepsilon_G$ | Gas holdup (volume fraction) | — |

**Reference:** Hughmark, G. A. (1967). *Ind. Eng. Chem. Process Des. Dev.*, 6(2), 218–220.
""")

# ── Sauter Mean Bubble Diameter ──────────────────────────────────────────
st.header("Sauter Mean Bubble Diameter (Calderbank)")
st.latex(r"d_{32} = 4.15 \; \frac{\sigma^{0.6}}{(P/V)^{0.4} \; \rho_L^{0.2}} + 0.0009")
st.markdown(r"""
Calderbank (1958) correlation for the Sauter mean bubble diameter in a
mechanically agitated aerated vessel.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\sigma$ | Surface tension | N/m |
| $P/V$ | Power per unit volume | W/m³ |
| $\rho_L$ | Liquid density | kg/m³ |
| $d_{32}$ | Sauter mean bubble diameter | m |

The first term captures breakup by turbulence; the constant 0.0009 m
represents the minimum stable bubble size.

**Reference:** Calderbank, P. H. (1958). Physical rate processes in industrial
fermentation, Part I: The interfacial area in gas-liquid contacting with
mechanical agitation. *Trans. Inst. Chem. Engrs.*, 36, 443–463.
""")

# ── Gas Flooding Speed ───────────────────────────────────────────────────
st.header("Gas Flooding Speed")
st.latex(r"Fl_G = \frac{Q_{gas}}{N \, D^3}")
st.latex(r"N_{flood} = \frac{Q_{gas}}{Fl_{crit} \, D^3}, \quad Fl_{crit} \approx 0.035")
st.markdown(r"""
The flooding number $Fl_G$ is the ratio of volumetric gas flow rate to
the impeller's volumetric throughput.  When $Fl_G$ exceeds a critical
value ($Fl_{crit} \approx 0.035$), the impeller can no longer
disperse the gas effectively — the impeller is said to be **flooded**.

| Symbol | Description | Units |
|--------|-------------|-------|
| $Q_{gas}$ | Volumetric gas flow rate | m³/s |
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |
| $Fl_{crit}$ | Critical flooding number | — |
| $N_{flood}$ | Minimum speed for gas dispersion | rev/s |

Operating at $N / N_{flood} > 1$ ensures the gas is dispersed.
Recommended: $N / N_{flood} > 1.3$ for reliable dispersion.

**Reference:** Nienow, A. W., Wisdom, D. J. and Middleton, J. C. (1977).
The effect of scale and geometry on flooding, recirculation and power in
gassed stirred vessels. *Proceedings 2nd European Conference on Mixing*, F1-1.
""")

# ══════════════════════════════════════════════════════════════════════════
# Liquid-Liquid Mass Transfer
# ══════════════════════════════════════════════════════════════════════════
st.divider()
st.title("📐 Liquid-Liquid Equations")

# ── Miscibility Screening ────────────────────────────────────────────────
st.header("Miscibility Screening")
st.markdown(r"""
The app uses a **two-tier** approach to assess whether two liquids are
immiscible (and therefore valid for L-L dispersion calculations):

**1. Known-pair lookup** — For built-in solvents, experimental miscibility
data from standard references (Perry's, CRC Handbook, Merck Index) is used
directly.  Pairs are classified as *miscible*, *partially miscible*, or
*immiscible*.

**2. Hansen Solubility Parameters (fallback)** — For custom fluids or
unknown pairs, the Hansen distance $R_a$ is computed:
""")
st.latex(r"R_a = \sqrt{4\,(\delta_{d1}-\delta_{d2})^2 + (\delta_{p1}-\delta_{p2})^2 + (\delta_{h1}-\delta_{h2})^2}")
st.markdown(r"""
Each solvent is described by three parameters: **dispersion** ($\delta_d$),
**polar** ($\delta_p$), and **hydrogen-bonding** ($\delta_h$), in MPa$^{0.5}$.
The factor of 4 on the dispersion term is empirical.

| $R_a$ (MPa$^{0.5}$) | Interpretation |
|----------------------|----------------|
| < 15 | Likely **miscible** |
| 15 – 25 | Partially miscible / borderline |
| > 25 | Likely **immiscible** |

> **Limitations:** The Hansen approach was originally developed for
> polymer–solvent systems.  It can give misleading results for aqueous
> systems because water's extreme $\delta_h \approx 42$ MPa$^{0.5}$ inflates
> the distance to almost all organics.  The known-pair lookup avoids this
> issue for built-in solvents.

**Reference:** Hansen, C. M. (2007). *Hansen Solubility Parameters:
A User's Handbook*, 2nd ed. CRC Press.
""")

# ── Weber Number ─────────────────────────────────────────────────────────
st.header("Impeller Weber Number (Liquid-Liquid)")
st.latex(r"We = \frac{\rho_c \, N^2 \, D^3}{\sigma_{LL}}")
st.markdown(r"""
The Weber number compares inertial (disrupting) forces to interfacial
tension (stabilising) forces for an immiscible liquid-liquid system.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\rho_c$ | Continuous-phase density | kg/m³ |
| $N$ | Impeller speed | rev/s |
| $D$ | Impeller diameter | m |
| $\sigma_{LL}$ | Interfacial tension between two liquids | N/m |

Higher We → finer drops; lower We → coarser drops or phase separation.
""")

# ── Sauter Mean Drop Diameter ────────────────────────────────────────────
st.header("Sauter Mean Drop Diameter (Hinze-Kolmogorov)")
st.latex(r"\frac{d_{32}}{D} = C_1 \, We^{-0.6} (1 + C_2 \, \phi_d)")
st.markdown(r"""
Hinze (1955) / Chen & Middleman (1967) correlation for the equilibrium
Sauter mean drop diameter in a turbulent stirred vessel.

| Symbol | Description | Typical value |
|--------|-------------|---------------|
| $C_1$ | Constant | 0.053 |
| $C_2$ | Holdup correction | 3.0 |
| $\phi_d$ | Dispersed-phase volume fraction | 0 – 0.5 |
| $D$ | Impeller diameter | m |

The factor $(1 + C_2 \phi_d)$ accounts for the increase in drop size due
to **coalescence** at higher dispersed-phase fractions.

**Reference:** Hinze, J. O. (1955). *AIChE J.*, 1(3), 289–295;
Chen, H. T. and Middleman, S. (1967). *AIChE J.*, 13(5), 989–995.
""")

# ── Minimum Dispersion Speed ─────────────────────────────────────────────
st.header("Minimum Dispersion Speed (Skelland & Seksaria)")
st.latex(r"N_{min} = C \left( \frac{\sigma_{LL}}{\rho_c \, D^3} \right)^{1/2} (1 + 2.5 \, \phi_d)")
st.markdown(r"""
Minimum impeller speed to maintain a liquid-liquid dispersion (Skelland &
Seksaria, 1978).  Below this speed, the dispersed phase coalesces and
separates.

| Symbol | Description | Units |
|--------|-------------|-------|
| $C$ | Geometry constant (≈ 1.03 for Rushton) | — |
| $\sigma_{LL}$ | Interfacial tension | N/m |
| $\rho_c$ | Continuous-phase density | kg/m³ |
| $D$ | Impeller diameter | m |
| $\phi_d$ | Dispersed-phase volume fraction | — |

**Reference:** Skelland, A. H. P. and Seksaria, R. (1978). Minimum impeller
speeds for liquid-liquid dispersion in baffled vessels. *Ind. Eng. Chem.
Process Des. Dev.*, 17(1), 56–61.
""")

# ── LL Mass Transfer ─────────────────────────────────────────────────────
st.header("Liquid-Liquid Mass-Transfer Coefficient")
st.latex(r"v_{slip} = (\varepsilon \, d_{32})^{1/3}")
st.latex(r"Re_d = \frac{\rho_c \, v_{slip} \, d_{32}}{\mu_c}")
st.latex(r"Sh = 2 + 0.6 \, Re_d^{1/2} \, Sc^{1/3}")
st.latex(r"k_{LL} = \frac{Sh \cdot D_{mol}}{d_{32}}")
st.markdown(r"""
Ranz-Marshall type correlation applied to liquid-liquid systems.
The slip velocity is estimated from the local turbulent energy dissipation
rate and the drop diameter.

| Symbol | Description | Units |
|--------|-------------|-------|
| $\varepsilon$ | Mass-specific energy dissipation rate | W/kg |
| $d_{32}$ | Sauter mean drop diameter | m |
| $v_{slip}$ | Turbulent slip velocity | m/s |
| $Re_d$ | Drop Reynolds number | — |
| $Sc = \mu_c / (\rho_c \, D_{mol})$ | Schmidt number | — |
| $D_{mol}$ | Molecular diffusivity of solute | m²/s |
| $k_{LL}$ | Liquid-liquid mass-transfer coefficient | m/s |

The **volumetric** coefficient is obtained by multiplying with the specific
interfacial area:

$$a_{LL} = \frac{6 \, \phi_d}{d_{32}}, \quad k_L a_{LL} = k_{LL} \cdot a_{LL}$$

**Reference:** Calderbank, P. H. and Moo-Young, M. B. (1961). The continuous
phase heat and mass-transfer properties of dispersions. *Chem. Eng. Sci.*,
16(1–2), 39–54.
""")

# ── Phase Separation ─────────────────────────────────────────────────────
st.header("Phase Separation Time Estimate")
st.latex(r"v_{drop} = \frac{\Delta\rho \, g \, d_{32}^2}{18 \, \mu_c}")
st.latex(r"t_{sep} \approx \frac{H}{v_{drop}}")
st.markdown(r"""
Stokes-regime settling velocity of the mean drop, combined with the
liquid height, provides a rough estimate of how quickly the dispersion
will separate once agitation stops.

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
""")
