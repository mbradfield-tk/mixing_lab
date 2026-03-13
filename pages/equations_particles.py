"""Equations – Solid-Particle Hydrodynamics."""

import streamlit as st

st.title("📐 Solid-Particle Equations")

# ── Terminal Settling Velocity ───────────────────────────────────────────
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

**Reference:** Clift, R., Grace, J.R. and Weber, M.E. (1978).
*Bubbles, Drops, and Particles.* Academic Press.
""")

st.header("Particle Reynolds Number")
st.latex(r"Re_p = \frac{\rho_L \, v_t \, d_p}{\mu}")
st.markdown("Classifies the flow regime around a settling particle.")

st.header("Zwietering Just-Suspended Speed")
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

st.header("Solid-Liquid Mass Transfer (Ranz-Marshall)")
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

st.header("Suspension Quality Assessment")
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
