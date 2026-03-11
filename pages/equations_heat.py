"""Equations – Heat Balance."""

import streamlit as st

st.title("📐 Heat Balance Equations")

# ── Heat generation ──────────────────────────────────────────────────────
st.header("H1. Reaction Rate (Molar)")

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
""")

st.header("H2. Heat Generation Rate")
st.latex(r"Q_{\text{gen}} = \lvert \Delta H_{rxn} \rvert \times r_{\text{total}}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\Delta H_{rxn}$ | Enthalpy of reaction (negative = exothermic) | kJ/mol |
| $r_{\\text{total}}$ | Total molar reaction rate | mol/s |
| $Q_{\\text{gen}}$ | Rate of heat release | W |
""")

# ── Heat removal ─────────────────────────────────────────────────────────
st.header("H3. Heat Removal Capacity (Jacket)")
st.latex(r"Q_{\text{cool}} = U \, A \, \Delta T")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $U$ | Overall heat-transfer coefficient | W/(m²·K) |
| $A$ | Heat-transfer area (jacket + bottom) | m² |
| $\\Delta T$ | $T_{\\text{process}} - T_{\\text{coolant}}$ | °C or K |
| $Q_{\\text{cool}}$ | Maximum steady-state heat removal | W |
""")

# ── Jacket area estimation ───────────────────────────────────────────────
st.header("H4. Jacket Heat-Transfer Area Estimation")
st.latex(r"A_{\text{cyl}} = \pi \, D_T \, H")
st.latex(r"A_{\text{total}} = A_{\text{cyl}} + f_{dish} \cdot \frac{\pi}{4} D_T^2")
st.markdown("""
The cylindrical wall area plus a bottom-dish correction factor:

| Bottom Dish Type | $f_{dish}$ |
|-----------------|------------|
| Flat | 1.00 |
| 2:1 Elliptical | 1.09 |
| Torispherical / DIN | 1.06 |
| Conical (~60°) | 1.20 |
""")

# ── U estimation ─────────────────────────────────────────────────────────
st.header("H5. Overall Heat-Transfer Coefficient Estimation")

st.subheader("H5a. Simple Material-Based Estimate (Fallback)")
st.markdown("""
Typical ranges for jacketed stirred tanks (mid-range values used when
reactor-specific $U$ is not provided):

| Wall Material | $U$ range (W/m²·K) |
|--------------|-------------------|
| Glass-lined | 100 – 250 |
| Stainless steel | 200 – 500 |
| Hastelloy | 200 – 450 |
| Carbon steel | 150 – 350 |

A simple agitation correction interpolates within the range: higher impeller
speed increases the internal heat-transfer coefficient $h_i$, which raises $U$.
""")

st.subheader("H5b. Individual Resistances Approach")
st.latex(r"\frac{1}{U} = \frac{1}{h_i} + \frac{x_w}{k_w} + \frac{1}{h_o} + R_f")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $h_i$ | Process-side (internal) heat-transfer coefficient | W/(m²·K) |
| $x_w$ | Wall thickness | m |
| $k_w$ | Wall thermal conductivity | W/(m·K) |
| $h_o$ | Utility-side (jacket) heat-transfer coefficient | W/(m²·K) |
| $R_f$ | Fouling resistance | m²·K/W |
""")

st.subheader("H5c. Process-Side Coefficient (Nusselt Correlation)")
st.markdown("""
For turbulent flow in a baffled, jacketed stirred vessel (Seider-Tate / DIN 28131 form):
""")
st.latex(r"\mathrm{Nu} = C \cdot \mathrm{Re}^{2/3} \cdot \mathrm{Pr}^{1/3} \cdot \left(\frac{\mu}{\mu_w}\right)^{0.14}")
st.latex(r"h_i = \frac{\mathrm{Nu} \cdot k_f}{D_T}")
st.markdown("""
where:
""")
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
""")

st.subheader("H5d. Wall Thermal Conductivity")
st.markdown("""
| Material | $k_w$ (W/(m·K)) |
|----------|-----------------|
| Stainless steel (316/304) | 16 |
| Hastelloy C-276 | 12 |
| Carbon steel | 50 |
| Glass lining | 1.0 |
| Titanium | 22 |
| Inconel | 15 |

For **glass-lined** vessels, the total wall resistance includes both the
glass lining (~1.5 mm, $k$ = 1.0 W/(m·K)) and the steel shell:

$$R_{wall} = \\frac{x_{glass}}{k_{glass}} + \\frac{x_{steel}}{k_{steel}}$$
""")

st.subheader("H5e. Jacket-Side Coefficient")
st.markdown("""
When jacket geometry and coolant flow data are not available, typical
values are used:

| Jacket Type | Typical $h_o$ (W/(m²·K)) |
|-------------|--------------------------|
| Simple jacket (water/glycol) | 1500 |
| Half-pipe coil | 2500 |
| Dimple jacket | 1200 |

These assume turbulent flow of aqueous-based coolants in the jacket.
""")

st.subheader("H5f. Fouling Resistance")
st.markdown("""
| Condition | $R_f$ (m²·K/W) |
|-----------|----------------|
| Clean (fresh equipment) | 0.0001 |
| Moderate (pharmaceutical) | 0.0002 |
| Heavy fouling | 0.001 |

Default value: **0.0002 m²·K/W** (moderate pharmaceutical process).
""")

# ── Heat balance assessment ──────────────────────────────────────────────
st.header("H6. Heat Balance Assessment")
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
