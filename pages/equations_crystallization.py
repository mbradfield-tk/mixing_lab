"""Equations – Crystallization."""

import streamlit as st

st.title("📐 Crystallization Equations")

# ── Supersaturation ──────────────────────────────────────────────────────
st.header("Supersaturation")

st.subheader("Relative Supersaturation")
st.latex(r"\sigma = \frac{c - c_{sat}}{c_{sat}}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $c$ | Actual solute concentration | g/L (or mol/L) |
| $c_{sat}$ | Equilibrium saturation concentration at temperature $T$ | g/L (or mol/L) |
| $\\sigma$ | Relative supersaturation (dimensionless) | — |

The driving force for both nucleation and crystal growth.

**Reference:** Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). *Handbook of Industrial Crystallization*, 3rd ed. Cambridge University Press, Ch. 2.
""")

st.subheader("Supersaturation Ratio")
st.latex(r"S = \frac{c}{c_{sat}} = 1 + \sigma")
st.markdown("""
An alternative definition commonly used in nucleation theory.
$S = 1$ at equilibrium; $S > 1$ indicates a supersaturated solution.
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
""")

# ── Metastable Zone Width ────────────────────────────────────────────────
st.header("Metastable Zone Width (MSZW)")
st.latex(r"\text{MSZW} = T_{sat} - T_{nuc}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $T_{sat}$ | Saturation (equilibrium) temperature | °C |
| $T_{nuc}$ | Temperature at which nucleation is first detected | °C |
| MSZW | Metastable zone width | °C |

The MSZW defines the maximum undercooling (or supersaturation) the system
can sustain before spontaneous primary nucleation occurs.

- **Narrow MSZW** → nucleation is easy to trigger → higher mixing sensitivity.
- **Wide MSZW** → more forgiving; local $\\sigma$ spikes are less likely to exceed the metastable limit.

Measured by turbidimetry, FBRM, or visual observation at a controlled cooling rate.

**Reference:** Myerson *et al.* (2019), Ch. 5.
""")

# ── Primary Nucleation ───────────────────────────────────────────────────
st.header("Primary Nucleation Rate")

st.subheader("Classical Nucleation Theory (CNT)")
st.latex(r"J = A \exp\!\left( -\frac{16 \pi \gamma^3 v_m^2}{3 k_B^3 T^3 (\ln S)^2} \right)")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $J$ | Nucleation rate (number of nuclei per unit volume per time) | #/(m³·s) |
| $A$ | Pre-exponential kinetic factor | #/(m³·s) |
| $\\gamma$ | Crystal–solution interfacial energy | J/m² |
| $v_m$ | Molecular volume of the solute in the crystal | m³ |
| $k_B$ | Boltzmann constant ($1.381 \\times 10^{-23}$) | J/K |
| $T$ | Absolute temperature | K |
| $S$ | Supersaturation ratio $c/c_{sat}$ | — |

CNT predicts an extremely steep dependence of $J$ on $\\sigma$ — a
small increase in local supersaturation (e.g., from poor mixing in a
feed plume) can increase the nucleation rate by orders of magnitude.

**Reference:** Myerson *et al.* (2019), Ch. 3; Mersmann, A. (2001). *Crystallization Technology Handbook*, Ch. 2.
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
far more sensitive to local $\\sigma$ than growth**.  This is the
fundamental reason mixing affects crystal size distribution (CSD).

**Reference:** Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical Reactions.* Wiley.
""")

# ── Induction Time ───────────────────────────────────────────────────────
st.header("Induction Time")
st.latex(r"t_{ind} \propto \frac{1}{J} \propto \exp\!\left( +\frac{16 \pi \gamma^3 v_m^2}{3 k_B^3 T^3 (\ln S)^2} \right)")
st.markdown("""
The induction time $t_{ind}$ is the experimentally observed delay between
achieving supersaturation and detecting the first nuclei (by turbidity
change, FBRM counts, etc.).

- $t_{ind}$ is inversely related to the nucleation rate $J$.
- Shorter $t_{ind}$ → faster nucleation → more mixing-sensitive.
- Measured at a specific $\\sigma$ and temperature.

In the Crystallization Sensitivity Protocol, $t_{ind}$ serves as the
crystallization analogue of $t_{rxn}$ used in the Reaction Sensitivity
Protocol.

**Reference:** Myerson *et al.* (2019), Ch. 5; Kashchiev, D. & van Rosmalen, G.M. (2003). *Cryst. Res. Technol.* 38, 555–574.
""")

# ── Crystal Growth ───────────────────────────────────────────────────────
st.header("Crystal Growth Rate")

st.subheader("Power-Law Growth Model")
st.latex(r"G = k_g \, \sigma^{\,g}")
st.markdown("""
| Symbol | Description | Typical range |
|--------|-------------|---------------|
| $G$ | Linear growth rate (face velocity) | m/s |
| $k_g$ | Growth rate constant | m/s (system-specific) |
| $g$ | Growth order | 1–2 |

This is the rate at which a crystal face advances.  For size-independent
growth (McCabe's $\\Delta L$ law), $G$ is the same for all crystals at a
given $\\sigma$.

**Reference:** Myerson *et al.* (2019), Ch. 6.
""")

st.subheader("Growth Time (Characteristic)")
st.latex(r"t_G = \frac{L_{target}}{G} = \frac{L_{target}}{k_g \, \sigma^{\,g}}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $L_{target}$ | Target crystal size | m (or µm) |
| $t_G$ | Time to grow a crystal to the target size | s |

This characteristic growth time is compared against macromixing time
($\\theta_{95}$) via the macro Damköhler number.

**Reference:** Green, D.A. (2019). Ch. 10 in Myerson *et al.* (2019).
""")

# ── Damköhler Numbers for Crystallization ────────────────────────────────
st.header("Damköhler Numbers for Crystallization")
st.markdown("""
The Damköhler number framework from Baldyga & Bourne (1999) is adapted for
crystallization by replacing the reaction time $t_{rxn}$ with crystallization
time scales:
""")

st.subheader("Micro Damköhler Number (Nucleation)")
st.latex(r"Da_{micro} = \frac{t_E}{t_{ind}}")
st.markdown("""
Compares micromixing time ($t_E$) to the induction time.  When $Da_{micro} > 1$,
nucleation occurs within the unmixed feed plume at locally high $\\sigma$.

| Symbol | Description |
|--------|-------------|
| $t_E$ | Engulfment micromixing time $= 17.3 \\sqrt{\\nu / \\varepsilon}$ |
| $t_{ind}$ | Induction time at target $\\sigma$ |
""")

st.subheader("Macro Damköhler Number (Growth)")
st.latex(r"Da_{macro} = \frac{\theta_{95}}{t_G}")
st.markdown("""
Compares blend time ($\\theta_{95}$) to the growth time.  When $Da_{macro} > 1$,
supersaturation is non-uniform across the vessel on the timescale of crystal
growth, leading to growth rate dispersion.

| Symbol | Description |
|--------|-------------|
| $\\theta_{95}$ | 95% macro-blend time $= 5.2\\,V / (N_Q\\,N\\,D^3)$ |
| $t_G$ | Growth time to target crystal size |
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

**Reference:** Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical Reactions.* Wiley;
Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). *Handbook of Industrial Mixing.* Wiley.
""")

# ── Nucleation vs Growth Competition ─────────────────────────────────────
st.header("Nucleation–Growth Competition")
st.latex(r"B \propto \sigma^{\,n} \qquad (n \sim 2{-}10)")
st.latex(r"G \propto \sigma^{\,g} \qquad (g \sim 1{-}2)")
st.markdown("""
Because $n \\gg g$, the nucleation rate increases much faster with local
supersaturation than the growth rate.  When mixing is poor:

1. Local $\\sigma$ in the feed plume is much higher than the bulk.
2. Nucleation accelerates disproportionately → **excess fines**.
3. More nuclei share the available solute → **smaller mean crystal size**.
4. Broader CSD and potentially incorrect polymorph (Ostwald's rule of stages).

This is the crystallization analogue of selectivity loss in competing chemical
reactions.

**Reference:** Baldyga, J. & Bourne, J.R. (1999); Green, D.A. (2019), Ch. 10 in Myerson *et al.* (2019).
""")

# ── Ostwald's Rule of Stages ────────────────────────────────────────────
st.header("Ostwald's Rule of Stages")
st.markdown("""
Ostwald's rule states that the **least stable (metastable) polymorph**
nucleates first, because it has the lower interfacial energy $\\gamma$
and hence the lower nucleation barrier.

At high local $\\sigma$ (e.g., in a poorly mixed feed plume), the
metastable form is preferentially nucleated.  Whether it persists or
converts to the stable form depends on solvent-mediated transformation
kinetics.

**Mixing implication:** Poor mixing → local $\\sigma$ spikes → metastable
polymorph nucleation → potential quality / regulatory failure if the
wrong form persists.

**Reference:** Myerson *et al.* (2019), Ch. 3.
""")

# ── Crystallization Enthalpy ─────────────────────────────────────────────
st.header("Crystallization Enthalpy")
st.latex(r"Q_{cryst} = |\Delta H_{cryst}| \times \dot{n}_{cryst}")
st.markdown("""
| Symbol | Description | Units |
|--------|-------------|-------|
| $\\Delta H_{cryst}$ | Enthalpy of crystallization (usually exothermic, < 0) | kJ/mol |
| $\\dot{n}_{cryst}$ | Molar crystallization rate | mol/s |
| $Q_{cryst}$ | Heat generation from crystallization | W |

For most organic pharmaceutical compounds, $|\\Delta H_{cryst}|$ is in the range
of 10–50 kJ/mol.

| $|\\Delta H_{cryst}|$ | Risk level |
|-----------------------|-----------|
| < 20 kJ/mol | 🟢 Low — heat effects unlikely to limit process |
| 20–40 kJ/mol | 🟡 Moderate — check at scale |
| > 40 kJ/mol | 🔴 High — heat removal may limit cooling/feed rates |

**Reference:** Myerson *et al.* (2019), Ch. 5; Mersmann (2001).
""")

st.subheader("Van 't Hoff Estimation of ΔH_cryst")
st.latex(r"\frac{d \ln c_{sat}}{d(1/T)} = -\frac{\Delta H_{cryst}}{R}")
st.markdown("""
When $\\Delta H_{cryst}$ is not measured directly, it can be estimated from the
slope of $\\ln c_{sat}$ vs. $1/T$ (a van 't Hoff plot of solubility data).

| Symbol | Description |
|--------|-------------|
| $R$ | Gas constant (8.314 J/(mol·K)) |
| $c_{sat}$ | Saturation concentration at temperature $T$ |

**Reference:** Myerson *et al.* (2019), Ch. 2.
""")

# ── Seeding ──────────────────────────────────────────────────────────────
st.header("Seeding — Effect on Mixing Sensitivity")
st.markdown("""
Seed crystals provide pre-existing surface area $A_{seed}$ that consumes
supersaturation via secondary nucleation and growth, reducing the reliance
on primary nucleation.

$$\\frac{d\\sigma}{dt} \\propto -k_g \\, \\sigma^g \\, A_{seed}$$

With adequate seed loading:
- The system operates at **lower effective $\\sigma$** → less sensitive to mixing.
- Growth dominates over nucleation → narrower CSD.
- **Rule of thumb:** A well-seeded process is 2–5× less sensitive to mixing
  than an unseeded one at the same target $\\sigma$.

**Reference:** Black, S.N. (2019). Ch. 13 in Myerson *et al.* (2019).
""")

# ── Mixing-Sensitive Crystallization Types ───────────────────────────────
st.header("Mixing Sensitivity by Crystallization Type")
st.markdown("""
| Type | How $\\sigma$ is generated | Dominant mixing concern | Inherent risk |
|------|--------------------------|------------------------|---------------|
| **Cooling** | Temperature reduction | Macromixing — uniformity of $T$ field | 🟢 Low |
| **Anti-solvent** | Addition of a poor solvent | Meso- & micromixing — feed plume $\\sigma$ spike | 🔴 High |
| **Reactive** | Chemical reaction produces insoluble product | Micromixing — instantaneous local $\\sigma$ | 🔴 High |
| **Evaporative** | Solvent removal at boiling surface | Macromixing — boiling zone $\\sigma$ gradient | 🟡 Moderate |
| **pH-shift** | Acid/base addition changes solubility | Meso- & micromixing — local pH plume | 🔴 High |

Anti-solvent, reactive, and pH-shift crystallizations are inherently **feed-sensitive**:
supersaturation is generated at the feed point, creating locally very high $\\sigma$ in
the feed plume before mixing disperses it.

**Reference:** Green, D.A. (2019), Ch. 10 in Myerson *et al.* (2019).
""")

# ── References ───────────────────────────────────────────────────────────
st.header("References")
st.markdown("""
1. Myerson, A.S., Erdemir, D. & Lee, A.Y. (2019). *Handbook of Industrial
   Crystallization*, 3rd ed. Cambridge University Press.
2. Green, D.A. (2019). Ch. 10: Understanding and Modeling Crystallizer Mixing
   and Suspension Flow. In Myerson *et al.* (2019).
3. Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical
   Reactions.* Wiley.
4. Karpinski, P.H. & Baldyga, J. (2019). Ch. 8: Precipitation Processes.
   In Myerson *et al.* (2019).
5. Black, S.N. (2019). Ch. 13: Crystallization in the Pharmaceutical
   Industry. In Myerson *et al.* (2019).
6. Mersmann, A. (2001). *Crystallization Technology Handbook*, 2nd ed.
   Marcel Dekker.
7. Kashchiev, D. & van Rosmalen, G.M. (2003). Review: Nucleation in solutions
   revisited. *Cryst. Res. Technol.* 38(7–8), 555–574.
8. Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). *Handbook of
   Industrial Mixing.* Wiley.
""")
