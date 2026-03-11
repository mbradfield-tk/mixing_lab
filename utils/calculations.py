"""
Hydrodynamic and mixing calculations for stirred-tank and continuous-flow reactors.

All equations referenced in the Equations Summary page of the Streamlit app.
"""

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Fundamental hydrodynamic parameters
# ---------------------------------------------------------------------------

def reynolds_number(N: float, D: float, rho: float, mu: float) -> float:
    """Impeller Reynolds number  Re = ρ N D² / μ"""
    return rho * N * D**2 / mu


def power_number_correlation(Re: float, Np_turb: float = 5.0) -> float:
    """
    Simplified Power-number model (turbulent plateau).
    Np ≈ Np_turb for Re > ~10 000; laminar correction for low Re.
    """
    if Re < 10:
        return 70 / Re          # laminar
    elif Re < 10000:
        return Np_turb * (Re / 10000)**0.18  # transitional (approx)
    return Np_turb              # turbulent


def impeller_power(Np: float, rho: float, N: float, D: float) -> float:
    """P = Np ρ N³ D⁵"""
    return Np * rho * N**3 * D**5


def power_per_volume(P: float, V: float) -> float:
    """ε = P / V  (W m⁻³)"""
    return P / V


def tip_speed(N: float, D: float) -> float:
    """u_tip = π N D"""
    return np.pi * N * D


def pumping_number_default() -> float:
    """Typical Nq for a pitched-blade turbine (down-pumping)."""
    return 0.79


def pumping_rate(Nq: float, N: float, D: float) -> float:
    """Q = Nq N D³"""
    return Nq * N * D**3


# ---------------------------------------------------------------------------
# Mixing times
# ---------------------------------------------------------------------------

def blend_time_turbulent(Nq: float, V: float, D: float, N: float) -> float:
    """
    Macro-blend (95 %) time using the circulation-model approach.
    θ_95 ≈ 5.2 V / (Nq N D³)   (Grenville correlation for turbulent flow)
    """
    Q = pumping_rate(Nq, N, D)
    if Q == 0:
        return np.inf
    return 5.2 * V / Q


def micromixing_time_engulfment(epsilon: float, nu: float) -> float:
    """
    Engulfment micro-mixing time (Baldyga & Bourne).
    t_E = 17.3 (ν / ε)^0.5
    """
    if epsilon <= 0:
        return np.inf
    return 17.3 * np.sqrt(nu / epsilon)


def kolmogorov_length(nu: float, epsilon: float) -> float:
    """η = (ν³ / ε)^(1/4)"""
    if epsilon <= 0:
        return np.inf
    return (nu**3 / epsilon)**0.25


def batchelor_length(nu: float, epsilon: float, D_mol: float) -> float:
    """λ_B = η · Sc^{-1/2},  Sc = ν / D_mol"""
    eta = kolmogorov_length(nu, epsilon)
    Sc = nu / D_mol if D_mol > 0 else 1e12
    return eta / np.sqrt(Sc)


# ---------------------------------------------------------------------------
# Turbulent energy-dissipation distribution
# ---------------------------------------------------------------------------

def epsilon_max_estimate(P: float, rho: float, D: float, N: float) -> float:
    """
    Local maximum dissipation rate near the impeller (order-of-magnitude).
    ε_max ≈ 0.5 Np N³ D² (Kresta & Wood)
    Uses simplified form: ε_max ≈ C · P / (ρ D³), C~1-5
    """
    if D == 0:
        return 0.0
    return 3.0 * P / (rho * D**3)


# ---------------------------------------------------------------------------
# kLa – volumetric mass-transfer coefficient (gas-liquid)
# ---------------------------------------------------------------------------

def kla_vant_riet(P_V: float, v_s: float, coalescing: bool = True) -> float:
    """
    Van 't Riet (1979) correlation for kLa in aerated stirred tanks.

    kLa = C₁ · (P/V)^C₂ · vₛ^C₃

    Coalescing (pure liquids):     C₁=0.026, C₂=0.4, C₃=0.5
    Non-coalescing (electrolytes): C₁=0.002, C₂=0.7, C₃=0.2

    Parameters
    ----------
    P_V : float   – gassed power per unit volume  (W/m³)
    v_s : float   – superficial gas velocity       (m/s)
    coalescing : bool – True for pure liquids, False for electrolytes

    Returns  kLa in 1/s.
    Reference: van 't Riet, K. (1979). Ind. Eng. Chem. Process Des. Dev., 18(3), 357-364.
    """
    if v_s <= 0 or P_V <= 0:
        return 0.0
    if coalescing:
        return 0.026 * P_V**0.4 * v_s**0.5
    return 0.002 * P_V**0.7 * v_s**0.2


def kla_surface(epsilon: float, nu: float, D_mol: float,
                D_tank: float, V: float) -> float:
    """
    Headspace-only (free-surface) volumetric mass-transfer coefficient.

    Uses the Lamont & Scott (1970) small-eddy surface-renewal model for
    the liquid-side mass-transfer coefficient:

        k_L = 0.4 · D_mol^{1/2} · (ε / ν)^{1/4}

    combined with the specific interfacial area of the flat free surface:

        a = A_surface / V = (π/4 · D_tank²) / V

    Parameters
    ----------
    epsilon : float  – mean energy dissipation rate  (W/m³ or W/kg, see note)
    nu      : float  – kinematic viscosity           (m²/s)
    D_mol   : float  – molecular diffusivity of dissolved gas  (m²/s)
    D_tank  : float  – tank diameter                 (m)
    V       : float  – liquid volume                 (m³)

    Returns  kLa_surface in 1/s.

    Note: ε should be in W/kg (= m²/s³) for dimensional consistency.
    The function converts internally if needed.

    Reference: Lamont, J. C. and Scott, D. S. (1970). An eddy cell model
    of mass transfer into the surface of a turbulent liquid.
    AIChE J., 16(4), 513–519.
    """
    if epsilon <= 0 or nu <= 0 or D_mol <= 0 or D_tank <= 0 or V <= 0:
        return 0.0
    kL = 0.4 * np.sqrt(D_mol) * (epsilon / nu) ** 0.25   # m/s
    A_surface = np.pi / 4 * D_tank**2                     # m²
    a = A_surface / V                                      # 1/m
    return kL * a                                          # 1/s


# ---------------------------------------------------------------------------
# Local (impeller-zone) micromixing time
# ---------------------------------------------------------------------------

def micromixing_time_local(epsilon_max: float, nu: float) -> float:
    """
    Engulfment micro-mixing time evaluated at the local maximum
    energy-dissipation rate near the impeller.

    t_E_local = 17.3 · (ν / ε_max)^0.5

    This gives a shorter (optimistic) micromixing time representative of
    the impeller discharge zone where reagents are often fed.
    """
    if epsilon_max <= 0:
        return np.inf
    return 17.3 * np.sqrt(nu / epsilon_max)


# ---------------------------------------------------------------------------
# Shear rate and shear stress
# ---------------------------------------------------------------------------

def average_shear_rate(P: float, mu: float, V: float) -> float:
    """
    Root-mean-square (Camp-Stein) average shear rate.

    γ̇_avg = √(P / (μ · V))

    Derived from equating viscous dissipation  μ·γ̇² = ε  (energy balance).
    Reference: Camp, T.R. and Stein, P.C. (1943). J. Boston Soc. Civil Eng., 30, 219-237.
    """
    if mu <= 0 or V <= 0:
        return 0.0
    return np.sqrt(P / (mu * V))


def maximum_shear_rate(epsilon_max: float, nu: float) -> float:
    """
    Maximum shear rate near the impeller.

    γ̇_max = √(ε_max / ν)
    """
    if nu <= 0 or epsilon_max <= 0:
        return 0.0
    return np.sqrt(epsilon_max / nu)


def shear_stress(mu: float, gamma_dot: float) -> float:
    """Newtonian shear stress  τ = μ · γ̇  (Pa)."""
    return mu * gamma_dot


# ---------------------------------------------------------------------------
# Solid-particle hydrodynamics
# ---------------------------------------------------------------------------

def settling_velocity(d_p: float, rho_p: float, rho_L: float,
                      mu: float, phi: float = 1.0) -> float:
    """Terminal settling velocity of a single particle.

    Uses the Stokes regime as the base estimate and then applies the
    Schiller-Naumann drag correction for intermediate Re_p (valid up
    to Re_p ≈ 1000).

    Parameters
    ----------
    d_p   : float  – particle diameter (m)
    rho_p : float  – particle density (kg/m³)
    rho_L : float  – liquid density (kg/m³)
    mu    : float  – dynamic viscosity of liquid (Pa·s)
    phi   : float  – sphericity / shape factor (0–1, 1 = sphere)

    Returns
    -------
    v_t : float – terminal settling velocity (m/s)
    """
    if d_p <= 0 or mu <= 0 or rho_L <= 0 or phi <= 0:
        return 0.0
    g = 9.81
    delta_rho = abs(rho_p - rho_L)
    # Stokes velocity
    v_stokes = (d_p**2 * g * delta_rho) / (18 * mu)
    # Shape factor correction (Haider & Levenspiel)
    v_stokes *= phi
    # Iterative Schiller-Naumann correction
    v_t = v_stokes
    for _ in range(20):
        Re_p = rho_L * v_t * d_p / mu
        if Re_p < 0.1:
            break
        Cd_corr = (24 / Re_p) * (1 + 0.15 * Re_p**0.687)
        v_t_new = np.sqrt(4 * g * d_p * delta_rho / (3 * Cd_corr * rho_L))
        v_t_new *= phi  # shape correction
        if abs(v_t_new - v_t) / max(v_t, 1e-30) < 1e-6:
            v_t = v_t_new
            break
        v_t = v_t_new
    return v_t


def particle_reynolds(d_p: float, v_t: float, rho_L: float,
                      mu: float) -> float:
    """Particle Reynolds number  Re_p = ρ_L · v_t · d_p / μ."""
    if mu <= 0:
        return 0.0
    return rho_L * v_t * d_p / mu


def zwietering_njs(S: float, nu: float, d_p: float, delta_rho: float,
                   rho_L: float, X: float, D_imp: float,
                   g: float = 9.81) -> float:
    """Zwietering (1958) just-suspended speed.

    N_js = S · ν^0.1 · d_p^0.2 · (g Δρ/ρ_L)^0.45 · X^0.13 · D^{-0.85}

    Parameters
    ----------
    S       : float  – Zwietering constant (geometry-dependent, typ. 1–10)
    nu      : float  – kinematic viscosity (m²/s)
    d_p     : float  – mass-mean particle diameter (m)
    delta_rho : float – |ρ_p − ρ_L| (kg/m³)
    rho_L   : float  – liquid density (kg/m³)
    X       : float  – solids loading (wt-% solids)
    D_imp   : float  – impeller diameter (m)
    g       : float  – gravitational acceleration (m/s²)

    Returns
    -------
    N_js : float – just-suspended impeller speed (rev/s)
    """
    if D_imp <= 0 or rho_L <= 0 or X <= 0 or d_p <= 0:
        return 0.0
    return (S * nu**0.1 * d_p**0.2
            * (g * delta_rho / rho_L)**0.45
            * X**0.13 * D_imp**(-0.85))


def solid_liquid_mass_transfer(d_p: float, v_slip: float, rho_L: float,
                               mu: float, D_mol: float) -> float:
    """Solid-liquid mass transfer coefficient via Ranz-Marshall.

    Sh = 2 + 0.6 · Re_p^{0.5} · Sc^{1/3}
    k_SL = Sh · D_mol / d_p

    Parameters
    ----------
    d_p    : float  – particle diameter (m)
    v_slip : float  – slip velocity (≈ settling velocity) (m/s)
    rho_L  : float  – liquid density (kg/m³)
    mu     : float  – dynamic viscosity (Pa·s)
    D_mol  : float  – molecular diffusivity of solute (m²/s)

    Returns
    -------
    k_SL : float – mass transfer coefficient (m/s)
    """
    if d_p <= 0 or mu <= 0 or D_mol <= 0:
        return 0.0
    Re_p = rho_L * v_slip * d_p / mu
    nu_val = mu / rho_L
    Sc = nu_val / D_mol if D_mol > 0 else 1e12
    Sh = 2.0 + 0.6 * np.sqrt(max(Re_p, 0)) * Sc**(1.0/3.0)
    return Sh * D_mol / d_p


def particle_suspension_criterion(N: float, N_js: float) -> str:
    """Qualitative suspension assessment based on N / N_js ratio."""
    if N_js <= 0:
        return "N/A"
    ratio = N / N_js
    if ratio < 0.7:
        return f"Poorly suspended (N/Njs={ratio:.2f})"
    elif ratio < 1.0:
        return f"Partially suspended (N/Njs={ratio:.2f})"
    elif ratio < 1.3:
        return f"Just suspended (N/Njs={ratio:.2f})"
    else:
        return f"Fully suspended (N/Njs={ratio:.2f})"


# ---------------------------------------------------------------------------
# Damköhler numbers
# ---------------------------------------------------------------------------

def damkohler_macro(t_blend: float, t_rxn: float) -> float:
    """Da_macro = θ_blend / t_rxn"""
    if t_rxn == 0:
        return np.inf
    return t_blend / t_rxn


def damkohler_micro(t_micro: float, t_rxn: float) -> float:
    """Da_micro = t_E / t_rxn"""
    if t_rxn == 0:
        return np.inf
    return t_micro / t_rxn


def damkohler_gl(kLa: float, t_rxn: float) -> float:
    """Gas-liquid Damköhler number  Da_GL = 1 / (kLa · t_rxn).

    Compares the characteristic mass-transfer time (1/kLa) to the
    reaction time.  Interpretation mirrors Da_macro / Da_micro:

      Da_GL << 1  →  mass transfer is fast, reaction-limited
      Da_GL ~ 1   →  comparable time scales
      Da_GL >> 1  →  mass-transfer-limited

    Returns 0 when kLa is zero (no gas-liquid transfer present).
    """
    if kLa <= 0:
        return 0.0
    if t_rxn <= 0:
        return np.inf
    return 1.0 / (kLa * t_rxn)


def mixing_sensitivity_assessment(Da_macro: float, Da_micro: float,
                                   Da_GL: float = 0.0) -> str:
    """
    Qualitative assessment based on Damköhler numbers.

    Da << 1  →  reaction-limited (insensitive to mixing)
    Da ~ 1   →  comparable time scales (potentially sensitive)
    Da >> 1  →  mixing-limited (sensitive)
    """
    labels = []
    for name, Da in [("Macro", Da_macro), ("Micro", Da_micro)]:
        if Da < 0.01:
            labels.append(f"{name}: Reaction-limited (Da={Da:.3g})")
        elif Da < 0.1:
            labels.append(f"{name}: Likely insensitive (Da={Da:.3g})")
        elif Da < 1:
            labels.append(f"{name}: Potentially sensitive (Da={Da:.3g})")
        elif Da < 10:
            labels.append(f"{name}: Mixing-sensitive (Da={Da:.3g})")
        else:
            labels.append(f"{name}: Strongly mixing-limited (Da={Da:.3g})")
    # G-L assessment (only when kLa is present)
    if Da_GL > 0:
        if Da_GL < 0.01:
            labels.append(f"G-L: Transfer-fast (Da_GL={Da_GL:.3g})")
        elif Da_GL < 0.1:
            labels.append(f"G-L: Likely insensitive (Da_GL={Da_GL:.3g})")
        elif Da_GL < 1:
            labels.append(f"G-L: Potentially transfer-limited (Da_GL={Da_GL:.3g})")
        elif Da_GL < 10:
            labels.append(f"G-L: Transfer-limited (Da_GL={Da_GL:.3g})")
        else:
            labels.append(f"G-L: Strongly transfer-limited (Da_GL={Da_GL:.3g})")
    return " | ".join(labels)


# ---------------------------------------------------------------------------
# Reaction time helpers
# ---------------------------------------------------------------------------

def half_life_first_order(k: float) -> float:
    """t_1/2 = ln2 / k"""
    if k <= 0:
        return np.inf
    return np.log(2) / k


def reaction_time_second_order(k: float, C0: float) -> float:
    """Characteristic time = 1 / (k C0)"""
    if k * C0 <= 0:
        return np.inf
    return 1.0 / (k * C0)


# ---------------------------------------------------------------------------
# Scale-up helpers
# ---------------------------------------------------------------------------

def scale_up_constant_tip_speed(N_small, D_small, D_large):
    """N_large such that tip speed is preserved."""
    return N_small * D_small / D_large


def scale_up_constant_P_V(N_small, D_small, D_large, Np_small=5.0, Np_large=5.0):
    """N_large such that P/V is preserved (geometric similarity assumed)."""
    # P/V ∝ Np N³ D⁵ / D³  ∝ Np N³ D²
    ratio = (Np_small / Np_large) * (D_small / D_large)**2
    return N_small * ratio**(1.0/3.0)


def scale_up_constant_Re(N_small, D_small, D_large):
    """N_large for constant Re (same fluid)."""
    return N_small * (D_small / D_large)**2


# ---------------------------------------------------------------------------
# Convenience: compute full parameter set for a reactor
# ---------------------------------------------------------------------------

def compute_reactor_hydro(
    N: float,          # rev/s
    D_imp: float,      # impeller diameter, m
    D_tank: float,     # tank diameter, m
    H: float,          # liquid height, m
    rho: float,        # kg/m³
    mu: float,         # Pa·s
    Np: float = None,  # power number (if known)
    Nq: float = None,  # pumping number (if known)
    v_s: float = 0.0,  # superficial gas velocity, m/s (0 = ungassed)
    coalescing: bool = True,  # True = pure liquid, False = electrolyte
    D_mol: float = 2.3e-9,   # molecular diffusivity of dissolved gas (m²/s)
) -> dict:
    """Return a dictionary of all computed hydrodynamic parameters."""
    V = np.pi / 4 * D_tank**2 * H
    nu = mu / rho
    Re = reynolds_number(N, D_imp, rho, mu)
    # Treat NaN the same as None (missing)
    if Np is None or (isinstance(Np, float) and np.isnan(Np)):
        Np = power_number_correlation(Re)
    if Nq is None or (isinstance(Nq, float) and np.isnan(Nq)):
        Nq = pumping_number_default()
    P = impeller_power(Np, rho, N, D_imp)
    eps = power_per_volume(P, V)            # W/m³
    u_tip = tip_speed(N, D_imp)
    Q = pumping_rate(Nq, N, D_imp)
    t_blend = blend_time_turbulent(Nq, V, D_imp, N)
    t_micro = micromixing_time_engulfment(eps, nu)
    eta = kolmogorov_length(nu, eps)
    eps_max = epsilon_max_estimate(P, rho, D_imp, N)  # W/m³
    eps_max_kg = eps_max / rho                        # W/kg

    # New parameters
    t_micro_local = micromixing_time_local(eps_max, nu)
    gamma_avg = average_shear_rate(P, mu, V)
    gamma_max = maximum_shear_rate(eps_max_kg, nu)
    tau_avg = shear_stress(mu, gamma_avg)
    kla = kla_vant_riet(eps, v_s, coalescing=coalescing)
    eps_kg = eps / rho  # W/kg for surface kLa (dimensional consistency)
    kla_surf = kla_surface(eps_kg, nu, D_mol, D_tank, V)

    return {
        "Volume (L)": V * 1000,
        "Re": Re,
        "Np": Np,
        "Power (W)": P,
        "P/V (W/m³)": eps,
        "P/V (W/L)": eps / 1000,
        "Tip speed (m/s)": u_tip,
        "Pumping rate (m³/s)": Q,
        "Blend time 95% (s)": t_blend,
        "Micromix time t_E (s)": t_micro,
        "Micromix time t_E_local (s)": t_micro_local,
        "Kolmogorov η (µm)": eta * 1e6,
        "ε_max (W/kg)": eps_max_kg,
        "Avg shear rate (1/s)": gamma_avg,
        "Max shear rate (1/s)": gamma_max,
        "Avg shear stress (Pa)": tau_avg,
        "kLa (1/s)": kla,
        "kLa_surface (1/s)": kla_surf,
        "ν (m²/s)": nu,
    }


def compute_damkohler_numbers(t_blend, t_micro, t_rxn,
                               kLa=0.0, kLa_surface=0.0):
    """Return Damköhler numbers (macro, micro, G-L) and assessment string."""
    Da_macro = damkohler_macro(t_blend, t_rxn)
    Da_micro = damkohler_micro(t_micro, t_rxn)
    # Use the larger of sparged or surface kLa for Da_GL
    kLa_eff = max(kLa, kLa_surface)
    Da_gl = damkohler_gl(kLa_eff, t_rxn)
    assessment = mixing_sensitivity_assessment(Da_macro, Da_micro, Da_gl)
    return {
        "Da_macro": Da_macro,
        "Da_micro": Da_micro,
        "Da_GL": Da_gl,
        "Assessment": assessment,
    }


# ---------------------------------------------------------------------------
# Heat generation & heat transfer
# ---------------------------------------------------------------------------

def reaction_rate_mol_per_s(order: str, k: float, C0: float,
                            V_L: float) -> float:
    """Instantaneous molar reaction rate (mol/s) at initial concentration.

    Parameters
    ----------
    order : str   – kinetic order ('1', 'pseudo-1', '2', 'pseudo-2')
    k     : float – rate constant (units consistent with order)
    C0    : float – initial concentration (mol/L)
    V_L   : float – liquid volume (litres)

    Returns
    -------
    r_total : float – total molar rate of consumption (mol/s)
    """
    if k <= 0 or C0 <= 0 or V_L <= 0:
        return 0.0
    if order in ("1", "pseudo-1"):
        r = k * C0  # mol/(L·s)
    elif order in ("2", "pseudo-2"):
        r = k * C0**2  # mol/(L·s)
    else:
        return 0.0
    return r * V_L  # mol/s  (r is per litre)


def heat_generation_rate(delta_H_kJ_mol: float, r_mol_per_s: float) -> float:
    """Rate of heat generation (W) from reaction.

    Q_rxn = |ΔH_rxn| × r

    Parameters
    ----------
    delta_H_kJ_mol : float – enthalpy of reaction (kJ/mol), negative = exothermic
    r_mol_per_s    : float – molar reaction rate (mol/s)

    Returns  Q_gen in Watts (positive = heat released).
    """
    return abs(delta_H_kJ_mol) * 1000.0 * r_mol_per_s  # kJ→J


def estimate_jacket_area(D_tank: float, H: float,
                         bottom_dish: str = "") -> float:
    """Estimate jacketed heat-transfer area (m²) from vessel geometry.

    Includes cylindrical wall + bottom dish.  The bottom dish area is
    approximated as:
      - 2:1 Elliptical:  A_dish ≈ 1.09 × π/4 × D²
      - Torispherical:   A_dish ≈ 1.06 × π/4 × D²
      - Conical:         A_dish ≈ 1.20 × π/4 × D²  (approx 60° cone)
      - Default:         A_dish ≈ π/4 × D²  (flat bottom)
    """
    if D_tank <= 0 or H <= 0:
        return 0.0
    A_cyl = np.pi * D_tank * H  # cylindrical wall
    A_flat = np.pi / 4 * D_tank**2
    dish = str(bottom_dish).lower() if bottom_dish else ""
    if "ellip" in dish:
        A_bottom = 1.09 * A_flat
    elif "torisph" in dish or "din" in dish:
        A_bottom = 1.06 * A_flat
    elif "conic" in dish:
        A_bottom = 1.20 * A_flat
    else:
        A_bottom = A_flat
    return A_cyl + A_bottom


def estimate_U(material: str = "", N_rps: float = 0.0) -> float:
    """Estimate overall heat-transfer coefficient U (W/m²·K).

    Typical ranges for jacketed stirred tanks:
      - Glass-lined:  100 – 250 W/m²·K
      - Stainless steel / Hastelloy: 200 – 500 W/m²·K
      - Carbon steel: 150 – 350 W/m²·K

    Returns a mid-range estimate based on wall material.  When agitation
    speed is provided, a simple correction nudges U toward the higher end
    of the range (better internal h_i at higher RPM).
    """
    mat = str(material).lower() if material else ""
    if "glass" in mat:
        U_lo, U_hi = 100.0, 250.0
    elif "hastel" in mat:
        U_lo, U_hi = 200.0, 450.0
    elif "carbon" in mat:
        U_lo, U_hi = 150.0, 350.0
    else:
        # Default to stainless steel range
        U_lo, U_hi = 200.0, 500.0
    # Simple agitation boost: interpolate from lo→hi as RPM ramps up
    # Assume "typical" max is ~3 rps; above that stays at U_hi
    frac = min(N_rps / 3.0, 1.0) if N_rps > 0 else 0.5
    return U_lo + frac * (U_hi - U_lo)


def heat_removal_capacity(U: float, A: float, dT: float) -> float:
    """Maximum steady-state heat-removal rate (W) through the jacket.

    Q_cool = U × A × ΔT

    Parameters
    ----------
    U  : float – overall heat-transfer coefficient (W/m²·K)
    A  : float – heat-transfer area (m²)
    dT : float – temperature driving force T_rxn − T_coolant (K or °C)

    Returns  Q_cool in Watts.
    """
    if U <= 0 or A <= 0 or dT <= 0:
        return 0.0
    return U * A * dT


def heat_balance_assessment(Q_gen: float, Q_cool: float) -> str:
    """Qualitative assessment comparing heat generation to cooling capacity."""
    if Q_gen <= 0:
        return "No heat generation"
    if Q_cool <= 0:
        return "⚠️ No cooling capacity estimated"
    ratio = Q_gen / Q_cool
    if ratio < 0.25:
        return f"Easily manageable (Q_gen/Q_cool = {ratio:.2f})"
    elif ratio < 0.5:
        return f"Comfortable margin (Q_gen/Q_cool = {ratio:.2f})"
    elif ratio < 0.75:
        return f"Moderate – monitor closely (Q_gen/Q_cool = {ratio:.2f})"
    elif ratio < 1.0:
        return f"⚠️ Tight – limited safety margin (Q_gen/Q_cool = {ratio:.2f})"
    else:
        return f"🔴 Insufficient cooling (Q_gen/Q_cool = {ratio:.2f})"
