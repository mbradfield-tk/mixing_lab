"""Mixing time calculations: blend time, micromixing, mesomixing, length scales."""

import numpy as np

from .hydrodynamics import pumping_rate

# Named constants for literature correlations
GRENVILLE_CONSTANT = 5.2        # Grenville (1992) blend-time coefficient
ENGULFMENT_CONSTANT = 17.3      # Baldyga & Bourne engulfment model constant
EPSILON_MAX_COEFF = 3.0         # Kresta & Wood ε_max / (P/ρD³) coefficient


def blend_time_turbulent(Nq: float, V: float, D: float, N: float) -> float:
    """
    Macro-blend (95 %) time using the circulation-model approach.
    θ_95 ≈ 5.2 V / (Nq N D³)   (Grenville correlation for turbulent flow)
    """
    Q = pumping_rate(Nq, N, D)
    if Q == 0:
        return np.inf
    return GRENVILLE_CONSTANT * V / Q


def micromixing_time_engulfment(epsilon: float, nu: float) -> float:
    """
    Engulfment micro-mixing time (Baldyga & Bourne).
    t_E = 17.3 (ν / ε)^0.5

    epsilon must be in W/kg (= m²/s³), NOT W/m³.
    """
    if epsilon <= 0:
        return np.inf
    return ENGULFMENT_CONSTANT * np.sqrt(nu / epsilon)


def micromixing_time_local(epsilon_max: float, nu: float) -> float:
    """
    Engulfment micro-mixing time evaluated at the local maximum
    energy-dissipation rate near the impeller.
    t_E_local = 17.3 · (ν / ε_max)^0.5
    """
    if epsilon_max <= 0:
        return np.inf
    return ENGULFMENT_CONSTANT * np.sqrt(nu / epsilon_max)


def mesomixing_time(epsilon: float, d_feed: float) -> float:
    """Mesomixing (turbulent dispersion) time constant.
    t_meso = 2 · (d_feed² / ε)^(1/3)
    """
    if epsilon <= 0 or d_feed <= 0:
        return np.inf
    return 2.0 * (d_feed**2 / epsilon) ** (1.0 / 3.0)


def kolmogorov_length(nu: float, epsilon: float) -> float:
    """η = (ν³ / ε)^(1/4).  epsilon must be in W/kg (= m²/s³)."""
    if epsilon <= 0:
        return np.inf
    return (nu**3 / epsilon)**0.25


def batchelor_length(nu: float, epsilon: float, D_mol: float) -> float:
    """λ_B = η · Sc^{-1/2},  Sc = ν / D_mol.  epsilon must be in W/kg."""
    eta = kolmogorov_length(nu, epsilon)
    Sc = nu / D_mol if D_mol > 0 else 1e12
    return eta / np.sqrt(Sc)


def epsilon_max_estimate(P: float, rho: float, D: float, N: float) -> float:
    """
    Local maximum dissipation rate near the impeller (order-of-magnitude).
    ε_max ≈ C · P / (ρ D³), C~3 (Kresta & Wood)
    """
    if D == 0:
        return 0.0
    return EPSILON_MAX_COEFF * P / (rho * D**3)


def average_shear_rate(P: float, mu: float, V: float) -> float:
    """Root-mean-square (Camp-Stein) average shear rate.  γ̇_avg = √(P / (μ · V))"""
    if mu <= 0 or V <= 0:
        return 0.0
    return np.sqrt(P / (mu * V))


def maximum_shear_rate(epsilon_max: float, nu: float) -> float:
    """Maximum shear rate near the impeller.  γ̇_max = √(ε_max / ν)"""
    if nu <= 0 or epsilon_max <= 0:
        return 0.0
    return np.sqrt(epsilon_max / nu)


def shear_stress(mu: float, gamma_dot: float) -> float:
    """Newtonian shear stress  τ = μ · γ̇  (Pa)."""
    return mu * gamma_dot
