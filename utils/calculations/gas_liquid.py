"""Gas-liquid mass transfer: kLa, holdup, bubble diameter, flooding."""

import numpy as np

# Named constants for correlations
HUGHMARK_CONSTANT = 0.505           # Hughmark (1967) gas holdup coefficient
CALDERBANK_D32_CONSTANT = 4.15      # Calderbank (1958) bubble diameter coefficient
CALDERBANK_D32_OFFSET = 0.0009      # Calderbank (1958) bubble diameter offset (m)


def kla_vant_riet(P_V: float, v_s: float, coalescing: bool = True) -> float:
    """Van 't Riet (1979) correlation for kLa in aerated stirred tanks."""
    if v_s <= 0 or P_V <= 0:
        return 0.0
    if coalescing:
        return 0.026 * P_V**0.4 * v_s**0.5
    return 0.002 * P_V**0.7 * v_s**0.2


def kla_surface(epsilon: float, nu: float, D_mol: float,
                D_tank: float, V: float) -> float:
    """Headspace-only (free-surface) kLa – Lamont & Scott (1970)."""
    if epsilon <= 0 or nu <= 0 or D_mol <= 0 or D_tank <= 0 or V <= 0:
        return 0.0
    kL = 0.4 * np.sqrt(D_mol) * (epsilon / nu) ** 0.25
    A_surface = np.pi / 4 * D_tank**2
    a = A_surface / V
    return kL * a


def gas_holdup_hughmark(v_s: float, P_V: float, mu: float,
                        sigma: float, rho: float) -> float:
    """Gas holdup (volume fraction) — simplified Hughmark (1967)."""
    if v_s <= 0 or P_V <= 0 or sigma <= 0:
        return 0.0
    eps_G = HUGHMARK_CONSTANT * v_s**0.47 * P_V**0.4 * (mu / sigma)**0.08
    return min(eps_G, 0.95)


def sauter_bubble_diameter(P_V: float, v_s: float, sigma: float,
                           rho: float) -> float:
    """Sauter mean bubble diameter d₃₂ — Calderbank (1958)."""
    if P_V <= 0 or sigma <= 0 or rho <= 0:
        return 0.0
    d32 = CALDERBANK_D32_CONSTANT * sigma**0.6 / (P_V**0.4 * rho**0.2) + CALDERBANK_D32_OFFSET
    return max(d32, 1e-6)


def gas_flooding_speed(Nq: float, D_imp: float, Q_gas: float) -> float:
    """Minimum impeller speed for complete gas dispersion."""
    if Q_gas <= 0 or D_imp <= 0:
        return 0.0
    Fl_crit = 0.035
    return Q_gas / (Fl_crit * D_imp**3)


def gas_flow_rate_from_vs(v_s: float, D_tank: float) -> float:
    """Convert superficial gas velocity to volumetric flow rate."""
    if v_s <= 0 or D_tank <= 0:
        return 0.0
    return v_s * np.pi / 4 * D_tank**2
