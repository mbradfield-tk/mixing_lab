"""Scale-up helpers for stirred-tank reactors."""


def scale_up_constant_tip_speed(N_small, D_small, D_large):
    """N_large such that tip speed is preserved."""
    return N_small * D_small / D_large


def scale_up_constant_P_V(N_small, D_small, D_large, Np_small=5.0, Np_large=5.0):
    """N_large such that P/V is preserved (geometric similarity assumed)."""
    ratio = (Np_small / Np_large) * (D_small / D_large)**2
    return N_small * ratio**(1.0/3.0)


def scale_up_constant_Re(N_small, D_small, D_large):
    """N_large for constant Re (same fluid)."""
    return N_small * (D_small / D_large)**2
