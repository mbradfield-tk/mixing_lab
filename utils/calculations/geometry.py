"""Bottom-dish geometry helpers and liquid height estimation."""

import numpy as np


def dish_geometry(D_tank: float, dish_type: str = "") -> tuple[float, float]:
    """Return (V_dish_m3, h_dish_m) for a vessel bottom dish."""
    if D_tank <= 0:
        return 0.0, 0.0

    dish = str(dish_type).lower().strip() if dish_type else ""

    if "conic" in dish:
        h_dish = D_tank / 2
        V_dish = np.pi / 12 * D_tank**2 * h_dish
    elif "torisph" in dish or "din" in dish or "dished" in dish:
        h_dish = 0.1935 * D_tank
        V_dish = 0.0847 * D_tank**3
    else:
        h_dish = D_tank / 4
        V_dish = np.pi * D_tank**3 / 24

    return V_dish, h_dish


def liquid_height_from_volume(V_L_litres: float, D_tank: float,
                              H_max: float, dish_type: str = "") -> float:
    """Compute liquid height (m) from fill volume, accounting for bottom dish."""
    if D_tank <= 0 or V_L_litres <= 0:
        return 0.0

    V_L_m3 = V_L_litres / 1000.0
    V_dish, h_dish = dish_geometry(D_tank, dish_type)
    A_cs = np.pi / 4 * D_tank**2

    if V_L_m3 <= V_dish and V_dish > 0:
        frac = V_L_m3 / V_dish
        return frac * h_dish
    else:
        H_cyl = (V_L_m3 - V_dish) / A_cs if A_cs > 0 else 0.0
        H_total = h_dish + H_cyl
        return min(H_total, H_max) if H_max > 0 else H_total
