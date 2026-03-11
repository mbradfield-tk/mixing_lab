"""
Temperature-dependent physical properties for common pharmaceutical solvents.
=============================================================================

Provides density (ρ), dynamic viscosity (μ), surface tension (σ), and
molecular diffusivity (D_mol) as functions of temperature for solvents
routinely used in small-molecule organic synthesis and pharmaceutical
development.

Correlations
------------
* **Density** – linear:  ρ(T) = ρ₂₅ + dρ/dT · (T − 25)  [kg/m³]
* **Viscosity** – Arrhenius:  μ(T) = μ₂₅ · exp[B·(1/T_K − 1/298.15)]
  where B = E_a/R  [K] is the activation-energy parameter.
* **Surface tension** – linear:  σ(T) = σ₂₅ + dσ/dT · (T − 25)  [N/m]
* **Diffusivity** – Stokes-Einstein scaling:
       D(T) = D₂₅ · (T_K / 298.15) · (μ₂₅ / μ(T))

All correlations are anchored at 25 °C so that the known reference values
are recovered exactly at that temperature.

Data sources: Perry's Chemical Engineers' Handbook 9th ed., CRC Handbook
of Chemistry and Physics, Yaws' Handbook, DIPPR, and primary literature.

Usage
-----
>>> from utils.solvent_properties import get_properties, list_solvents
>>> props = get_properties("Water", 37.0)
>>> props["rho_kg_m3"], props["mu_Pa_s"]
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

R_GAS = 8.314  # J/(mol·K)

# ---------------------------------------------------------------------------
# Data class for a single solvent
# ---------------------------------------------------------------------------

@dataclass
class SolventData:
    """Container for temperature-dependent correlation parameters."""
    name: str                          # Display name
    cas: str                           # CAS number
    mw: float                          # Molecular weight (g/mol)

    # Phase boundaries
    mp_C: float                        # Melting point (°C)
    bp_C: float                        # Normal boiling point (°C)

    # Density at 25 °C and linear slope
    rho_25: float                      # kg/m³ at 25 °C
    drho_dT: float                     # kg/m³/°C  (negative = typical)

    # Viscosity at 25 °C and Arrhenius activation energy
    mu_25: float                       # Pa·s at 25 °C
    Ea_mu: float                       # Activation energy for viscous flow (J/mol)

    # Surface tension at 25 °C and linear slope
    sig_25: float                      # N/m at 25 °C
    dsig_dT: float                     # N/m/°C  (negative = typical)

    # Reference diffusivity at 25 °C [m²/s]
    D_ref_25: float

    # Optional
    notes: str = ""


# ---------------------------------------------------------------------------
# Solvent database
# ---------------------------------------------------------------------------

SOLVENT_DB: Dict[str, SolventData] = {}

def _add(s: SolventData):
    SOLVENT_DB[s.name] = s

# ─── Water ────────────────────────────────────────────────────────────────
_add(SolventData(
    name="Water", cas="7732-18-5", mw=18.015,
    mp_C=0.0, bp_C=100.0,
    rho_25=997.0, drho_dT=-0.26,
    mu_25=8.90e-4, Ea_mu=15500.0,
    sig_25=0.0720, dsig_dT=-1.50e-4,
    D_ref_25=2.3e-9,
))

# ─── Methanol ────────────────────────────────────────────────────────────
_add(SolventData(
    name="Methanol", cas="67-56-1", mw=32.04,
    mp_C=-97.6, bp_C=64.7,
    rho_25=787.0, drho_dT=-0.94,
    mu_25=5.44e-4, Ea_mu=10800.0,
    sig_25=0.0223, dsig_dT=-7.7e-5,
    D_ref_25=1.6e-9,
))

# ─── Ethanol ─────────────────────────────────────────────────────────────
_add(SolventData(
    name="Ethanol", cas="64-17-5", mw=46.07,
    mp_C=-114.1, bp_C=78.4,
    rho_25=789.0, drho_dT=-0.85,
    mu_25=1.09e-3, Ea_mu=14000.0,
    sig_25=0.0220, dsig_dT=-8.3e-5,
    D_ref_25=1.2e-9,
))

# ─── Isopropanol (IPA) ──────────────────────────────────────────────────
_add(SolventData(
    name="Isopropanol (IPA)", cas="67-63-0", mw=60.10,
    mp_C=-89.5, bp_C=82.6,
    rho_25=786.0, drho_dT=-0.87,
    mu_25=2.04e-3, Ea_mu=18000.0,
    sig_25=0.0210, dsig_dT=-8.0e-5,
    D_ref_25=1.0e-9,
))

# ─── Acetone ─────────────────────────────────────────────────────────────
_add(SolventData(
    name="Acetone", cas="67-64-1", mw=58.08,
    mp_C=-94.7, bp_C=56.1,
    rho_25=784.0, drho_dT=-1.19,
    mu_25=3.06e-4, Ea_mu=7100.0,
    sig_25=0.0234, dsig_dT=-1.12e-4,
    D_ref_25=2.4e-9,
))

# ─── MEK (Methyl Ethyl Ketone / 2-Butanone) ─────────────────────────────
_add(SolventData(
    name="MEK (2-Butanone)", cas="78-93-3", mw=72.11,
    mp_C=-86.7, bp_C=79.6,
    rho_25=800.0, drho_dT=-1.05,
    mu_25=4.05e-4, Ea_mu=8000.0,
    sig_25=0.0243, dsig_dT=-1.00e-4,
    D_ref_25=1.8e-9,
))

# ─── Acetonitrile ────────────────────────────────────────────────────────
_add(SolventData(
    name="Acetonitrile", cas="75-05-8", mw=41.05,
    mp_C=-43.8, bp_C=82.0,
    rho_25=786.0, drho_dT=-1.00,
    mu_25=3.69e-4, Ea_mu=7500.0,
    sig_25=0.0290, dsig_dT=-1.04e-4,
    D_ref_25=2.4e-9,
))

# ─── DCM (Dichloromethane) ───────────────────────────────────────────────
_add(SolventData(
    name="DCM (Dichloromethane)", cas="75-09-2", mw=84.93,
    mp_C=-96.7, bp_C=39.6,
    rho_25=1326.0, drho_dT=-1.73,
    mu_25=4.13e-4, Ea_mu=6500.0,
    sig_25=0.0280, dsig_dT=-1.18e-4,
    D_ref_25=2.1e-9,
))

# ─── Chloroform ──────────────────────────────────────────────────────────
_add(SolventData(
    name="Chloroform", cas="67-66-3", mw=119.38,
    mp_C=-63.5, bp_C=61.2,
    rho_25=1480.0, drho_dT=-1.73,
    mu_25=5.36e-4, Ea_mu=7000.0,
    sig_25=0.0271, dsig_dT=-1.12e-4,
    D_ref_25=2.0e-9,
))

# ─── Ethyl Acetate ───────────────────────────────────────────────────────
_add(SolventData(
    name="Ethyl Acetate", cas="141-78-6", mw=88.11,
    mp_C=-83.6, bp_C=77.1,
    rho_25=902.0, drho_dT=-1.17,
    mu_25=4.26e-4, Ea_mu=7500.0,
    sig_25=0.0238, dsig_dT=-1.10e-4,
    D_ref_25=2.2e-9,
))

# ─── THF (Tetrahydrofuran) ──────────────────────────────────────────────
_add(SolventData(
    name="THF", cas="109-99-9", mw=72.11,
    mp_C=-108.4, bp_C=66.0,
    rho_25=889.0, drho_dT=-1.05,
    mu_25=4.63e-4, Ea_mu=7200.0,
    sig_25=0.0268, dsig_dT=-9.5e-5,
    D_ref_25=2.0e-9,
))

# ─── Toluene ─────────────────────────────────────────────────────────────
_add(SolventData(
    name="Toluene", cas="108-88-3", mw=92.14,
    mp_C=-95.0, bp_C=110.6,
    rho_25=867.0, drho_dT=-0.87,
    mu_25=5.54e-4, Ea_mu=8500.0,
    sig_25=0.0280, dsig_dT=-1.04e-4,
    D_ref_25=2.0e-9,
))

# ─── DMF (Dimethylformamide) ────────────────────────────────────────────
_add(SolventData(
    name="DMF", cas="68-12-2", mw=73.09,
    mp_C=-60.5, bp_C=153.0,
    rho_25=944.0, drho_dT=-0.87,
    mu_25=8.02e-4, Ea_mu=10000.0,
    sig_25=0.0370, dsig_dT=-1.12e-4,
    D_ref_25=1.5e-9,
))

# ─── DMSO (Dimethyl sulfoxide) ──────────────────────────────────────────
_add(SolventData(
    name="DMSO", cas="67-68-5", mw=78.13,
    mp_C=18.5, bp_C=189.0,
    rho_25=1100.0, drho_dT=-0.73,
    mu_25=1.99e-3, Ea_mu=14500.0,
    sig_25=0.0436, dsig_dT=-1.10e-4,
    D_ref_25=0.9e-9,
))

# ─── Heptane ─────────────────────────────────────────────────────────────
_add(SolventData(
    name="Heptane", cas="142-82-5", mw=100.20,
    mp_C=-90.6, bp_C=98.4,
    rho_25=684.0, drho_dT=-0.81,
    mu_25=3.87e-4, Ea_mu=7500.0,
    sig_25=0.0200, dsig_dT=-8.8e-5,
    D_ref_25=2.5e-9,
))

# ─── Hexane ──────────────────────────────────────────────────────────────
_add(SolventData(
    name="Hexane", cas="110-54-3", mw=86.18,
    mp_C=-95.3, bp_C=68.7,
    rho_25=655.0, drho_dT=-0.90,
    mu_25=2.94e-4, Ea_mu=6500.0,
    sig_25=0.0179, dsig_dT=-9.6e-5,
    D_ref_25=2.7e-9,
))

# ─── MTBE (Methyl tert-butyl ether) ─────────────────────────────────────
_add(SolventData(
    name="MTBE", cas="1634-04-4", mw=88.15,
    mp_C=-108.6, bp_C=55.2,
    rho_25=740.0, drho_dT=-1.09,
    mu_25=3.40e-4, Ea_mu=7000.0,
    sig_25=0.0190, dsig_dT=-9.8e-5,
    D_ref_25=2.3e-9,
))

# ─── Acetic Acid ────────────────────────────────────────────────────────
_add(SolventData(
    name="Acetic Acid", cas="64-19-7", mw=60.05,
    mp_C=16.6, bp_C=117.9,
    rho_25=1049.0, drho_dT=-0.82,
    mu_25=1.13e-3, Ea_mu=10500.0,
    sig_25=0.0271, dsig_dT=-8.4e-5,
    D_ref_25=1.3e-9,
))

# ─── NMP (N-Methyl-2-pyrrolidone) ───────────────────────────────────────
_add(SolventData(
    name="NMP", cas="872-50-4", mw=99.13,
    mp_C=-24.4, bp_C=202.0,
    rho_25=1028.0, drho_dT=-0.70,
    mu_25=1.67e-3, Ea_mu=12000.0,
    sig_25=0.0410, dsig_dT=-9.0e-5,
    D_ref_25=1.1e-9,
))

# ─── 2-MeTHF (2-Methyltetrahydrofuran) ─────────────────────────────────
_add(SolventData(
    name="2-MeTHF", cas="96-47-9", mw=86.13,
    mp_C=-136.0, bp_C=80.3,
    rho_25=855.0, drho_dT=-1.00,
    mu_25=4.60e-4, Ea_mu=7500.0,
    sig_25=0.0245, dsig_dT=-9.5e-5,
    D_ref_25=1.8e-9,
))

# ─── 1,4-Dioxane ────────────────────────────────────────────────────────
_add(SolventData(
    name="1,4-Dioxane", cas="123-91-1", mw=88.11,
    mp_C=11.8, bp_C=101.1,
    rho_25=1033.0, drho_dT=-0.88,
    mu_25=1.18e-3, Ea_mu=12000.0,
    sig_25=0.0330, dsig_dT=-1.02e-4,
    D_ref_25=1.7e-9,
))

# ─── Diethyl Ether ──────────────────────────────────────────────────────
_add(SolventData(
    name="Diethyl Ether", cas="60-29-7", mw=74.12,
    mp_C=-116.3, bp_C=34.6,
    rho_25=713.0, drho_dT=-1.25,
    mu_25=2.22e-4, Ea_mu=5500.0,
    sig_25=0.0170, dsig_dT=-1.12e-4,
    D_ref_25=2.6e-9,
))


# ---------------------------------------------------------------------------
# Property computation
# ---------------------------------------------------------------------------

def _liquid_range(solvent: SolventData) -> tuple[float, float]:
    """Return (T_min_C, T_max_C) for the liquid phase."""
    return (solvent.mp_C, solvent.bp_C)


def density(T_C: float, solvent: SolventData) -> float:
    """Density [kg/m³] at temperature T_C [°C]."""
    return solvent.rho_25 + solvent.drho_dT * (T_C - 25.0)


def viscosity(T_C: float, solvent: SolventData) -> float:
    """Dynamic viscosity [Pa·s] at temperature T_C [°C].

    Arrhenius form: μ(T) = μ₂₅ · exp[B·(1/T − 1/T_ref)]
    where B = E_a / R.
    """
    T_K = T_C + 273.15
    T_ref_K = 298.15
    B = solvent.Ea_mu / R_GAS
    return solvent.mu_25 * math.exp(B * (1.0 / T_K - 1.0 / T_ref_K))


def surface_tension(T_C: float, solvent: SolventData) -> float:
    """Surface tension [N/m] at temperature T_C [°C]."""
    return max(solvent.sig_25 + solvent.dsig_dT * (T_C - 25.0), 0.0)


def diffusivity(T_C: float, solvent: SolventData) -> float:
    """Molecular diffusivity [m²/s] (Stokes-Einstein scaling from 25 °C ref).

    D(T) = D₂₅ · (T_K / 298.15) · (μ₂₅ / μ(T))
    """
    T_K = T_C + 273.15
    mu_T = viscosity(T_C, solvent)
    if mu_T <= 0:
        return solvent.D_ref_25
    return solvent.D_ref_25 * (T_K / 298.15) * (solvent.mu_25 / mu_T)


def get_properties(solvent_name: str, T_C: float) -> dict:
    """Return a dict of all physical properties at the given temperature.

    Keys: rho_kg_m3, mu_Pa_s, D_mol_m2_s, surface_tension_N_m,
          T_C, solvent, bp_C, mp_C, in_range
    """
    if solvent_name not in SOLVENT_DB:
        raise KeyError(f"Unknown solvent: {solvent_name!r}.  "
                       f"Available: {sorted(SOLVENT_DB.keys())}")
    s = SOLVENT_DB[solvent_name]
    T_lo, T_hi = _liquid_range(s)
    in_range = (T_lo <= T_C <= T_hi)
    return {
        "solvent": solvent_name,
        "T_C": T_C,
        "in_range": in_range,
        "rho_kg_m3": density(T_C, s),
        "mu_Pa_s": viscosity(T_C, s),
        "D_mol_m2_s": diffusivity(T_C, s),
        "surface_tension_N_m": surface_tension(T_C, s),
        "bp_C": s.bp_C,
        "mp_C": s.mp_C,
        "mw": s.mw,
        "cas": s.cas,
    }


def list_solvents() -> list[str]:
    """Return sorted list of available solvent names."""
    return sorted(SOLVENT_DB.keys())


def solvent_info_table() -> list[dict]:
    """Return a list-of-dicts summary for display (one row per solvent)."""
    rows = []
    for name, s in sorted(SOLVENT_DB.items()):
        rows.append({
            "Solvent": s.name,
            "CAS": s.cas,
            "MW": s.mw,
            "m.p. (°C)": s.mp_C,
            "b.p. (°C)": s.bp_C,
            "ρ₂₅ (kg/m³)": f"{s.rho_25:.1f}",
            "μ₂₅ (mPa·s)": f"{s.mu_25 * 1000:.3f}",
            "σ₂₅ (mN/m)": f"{s.sig_25 * 1000:.1f}",
            "D₂₅ (m²/s)": f"{s.D_ref_25:.2e}",
        })
    return rows
