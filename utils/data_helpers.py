"""Shared data-access helpers used across multiple pages."""

import pathlib
import numpy as np
import pandas as pd
import streamlit as st

from utils.solvent_properties import SOLVENT_DB, get_properties, is_known_solvent

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
_IMG_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "reactors"
_IMG_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp"}


def build_search_names(df: pd.DataFrame) -> pd.DataFrame:
    """Add a 'search_name' column combining reactor_name, probes, and impeller_type.

    Format: "Reactor Name (probes; impeller_type)"
    If probes is empty, only impeller_type is shown; if both empty, no parenthetical.
    """
    if "reactor_name" not in df.columns:
        return df

    def _make_search_name(row):
        name = str(row.get("reactor_name", "")) if pd.notna(row.get("reactor_name")) else ""
        probes = str(row.get("probes", "")).strip() if pd.notna(row.get("probes")) else ""
        imp_type = str(row.get("impeller_type", "")).strip() if pd.notna(row.get("impeller_type")) else ""
        parts = [p for p in [probes, imp_type] if p]
        if parts:
            return f"{name} ({'; '.join(parts)})"
        return name

    df = df.copy()
    df["search_name"] = df.apply(_make_search_name, axis=1)
    return df


def reactor_search_name(reactors_df: pd.DataFrame, reactor_name: str) -> str:
    """Look up the search_name for a given reactor_name."""
    if "search_name" not in reactors_df.columns:
        return reactor_name
    row = reactors_df[reactors_df["reactor_name"] == reactor_name]
    if row.empty:
        return reactor_name
    return str(row.iloc[0]["search_name"])


def load_db(key: str, filename: str, columns: list[str] | None = None) -> pd.DataFrame:
    """Load a CSV from data/ into session state (once), returning the DataFrame."""
    if key not in st.session_state:
        p = DATA_DIR / filename
        if p.exists():
            try:
                df = pd.read_csv(p)
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
                st.error(f"Could not read '{filename}': {exc}. Starting with an empty table.")
                df = pd.DataFrame(columns=columns or [])
        else:
            df = pd.DataFrame(columns=columns or [])
        # Add search_name for reactor databases
        if "reactor_name" in df.columns:
            df = build_search_names(df)
        st.session_state[key] = df
    return st.session_state[key]


def safe_float(val, default=0.0) -> float:
    """Return float(val) if non-NaN, else default."""
    try:
        v = float(val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


def all_fluid_names(custom_fluids: pd.DataFrame) -> list[str]:
    """Build combined fluid list: built-in solvents + custom fluids."""
    solvent_names = sorted(SOLVENT_DB.keys())
    custom_names = custom_fluids["fluid_name"].tolist() if not custom_fluids.empty else []
    return solvent_names + custom_names


def get_fluid_props(fluid_name: str, custom_fluids: pd.DataFrame,
                    T_C: float = 25.0, P_atm: float = 1.0) -> dict:
    """Resolve fluid properties from built-in solvents or custom fluids.

    Returns a dict with keys: rho, mu, D_mol, sigma, Cp, k_fluid,
    is_solvent, in_range, mp_C, bp_at_P_C (last three only for solvents).
    """
    if is_known_solvent(fluid_name):
        fprops = get_properties(fluid_name, T_C, P_atm)
        return {
            "rho": fprops["rho_kg_m3"],
            "mu": fprops["mu_Pa_s"],
            "D_mol": fprops["D_mol_m2_s"],
            "sigma": fprops.get("surface_tension_N_m", 0.072),
            "Cp": fprops.get("Cp_J_per_kgK", 0.0),
            "k_fluid": fprops.get("k_W_per_mK", 0.0),
            "is_solvent": True,
            "in_range": fprops["in_range"],
            "mp_C": fprops.get("mp_C", 0.0),
            "bp_at_P_C": fprops.get("bp_at_P_C", 0.0),
        }
    else:
        cust = custom_fluids[custom_fluids["fluid_name"] == fluid_name]
        if cust.empty:
            st.error(f"Fluid '{fluid_name}' not found in built-in or custom databases.")
            st.stop()
        row = cust.iloc[0]
        return {
            "rho": float(row["rho_kg_m3"]),
            "mu": float(row["mu_Pa_s"]),
            "D_mol": float(row["D_mol_m2_s"]),
            "sigma": float(row.get("surface_tension_N_m", 0.072)),
            "Cp": float(row.get("Cp_J_per_kgK", 0.0)),
            "k_fluid": float(row.get("k_W_per_mK", 0.0)),
            "is_solvent": False,
            "in_range": True,
            "mp_C": 0.0,
            "bp_at_P_C": 0.0,
        }


def safe_iloc(df: pd.DataFrame, column: str, value, label: str = "record"):
    """Filter df by column == value and return the first row, with guard."""
    filtered = df[df[column] == value]
    if filtered.empty:
        st.error(f"{label} '{value}' not found in database. It may have been deleted.")
        st.stop()
    return filtered.iloc[0]


@st.cache_data(show_spinner=False)
def _build_reactor_image_index(dir_mtime: float) -> dict[str, str]:
    """Map image stem -> file path for all images in the reactors dir.

    Cached on the directory mtime so the (potentially network-backed)
    directory is only scanned when its contents change, not every rerun.
    """
    index: dict[str, str] = {}
    if not _IMG_DIR.exists():
        return index
    for p in _IMG_DIR.iterdir():
        if p.is_file() and p.suffix.lower() in _IMG_SUFFIXES:
            index[p.stem] = str(p)
    return index


def _reactor_image_index() -> dict[str, pathlib.Path]:
    """Return the reactor image index, refreshed when the directory changes."""
    try:
        mtime = _IMG_DIR.stat().st_mtime
    except OSError:
        return {}
    return {stem: pathlib.Path(path)
            for stem, path in _build_reactor_image_index(mtime).items()}


def find_reactor_image(reactors_df: pd.DataFrame, reactor_name: str,
                       suffix: str) -> pathlib.Path | None:
    """Find a reactor image (iso, side, etc.) by reactor_id prefix."""
    row = reactors_df[reactors_df["reactor_name"] == reactor_name]
    if row.empty:
        return None
    prefix = str(row.iloc[0].get("reactor_id", ""))
    if not prefix or prefix == "nan":
        return None
    return _reactor_image_index().get(prefix + "_" + suffix)
