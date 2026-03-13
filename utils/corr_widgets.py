"""Streamlit widgets for correlation mode selection.

Provides:
- ``render_correlation_matrix`` — single-select (radio) per parameter (Page 7).
- ``render_correlation_matrix_multi`` — multi-select per parameter (Page 5).
"""

from __future__ import annotations

import streamlit as st
from utils.rom_registry import (
    SUPPORTED_PARAMS,
    PARAM_DISPLAY,
    available_param_modes,
    has_any_alt_correlations,
)

ALL_MODES = ["Literature", "ROM", "Experimental"]

# Colours used when overlaying envelopes per mode
MODE_COLORS: dict[str, str] = {
    "Literature": "#1f77b4",
    "ROM": "#2ca02c",
    "Experimental": "#ff7f0e",
}


def render_correlation_matrix(
    reactor_name: str,
    key_prefix: str = "corr_matrix",
) -> dict[str, str]:
    """Single-select per parameter (radio buttons). Used on Page 7.

    Returns dict[str, str] mapping every supported param to its chosen mode.
    """
    param_modes = available_param_modes(reactor_name)
    choices: dict[str, str] = {p: "Literature" for p in SUPPORTED_PARAMS}
    alt_params = [p for p in SUPPORTED_PARAMS if len(param_modes[p]) > 1]

    if not alt_params:
        st.info("No ROM or Experimental correlations registered for this reactor.")
        return choices

    for p in alt_params:
        avail = param_modes[p]
        label = PARAM_DISPLAY.get(p, p)
        chosen = st.radio(
            label,
            options=avail,
            horizontal=True,
            key=f"{key_prefix}_{p}",
        )
        choices[p] = chosen

    return choices


def render_correlation_matrix_multi(
    reactor_name: str,
    key_prefix: str = "corr_matrix_m",
) -> dict[str, list[str]]:
    """Multi-select per parameter (checkboxes). Used on Page 5.

    Returns dict[str, list[str]] mapping every supported param to a
    *list* of selected modes.  Literature is always included.
    Parameters with no alt correlations default to ``["Literature"]``.
    """
    param_modes = available_param_modes(reactor_name)
    choices: dict[str, list[str]] = {p: ["Literature"] for p in SUPPORTED_PARAMS}
    alt_params = [p for p in SUPPORTED_PARAMS if len(param_modes[p]) > 1]

    if not alt_params:
        st.info("No ROM or Experimental correlations registered for this reactor.")
        return choices

    # Header
    cols = st.columns([3] + [1] * len(ALL_MODES))
    cols[0].markdown("**Parameter**")
    for i, m in enumerate(ALL_MODES):
        cols[i + 1].markdown(f"**{m}**")

    for p in alt_params:
        avail = set(param_modes[p])
        label = PARAM_DISPLAY.get(p, p)
        cols = st.columns([3] + [1] * len(ALL_MODES))
        cols[0].markdown(label)
        selected: list[str] = []
        for i, m in enumerate(ALL_MODES):
            disabled = m not in avail
            default = (m == "Literature")  # Literature checked by default
            checked = cols[i + 1].checkbox(
                m,
                value=default,
                disabled=disabled,
                key=f"{key_prefix}_{p}_{m}",
                label_visibility="collapsed",
            )
            if checked and not disabled:
                selected.append(m)
        # Ensure at least Literature is selected
        if not selected:
            selected = ["Literature"]
        choices[p] = selected

    return choices


def priority_mode_dict(
    selections: dict[str, list[str]],
) -> dict[str, str]:
    """Collapse multi-selections to a single mode per param using priority.

    Priority: Experimental > ROM > Literature.
    """
    result: dict[str, str] = {}
    for p, modes in selections.items():
        if "Experimental" in modes:
            result[p] = "Experimental"
        elif "ROM" in modes:
            result[p] = "ROM"
        else:
            result[p] = "Literature"
    return result


def active_modes_set(
    selections: dict[str, list[str]],
) -> list[str]:
    """Return deduplicated, deterministically-ordered list of all selected modes."""
    seen: set[str] = set()
    for modes in selections.values():
        seen.update(modes)
    return [m for m in ALL_MODES if m in seen]


def build_mode_dict_for(
    mode_label: str,
    selections: dict[str, list[str]],
) -> dict[str, str]:
    """Build a per-param mode dict for a given *mode_label*.

    For each param, uses *mode_label* if that mode was selected,
    otherwise falls back to Literature.
    """
    return {
        p: mode_label if mode_label in modes else "Literature"
        for p, modes in selections.items()
    }
