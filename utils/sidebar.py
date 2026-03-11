"""Shared sidebar branding – call ``render_sidebar()`` at the top of every page."""

from pathlib import Path
import streamlit as st

_LOGO = Path(__file__).resolve().parent.parent / "images" / "general" / "logo.png"


def render_sidebar() -> None:
    """Display the Mixing Lab logo above the navigation in the sidebar."""
    if _LOGO.exists():
        st.logo(str(_LOGO))
