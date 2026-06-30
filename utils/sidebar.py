"""Shared sidebar branding – call ``render_sidebar()`` at the top of every page."""

import base64
from functools import lru_cache
from pathlib import Path
import streamlit as st

_LOGO = Path(__file__).resolve().parent.parent / "images" / "general" / "logo.png"


@lru_cache(maxsize=1)
def _logo_data_uri() -> str | None:
    """Return the logo encoded as a base64 ``data:`` URI (cached)."""
    if not _LOGO.exists():
        return None
    encoded = base64.b64encode(_LOGO.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def render_sidebar() -> None:
    """Display the Mixing Lab logo above the navigation in the sidebar."""
    logo = _logo_data_uri()
    if logo:
        st.logo(logo)
