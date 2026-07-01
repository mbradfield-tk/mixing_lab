"""Shared sidebar branding – call ``render_sidebar()`` at the top of every page."""

import base64
import io
from functools import lru_cache
from pathlib import Path
import streamlit as st

_LOGO = Path(__file__).resolve().parent.parent / "images" / "general" / "logo.png"

# Displayed height of the sidebar logo, in pixels. Change this to resize.
LOGO_HEIGHT_PX = 140

# Make a solid (near-white) logo background transparent so it blends with the
# app theme. Set to ``False`` to serve the original image unchanged.
MAKE_BG_TRANSPARENT = True
# Colour tolerance for flood-filling the background from the image corners.
_BG_FILL_THRESHOLD = 40


@lru_cache(maxsize=4)
def _transparent_logo_bytes(path_str: str, mtime: float) -> bytes | None:
    """Return PNG bytes with the solid border background flood-filled to alpha 0.

    Only pixels connected to the image edges within ``_BG_FILL_THRESHOLD`` of
    the corner colour are removed, so interior detail (including whites) is
    preserved. Falls back to the raw bytes if Pillow is unavailable.
    """
    path = Path(path_str)
    if not path.exists():
        return None
    raw = path.read_bytes()
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return raw

    im = Image.open(io.BytesIO(raw)).convert("RGBA")
    w, h = im.size
    rgb = im.convert("RGB")
    sentinel = (0, 255, 0)  # marks background region after flood fill
    for corner in ((0, 0), (w - 1, 0), (0, h - 1), (w - 1, h - 1)):
        ImageDraw.floodfill(rgb, corner, sentinel, thresh=_BG_FILL_THRESHOLD)

    px_rgb = rgb.load()
    px_out = im.load()
    for y in range(h):
        for x in range(w):
            if px_rgb[x, y] == sentinel:
                r, g, b, _ = px_out[x, y]
                px_out[x, y] = (r, g, b, 0)

    buf = io.BytesIO()
    im.save(buf, format="PNG")
    return buf.getvalue()


@lru_cache(maxsize=4)
def _logo_data_uri(path_str: str, mtime: float) -> str | None:
    """Return the logo encoded as a base64 ``data:`` URI.

    ``mtime`` is part of the cache key so the cache invalidates automatically
    whenever the image file is modified on disk.
    """
    path = Path(path_str)
    if not path.exists():
        return None
    if MAKE_BG_TRANSPARENT:
        data = _transparent_logo_bytes(path_str, mtime)
    else:
        data = path.read_bytes()
    if not data:
        return None
    encoded = base64.b64encode(data).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def render_sidebar() -> None:
    """Display the Mixing Lab logo above the navigation in the sidebar."""
    if not _LOGO.exists():
        return
    logo = _logo_data_uri(str(_LOGO), _LOGO.stat().st_mtime)
    if not logo:
        return

    # ``st.logo`` is the only element Streamlit renders *above* the navigation
    # menu, so use it for placement, then enlarge it with targeted CSS.
    st.logo(logo, size="large")
    st.markdown(
        f"""
        <style>
        [data-testid="stSidebarHeader"] {{
            height: auto !important;
            min-height: 0 !important;
            max-height: none !important;
            display: flex !important;
            align-items: center;
            justify-content: flex-start;
            overflow: visible !important;
            padding-top: 1rem;
            padding-bottom: 0.5rem;
        }}
        [data-testid="stSidebarLogo"] {{
            height: {LOGO_HEIGHT_PX}px !important;
            max-height: none !important;
            width: auto !important;
            object-fit: contain;
        }}
        </style>
        """,
        unsafe_allow_html=True,
    )
