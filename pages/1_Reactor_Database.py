"""Page 1 – Reactor Database: browse, add, edit, import/export reactor geometries."""

import streamlit as st

import pandas as pd
import pathlib
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Arc
import plotly.graph_objects as go

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
REACTOR_CSV = DATA_DIR / "reactors.csv"
CFD_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "CFD"
REACTOR_IMG_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "reactors"
IMG_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp"}

# Columns the app expects (superset – legacy + enriched)
CORE_COLS = [
    "reactor_name", "owner", "type", "scale",
    "D_tank_m", "H_m", "D_imp_m", "impeller_type", "Np", "Nq",
    "N_rpm_min", "N_rpm_max", "N_rps",
    "V_L_min", "V_L_max", "V_L",
    "material", "baffles",
    "bottom_dish", "top_dish",
    "impeller_count",
    "imp1_clearance_m", "imp1_height_m",
    "D_imp2_m", "Np2", "imp2_clearance_m", "imp2_height_m",
    "D_imp3_m", "Np3", "imp3_clearance_m", "imp3_height_m",
    "Zwietering_S", "GMB_z",
    "wall_thickness_mm", "OD_m", "knuckle_radius_m",
    "notes",
]


def _load_reactors() -> pd.DataFrame:
    if REACTOR_CSV.exists():
        df = pd.read_csv(REACTOR_CSV)
        # Ensure all expected columns exist (for older CSVs)
        for c in CORE_COLS:
            if c not in df.columns:
                df[c] = np.nan
        return df
    return pd.DataFrame(columns=CORE_COLS)


def _save_reactors(df: pd.DataFrame):
    df.to_csv(REACTOR_CSV, index=False)


# ── Load into session state ──────────────────────────────────────────────
if "reactor_db" not in st.session_state:
    st.session_state.reactor_db = _load_reactors()
else:
    # Ensure any columns added after the initial load are present
    for c in CORE_COLS:
        if c not in st.session_state.reactor_db.columns:
            st.session_state.reactor_db[c] = np.nan

st.title("⚗️ Reactor Database")

_is_admin = st.session_state.get("admin_authenticated", False)
_ADMIN_HINT = "Log in via Admin Tools to enable editing."

tab_browse, tab_add, tab_import = st.tabs(["Browse & Edit", "Add Reactor", "Import / Export"])

# ── Browse & Edit ─────────────────────────────────────────────────────────
with tab_browse:
    st.markdown("Filter and edit the reactor database. Changes are saved when you click **Save changes**.")

    col_f0, col_f1, col_f2, col_f3 = st.columns(4)
    with col_f0:
        filt_owner = st.multiselect("Owner", options=sorted(st.session_state.reactor_db["owner"].dropna().unique().tolist()), default=None)
    with col_f1:
        filt_scale = st.multiselect("Scale", options=st.session_state.reactor_db["scale"].dropna().unique().tolist(), default=None)
    with col_f2:
        filt_type = st.multiselect("Type", options=st.session_state.reactor_db["type"].dropna().unique().tolist(), default=None)
    with col_f3:
        filt_impeller = st.multiselect("Impeller", options=st.session_state.reactor_db["impeller_type"].dropna().unique().tolist(), default=None)

    df_display = st.session_state.reactor_db.copy()
    if filt_owner:
        df_display = df_display[df_display["owner"].isin(filt_owner)]
    if filt_scale:
        df_display = df_display[df_display["scale"].isin(filt_scale)]
    if filt_type:
        df_display = df_display[df_display["type"].isin(filt_type)]
    if filt_impeller:
        df_display = df_display[df_display["impeller_type"].isin(filt_impeller)]

    _filters_active = bool(filt_owner or filt_scale or filt_type or filt_impeller)

    edited_df = st.data_editor(
        df_display,
        num_rows="dynamic",
        width='content',
        key="reactor_editor",
    )

    col_s1, col_s2 = st.columns([1, 5])
    with col_s1:
        if st.button("💾 Save changes", key="save_reactor_browse",
                     disabled=not _is_admin, help=None if _is_admin else _ADMIN_HINT):
            if _filters_active:
                # Merge edits back: keep unfiltered rows, replace filtered rows with edits
                full = st.session_state.reactor_db.copy()
                full = full.drop(index=df_display.index, errors='ignore')
                st.session_state.reactor_db = pd.concat([full, edited_df], ignore_index=True)
            else:
                st.session_state.reactor_db = edited_df.copy()
            _save_reactors(st.session_state.reactor_db)
            st.success("Reactor database saved.")

    # ── Reactor Viewer ────────────────────────────────────────────────────
    st.divider()
    st.subheader("🔍 Reactor Viewer")

    reactor_names = df_display["reactor_name"].dropna().tolist()
    if reactor_names:
        selected_reactor = st.selectbox(
            "Select a reactor",
            options=["(none)"] + reactor_names,
            key="img_reactor_select",
        )

        if selected_reactor != "(none)":
            # Build a file-matching prefix: "Cambrex – R-802" → "Cambrex_R-802"
            prefix = selected_reactor.replace(" – ", "_").replace(" - ", "_")

            def _find_images(folder: pathlib.Path) -> list[pathlib.Path]:
                if not folder.exists():
                    return []
                return sorted([
                    p for p in folder.iterdir()
                    if p.is_file()
                    and p.suffix.lower() in IMG_SUFFIXES
                    and p.name.startswith(prefix)
                ])

            reactor_imgs = _find_images(REACTOR_IMG_DIR)
            cfd_imgs = _find_images(CFD_DIR)

            # ── Action buttons ────────────────────────────────────────────
            _btn_col1, _btn_col2, _btn_col3, _ = st.columns([1, 1, 1, 2])
            with _btn_col1:
                if st.button("📐 Draw Reactor", key="draw_reactor_btn"):
                    st.session_state["_show_reactor_schematic"] = True
            with _btn_col2:
                _render_mode = st.radio(
                    "Mode", ["2D", "3D"], horizontal=True,
                    key="reactor_render_mode", label_visibility="collapsed",
                )
            with _btn_col3:
                if st.button("🖼️ Show Images", key="show_reactor_imgs_btn"):
                    st.session_state["_show_reactor_images"] = True

            # ── Reactor Schematic ─────────────────────────────────────────
            if st.session_state.get("_show_reactor_schematic"):
                _r = df_display[df_display["reactor_name"] == selected_reactor].iloc[0]
                _D = float(_r["D_tank_m"]) if pd.notna(_r.get("D_tank_m")) else None
                _H = float(_r["H_m"]) if pd.notna(_r.get("H_m")) else None

                if _D is None or _H is None or _D <= 0 or _H <= 0:
                    st.warning("Insufficient geometry data (D_tank_m, H_m) to draw a schematic.")

                # ── Shared geometry helpers ────────────────────────────
                elif True:  # scope block
                    _R = _D / 2.0  # vessel radius
                    _bottom = str(_r.get("bottom_dish", "")).strip() if pd.notna(_r.get("bottom_dish")) else ""
                    _top = str(_r.get("top_dish", "")).strip() if pd.notna(_r.get("top_dish")) else ""
                    _n_imp = int(_r["impeller_count"]) if pd.notna(_r.get("impeller_count")) else 1

                    # Dish depth heuristics (fraction of R)
                    def _dish_depth(dish_type: str, radius: float) -> float:
                        dt = dish_type.lower()
                        if "elliptical" in dt or "2:1" in dt:
                            return radius / 2.0
                        elif "torispherical" in dt or "din" in dt:
                            return radius * 0.25
                        elif "conical" in dt:
                            return radius * 0.6
                        elif "dished" in dt:
                            return radius * 0.2
                        elif "flat" in dt:
                            return 0.0
                        return 0.0

                    _bot_depth = _dish_depth(_bottom, _R)
                    _top_depth = _dish_depth(_top, _R)

                    # Impeller geometry (shared between 2D & 3D)
                    _imp_data = [
                        ("D_imp_m", "imp1_clearance_m", "imp1_height_m"),
                        ("D_imp2_m", "imp2_clearance_m", "imp2_height_m"),
                        ("D_imp3_m", "imp3_clearance_m", "imp3_height_m"),
                    ]
                    _imp_colors_hex = ["#1976D2", "#F57C00", "#388E3C"]

                    _impellers = []  # list of (d_imp, clr_y, h_imp, color)
                    for i in range(_n_imp):
                        d_col, c_col, h_col = _imp_data[i]
                        _d_imp = float(_r[d_col]) if pd.notna(_r.get(d_col)) else None
                        _clr = float(_r[c_col]) if pd.notna(_r.get(c_col)) else None
                        _h_imp = float(_r[h_col]) if pd.notna(_r.get(h_col)) else None
                        if _d_imp is None or _d_imp <= 0:
                            continue
                        if _clr is None or _clr <= 0:
                            _clr = _H * (i + 1) / (_n_imp + 1)
                        if _h_imp is None or _h_imp <= 0:
                            _h_imp = _d_imp * 0.15
                        _impellers.append((_d_imp, _clr, _h_imp, _imp_colors_hex[i % len(_imp_colors_hex)]))

                    _lowest_imp_y = min((c[1] for c in _impellers), default=None)

                    if _render_mode == "2D":
                        # ── 2-D matplotlib schematic ──────────────────
                        _total_h = _bot_depth + _H + _top_depth
                        _margin = max(_total_h * 0.1, 0.02)

                        fig, ax = plt.subplots(1, 1, figsize=(2.5, 4.5))
                        ax.set_aspect("equal")
                        ax.set_xlim(-_R - _margin, _R + _margin)
                        ax.set_ylim(-_bot_depth - _margin, _H + _top_depth + _margin)
                        ax.set_axis_off()

                        _wall_lw = 2.0
                        _wall_color = "#333333"

                        # Straight walls
                        ax.plot([-_R, -_R], [0, _H], color=_wall_color, lw=_wall_lw)
                        ax.plot([_R, _R], [0, _H], color=_wall_color, lw=_wall_lw)

                        # Bottom dish
                        if _bot_depth > 0:
                            ax.add_patch(Arc((0, 0), _D, _bot_depth * 2,
                                             theta1=180, theta2=360,
                                             color=_wall_color, lw=_wall_lw))
                        else:
                            ax.plot([-_R, _R], [0, 0], color=_wall_color, lw=_wall_lw)

                        # Top dish
                        if _top_depth > 0:
                            ax.add_patch(Arc((0, _H), _D, _top_depth * 2,
                                             theta1=0, theta2=180,
                                             color=_wall_color, lw=_wall_lw))
                        else:
                            ax.plot([-_R, _R], [_H, _H], color=_wall_color, lw=_wall_lw)

                        # Impellers
                        for idx_imp, (_d_imp, _clr, _h_imp, _color) in enumerate(_impellers):
                            _r_imp = _d_imp / 2.0
                            rect = patches.FancyBboxPatch(
                                (-_r_imp, _clr - _h_imp / 2.0),
                                _d_imp, _h_imp,
                                boxstyle="round,pad=0.002",
                                facecolor=_color, edgecolor=_color,
                                alpha=0.7, lw=1.5, zorder=4,
                            )
                            ax.add_patch(rect)
                            ax.text(_r_imp + _R * 0.08, _clr,
                                    f"Imp {idx_imp+1}\n⌀{_d_imp*1000:.0f} mm",
                                    fontsize=7, va="center", color=_color)

                        # Shaft
                        _shaft_top = _H + _top_depth * 0.9
                        _shaft_bot = min(_lowest_imp_y - _H * 0.05, 0.0) if _lowest_imp_y is not None else 0.0
                        ax.plot([0, 0], [_shaft_bot, _shaft_top],
                                color="#555555", lw=1.5, zorder=3)

                        # Dimension annotations
                        _dim_color = "#888888"
                        _dim_fs = 7
                        _arr_y = -_bot_depth * 0.6 if _bot_depth > 0 else -_margin * 0.6
                        ax.annotate("", xy=(_R, _arr_y), xytext=(-_R, _arr_y),
                                    arrowprops=dict(arrowstyle="<->", color=_dim_color, lw=1))
                        ax.text(0, _arr_y - _margin * 0.25,
                                f"⌀ {_D*1000:.0f} mm",
                                ha="center", va="top", fontsize=_dim_fs, color=_dim_color)

                        _hx = _R + _margin * 0.4
                        ax.annotate("", xy=(_hx, _H), xytext=(_hx, 0),
                                    arrowprops=dict(arrowstyle="<->", color=_dim_color, lw=1))
                        ax.text(_hx + _margin * 0.15, _H / 2,
                                f"H {_H*1000:.0f} mm",
                                ha="left", va="center", fontsize=_dim_fs,
                                color=_dim_color, rotation=90)

                        ax.set_title(selected_reactor, fontsize=10, fontweight="bold", pad=10)
                        fig.tight_layout()
                        _col_fig, _ = st.columns([1, 2])
                        with _col_fig:
                            st.pyplot(fig)
                        plt.close(fig)

                    else:
                        # ── 3-D Plotly interactive schematic ──────────
                        _n_circ = 60  # circumferential resolution
                        _theta = np.linspace(0, 2 * np.pi, _n_circ)
                        _cos = np.cos(_theta)
                        _sin = np.sin(_theta)

                        _vessel_color = "rgba(180,180,180,0.25)"
                        _traces: list = []

                        # Helper: build a dish surface mesh
                        def _dish_surface(radius, depth, dish_type, z_base, flip_up):
                            """Return (x, y, z) arrays for a dish surface.
                            flip_up=False → bowl goes below z_base; True → dome above."""
                            _n_d = 20
                            dt = dish_type.lower()
                            t_param = np.linspace(0, 1, _n_d)  # 0=rim, 1=centre

                            if "conical" in dt:
                                # Linear cone: r shrinks linearly, z deepens linearly
                                r_profile = radius * (1 - t_param)
                                z_profile = depth * t_param
                            elif "elliptical" in dt or "2:1" in dt:
                                # True 2:1 ellipse: r² / R² + z² / depth² = 1
                                phi = np.linspace(0, np.pi / 2, _n_d)
                                r_profile = radius * np.cos(phi)
                                z_profile = depth * np.sin(phi)
                            elif "torispherical" in dt or "din" in dt:
                                # Approximate torispherical as a flatter ellipse
                                phi = np.linspace(0, np.pi / 2, _n_d)
                                r_profile = radius * np.cos(phi)
                                z_profile = depth * np.sin(phi)
                            elif "dished" in dt:
                                # Shallow spherical segment
                                phi = np.linspace(0, np.pi / 2, _n_d)
                                r_profile = radius * np.cos(phi)
                                z_profile = depth * np.sin(phi)
                            else:
                                # Default ellipsoidal
                                phi = np.linspace(0, np.pi / 2, _n_d)
                                r_profile = radius * np.cos(phi)
                                z_profile = depth * np.sin(phi)

                            xs = np.outer(r_profile, _cos)
                            ys = np.outer(r_profile, _sin)
                            if flip_up:
                                zs = z_base + np.outer(z_profile, np.ones(_n_circ))
                            else:
                                zs = z_base - np.outer(z_profile, np.ones(_n_circ))
                            return xs, ys, zs

                        # -- Cylindrical wall (shaded surface)
                        _n_z_wall = 20
                        _zz_wall = np.linspace(0, _H, _n_z_wall)
                        _x_wall = np.outer(np.ones(_n_z_wall), _R * _cos)
                        _y_wall = np.outer(np.ones(_n_z_wall), _R * _sin)
                        _z_wall = np.outer(_zz_wall, np.ones(_n_circ))
                        _traces.append(go.Surface(
                            x=_x_wall, y=_y_wall, z=_z_wall,
                            colorscale=[[0, _vessel_color], [1, _vessel_color]],
                            showscale=False, opacity=0.3,
                            hoverinfo="skip",
                        ))

                        # -- Bottom dish (surface)
                        if _bot_depth > 0:
                            _xb, _yb, _zb = _dish_surface(_R, _bot_depth, _bottom, 0.0, flip_up=False)
                            _traces.append(go.Surface(
                                x=_xb, y=_yb, z=_zb,
                                colorscale=[[0, _vessel_color], [1, _vessel_color]],
                                showscale=False, opacity=0.3,
                                hoverinfo="skip",
                            ))
                        else:
                            _traces.append(go.Scatter3d(
                                x=(_R * _cos).tolist(), y=(_R * _sin).tolist(),
                                z=[0] * _n_circ,
                                mode="lines", line=dict(color="grey", width=2),
                                showlegend=False, hoverinfo="skip",
                            ))

                        # -- Top dish (surface)
                        if _top_depth > 0:
                            _xt, _yt, _zt = _dish_surface(_R, _top_depth, _top, _H, flip_up=True)
                            _traces.append(go.Surface(
                                x=_xt, y=_yt, z=_zt,
                                colorscale=[[0, _vessel_color], [1, _vessel_color]],
                                showscale=False, opacity=0.3,
                                hoverinfo="skip",
                            ))
                        else:
                            _traces.append(go.Scatter3d(
                                x=(_R * _cos).tolist(), y=(_R * _sin).tolist(),
                                z=[_H] * _n_circ,
                                mode="lines", line=dict(color="grey", width=2),
                                showlegend=False, hoverinfo="skip",
                            ))

                        # -- Shaft
                        _shaft_top = _H + _top_depth * 0.9
                        _shaft_bot = min(_lowest_imp_y - _H * 0.05, 0.0) if _lowest_imp_y is not None else 0.0
                        _shaft_r = _R * 0.03
                        for t in np.linspace(0, 2 * np.pi, 8, endpoint=False):
                            _traces.append(go.Scatter3d(
                                x=[_shaft_r * np.cos(t)] * 2,
                                y=[_shaft_r * np.sin(t)] * 2,
                                z=[_shaft_bot, _shaft_top],
                                mode="lines", line=dict(color="#555555", width=3),
                                showlegend=False, hoverinfo="skip",
                            ))

                        # -- Impellers (3-D cylinders)
                        for idx_imp, (_d_imp, _clr, _h_imp, _color) in enumerate(_impellers):
                            _r_imp = _d_imp / 2.0
                            _z_bot = _clr - _h_imp / 2.0
                            _z_top = _clr + _h_imp / 2.0

                            # Top & bottom disk surfaces
                            _n_r = 10
                            _rr = np.linspace(0, _r_imp, _n_r)
                            for _z_disk in (_z_bot, _z_top):
                                _xi = np.outer(_rr, _cos)
                                _yi = np.outer(_rr, _sin)
                                _zi = np.full_like(_xi, _z_disk)
                                _traces.append(go.Surface(
                                    x=_xi, y=_yi, z=_zi,
                                    colorscale=[[0, _color], [1, _color]],
                                    showscale=False, opacity=0.7,
                                    hoverinfo="skip",
                                ))

                            # Side surface (cylinder wall)
                            _n_z_cyl = 6
                            _zz_cyl = np.linspace(_z_bot, _z_top, _n_z_cyl)
                            _x_cyl = np.outer(np.ones(_n_z_cyl), _r_imp * _cos)
                            _y_cyl = np.outer(np.ones(_n_z_cyl), _r_imp * _sin)
                            _z_cyl = np.outer(_zz_cyl, np.ones(_n_circ))
                            _traces.append(go.Surface(
                                x=_x_cyl, y=_y_cyl, z=_z_cyl,
                                colorscale=[[0, _color], [1, _color]],
                                showscale=False, opacity=0.7,
                                hoverinfo="skip",
                                name=f"Imp {idx_imp+1} (⌀{_d_imp*1000:.0f} mm)",
                            ))

                        # -- Layout
                        _total_z = _bot_depth + _H + _top_depth
                        _max_dim = max(_D, _total_z)
                        _aspect_x = _D / _max_dim
                        _aspect_y = _D / _max_dim
                        _aspect_z = _total_z / _max_dim
                        _pad_z = _total_z * 0.08
                        _pad_xy = _R * 0.15
                        _layout = go.Layout(
                            title=dict(text=selected_reactor, x=0.5),
                            scene=dict(
                                xaxis=dict(visible=False, range=[-_R - _pad_xy, _R + _pad_xy]),
                                yaxis=dict(visible=False, range=[-_R - _pad_xy, _R + _pad_xy]),
                                zaxis=dict(visible=False, range=[
                                    -_bot_depth - _pad_z,
                                    _H + _top_depth + _pad_z,
                                ]),
                                aspectmode="manual",
                                aspectratio=dict(x=_aspect_x, y=_aspect_y, z=_aspect_z),
                                camera=dict(eye=dict(x=1.6, y=1.6, z=0.9)),
                            ),
                            margin=dict(l=0, r=0, t=40, b=0),
                            height=550,
                        )
                        _fig3d = go.Figure(data=_traces, layout=_layout)
                        st.plotly_chart(_fig3d, use_container_width=True)

            # ── Reactor Images ────────────────────────────────────────────
            if st.session_state.get("_show_reactor_images"):
                if not reactor_imgs and not cfd_imgs:
                    st.info(
                        f"No images found for **{selected_reactor}**. "
                        "Place files in `images/reactors/` or `images/CFD/` "
                        "named `{Owner}_{Reactor}_*.png`."
                    )

                if reactor_imgs:
                    st.markdown("#### Reactor Photos")
                    cols = st.columns(min(len(reactor_imgs), 4))
                    for idx, img_path in enumerate(reactor_imgs):
                        with cols[idx % len(cols)]:
                            label = img_path.stem.removeprefix(prefix).lstrip("_") or img_path.stem
                            st.image(str(img_path), caption=label, width='content')

                if cfd_imgs:
                    st.markdown("#### CFD Results")

                    groups: dict[str, list[pathlib.Path]] = {}
                    for m in cfd_imgs:
                        stem = m.stem
                        tag = stem.split("_CFD_")[-1] if "_CFD_" in stem else stem
                        cat = tag.rstrip("0123456789") or tag
                        groups.setdefault(cat, []).append(m)

                    for cat, imgs in groups.items():
                        st.markdown(f"**{cat.replace('_', ' ').title()}**")
                        cols = st.columns(min(len(imgs), 3))
                        for idx, img_path in enumerate(imgs):
                            with cols[idx % len(cols)]:
                                st.image(
                                    str(img_path),
                                    caption=img_path.stem.split("_CFD_")[-1],
                                    width='content',
                                )

    else:
        st.info("No reactors in the filtered view.")

# ── Add Reactor ───────────────────────────────────────────────────────────
with tab_add:
    st.markdown("Add a new reactor to the database.")

    # Impeller count selector OUTSIDE the form so it triggers an immediate rerun
    impeller_count = st.number_input(
        "How many impellers?", min_value=1, max_value=3, value=1,
        key="add_reactor_imp_count",
    )

    with st.form("add_reactor_form"):

        # ── Row 1: Identity & Vessel ──────────────────────────────────────
        st.markdown("##### Identity & Vessel Geometry")
        c1, c2, c3, c4 = st.columns(4)
        with c1:
            name = st.text_input("Reactor name *")
            owner = st.text_input("Owner / site")
            rtype = st.selectbox("Type", ["Batch", "Continuous", "Semi-batch", "Fed-batch"])
            scale = st.selectbox("Scale", ["Lab", "Pilot", "Manufacturing"])
        with c2:
            D_tank = st.number_input("Tank ID (m)", min_value=0.0, value=0.10, format="%.4f")
            H = st.number_input("Height tan-tan (m)", min_value=0.0, value=0.13, format="%.4f")
            OD = st.number_input("Outside diameter (m)", min_value=0.0, value=0.0, format="%.4f")
            wall_thickness = st.number_input("Wall thickness (mm)", min_value=0.0, value=0.0, format="%.2f")
        with c3:
            V_L_min = st.number_input("Volume min (L)", min_value=0.0, value=0.0, format="%.2f")
            V_L_max = st.number_input("Volume max (L)", min_value=0.0, value=1.0, format="%.2f")
            bottom_dish = st.text_input("Bottom dish type", "")
            top_dish = st.text_input("Top dish type", "")
        with c4:
            material = st.text_input("Material of construction", "")
            baffles = st.selectbox("Baffles", ["", "Yes", "No"])
            knuckle_radius = st.number_input("Knuckle radius (m)", min_value=0.0, value=0.0, format="%.4f")

        # ── Row 2: Agitation ──────────────────────────────────────────────
        st.markdown("##### Agitation")
        a1, a2, a3 = st.columns(3)
        with a1:
            N_rpm_min = st.number_input("RPM min", min_value=0.0, value=0.0, format="%.0f")
        with a2:
            N_rpm_max = st.number_input("RPM max", min_value=0.0, value=400.0, format="%.0f")
        with a3:
            N_rps = st.number_input("Default speed (rev/s)", min_value=0.0, value=0.0, format="%.2f",
                                    help="Leave at 0 to auto-compute from RPM midpoint.")

        # ── Impeller 1 (always shown) ────────────────────────────────────
        st.markdown("##### Impeller 1 (primary)")
        i1a, i1b, i1c, i1d = st.columns(4)
        with i1a:
            D_imp = st.number_input("Diameter (m)", min_value=0.0, value=0.05, format="%.4f", key="imp1_d")
            impeller_type = st.selectbox("Type", [
                "Pitched-blade turbine", "Rushton turbine", "Retreat-curve impeller",
                "Anchor", "Magnetic stir bar", "A310 hydrofoil", "A320 hydrofoil", "Other",
            ], key="imp1_type")
        with i1b:
            Np_val = st.number_input("Np", min_value=0.0, value=0.0, format="%.2f", key="imp1_np")
            Nq_val = st.number_input("Nq", min_value=0.0, value=0.0, format="%.2f", key="imp1_nq")
        with i1c:
            imp1_clearance = st.number_input("Clearance (m)", min_value=0.0, value=0.0, format="%.4f", key="imp1_clr")
        with i1d:
            imp1_height = st.number_input("Height (m)", min_value=0.0, value=0.0, format="%.4f", key="imp1_h")

        # ── Impeller 2 (shown when count ≥ 2) ───────────────────────────
        D_imp2 = Np2_val = imp2_clearance = imp2_height = 0.0
        if impeller_count >= 2:
            st.markdown("##### Impeller 2")
            i2a, i2b, i2c, i2d = st.columns(4)
            with i2a:
                D_imp2 = st.number_input("Diameter (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_d")
            with i2b:
                Np2_val = st.number_input("Np", min_value=0.0, value=0.0, format="%.2f", key="imp2_np")
            with i2c:
                imp2_clearance = st.number_input("Clearance (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_clr")
            with i2d:
                imp2_height = st.number_input("Height (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_h")

        # ── Impeller 3 (shown when count ≥ 3) ───────────────────────────
        D_imp3 = Np3_val = imp3_clearance = imp3_height = 0.0
        if impeller_count >= 3:
            st.markdown("##### Impeller 3")
            i3a, i3b, i3c, i3d = st.columns(4)
            with i3a:
                D_imp3 = st.number_input("Diameter (m)", min_value=0.0, value=0.0, format="%.4f", key="imp3_d")
            with i3b:
                Np3_val = st.number_input("Np", min_value=0.0, value=0.0, format="%.2f", key="imp3_np")
            with i3c:
                imp3_clearance = st.number_input("Clearance (m)", min_value=0.0, value=0.0, format="%.4f", key="imp3_clr")
            with i3d:
                imp3_height = st.number_input("Height (m)", min_value=0.0, value=0.0, format="%.4f", key="imp3_h")

        # ── Row 6: Scale-up parameters ───────────────────────────────────
        st.markdown("##### Scale-up Parameters")
        s1, s2 = st.columns(2)
        with s1:
            zwietering_s = st.number_input("Zwietering S parameter", min_value=0.0, value=0.0, format="%.2f")
        with s2:
            gmb_z = st.number_input("GMB z parameter", min_value=0.0, value=0.0, format="%.4f")

        notes = st.text_area("Notes")
        submitted = st.form_submit_button("Add reactor")

        def _or_nan(v):
            """Return NaN for zero / empty values."""
            return v if v else np.nan

        if submitted and name:
            new_row = pd.DataFrame([{
                "reactor_name": name,
                "owner": owner,
                "type": rtype,
                "scale": scale,
                "D_tank_m": _or_nan(D_tank),
                "H_m": _or_nan(H),
                "D_imp_m": _or_nan(D_imp),
                "impeller_type": impeller_type,
                "Np": _or_nan(Np_val),
                "Nq": _or_nan(Nq_val),
                "N_rpm_min": _or_nan(N_rpm_min),
                "N_rpm_max": _or_nan(N_rpm_max),
                "N_rps": N_rps if N_rps > 0 else ((N_rpm_min + N_rpm_max) / 2 / 60 if N_rpm_max > 0 else np.nan),
                "V_L_min": _or_nan(V_L_min),
                "V_L_max": _or_nan(V_L_max),
                "V_L": V_L_max if V_L_max > 0 else np.nan,
                "material": material,
                "baffles": baffles,
                "bottom_dish": bottom_dish,
                "top_dish": top_dish,
                "impeller_count": impeller_count,
                "imp1_clearance_m": _or_nan(imp1_clearance),
                "imp1_height_m": _or_nan(imp1_height),
                "D_imp2_m": _or_nan(D_imp2),
                "Np2": _or_nan(Np2_val),
                "imp2_clearance_m": _or_nan(imp2_clearance),
                "imp2_height_m": _or_nan(imp2_height),
                "D_imp3_m": _or_nan(D_imp3),
                "Np3": _or_nan(Np3_val),
                "imp3_clearance_m": _or_nan(imp3_clearance),
                "imp3_height_m": _or_nan(imp3_height),
                "Zwietering_S": _or_nan(zwietering_s),
                "GMB_z": _or_nan(gmb_z),
                "wall_thickness_mm": _or_nan(wall_thickness),
                "OD_m": _or_nan(OD),
                "knuckle_radius_m": _or_nan(knuckle_radius),
                "notes": notes,
            }])
            st.session_state.reactor_db = pd.concat(
                [st.session_state.reactor_db, new_row], ignore_index=True)
            _save_reactors(st.session_state.reactor_db)
            st.success(f"Added **{name}** to the reactor database.")
        elif submitted:
            st.warning("Please enter a reactor name.")

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.markdown("### Export")
    st.download_button(
        "⬇️ Download reactor database (CSV)",
        data=st.session_state.reactor_db.to_csv(index=False).encode("utf-8"),
        file_name="reactors_export.csv",
        mime="text/csv",
    )
    st.markdown("### Import")
    uploaded = st.file_uploader("Upload a CSV file", type=["csv"], key="reactor_upload")
    if uploaded is not None:
        try:
            new_df = pd.read_csv(uploaded)
            st.dataframe(new_df.head())
            mode = st.radio("Import mode", ["Replace existing database", "Append to existing database"], key="reactor_import_mode")
            if st.button("Confirm import", key="reactor_import_confirm",
                         disabled=not _is_admin, help=None if _is_admin else _ADMIN_HINT):
                if mode.startswith("Replace"):
                    st.session_state.reactor_db = new_df
                else:
                    st.session_state.reactor_db = pd.concat(
                        [st.session_state.reactor_db, new_df], ignore_index=True)
                _save_reactors(st.session_state.reactor_db)
                st.success("Reactor database updated.")
        except Exception as e:
            st.error(f"Error reading CSV: {e}")
