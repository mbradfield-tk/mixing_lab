"""Page 1 – Vessel Database: browse, add, edit, import/export reactor geometries."""

import streamlit as st

import pandas as pd
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Arc
import plotly.graph_objects as go

from utils.data_helpers import (
    build_search_names, reactor_search_name,
    find_reactor_model_3d, render_reactor_3d,
)
from utils.validation import TEMP_MIN_C, TEMP_MAX_C, name_exists

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
REACTOR_CSV = DATA_DIR / "reactors.csv"
MEASUREMENTS_CSV = DATA_DIR / "measured" / "measurements_all.csv"
CFD_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "CFD"
REACTOR_IMG_DIR = pathlib.Path(__file__).resolve().parent.parent / "images" / "reactors"
IMG_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".bmp", ".webp"}

# Load HTM names for heat_transfer_medium dropdown
_HTM_CSV = DATA_DIR / "HTM.csv"


@st.cache_data(show_spinner=False)
def _load_htm_names(csv_path: str, mtime: float) -> list[str]:
    """Load HTM names once; cache keyed on path + file mtime."""
    return pd.read_csv(csv_path)["htm_name"].tolist()


_HTM_NAMES: list[str] = []
if _HTM_CSV.exists():
    _HTM_NAMES = _load_htm_names(str(_HTM_CSV), _HTM_CSV.stat().st_mtime)

# Columns the app expects (superset – legacy + enriched)
CORE_COLS = [
    "reactor_id",
    "reactor_name", "owner", "tag", "location", "manufacturer", "manufacturer_model",
    "type", "scale",
    "D_tank_m", "H_m", "D_imp_m", "impeller_type", "impeller_flow", "impeller_model", "Np", "Nq",
    "N_rpm_min", "N_rpm_max", "N_rps",
    "V_L_min", "V_L_max", "V_L",
    "T_max_C", "P_max_atm",
    "shell_material", "lining", "lining_material", "baffles",
    "bottom_dish", "top_dish",
    "impeller_count",
    "imp1_clearance_m", "imp1_height_m",
    "D_imp2_m", "impeller_type2", "impeller_flow2", "impeller_model2", "Np2", "imp2_clearance_m", "imp2_height_m",
    "D_imp3_m", "impeller_type3", "impeller_flow3", "impeller_model3", "Np3", "imp3_clearance_m", "imp3_height_m",
    "Zwietering_S", "GMB_z",
    "wall_thickness_mm", "OD_m", "knuckle_radius_m",
    "notes",
    "instrumentation", "discharge_location", "insulated",
    "gas_addition", "gas_feed_control",
    "no_ports", "motor_power_kW", "aux_units",
    "cip", "heating_cooling", "heat_transfer_medium", "heat_exchanger",
    "probes",
]


# Columns that must be treated as text (not numeric) in the data editor
_STR_COLS = {
    "reactor_id",
    "reactor_name", "owner", "tag", "location", "manufacturer", "manufacturer_model",
    "type", "scale",
    "impeller_type", "impeller_type2", "impeller_type3",
    "impeller_flow", "impeller_flow2", "impeller_flow3",
    "impeller_model", "impeller_model2", "impeller_model3",
    "shell_material", "lining", "lining_material", "baffles",
    "bottom_dish", "top_dish",
    "notes",
    "instrumentation", "discharge_location", "insulated",
    "gas_addition", "gas_feed_control", "aux_units",
    "cip", "heating_cooling", "heat_transfer_medium", "heat_exchanger",
    "probes",
}


def _load_reactors() -> pd.DataFrame:
    if REACTOR_CSV.exists():
        df = pd.read_csv(REACTOR_CSV)
        # Ensure all expected columns exist (for older CSVs)
        for c in CORE_COLS:
            if c not in df.columns:
                df[c] = np.nan
        # Force string columns to object dtype so data_editor allows text input
        for c in _STR_COLS:
            if c in df.columns:
                df[c] = df[c].astype(object)
        return df
    return pd.DataFrame(columns=CORE_COLS)


def _save_reactors(df: pd.DataFrame):
    df.to_csv(REACTOR_CSV, index=False)


def _next_reactor_id(df: pd.DataFrame) -> str:
    """Generate the next sequential reactor ID (RX-001, RX-002, …)."""
    existing = df["reactor_id"].dropna().tolist()
    nums = []
    for rid in existing:
        rid = str(rid)
        if rid.startswith("RX-") and rid[3:].isdigit():
            nums.append(int(rid[3:]))
    return f"RX-{max(nums, default=0) + 1:03d}"


def _assign_missing_reactor_ids(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[tuple[str, str]]]:
    """Ensure every row has a unique RX-NNN reactor_id.

    Existing IDs are never modified. Rows with a missing or blank reactor_id
    receive the next available sequential ID.

    Returns
    -------
    df : updated copy of the dataframe
    assigned : list of (reactor_name, new_id) for each row that was assigned
    """
    df = df.copy()
    if "reactor_id" not in df.columns:
        df.insert(0, "reactor_id", np.nan)

    # Collect numbers already in use
    used: set[int] = set()
    for rid in df["reactor_id"].dropna():
        s = str(rid).strip()
        if s.upper().startswith("RX-") and s[3:].isdigit():
            used.add(int(s[3:]))

    def _next(used_set: set[int]) -> str:
        n = max(used_set, default=0) + 1
        used_set.add(n)
        return f"RX-{n:03d}"

    assigned: list[tuple[str, str]] = []
    df = df.reset_index(drop=True)
    for idx, row in df.iterrows():
        rid = row.get("reactor_id")
        if pd.isna(rid) or str(rid).strip() == "":
            new_id = _next(used)
            df.at[idx, "reactor_id"] = new_id
            assigned.append((str(row.get("reactor_name", f"row {idx}")), new_id))
    return df, assigned


# ── Load into session state ──────────────────────────────────────────────
if "reactor_db" not in st.session_state:
    _loaded = _load_reactors()
    _loaded, _startup_assigned = _assign_missing_reactor_ids(_loaded)
    if _startup_assigned:
        _save_reactors(_loaded)
        st.session_state["_startup_assigned_ids"] = _startup_assigned
    st.session_state.reactor_db = build_search_names(_loaded)
else:
    # Ensure any columns added after the initial load are present
    for c in CORE_COLS:
        if c not in st.session_state.reactor_db.columns:
            st.session_state.reactor_db[c] = np.nan
    # Keep string columns as object dtype so data_editor allows text input
    for c in _STR_COLS:
        if c in st.session_state.reactor_db.columns:
            st.session_state.reactor_db[c] = st.session_state.reactor_db[c].astype(object)
    # Rebuild search_name in case probes or impeller_type changed
    st.session_state.reactor_db = build_search_names(st.session_state.reactor_db)

st.title("⚗️ Vessel Database")

# Warn if any reactor IDs were auto-assigned on this load
_startup_assigned = st.session_state.pop("_startup_assigned_ids", [])
if _startup_assigned:
    _names_list = "\n".join(
        f"- **{_n}** → `{_id}`" for _n, _id in _startup_assigned
    )
    st.warning(
        f"⚠️ {len(_startup_assigned)} reactor(s) in `reactors.csv` were missing a "
        f"`reactor_id` and have been assigned new IDs (CSV updated automatically):\n\n{_names_list}"
    )

_is_admin = st.session_state.get("admin_authenticated", False)
_ADMIN_HINT = "Log in via Admin Tools to enable editing."

tab_browse, tab_add, tab_import = st.tabs(["Browse & Edit", "Add Reactor", "Import / Export"])

# ── Browse & Edit ─────────────────────────────────────────────────────────
with tab_browse:
    st.markdown("Filter and edit the vessel database. Changes are saved when you click **Save changes**.")

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

    # Column visibility selector
    _all_cols = df_display.columns.tolist()
    _default_cols = [c for c in ["reactor_id", "reactor_name", "owner", "tag", "location", "type", "scale",
                                  "D_tank_m", "H_m", "D_imp_m", "impeller_type",
                                  "N_rpm_min", "N_rpm_max", "V_L_min", "V_L_max",
                                  "shell_material"] if c in _all_cols]
    with st.expander("Column visibility", expanded=False):
        _preset_col1, _preset_col2 = st.columns(2)
        with _preset_col1:
            if st.button("Select all columns", key="cols_all"):
                st.session_state["browse_columns"] = _all_cols
                st.rerun()
        with _preset_col2:
            if st.button("Reset to defaults", key="cols_default"):
                st.session_state["browse_columns"] = _default_cols
                st.rerun()
        visible_cols = st.multiselect(
            "Columns to display",
            options=_all_cols,
            default=st.session_state.get("browse_columns", _default_cols),
            key="browse_columns",
        )

    # Apply column selection for display; keep full df for saving
    _visible_cols = visible_cols if visible_cols else _all_cols
    df_editor = df_display[_visible_cols]

    edited_df = st.data_editor(
        df_editor,
        num_rows="dynamic",
        width='content',
        key="reactor_editor",
    )

    col_s1, col_s2 = st.columns([1, 5])
    with col_s1:
        if st.button("💾 Save changes", key="save_reactor_browse",
                     disabled=not _is_admin, help=None if _is_admin else _ADMIN_HINT):
            # Merge visible-column edits back into the full dataframe
            _save_df = df_display.copy()
            # Update only the columns that were visible/editable
            for col in edited_df.columns:
                if col in _save_df.columns:
                    _save_df[col] = edited_df[col].values

            if _filters_active:
                full = st.session_state.reactor_db.copy()
                full = full.drop(index=df_display.index, errors='ignore')
                st.session_state.reactor_db = pd.concat([full, _save_df], ignore_index=True)
            else:
                st.session_state.reactor_db = _save_df
            _save_reactors(st.session_state.reactor_db)
            st.toast("Vessel database saved.")
            st.rerun()

    # ── Reactor Viewer ────────────────────────────────────────────────────
    st.divider()

    reactor_names = df_display["reactor_name"].dropna().tolist()
    if reactor_names:
        selected_reactor = st.selectbox(
            "Select a reactor",
            options=["(none)"] + reactor_names,
            format_func=lambda n: reactor_search_name(df_display, n),
            key="img_reactor_select",
        )

        if selected_reactor != "(none)":
            st.subheader(f"🔍 {selected_reactor}")
            # Use reactor_id as file-matching prefix (e.g. "RX-012")
            _sel_row = df_display[df_display["reactor_name"] == selected_reactor].iloc[0]
            prefix = str(_sel_row.get("reactor_id", "")) if pd.notna(_sel_row.get("reactor_id")) else ""

            def _find_images(folder: pathlib.Path) -> list[pathlib.Path]:
                if not folder.exists():
                    return []
                return sorted([
                    p for p in folder.iterdir()
                    if p.is_file()
                    and p.suffix.lower() in IMG_SUFFIXES
                    and (p.stem == prefix or p.name.startswith(prefix + "_"))
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
                if st.button("🖼️ CFD", key="show_reactor_cfd_btn"):
                    st.session_state["_show_reactor_cfd"] = True

            # ── Reactor Schematic ─────────────────────────────────────────
            if st.session_state.get("_show_reactor_schematic"):
                _r = _sel_row
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

                    # ── Dish geometry heuristics ──────────────────────────
                    # Approximate the head depth & profile shape from its
                    # description. Depths follow standard pressure-vessel heads,
                    # expressed as a fraction of the tank radius R (= D/2):
                    #   Hemispherical / round       depth = R       (0.50 D)
                    #   DIN 28013 Korbbogen         depth ≈ 0.51 R  (0.255 D)
                    #   Semi-ellipsoidal / 2:1      depth = 0.50 R  (0.25 D)
                    #   DIN 28011 Klöpper           depth ≈ 0.39 R  (0.194 D)
                    #   ASME F&D / torispherical    depth ≈ 0.34 R  (0.17 D)
                    #   Conical                     depth = R       (~45° cone)
                    #   Shallow / standard dished   depth ≈ 0.20 R
                    #   Flat                        depth = 0
                    def _dish_depth(dish_type: str, radius: float) -> float:
                        dt = dish_type.lower().strip()
                        if not dt or "flat" in dt or "none" in dt:
                            return 0.0
                        if "cone" in dt or "conical" in dt:
                            return radius
                        if "hemi" in dt or "round" in dt:
                            return radius
                        if "korbbogen" in dt or "28013" in dt:
                            return radius * 0.51
                        if ("klopper" in dt or "kloepper" in dt or "klöpper" in dt
                                or "28011" in dt or ("din" in dt and "tori" in dt)):
                            return radius * 0.39
                        if ("tori" in dt or "asme" in dt or "f&d" in dt
                                or "f & d" in dt or "flanged" in dt):
                            return radius * 0.34
                        if "ellip" in dt or "2:1" in dt:
                            return radius * 0.50
                        if "dish" in dt:
                            return radius * 0.20
                        return radius * 0.20  # unknown but non-flat → shallow dish

                    def _dish_shape(dish_type: str) -> str:
                        """Classify the dish profile: 'flat', 'cone' or 'curved'."""
                        dt = dish_type.lower().strip()
                        if not dt or "flat" in dt or "none" in dt:
                            return "flat"
                        if "cone" in dt or "conical" in dt:
                            return "cone"
                        return "curved"

                    _bot_depth = _dish_depth(_bottom, _R)
                    _top_depth = _dish_depth(_top, _R)
                    _bot_shape = _dish_shape(_bottom)
                    _top_shape = _dish_shape(_top)

                    # Impeller geometry (shared between 2D & 3D)
                    _imp_data = [
                        ("D_imp_m", "imp1_clearance_m", "imp1_height_m", "impeller_type"),
                        ("D_imp2_m", "imp2_clearance_m", "imp2_height_m", "impeller_type2"),
                        ("D_imp3_m", "imp3_clearance_m", "imp3_height_m", "impeller_type3"),
                    ]
                    _imp_colors_hex = ["#1976D2", "#F57C00", "#388E3C"]

                    # Clearance (off-bottom) is measured from the LOWEST interior
                    # point of the vessel – i.e. the bottom of the dish at
                    # y = -_bot_depth – not from the tangent line at y = 0.
                    # Each tuple stores the absolute centre height ``cy``.
                    _impellers = []  # (d_imp, cy, h_imp, color, itype)
                    for i in range(_n_imp):
                        d_col, c_col, h_col, t_col = _imp_data[i]
                        _d_imp = float(_r[d_col]) if pd.notna(_r.get(d_col)) else None
                        _clr = float(_r[c_col]) if pd.notna(_r.get(c_col)) else None
                        _h_imp = float(_r[h_col]) if pd.notna(_r.get(h_col)) else None
                        _itype = str(_r.get(t_col, "")) if pd.notna(_r.get(t_col)) else ""
                        if _d_imp is None or _d_imp <= 0:
                            continue
                        if _clr is None or _clr <= 0:
                            # Default: spread over the full interior height
                            # (dish bottom → top tangent).
                            _clr = (_bot_depth + _H) * (i + 1) / (_n_imp + 1)
                        if _h_imp is None or _h_imp <= 0:
                            _h_imp = _d_imp * 0.15
                        _cy = -_bot_depth + _clr  # absolute centre height
                        _impellers.append((_d_imp, _cy, _h_imp, _imp_colors_hex[i % len(_imp_colors_hex)], _itype))

                    _lowest_imp_y = min((c[1] for c in _impellers), default=None)

                    # ── Liquid fill level (shared between 2D & 3D) ────────
                    # Interior radius as a function of absolute height z, where
                    # z = 0 is the bottom tangent line and z = -_bot_depth is
                    # the lowest interior point of the bottom dish.
                    def _radius_at(z: float) -> float:
                        if z >= 0.0:
                            return _R  # straight wall (and beyond)
                        if _bot_depth <= 0:
                            return _R  # flat bottom
                        d = -z  # depth below the bottom tangent line
                        if d >= _bot_depth:
                            return 0.0
                        if _bot_shape == "cone":
                            return _R * (z + _bot_depth) / _bot_depth
                        # curved heads approximated as a half-ellipse
                        return _R * float(np.sqrt(max(0.0, 1.0 - (d / _bot_depth) ** 2)))

                    # ── Impeller fit validation ───────────────────────────
                    # Warn if an impeller block intersects the vessel wall
                    # (radially wider than the local interior radius) or spills
                    # past the dish / top boundaries – usually a sign that the
                    # impeller diameter, height or clearance is mis-entered.
                    _imp_warnings: list[str] = []
                    _tol = _R * 1e-3
                    for _wi, (_d_w, _cy_w, _h_w, _c_w, _t_w) in enumerate(_impellers, 1):
                        _r_w = _d_w / 2.0
                        _zb_w = _cy_w - _h_w / 2.0
                        _zt_w = _cy_w + _h_w / 2.0
                        # Narrowest interior radius spanned by the impeller height
                        _local_r = min(_radius_at(z) for z in np.linspace(_zb_w, _zt_w, 12))
                        if _r_w > _R + _tol:
                            _imp_warnings.append(
                                f"Impeller {_wi}: diameter ⌀{_d_w*1000:.0f} mm exceeds the "
                                f"tank ID ⌀{_D*1000:.0f} mm."
                            )
                        elif _r_w > _local_r + _tol:
                            _imp_warnings.append(
                                f"Impeller {_wi}: the blade (⌀{_d_w*1000:.0f} mm) overlaps the "
                                f"dish wall at this height – check the diameter, height or clearance."
                            )
                        if _zb_w < -_bot_depth - _tol:
                            _imp_warnings.append(
                                f"Impeller {_wi}: the block extends below the bottom of the "
                                f"vessel – check the clearance or height."
                            )
                        if _zt_w > _H + _top_depth + _tol:
                            _imp_warnings.append(
                                f"Impeller {_wi}: the block extends above the top of the "
                                f"vessel – check the clearance or height."
                            )
                    if _imp_warnings:
                        st.error(
                            "⚠️ The impeller size or offset might be incorrect — the "
                            "drawn impeller overlaps the vessel walls:\n\n"
                            + "\n".join(f"- {_m}" for _m in _imp_warnings)
                        )

                    # Cumulative interior volume vs. height, less an estimated
                    # impeller metal displacement (bounding swept disc × solidity).
                    _IMP_SOLIDITY = 0.20
                    _z_grid = np.linspace(-_bot_depth, _H, 400)
                    _rad_grid = np.array([_radius_at(z) for z in _z_grid])
                    _area_grid = np.pi * _rad_grid ** 2
                    _dz = (_H + _bot_depth) / (len(_z_grid) - 1) if len(_z_grid) > 1 else 0.0
                    _cum_vessel = np.concatenate(
                        [[0.0], np.cumsum((_area_grid[:-1] + _area_grid[1:]) / 2.0 * _dz)]
                    )
                    _disp_below = np.zeros_like(_z_grid)
                    for (_d_i, _cy_i, _h_i, _c_i, _t_i) in _impellers:
                        _ri = _d_i / 2.0
                        _vdisp = np.pi * _ri ** 2 * _h_i * _IMP_SOLIDITY
                        _z0 = _cy_i - _h_i / 2.0
                        _z1 = _cy_i + _h_i / 2.0
                        if _z1 > _z0:
                            _disp_below += _vdisp * np.clip((_z_grid - _z0) / (_z1 - _z0), 0.0, 1.0)
                    _cap_grid = np.clip(_cum_vessel - _disp_below, 0.0, None)
                    _total_vol_L = float(_cap_grid[-1]) * 1000.0

                    _liquid_level = None  # absolute height of liquid surface (z)
                    _liquid_fill_L = None
                    _show_dims = st.checkbox(
                        "📏 Show dimensions", value=True, key="reactor_show_dims",
                        help="Display the tank diameter and height dimension "
                             "annotations on the 2D schematic.",
                    )
                    _show_liquid = st.checkbox(
                        "💧 Show liquid level", key="reactor_show_liquid",
                        help="Display the liquid surface for a given fill volume. "
                             "The level accounts for the dish and straight-wall "
                             "volumes and an estimated impeller displacement.",
                    )
                    if _show_liquid and _total_vol_L > 0:
                        _liquid_fill_L = st.number_input(
                            "Liquid fill volume (L)", min_value=0.0,
                            max_value=float(round(_total_vol_L, 2)),
                            value=float(round(_total_vol_L * 0.7, 2)),
                            step=float(max(round(_total_vol_L / 100.0, 2), 0.01)),
                            key="reactor_liquid_fill_L",
                            help=f"Working volume to full brim ≈ {_total_vol_L:,.1f} L",
                        )
                        # Invert the capacity curve to find the surface height.
                        _liquid_level = float(np.interp(
                            _liquid_fill_L / 1000.0, _cap_grid, _z_grid))
                        _fill_pct = (_liquid_fill_L / _total_vol_L * 100.0
                                     if _total_vol_L > 0 else 0.0)
                        st.caption(
                            f"Liquid surface at **{_liquid_level * 1000:+.0f} mm** "
                            f"relative to the bottom tangent line "
                            f"({_fill_pct:.0f}% of brim-full volume)."
                        )

                    if _render_mode == "2D":
                        # ── 2-D matplotlib schematic ──────────────────
                        _total_h = _bot_depth + _H + _top_depth

                        # Size-aware padding so dimension labels never crowd or
                        # overlap the vessel. Height dimension lives on the LEFT,
                        # impeller labels on the RIGHT, diameter BELOW the dish.
                        _ref = max(_R, _total_h)        # characteristic size
                        _gap = _ref * 0.10              # offset of dim lines from vessel
                        _left_pad = _R * 0.95           # room for height dimension (left)
                        _right_pad = _R * 1.05          # room for impeller labels (right)
                        _bot_pad = _gap + _ref * 0.16   # room for diameter dimension (below)
                        _top_pad = _ref * 0.06
                        if _liquid_level is not None:
                            _top_pad += _ref * 0.10     # room for the volume label

                        fig, ax = plt.subplots(1, 1, figsize=(3.4, 4.8))
                        ax.set_aspect("equal")
                        ax.set_xlim(-_R - _left_pad, _R + _right_pad)
                        ax.set_ylim(-_bot_depth - _bot_pad, _H + _top_depth + _top_pad)
                        ax.set_axis_off()

                        _wall_lw = 2.0
                        _wall_color = "#333333"

                        # Straight walls
                        ax.plot([-_R, -_R], [0, _H], color=_wall_color, lw=_wall_lw)
                        ax.plot([_R, _R], [0, _H], color=_wall_color, lw=_wall_lw)

                        # Bottom dish (shape-aware)
                        if _bot_depth <= 0 or _bot_shape == "flat":
                            ax.plot([-_R, _R], [0, 0], color=_wall_color, lw=_wall_lw)
                        elif _bot_shape == "cone":
                            ax.plot([-_R, 0], [0, -_bot_depth], color=_wall_color, lw=_wall_lw)
                            ax.plot([_R, 0], [0, -_bot_depth], color=_wall_color, lw=_wall_lw)
                        else:
                            ax.add_patch(Arc((0, 0), _D, _bot_depth * 2,
                                             theta1=180, theta2=360,
                                             color=_wall_color, lw=_wall_lw))

                        # Top dish (shape-aware)
                        if _top_depth <= 0 or _top_shape == "flat":
                            ax.plot([-_R, _R], [_H, _H], color=_wall_color, lw=_wall_lw)
                        elif _top_shape == "cone":
                            ax.plot([-_R, 0], [_H, _H + _top_depth], color=_wall_color, lw=_wall_lw)
                            ax.plot([_R, 0], [_H, _H + _top_depth], color=_wall_color, lw=_wall_lw)
                        else:
                            ax.add_patch(Arc((0, _H), _D, _top_depth * 2,
                                             theta1=0, theta2=180,
                                             color=_wall_color, lw=_wall_lw))

                        # Liquid fill (drawn behind the impellers)
                        if _liquid_level is not None:
                            _z_liq = np.linspace(-_bot_depth, _liquid_level, 80)
                            _r_liq = np.array([_radius_at(z) for z in _z_liq])
                            _xs_liq = np.concatenate([_r_liq, -_r_liq[::-1]])
                            _ys_liq = np.concatenate([_z_liq, _z_liq[::-1]])
                            ax.fill(_xs_liq, _ys_liq, color="#4FC3F7",
                                    alpha=0.30, lw=0, zorder=1)
                            _r_surf = _radius_at(_liquid_level)
                            ax.plot([-_r_surf, _r_surf], [_liquid_level, _liquid_level],
                                    color="#0288D1", lw=1.4, zorder=2)
                            # Volume label placed just above the top of the vessel
                            ax.text(0, _H + _top_depth + _ref * 0.02,
                                    f"{_liquid_fill_L:,.1f} L",
                                    ha="center", va="bottom", fontsize=7,
                                    color="#0277BD", zorder=2)

                        # Impellers (clearance measured from dish bottom)
                        for idx_imp, (_d_imp, _cy, _h_imp, _color, _itype) in enumerate(_impellers):
                            _r_imp = _d_imp / 2.0
                            if "chevron" in _itype.lower():
                                # Draw a downward chevron (V-band) instead of a
                                # plain rectangle to convey the blade geometry.
                                # Chevron impellers run in conical-bottom
                                # vessels, so the blade V follows the same angle
                                # as the cone (drop/run = _bot_depth / _R).
                                if _bot_shape == "cone" and _R > 0:
                                    _v_drop = _r_imp * (_bot_depth / _R)
                                else:
                                    _v_drop = max(_h_imp, _r_imp * 0.5)
                                _half_v = (_v_drop + _h_imp) / 2.0
                                _pts = [
                                    (-_r_imp, _cy + _half_v),
                                    (0.0,     _cy + _half_v - _v_drop),
                                    (_r_imp,  _cy + _half_v),
                                    (_r_imp,  _cy + _half_v - _h_imp),
                                    (0.0,     _cy - _half_v),
                                    (-_r_imp, _cy + _half_v - _h_imp),
                                ]
                                ax.add_patch(patches.Polygon(
                                    _pts, closed=True,
                                    facecolor=_color, edgecolor=_color,
                                    alpha=0.7, lw=1.5, zorder=4,
                                ))
                            else:
                                ax.add_patch(patches.FancyBboxPatch(
                                    (-_r_imp, _cy - _h_imp / 2.0),
                                    _d_imp, _h_imp,
                                    boxstyle="round,pad=0.002",
                                    facecolor=_color, edgecolor=_color,
                                    alpha=0.7, lw=1.5, zorder=4,
                                ))
                            # Leader line from the impeller tip to its label
                            if _show_dims:
                                ax.plot([_r_imp, _R + _right_pad * 0.12], [_cy, _cy],
                                        color=_color, lw=0.6, alpha=0.5, zorder=3)
                                ax.text(_R + _right_pad * 0.15, _cy,
                                        f"Imp {idx_imp+1}  ⌀{_d_imp*1000:.0f} mm",
                                        fontsize=7, va="center", ha="left", color=_color)
                        # Shaft
                        _shaft_top = _H + _top_depth * 0.9
                        _shaft_bot = min(_lowest_imp_y - _H * 0.05, 0.0) if _lowest_imp_y is not None else 0.0
                        ax.plot([0, 0], [_shaft_bot, _shaft_top],
                                color="#555555", lw=1.5, zorder=3)

                        # ── Dimension annotations ─────────────────────
                        if _show_dims:
                            _dim_color = "#888888"
                            _dim_fs = 7
                            _wit_lw = 0.6  # witness (extension) line weight

                            # Diameter: dimension line fully BELOW the bottom dish
                            _arr_y = -_bot_depth - _gap
                            for _sx in (-_R, _R):
                                ax.plot([_sx, _sx], [0, _arr_y],
                                        color=_dim_color, lw=_wit_lw, zorder=2)
                            ax.annotate("", xy=(_R, _arr_y), xytext=(-_R, _arr_y),
                                        arrowprops=dict(arrowstyle="<->", color=_dim_color, lw=1))
                            ax.text(0, _arr_y - _ref * 0.04,
                                    f"⌀ {_D*1000:.0f} mm",
                                    ha="center", va="top", fontsize=_dim_fs, color=_dim_color)

                            # Height: dimension line to the LEFT of the vessel
                            _hx = -_R - _left_pad * 0.5
                            for _sy in (0.0, _H):
                                ax.plot([-_R, _hx], [_sy, _sy],
                                        color=_dim_color, lw=_wit_lw, zorder=2)
                            ax.annotate("", xy=(_hx, _H), xytext=(_hx, 0),
                                        arrowprops=dict(arrowstyle="<->", color=_dim_color, lw=1))
                            ax.text(_hx - _R * 0.06, _H / 2,
                                    f"H {_H*1000:.0f} mm",
                                    ha="right", va="center", fontsize=_dim_fs,
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

                            if "cone" in dt or "conical" in dt:
                                # Linear cone: r shrinks linearly, z deepens linearly
                                r_profile = radius * (1 - t_param)
                                z_profile = depth * t_param
                            elif "hemi" in dt or "round" in dt:
                                # Hemisphere: r² + z² = R² (depth = radius)
                                phi = np.linspace(0, np.pi / 2, _n_d)
                                r_profile = radius * np.cos(phi)
                                z_profile = radius * np.sin(phi)
                            else:
                                # All other curved heads (2:1 ellipsoidal,
                                # torispherical, Klöpper, Korbbogen, dished, …)
                                # are approximated by an ellipse of the matching
                                # depth so 2-D and 3-D views stay consistent.
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

                        # -- Liquid surface (translucent disc at the fill level)
                        if _liquid_level is not None:
                            _r_surf = _radius_at(_liquid_level)
                            if _r_surf > 0:
                                _n_rs = 12
                                _rr_s = np.linspace(0, _r_surf, _n_rs)
                                _xs_s = np.outer(_rr_s, _cos)
                                _ys_s = np.outer(_rr_s, _sin)
                                _zs_s = np.full_like(_xs_s, _liquid_level)
                                _traces.append(go.Surface(
                                    x=_xs_s, y=_ys_s, z=_zs_s,
                                    colorscale=[[0, "#4FC3F7"], [1, "#4FC3F7"]],
                                    showscale=False, opacity=0.45,
                                    hoverinfo="skip",
                                    name=f"Liquid {_liquid_fill_L:,.1f} L",
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

                        # -- Impellers (3-D)
                        for idx_imp, (_d_imp, _clr, _h_imp, _color, _itype) in enumerate(_impellers):
                            _r_imp = _d_imp / 2.0

                            if "chevron" in _itype.lower():
                                # Chevron: two coaxial cones (rim high, apex low).
                                # The blade V follows the conical-bottom angle
                                # (drop/run = _bot_depth / _R).
                                if _bot_shape == "cone" and _R > 0:
                                    _v_drop = _r_imp * (_bot_depth / _R)
                                else:
                                    _v_drop = max(_h_imp, _r_imp * 0.5)
                                _half_v = (_v_drop + _h_imp) / 2.0
                                _n_r = 12
                                _rr = np.linspace(0, _r_imp, _n_r)
                                _frac = _rr / _r_imp if _r_imp > 0 else _rr
                                _z_top_prof = (_clr + _half_v - _v_drop) + _v_drop * _frac
                                _z_bot_prof = (_clr - _half_v) + _v_drop * _frac
                                for _zp in (_z_top_prof, _z_bot_prof):
                                    _xi = np.outer(_rr, _cos)
                                    _yi = np.outer(_rr, _sin)
                                    _zi = np.outer(_zp, np.ones(_n_circ))
                                    _traces.append(go.Surface(
                                        x=_xi, y=_yi, z=_zi,
                                        colorscale=[[0, _color], [1, _color]],
                                        showscale=False, opacity=0.7,
                                        hoverinfo="skip",
                                        name=f"Imp {idx_imp+1} chevron (⌀{_d_imp*1000:.0f} mm)",
                                    ))
                                # Outer rim band joining the two cones
                                _zz_rim = np.linspace(_clr + _half_v - _h_imp, _clr + _half_v, 4)
                                _x_rim = np.outer(np.ones(4), _r_imp * _cos)
                                _y_rim = np.outer(np.ones(4), _r_imp * _sin)
                                _z_rim = np.outer(_zz_rim, np.ones(_n_circ))
                                _traces.append(go.Surface(
                                    x=_x_rim, y=_y_rim, z=_z_rim,
                                    colorscale=[[0, _color], [1, _color]],
                                    showscale=False, opacity=0.7, hoverinfo="skip",
                                ))
                                continue

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
                        st.plotly_chart(_fig3d, width='stretch')

            # ── Navigable 3D Model (shown automatically when a GLB exists) ─
            _model_path = find_reactor_model_3d(df_display, selected_reactor)
            if _model_path is not None:
                st.markdown("#### Interactive 3D Model")
                with st.container(border=True):
                    render_reactor_3d(_model_path, height=420, auto_rotate=False)
                    st.caption("Drag to rotate · scroll to zoom · right-drag to pan")

            # ── Reactor Images (shown automatically) ─────────────────────
            if reactor_imgs:
                st.markdown("#### Reactor Photos")
                with st.container(border=True):
                    cols = st.columns(min(len(reactor_imgs), 4))
                    for idx, img_path in enumerate(reactor_imgs):
                        with cols[idx % len(cols)]:
                            label = img_path.stem.removeprefix(prefix).lstrip("_") or img_path.stem
                            st.image(str(img_path), caption=label, width='stretch')
            elif not cfd_imgs:
                st.info(
                    f"No images found for **{selected_reactor}** (ID: `{prefix}`). "
                    "Place files in `images/reactors/` or `images/CFD/` "
                    f"named `{prefix}_*.png`."
                )

            # ── CFD Images (shown on button click) ────────────────────────
            if st.session_state.get("_show_reactor_cfd"):
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
                                    width='stretch',
                                )
                else:
                    st.info(
                        f"No CFD images found for **{selected_reactor}** (ID: `{prefix}`). "
                        f"Place files in `images/CFD/` named `{prefix}_CFD_*.png`."
                    )

            # ── Measured Data ─────────────────────────────────────────────
            st.divider()
            st.markdown("#### 📊 Data")

            if MEASUREMENTS_CSV.exists():
                _meas_all = pd.read_csv(MEASUREMENTS_CSV)
                _meas_rx = _meas_all[_meas_all["reactor_id"] == prefix].copy()

                if _meas_rx.empty:
                    st.info(f"No measured data found for **{selected_reactor}** (`{prefix}`) in `measurements_all.csv`.")
                else:
                    # Ensure expected columns exist
                    if "type" not in _meas_rx.columns:
                        _meas_rx["type"] = "measured"
                    if "fluid" not in _meas_rx.columns:
                        _meas_rx["fluid"] = "unknown"

                    _meas_rx["value"] = pd.to_numeric(_meas_rx["value"], errors="coerce")
                    _meas_rx = _meas_rx.dropna(subset=["variable", "value"])

                    # Variable labels: "variable (units)"
                    _var_units = (
                        _meas_rx.groupby("variable")["units"]
                        .first()
                        .to_dict()
                    )
                    def _vlabel(v: str) -> str:
                        u = _var_units.get(v, "")
                        return f"{v} ({u})" if u else v

                    _set_vars = sorted(_meas_rx.loc[_meas_rx["type"] == "set", "variable"].unique().tolist())
                    _meas_vars = sorted(_meas_rx.loc[_meas_rx["type"] == "measured", "variable"].unique().tolist())

                    if not _meas_vars or not _set_vars:
                        st.warning(
                            "Need at least one `set` variable and one `measured` variable to plot. "
                            "Check the `type` column in `measurements_all.csv`."
                        )
                    else:
                        _dc1, _dc2, _dc3, _dc4 = st.columns(4)
                        with _dc1:
                            _y_var = st.selectbox(
                                "Y-axis (measured)",
                                options=_meas_vars,
                                format_func=_vlabel,
                                key="data_y_var",
                            )
                        with _dc2:
                            _x_default = _set_vars.index("stir speed") if "stir speed" in _set_vars else 0
                            _x_var = st.selectbox(
                                "X-axis (set)",
                                options=_set_vars,
                                index=_x_default,
                                format_func=_vlabel,
                                key="data_x_var",
                            )
                        with _dc3:
                            _series_options = [v for v in _set_vars if v != _x_var]
                            _series_default_label = "fill volume"
                            if _series_default_label in _series_options:
                                _series_default = _series_options.index(_series_default_label) + 1  # +1 for "(none)"
                            else:
                                _series_default = 0
                            _series_var = st.selectbox(
                                "Series (set)",
                                options=["(none)"] + _series_options,
                                index=_series_default,
                                key="data_series_var",
                            )
                        with _dc4:
                            _fluids = sorted(_meas_rx["fluid"].dropna().unique().tolist())
                            _tile_fluid = st.checkbox(
                                "Tile by fluid",
                                value=False,
                                key="data_tile_fluid",
                                disabled=len(_fluids) < 2,
                            )

                        # Pivot long → wide per (exp_id, meas_id, fluid)
                        _pivot_cols = ["exp_id", "meas_id", "fluid"]
                        _wide = _meas_rx.pivot_table(
                            index=_pivot_cols,
                            columns="variable",
                            values="value",
                            aggfunc="first",
                        ).reset_index()

                        # Check selected variables exist in pivoted data
                        _needed = [_x_var, _y_var]
                        if _series_var != "(none)":
                            _needed.append(_series_var)
                        _missing_vars = [v for v in _needed if v not in _wide.columns]
                        if _missing_vars:
                            st.warning(f"Variables not found in data: {', '.join(_missing_vars)}")
                        else:
                            _plot_df = _wide.dropna(subset=[_x_var, _y_var]).copy()

                            if _plot_df.empty:
                                st.info("No data points available for the selected variable combination.")
                            else:
                                _x_label = _vlabel(_x_var)
                                _y_label = _vlabel(_y_var)

                                if _series_var != "(none)":
                                    _plot_df[_series_var] = _plot_df[_series_var].astype(str)

                                # Aggregate duplicate x-y pairs: mean ± std
                                def _agg_series(df_in: pd.DataFrame) -> pd.DataFrame:
                                    agg = df_in.groupby(_x_var, as_index=False)[_y_var].agg(
                                        _y_mean="mean", _y_std="std", _y_count="count"
                                    )
                                    agg["_y_std"] = agg["_y_std"].fillna(0)
                                    return agg.sort_values(_x_var)

                                import plotly.colors as _pc
                                _color_seq = _pc.qualitative.Plotly
                                _color_state = [0]  # mutable counter

                                def _next_color() -> str:
                                    c = _color_seq[_color_state[0] % len(_color_seq)]
                                    _color_state[0] += 1
                                    return c

                                if _tile_fluid and len(_fluids) >= 2:
                                    from plotly.subplots import make_subplots
                                    _plot_fluids = sorted(_plot_df["fluid"].dropna().unique().tolist())
                                    _n_fluids = len(_plot_fluids)
                                    _fig_data = make_subplots(
                                        rows=1, cols=_n_fluids,
                                        subplot_titles=_plot_fluids,
                                        shared_yaxes=True,
                                    )
                                    _legend_shown: set[str] = set()
                                    _series_colors: dict[str, str] = {}
                                    for _fi, _fl in enumerate(_plot_fluids, 1):
                                        _sub = _plot_df[_plot_df["fluid"] == _fl]
                                        if _series_var != "(none)":
                                            for _sv in sorted(_sub[_series_var].unique()):
                                                _ss = _sub[_sub[_series_var] == _sv]
                                                _sname = f"{_series_var}={_sv}"
                                                if _sname not in _series_colors:
                                                    _series_colors[_sname] = _next_color()
                                                _sc = _series_colors[_sname]
                                                _agg = _agg_series(_ss)
                                                _fig_data.add_trace(
                                                    go.Scatter(
                                                        x=_agg[_x_var], y=_agg["_y_mean"],
                                                        error_y=dict(type="data", array=_agg["_y_std"].tolist(), visible=True),
                                                        mode="markers+lines",
                                                        marker=dict(color=_sc),
                                                        line=dict(color=_sc),
                                                        name=_sname,
                                                        legendgroup=_sname,
                                                        showlegend=_sname not in _legend_shown,
                                                    ),
                                                    row=1, col=_fi,
                                                )
                                                _legend_shown.add(_sname)
                                        else:
                                            _sc = _next_color()
                                            _agg = _agg_series(_sub)
                                            _fig_data.add_trace(
                                                go.Scatter(
                                                    x=_agg[_x_var], y=_agg["_y_mean"],
                                                    error_y=dict(type="data", array=_agg["_y_std"].tolist(), visible=True),
                                                    mode="markers+lines",
                                                    marker=dict(color=_sc),
                                                    line=dict(color=_sc),
                                                    name=_fl,
                                                    showlegend=True,
                                                ),
                                                row=1, col=_fi,
                                            )
                                    _fig_data.update_xaxes(title_text=_x_label)
                                    _fig_data.update_yaxes(title_text=_y_label, col=1)
                                    _fig_data.update_layout(height=450, title=f"{_y_label} vs {_x_label}")
                                    st.plotly_chart(_fig_data, width='stretch')
                                else:
                                    _fig_data = go.Figure()
                                    if _series_var != "(none)":
                                        for _sv in sorted(_plot_df[_series_var].unique()):
                                            _ss = _plot_df[_plot_df[_series_var] == _sv]
                                            _sc = _next_color()
                                            _agg = _agg_series(_ss)
                                            _fig_data.add_trace(go.Scatter(
                                                x=_agg[_x_var], y=_agg["_y_mean"],
                                                error_y=dict(type="data", array=_agg["_y_std"].tolist(), visible=True),
                                                mode="markers+lines",
                                                marker=dict(color=_sc),
                                                line=dict(color=_sc),
                                                name=f"{_series_var}={_sv}",
                                            ))
                                    else:
                                        _sc = _next_color()
                                        _agg = _agg_series(_plot_df)
                                        _fig_data.add_trace(go.Scatter(
                                            x=_agg[_x_var], y=_agg["_y_mean"],
                                            error_y=dict(type="data", array=_agg["_y_std"].tolist(), visible=True),
                                            mode="markers+lines",
                                            marker=dict(color=_sc),
                                            line=dict(color=_sc),
                                            name=_y_var,
                                        ))
                                    _fig_data.update_layout(
                                        xaxis_title=_x_label,
                                        yaxis_title=_y_label,
                                        title=f"{_y_label} vs {_x_label}",
                                        height=450,
                                    )
                                    st.plotly_chart(_fig_data, width='stretch')

                                # Data table
                                _table_cols = ["exp_id", "meas_id", "fluid"]
                                _table_cols += [c for c in _needed if c not in _table_cols]
                                _table_cols = [c for c in _table_cols if c in _plot_df.columns]
                                st.markdown("##### Plotted data")
                                st.dataframe(_plot_df[_table_cols], width='stretch', hide_index=True)
            else:
                st.info("`measurements_all.csv` not found in `data/measured/`.")

    else:
        st.info("No reactors in the filtered view.")

# ── Add Reactor ───────────────────────────────────────────────────────────
with tab_add:
    st.markdown("Add a new reactor to the database.")

    # ── Template selector (outside form so rerun populates fields) ────
    _all_reactor_names = st.session_state.reactor_db["reactor_name"].dropna().tolist()
    template_reactor = st.selectbox(
        "Use existing reactor as template",
        options=["(none)"] + _all_reactor_names,
        key="add_reactor_template",
    )

    _TYPE_OPTIONS = ["Batch", "Continuous", "Semi-batch", "Fed-batch"]
    _SCALE_OPTIONS = ["Lab", "Pilot", "Manufacturing"]
    _BAFFLE_OPTIONS = ["", "Yes", "No"]
    _YES_NO_OPTIONS = ["", "Yes", "No"]
    _DISCHARGE_OPTIONS = ["", "Top", "Bottom", "Side"]
    _GAS_ADDITION_OPTIONS = ["", "Headspace", "Sparged"]
    _GAS_FEED_CONTROL_OPTIONS = ["", "Continuous", "Constant Pressure"]
    _AUX_UNITS_OPTIONS = ["", "Sparger", "Mixer", "Condenser", "Baffle"]
    _HEAT_EXCHANGER_OPTIONS = ["", "Jacket", "Internal"]
    _HTM_OPTIONS = [""] + _HTM_NAMES
    _IMP_TYPE_OPTIONS = [
        "Pitched-blade turbine", "Rushton turbine", "Retreat-curve impeller",
        "Anchor", "Magnetic stir bar", "A310 hydrofoil", "A320 hydrofoil", "Other",
    ]
    _IMP_FLOW_OPTIONS = ["", "Axial", "Radial"]

    if template_reactor != "(none)" and st.session_state.get("_applied_template") != template_reactor:
        st.session_state["_applied_template"] = template_reactor
        _tpl = st.session_state.reactor_db[
            st.session_state.reactor_db["reactor_name"] == template_reactor
        ].iloc[0]

        def _tv(col, default=0.0):
            v = _tpl.get(col)
            return float(v) if pd.notna(v) else default

        def _ts(col, default=""):
            v = _tpl.get(col)
            return str(v).strip() if pd.notna(v) else default

        def _tsel(col, options, default_idx=0):
            v = _ts(col)
            return v if v in options else options[default_idx]

        # Push template values into widget session-state keys
        st.session_state["add_owner"] = _ts("owner")
        st.session_state["add_tag"] = _ts("tag")
        st.session_state["add_location"] = _ts("location")
        st.session_state["add_manufacturer"] = _ts("manufacturer")
        st.session_state["add_manufacturer_model"] = _ts("manufacturer_model")
        st.session_state["add_type"] = _tsel("type", _TYPE_OPTIONS)
        st.session_state["add_scale"] = _tsel("scale", _SCALE_OPTIONS)
        st.session_state["add_D_tank"] = _tv("D_tank_m", 0.10)
        st.session_state["add_H"] = _tv("H_m", 0.13)
        st.session_state["add_OD"] = _tv("OD_m")
        st.session_state["add_wall"] = _tv("wall_thickness_mm")
        st.session_state["add_V_min"] = _tv("V_L_min")
        st.session_state["add_V_max"] = _tv("V_L_max", 1.0)
        st.session_state["add_T_max"] = _tv("T_max_C")
        st.session_state["add_P_max"] = _tv("P_max_atm")
        st.session_state["add_bottom_dish"] = _ts("bottom_dish")
        st.session_state["add_top_dish"] = _ts("top_dish")
        st.session_state["add_shell_material"] = _ts("shell_material")
        st.session_state["add_lining"] = _tsel("lining", _YES_NO_OPTIONS)
        st.session_state["add_lining_material"] = _ts("lining_material")
        st.session_state["add_baffles"] = _tsel("baffles", _BAFFLE_OPTIONS)
        st.session_state["add_knuckle"] = _tv("knuckle_radius_m")
        st.session_state["add_rpm_min"] = _tv("N_rpm_min")
        st.session_state["add_rpm_max"] = _tv("N_rpm_max", 400.0)
        st.session_state["add_rps"] = _tv("N_rps")
        _imp_count = max(1, min(3, int(_tv("impeller_count", 1))))
        st.session_state["add_reactor_imp_count"] = _imp_count
        st.session_state["imp1_d"] = _tv("D_imp_m", 0.05)
        st.session_state["imp1_type"] = _tsel("impeller_type", _IMP_TYPE_OPTIONS)
        st.session_state["imp1_flow"] = _tsel("impeller_flow", _IMP_FLOW_OPTIONS)
        st.session_state["imp1_model"] = _ts("impeller_model")
        st.session_state["imp1_np"] = _tv("Np")
        st.session_state["imp1_nq"] = _tv("Nq")
        st.session_state["imp1_clr"] = _tv("imp1_clearance_m")
        st.session_state["imp1_h"] = _tv("imp1_height_m")
        st.session_state["imp2_d"] = _tv("D_imp2_m")
        st.session_state["imp2_type"] = _tsel("impeller_type2", _IMP_TYPE_OPTIONS)
        st.session_state["imp2_flow"] = _tsel("impeller_flow2", _IMP_FLOW_OPTIONS)
        st.session_state["imp2_model"] = _ts("impeller_model2")
        st.session_state["imp2_np"] = _tv("Np2")
        st.session_state["imp2_clr"] = _tv("imp2_clearance_m")
        st.session_state["imp2_h"] = _tv("imp2_height_m")
        st.session_state["imp3_d"] = _tv("D_imp3_m")
        st.session_state["imp3_type"] = _tsel("impeller_type3", _IMP_TYPE_OPTIONS)
        st.session_state["imp3_flow"] = _tsel("impeller_flow3", _IMP_FLOW_OPTIONS)
        st.session_state["imp3_model"] = _ts("impeller_model3")
        st.session_state["imp3_np"] = _tv("Np3")
        st.session_state["imp3_clr"] = _tv("imp3_clearance_m")
        st.session_state["imp3_h"] = _tv("imp3_height_m")
        st.session_state["add_zwietering"] = _tv("Zwietering_S")
        st.session_state["add_gmb"] = _tv("GMB_z")
        st.session_state["add_notes"] = _ts("notes")
        st.session_state["add_instrumentation"] = _ts("instrumentation")
        st.session_state["add_discharge_location"] = _tsel("discharge_location", _DISCHARGE_OPTIONS)
        st.session_state["add_insulated"] = _tsel("insulated", _YES_NO_OPTIONS)
        st.session_state["add_gas_addition"] = _tsel("gas_addition", _GAS_ADDITION_OPTIONS)
        st.session_state["add_gas_feed_control"] = _tsel("gas_feed_control", _GAS_FEED_CONTROL_OPTIONS)
        st.session_state["add_no_ports"] = _tv("no_ports")
        st.session_state["add_motor_power_kW"] = _tv("motor_power_kW")
        st.session_state["add_aux_units"] = _tsel("aux_units", _AUX_UNITS_OPTIONS)
        st.session_state["add_cip"] = _tsel("cip", _YES_NO_OPTIONS)
        st.session_state["add_heating_cooling"] = _tsel("heating_cooling", _YES_NO_OPTIONS)
        st.session_state["add_heat_transfer_medium"] = _tsel("heat_transfer_medium", _HTM_OPTIONS)
        st.session_state["add_heat_exchanger"] = _tsel("heat_exchanger", _HEAT_EXCHANGER_OPTIONS)
        st.session_state["add_probes"] = _ts("probes")
        st.rerun()
    elif template_reactor == "(none)":
        st.session_state.pop("_applied_template", None)

    # Impeller count selector OUTSIDE the form so it triggers an immediate rerun
    impeller_count = st.number_input(
        "How many impellers?", min_value=1, max_value=3, value=1,
        key="add_reactor_imp_count",
    )

    with st.form("add_reactor_form"):

        # ── Vessel Identity ───────────────────────────────────────────────
        st.markdown("##### Vessel Identity")
        _auto_id = _next_reactor_id(st.session_state.reactor_db)
        st.text_input("Reactor ID (auto-generated)", value=_auto_id, disabled=True, key="_display_reactor_id")
        id1, id2, id3, id4 = st.columns(4)
        with id1:
            owner = st.text_input("Owner / site *", max_chars=80, key="add_owner")
            tag = st.text_input("Tag", max_chars=80, key="add_tag")
        with id2:
            location = st.text_input("Location", key="add_location")
            rtype = st.selectbox("Type", _TYPE_OPTIONS, key="add_type")
        with id3:
            manufacturer = st.text_input("Manufacturer", key="add_manufacturer")
            scale = st.selectbox("Scale", _SCALE_OPTIONS, key="add_scale")
        with id4:
            manufacturer_model = st.text_input("Manufacturer model *", key="add_manufacturer_model")


        # ── Geometry ──────────────────────────────────────────────────────
        st.markdown("##### Geometry")
        g1, g2, g3, g4 = st.columns(4)
        with g1:
            D_tank = st.number_input("Tank ID (m)", min_value=0.0, value=0.10, format="%.4f", key="add_D_tank")
            H = st.number_input("Height tan-tan (m)", min_value=0.0, value=0.13, format="%.4f", key="add_H")
            OD = st.number_input("Outside diameter (m)", min_value=0.0, value=0.0, format="%.4f", key="add_OD")
        with g2:
            wall_thickness = st.number_input("Wall thickness (mm)", min_value=0.0, value=0.0, format="%.2f", key="add_wall")
            knuckle_radius = st.number_input("Knuckle radius (m)", min_value=0.0, value=0.0, format="%.4f", key="add_knuckle")
            baffles = st.selectbox("Baffles", _BAFFLE_OPTIONS, key="add_baffles")
        with g3:
            bottom_dish = st.text_input("Bottom dish type", "", key="add_bottom_dish")
            top_dish = st.text_input("Top dish type", "", key="add_top_dish")
        with g4:
            shell_material = st.text_input("Shell material", "", key="add_shell_material")
            lining = st.selectbox("Lining?", _YES_NO_OPTIONS, key="add_lining")
            lining_material = st.text_input("Lining material", "", key="add_lining_material")

        # ── Operating Ranges ──────────────────────────────────────────────
        st.markdown("##### Operating Ranges")
        or1, or2, or3, or4, or5 = st.columns(5)
        with or1:
            N_rpm_min = st.number_input("RPM min", min_value=0.0, value=0.0, format="%.0f", key="add_rpm_min")
        with or2:
            N_rpm_max = st.number_input("RPM max", min_value=0.0, value=400.0, format="%.0f", key="add_rpm_max")
        with or3:
            N_rps = st.number_input("Default speed (rev/s)", min_value=0.0, value=0.0, format="%.2f",
                                    key="add_rps", help="Leave at 0 to auto-compute from RPM midpoint.")
        with or4:
            V_L_min = st.number_input("Volume min (L)", min_value=0.0, value=0.0, format="%.2f", key="add_V_min")
        with or5:
            V_L_max = st.number_input("Volume max (L)", min_value=0.0, value=1.0, format="%.2f", key="add_V_max")
        or6, or7, _ = st.columns(3)
        with or6:
            T_max_C = st.number_input("Temperature max (°C)", min_value=TEMP_MIN_C, max_value=TEMP_MAX_C, value=0.0, format="%.1f", key="add_T_max")
        with or7:
            P_max_atm = st.number_input("Pressure max (atm)", min_value=0.0, value=0.0, format="%.2f", key="add_P_max")

        # ── Impeller 1 (always shown) ────────────────────────────────────
        st.markdown("##### Impeller 1 (primary)")
        i1a, i1b, i1c, i1d = st.columns(4)
        with i1a:
            D_imp = st.number_input("Diameter (m)", min_value=0.0, value=0.05, format="%.4f", key="imp1_d")
            impeller_type = st.selectbox("Type", _IMP_TYPE_OPTIONS, key="imp1_type")
            impeller_flow = st.selectbox("Flow", _IMP_FLOW_OPTIONS, key="imp1_flow")
            impeller_model_val = st.text_input("Impeller model", "", key="imp1_model")
        with i1b:
            Np_val = st.number_input("Np", min_value=0.0, value=0.0, format="%.2f", key="imp1_np")
            Nq_val = st.number_input("Nq", min_value=0.0, value=0.0, format="%.2f", key="imp1_nq")
        with i1c:
            imp1_clearance = st.number_input("Clearance (m)", min_value=0.0, value=0.0, format="%.4f", key="imp1_clr")
        with i1d:
            imp1_height = st.number_input("Height (m)", min_value=0.0, value=0.0, format="%.4f", key="imp1_h")

        # ── Impeller 2 (shown when count ≥ 2) ───────────────────────────
        D_imp2 = Np2_val = imp2_clearance = imp2_height = 0.0
        impeller_type2 = impeller_flow2 = impeller_model2 = ""
        if impeller_count >= 2:
            st.markdown("##### Impeller 2")
            i2a, i2b, i2c, i2d = st.columns(4)
            with i2a:
                D_imp2 = st.number_input("Diameter (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_d")
                impeller_type2 = st.selectbox("Type", _IMP_TYPE_OPTIONS, key="imp2_type")
                impeller_flow2 = st.selectbox("Flow", _IMP_FLOW_OPTIONS, key="imp2_flow")
                impeller_model2 = st.text_input("Impeller model", "", key="imp2_model")
            with i2b:
                Np2_val = st.number_input("Np", min_value=0.0, value=0.0, format="%.2f", key="imp2_np")
            with i2c:
                imp2_clearance = st.number_input("Clearance (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_clr")
            with i2d:
                imp2_height = st.number_input("Height (m)", min_value=0.0, value=0.0, format="%.4f", key="imp2_h")

        # ── Impeller 3 (shown when count ≥ 3) ───────────────────────────
        D_imp3 = Np3_val = imp3_clearance = imp3_height = 0.0
        impeller_type3 = impeller_flow3 = impeller_model3 = ""
        if impeller_count >= 3:
            st.markdown("##### Impeller 3")
            i3a, i3b, i3c, i3d = st.columns(4)
            with i3a:
                D_imp3 = st.number_input("Diameter (m)", min_value=0.0, value=0.0, format="%.4f", key="imp3_d")
                impeller_type3 = st.selectbox("Type", _IMP_TYPE_OPTIONS, key="imp3_type")
                impeller_flow3 = st.selectbox("Flow", _IMP_FLOW_OPTIONS, key="imp3_flow")
                impeller_model3 = st.text_input("Impeller model", "", key="imp3_model")
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
            zwietering_s = st.number_input("Zwietering S parameter", min_value=0.0, value=0.0, format="%.2f", key="add_zwietering")
        with s2:
            gmb_z = st.number_input("GMB z parameter", min_value=0.0, value=0.0, format="%.4f", key="add_gmb")

        # ── Row 7: Equipment & Process ────────────────────────────────────
        st.markdown("##### Equipment & Process")
        e1, e2, e3, e4 = st.columns(4)
        with e1:
            instrumentation = st.text_input("Instrumentation", "", key="add_instrumentation")
            discharge_location = st.selectbox("Discharge location", _DISCHARGE_OPTIONS, key="add_discharge_location")
            insulated = st.selectbox("Insulated?", _YES_NO_OPTIONS, key="add_insulated")
        with e2:
            gas_addition = st.selectbox("Gas addition type", _GAS_ADDITION_OPTIONS, key="add_gas_addition")
            gas_feed_control = st.selectbox("Gas feed control", _GAS_FEED_CONTROL_OPTIONS, key="add_gas_feed_control")
            aux_units = st.selectbox("Auxiliary units", _AUX_UNITS_OPTIONS, key="add_aux_units")
        with e3:
            no_ports = st.number_input("Number of ports", min_value=0.0, value=0.0, format="%.0f", key="add_no_ports")
            motor_power_kW = st.number_input("Motor power (kW)", min_value=0.0, value=0.0, format="%.2f", key="add_motor_power_kW")
            heat_exchanger = st.selectbox("Heat exchanger", _HEAT_EXCHANGER_OPTIONS, key="add_heat_exchanger")
        with e4:
            cip = st.selectbox("CIP?", _YES_NO_OPTIONS, key="add_cip")
            heating_cooling = st.selectbox("Heating/Cooling?", _YES_NO_OPTIONS, key="add_heating_cooling")
            heat_transfer_medium = st.selectbox("Heat transfer medium", _HTM_OPTIONS, key="add_heat_transfer_medium")
            
        probes = st.text_input("Probes (comma-separated)", key="add_probes",
                               help="e.g. Tr, EasyViewer, FBRM")
        notes = st.text_area("Notes", key="add_notes")
        submitted = st.form_submit_button("Add reactor")

        def _or_nan(v):
            """Return NaN for zero / empty values."""
            return v if v else np.nan

        # Auto-generate reactor_name from owner + tag
        name = f"{owner} \u2013 {tag}".strip(" \u2013 ") if (owner or tag) else ""

        # ── Cross-field validation ──────────────────────────────────────
        _reactor_errors = []
        if submitted and owner and manufacturer_model:
            if D_tank > 0 and D_imp > 0 and D_imp >= D_tank:
                _reactor_errors.append(
                    f"Impeller diameter ({D_imp:g} m) must be smaller than tank ID ({D_tank:g} m).")
            if H > 0 and imp1_clearance > 0 and imp1_clearance >= H:
                _reactor_errors.append(
                    f"Impeller 1 clearance ({imp1_clearance:g} m) must be less than height tan-tan ({H:g} m).")
            if N_rpm_min > 0 and N_rpm_max > 0 and N_rpm_min > N_rpm_max:
                _reactor_errors.append(
                    f"RPM min ({N_rpm_min:g}) must be ≤ RPM max ({N_rpm_max:g}).")
            if V_L_min > 0 and V_L_max > 0 and V_L_min > V_L_max:
                _reactor_errors.append(
                    f"Volume min ({V_L_min:g} L) must be ≤ volume max ({V_L_max:g} L).")
            if name_exists(st.session_state.reactor_db.get("reactor_name", []), name):
                _reactor_errors.append(f"A reactor named “{name}” already exists.")

        if submitted and owner and manufacturer_model and not _reactor_errors:
            new_row = pd.DataFrame([{
                "reactor_id": _auto_id,
                "reactor_name": name,
                "owner": owner,
                "tag": tag,
                "location": location,
                "manufacturer": manufacturer,
                "manufacturer_model": manufacturer_model,
                "type": rtype,
                "scale": scale,
                "D_tank_m": _or_nan(D_tank),
                "H_m": _or_nan(H),
                "D_imp_m": _or_nan(D_imp),
                "impeller_type": impeller_type,
                "impeller_flow": impeller_flow,
                "impeller_model": impeller_model_val,
                "Np": _or_nan(Np_val),
                "Nq": _or_nan(Nq_val),
                "N_rpm_min": _or_nan(N_rpm_min),
                "N_rpm_max": _or_nan(N_rpm_max),
                "N_rps": N_rps if N_rps > 0 else ((N_rpm_min + N_rpm_max) / 2 / 60 if N_rpm_max > 0 else np.nan),
                "V_L_min": _or_nan(V_L_min),
                "V_L_max": _or_nan(V_L_max),
                "V_L": V_L_max if V_L_max > 0 else np.nan,
                "T_max_C": _or_nan(T_max_C),
                "P_max_atm": _or_nan(P_max_atm),
                "shell_material": shell_material,
                "lining": lining,
                "lining_material": lining_material,
                "baffles": baffles,
                "bottom_dish": bottom_dish,
                "top_dish": top_dish,
                "impeller_count": impeller_count,
                "imp1_clearance_m": _or_nan(imp1_clearance),
                "imp1_height_m": _or_nan(imp1_height),
                "D_imp2_m": _or_nan(D_imp2),
                "impeller_type2": impeller_type2,
                "impeller_flow2": impeller_flow2,
                "impeller_model2": impeller_model2,
                "Np2": _or_nan(Np2_val),
                "imp2_clearance_m": _or_nan(imp2_clearance),
                "imp2_height_m": _or_nan(imp2_height),
                "D_imp3_m": _or_nan(D_imp3),
                "impeller_type3": impeller_type3,
                "impeller_flow3": impeller_flow3,
                "impeller_model3": impeller_model3,
                "Np3": _or_nan(Np3_val),
                "imp3_clearance_m": _or_nan(imp3_clearance),
                "imp3_height_m": _or_nan(imp3_height),
                "Zwietering_S": _or_nan(zwietering_s),
                "GMB_z": _or_nan(gmb_z),
                "wall_thickness_mm": _or_nan(wall_thickness),
                "OD_m": _or_nan(OD),
                "knuckle_radius_m": _or_nan(knuckle_radius),
                "notes": notes,
                "instrumentation": instrumentation,
                "discharge_location": discharge_location,
                "insulated": insulated,
                "gas_addition": gas_addition,
                "gas_feed_control": gas_feed_control,
                "no_ports": _or_nan(no_ports),
                "motor_power_kW": _or_nan(motor_power_kW),
                "aux_units": aux_units,
                "cip": cip,
                "heating_cooling": heating_cooling,
                "heat_transfer_medium": heat_transfer_medium,
                "heat_exchanger": heat_exchanger,
                "probes": probes,
            }])
            st.session_state.reactor_db = pd.concat(
                [st.session_state.reactor_db, new_row], ignore_index=True)
            _save_reactors(st.session_state.reactor_db)
            st.success(f"Added **{name}** to the vessel database.")
        elif submitted and _reactor_errors:
            for _e in _reactor_errors:
                st.error(_e)
        elif submitted:
            st.warning("Please enter an owner and a manufacturer model.")

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.markdown("### Export")
    st.download_button(
        "⬇️ Download vessel database (CSV)",
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
                    result_df = new_df
                else:
                    result_df = pd.concat(
                        [st.session_state.reactor_db, new_df], ignore_index=True)

                # Assign unique IDs to any rows that are missing one
                result_df, _newly_assigned = _assign_missing_reactor_ids(result_df)
                if _newly_assigned:
                    _names_list = "\n".join(
                        f"- **{name}** → `{rid}`" for name, rid in _newly_assigned
                    )
                    st.warning(
                        f"⚠️ {len(_newly_assigned)} reactor(s) were missing a "
                        f"`reactor_id` and have been assigned new IDs:\n\n{_names_list}"
                    )

                st.session_state.reactor_db = result_df
                _save_reactors(st.session_state.reactor_db)
                st.success("Vessel database updated.")
        except Exception as e:
            st.error(f"Error reading CSV: {e}")
