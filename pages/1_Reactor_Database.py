"""Page 1 – Reactor Database: browse, add, edit, import/export reactor geometries."""

import streamlit as st

import pandas as pd
import pathlib
import numpy as np
from PIL import Image

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

    edited_df = st.data_editor(
        df_display,
        num_rows="dynamic",
        use_container_width=False,
        key="reactor_editor",
    )

    col_s1, col_s2 = st.columns([1, 5])
    with col_s1:
        if st.button("💾 Save changes", key="save_reactor_browse"):
            st.session_state.reactor_db = edited_df.copy()
            _save_reactors(st.session_state.reactor_db)
            st.success("Reactor database saved.")

    # ── Image Viewer ────────────────────────────────────────────────────────
    st.divider()
    st.subheader("🖼️ Reactor Images")

    reactor_names = df_display["reactor_name"].dropna().tolist()
    if reactor_names:
        selected_reactor = st.selectbox(
            "Select a reactor to view images",
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

            if not reactor_imgs and not cfd_imgs:
                st.info(
                    f"No images found for **{selected_reactor}**. "
                    "Place files in `images/reactors/` or `images/CFD/` "
                    "named `{Owner}_{Reactor}_*.png`."
                )

            # ── Reactor photos ────────────────────────────────────────────
            if reactor_imgs:
                st.markdown("#### Reactor Photos")
                cols = st.columns(min(len(reactor_imgs), 4))
                for idx, img_path in enumerate(reactor_imgs):
                    with cols[idx % len(cols)]:
                        # Caption from the part after the reactor name
                        label = img_path.stem.removeprefix(prefix).lstrip("_") or img_path.stem
                        st.image(str(img_path), caption=label, use_container_width=False)

            # ── CFD images ────────────────────────────────────────────────
            if cfd_imgs:
                st.markdown("#### CFD Results")

                # Group by type (EDR, velocity, …)
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
                                use_container_width=False,
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
            if st.button("Confirm import", key="reactor_import_confirm"):
                if mode.startswith("Replace"):
                    st.session_state.reactor_db = new_df
                else:
                    st.session_state.reactor_db = pd.concat(
                        [st.session_state.reactor_db, new_df], ignore_index=True)
                _save_reactors(st.session_state.reactor_db)
                st.success("Reactor database updated.")
        except Exception as e:
            st.error(f"Error reading CSV: {e}")
