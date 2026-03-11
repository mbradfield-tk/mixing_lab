"""Page 4 – Particle Database: browse, define, import/export particle properties."""

import streamlit as st

import pandas as pd
import pathlib

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
PARTICLE_CSV = DATA_DIR / "particles.csv"

COLUMNS = [
    "particle_name", "rho_p_kg_m3", "d10_um", "d50_um", "d90_um",
    "shape_description", "shape_factor", "notes",
]

SHAPE_OPTIONS = [
    "spherical", "equant", "plate", "needle", "fibrous",
    "irregular", "rhombohedral", "dendritic",
]


def _load_particles() -> pd.DataFrame:
    if PARTICLE_CSV.exists():
        return pd.read_csv(PARTICLE_CSV)
    return pd.DataFrame(columns=COLUMNS)


def _save_particles(df: pd.DataFrame):
    df.to_csv(PARTICLE_CSV, index=False)


if "particle_db" not in st.session_state:
    st.session_state.particle_db = _load_particles()

st.title("❉ Particle Database")

tab_browse, tab_add, tab_import = st.tabs([
    "Browse & Edit", "Add Particle", "Import / Export",
])

# ── Browse & Edit ─────────────────────────────────────────────────────────
with tab_browse:
    st.markdown("Edit particle properties directly in the table.  Click **Save changes** to persist.")

    search_term = st.text_input("🔍 Search by particle name", key="particle_search",
                                placeholder="e.g. glass, API, lactose…")

    df_display = st.session_state.particle_db.copy()
    if search_term:
        mask = df_display["particle_name"].str.contains(search_term, case=False, na=False)
        df_display = df_display[mask]

    edited = st.data_editor(
        df_display,
        num_rows="dynamic",
        width="stretch",
        column_config={
            "rho_p_kg_m3": st.column_config.NumberColumn("ρ_p (kg/m³)", format="%.0f"),
            "d10_um": st.column_config.NumberColumn("D10 (µm)", format="%.1f"),
            "d50_um": st.column_config.NumberColumn("D50 (µm)", format="%.1f"),
            "d90_um": st.column_config.NumberColumn("D90 (µm)", format="%.1f"),
            "shape_description": st.column_config.SelectboxColumn(
                "Shape", options=SHAPE_OPTIONS, width="medium",
            ),
            "shape_factor": st.column_config.NumberColumn("φ (shape factor)", format="%.2f",
                                                           help="Sphericity 0–1, 1.0 = perfect sphere"),
        },
        key="particle_editor",
    )

    if st.button("💾 Save changes", key="save_particle"):
        if search_term:
            full = st.session_state.particle_db.copy()
            full.update(edited)
            st.session_state.particle_db = full
        else:
            st.session_state.particle_db = edited.copy()
        _save_particles(st.session_state.particle_db)
        st.success("Particle database saved.")

# ── Add Particle ──────────────────────────────────────────────────────────
with tab_add:
    with st.form("add_particle"):
        c1, c2 = st.columns(2)
        with c1:
            name = st.text_input("Particle name *")
            rho_p = st.number_input("Particle density ρ_p (kg/m³)", min_value=1.0,
                                     value=1400.0, format="%.0f")
            d10 = st.number_input("D10 (µm)", min_value=0.1, value=20.0, format="%.1f")
            d50 = st.number_input("D50 (µm)", min_value=0.1, value=75.0, format="%.1f")
        with c2:
            d90 = st.number_input("D90 (µm)", min_value=0.1, value=150.0, format="%.1f")
            shape = st.selectbox("Shape description", SHAPE_OPTIONS, index=0)
            phi = st.number_input("Shape factor φ (sphericity)", min_value=0.01,
                                   max_value=1.0, value=0.85, format="%.2f",
                                   help="1.0 = sphere, lower = more irregular")
            notes = st.text_input("Notes", "")
        submitted = st.form_submit_button("Add particle")
        if submitted and name:
            new = pd.DataFrame([{
                "particle_name": name,
                "rho_p_kg_m3": rho_p,
                "d10_um": d10,
                "d50_um": d50,
                "d90_um": d90,
                "shape_description": shape,
                "shape_factor": phi,
                "notes": notes,
            }])
            st.session_state.particle_db = pd.concat(
                [st.session_state.particle_db, new], ignore_index=True)
            _save_particles(st.session_state.particle_db)
            st.success(f"Added **{name}**.")
        elif submitted:
            st.warning("Enter a particle name.")

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.download_button(
        "⬇️ Download particle database (CSV)",
        data=st.session_state.particle_db.to_csv(index=False).encode("utf-8"),
        file_name="particles_export.csv",
        mime="text/csv",
    )
    uploaded = st.file_uploader("Upload CSV", type=["csv"], key="particle_upload")
    if uploaded:
        try:
            new_df = pd.read_csv(uploaded)
            st.dataframe(new_df.head())
            mode = st.radio("Import mode", ["Replace", "Append"], key="particle_import_mode")
            if st.button("Confirm import", key="particle_import_confirm"):
                if mode == "Replace":
                    st.session_state.particle_db = new_df
                else:
                    st.session_state.particle_db = pd.concat(
                        [st.session_state.particle_db, new_df], ignore_index=True)
                _save_particles(st.session_state.particle_db)
                st.success("Particle database updated.")
        except Exception as e:
            st.error(f"Error: {e}")
