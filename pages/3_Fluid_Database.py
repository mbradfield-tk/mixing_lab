"""Page 3 – Fluid System Database: browse, define, import/export fluid properties."""

import streamlit as st

import pandas as pd
import numpy as np
import pathlib

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, list_solvents, solvent_info_table,
    density, viscosity, surface_tension, diffusivity,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
FLUID_CSV = DATA_DIR / "fluids.csv"


def _load_fluids() -> pd.DataFrame:
    if FLUID_CSV.exists():
        return pd.read_csv(FLUID_CSV)
    return pd.DataFrame(columns=[
        "fluid_name", "rho_kg_m3", "mu_Pa_s", "D_mol_m2_s",
        "surface_tension_N_m", "notes",
    ])


def _save_fluids(df: pd.DataFrame):
    df.to_csv(FLUID_CSV, index=False)


if "fluid_db" not in st.session_state:
    st.session_state.fluid_db = _load_fluids()

st.title("💧 Fluid System Database")

tab_browse, tab_add, tab_solvent, tab_import = st.tabs([
    "Browse & Edit", "Add Fluid", "Solvent Properties (T)", "Import / Export",
])

# ── Browse & Edit ─────────────────────────────────────────────────────────
with tab_browse:
    st.markdown("Edit fluid properties directly in the table.  Click **Save changes** to persist.")

    search_term = st.text_input("🔍 Search by fluid name", key="fluid_search",
                                placeholder="e.g. water, methanol, toluene…")

    df_display = st.session_state.fluid_db.copy()
    if search_term:
        mask = df_display["fluid_name"].str.contains(search_term, case=False, na=False)
        df_display = df_display[mask]

    edited = st.data_editor(
        df_display,
        num_rows="dynamic",
        use_container_width=False,
        column_config={
            "rho_kg_m3": st.column_config.NumberColumn("ρ (kg/m³)", format="%.1f"),
            "mu_Pa_s": st.column_config.NumberColumn("μ (Pa·s)", format="%.6f"),
            "D_mol_m2_s": st.column_config.NumberColumn("D_mol (m²/s)", format="%.2e"),
            "surface_tension_N_m": st.column_config.NumberColumn("σ (N/m)", format="%.4f"),
        },
        key="fluid_editor",
    )

    if st.button("💾 Save changes", key="save_fluid"):
        if search_term:
            # Merge edits back into the full dataframe
            full = st.session_state.fluid_db.copy()
            full.update(edited)
            st.session_state.fluid_db = full
        else:
            st.session_state.fluid_db = edited.copy()
        _save_fluids(st.session_state.fluid_db)
        st.success("Fluid database saved.")

# ── Add Fluid ─────────────────────────────────────────────────────────────
with tab_add:
    with st.form("add_fluid"):
        c1, c2 = st.columns(2)
        with c1:
            name = st.text_input("Fluid name *")
            rho = st.number_input("Density ρ (kg/m³)", min_value=1.0, value=997.0, format="%.1f")
            mu = st.number_input("Dynamic viscosity μ (Pa·s)", min_value=1e-6, value=0.00089, format="%.6f")
        with c2:
            D_mol = st.number_input("Molecular diffusivity D (m²/s)", min_value=1e-12, value=2.3e-9, format="%.2e")
            sigma = st.number_input("Surface tension σ (N/m)", min_value=0.0, value=0.072, format="%.4f")
            notes = st.text_input("Notes", "")
        submitted = st.form_submit_button("Add fluid")
        if submitted and name:
            new = pd.DataFrame([{
                "fluid_name": name,
                "rho_kg_m3": rho,
                "mu_Pa_s": mu,
                "D_mol_m2_s": D_mol,
                "surface_tension_N_m": sigma,
                "notes": notes,
            }])
            st.session_state.fluid_db = pd.concat(
                [st.session_state.fluid_db, new], ignore_index=True)
            _save_fluids(st.session_state.fluid_db)
            st.success(f"Added **{name}**.")
        elif submitted:
            st.warning("Enter a fluid name.")

# ── Solvent Properties at Temperature ─────────────────────────────────────
with tab_solvent:
    st.markdown(
        "Compute physical properties for common pharmaceutical solvents at any "
        "liquid-phase temperature.  Correlations use literature-fitted density, "
        "Arrhenius viscosity, linear surface tension, and Stokes-Einstein diffusivity."
    )

    # --- Reference table of available solvents ---
    with st.expander("Available solvents (properties at 25 °C)", expanded=False):
        st.dataframe(pd.DataFrame(solvent_info_table()), use_container_width=False, hide_index=True)

    # --- Solvent selector + temperature ---
    s_col1, s_col2 = st.columns(2)
    with s_col1:
        solvent_name = st.selectbox("Solvent", list_solvents(), key="solv_sel")
    with s_col2:
        sd = SOLVENT_DB[solvent_name]
        T_C = st.number_input(
            f"Temperature (°C)  [liquid range: {sd.mp_C:.1f} – {sd.bp_C:.1f}]",
            min_value=-200.0, max_value=400.0,
            value=25.0, step=1.0, format="%.1f",
            key="solv_temp",
        )

    props = get_properties(solvent_name, T_C)

    if not props["in_range"]:
        st.warning(
            f"⚠️ {T_C:.1f} °C is outside the liquid range "
            f"({sd.mp_C:.0f} – {sd.bp_C:.0f} °C) for {solvent_name}.  "
            f"Values are extrapolated and may be unreliable."
        )

    # --- Property results ---
    st.subheader(f"{solvent_name} at {T_C:.1f} °C")
    pc1, pc2, pc3, pc4 = st.columns(4)
    pc1.metric("ρ (kg/m³)", f"{props['rho_kg_m3']:.2f}")
    pc2.metric("μ (Pa·s)", f"{props['mu_Pa_s']:.6f}")
    pc3.metric("σ (N/m)", f"{props['surface_tension_N_m']:.4f}")
    pc4.metric("D_mol (m²/s)", f"{props['D_mol_m2_s']:.3e}")

    st.caption(f"MW = {props['mw']:.2f} g/mol  ·  CAS {props['cas']}  ·  "
               f"b.p. = {props['bp_C']:.1f} °C  ·  m.p. = {props['mp_C']:.1f} °C")

    # --- Property-vs-temperature curves ---
    st.subheader("Property vs Temperature")
    T_lo = sd.mp_C
    T_hi = sd.bp_C
    T_arr = np.linspace(T_lo, T_hi, 200)

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=2, cols=2, subplot_titles=[
        "Density ρ (kg/m³)", "Viscosity μ (Pa·s)",
        "Surface tension σ (N/m)", "Diffusivity D (m²/s)",
    ], vertical_spacing=0.12, horizontal_spacing=0.10)

    rho_arr = [density(T, sd) for T in T_arr]
    mu_arr = [viscosity(T, sd) for T in T_arr]
    sig_arr = [surface_tension(T, sd) for T in T_arr]
    D_arr = [diffusivity(T, sd) for T in T_arr]

    for r, c, y_arr, name in [
        (1, 1, rho_arr, "ρ"), (1, 2, mu_arr, "μ"),
        (2, 1, sig_arr, "σ"), (2, 2, D_arr, "D"),
    ]:
        fig.add_trace(go.Scatter(
            x=T_arr, y=y_arr, mode="lines",
            line=dict(width=2), name=name, showlegend=False,
        ), row=r, col=c)
        # marker at selected temperature
        idx = np.argmin(np.abs(T_arr - T_C))
        fig.add_trace(go.Scatter(
            x=[T_C], y=[y_arr[idx]], mode="markers",
            marker=dict(size=10, color="red", symbol="circle"),
            name=f"{T_C:.0f} °C", showlegend=False,
        ), row=r, col=c)
        fig.update_xaxes(title_text="T (°C)", row=r, col=c)

    fig.update_layout(height=550, margin=dict(t=40, b=40))
    st.plotly_chart(fig, use_container_width=False)

    # --- Add to fluid database ---
    st.subheader("Add to Fluid Database")
    auto_name = f"{solvent_name} ({T_C:.0f} °C)"
    # Update the default name whenever solvent or temperature changes
    if st.session_state.get("_solv_prev_auto") != auto_name:
        st.session_state["solv_add_name"] = auto_name
        st.session_state["_solv_prev_auto"] = auto_name
    fluid_name = st.text_input("Fluid name for database entry", key="solv_add_name")

    if st.button("➕ Add to fluid database", key="solv_add_btn"):
        new_row = pd.DataFrame([{
            "fluid_name": fluid_name,
            "rho_kg_m3": round(props["rho_kg_m3"], 2),
            "mu_Pa_s": round(props["mu_Pa_s"], 8),
            "D_mol_m2_s": props["D_mol_m2_s"],
            "surface_tension_N_m": round(props["surface_tension_N_m"], 5),
            "notes": f"Auto-generated from {solvent_name} at {T_C:.1f} °C",
        }])
        st.session_state.fluid_db = pd.concat(
            [st.session_state.fluid_db, new_row], ignore_index=True)
        _save_fluids(st.session_state.fluid_db)
        st.success(f"Added **{fluid_name}** to the fluid database.")
        st.rerun()

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.download_button(
        "⬇️ Download fluid database (CSV)",
        data=st.session_state.fluid_db.to_csv(index=False).encode("utf-8"),
        file_name="fluids_export.csv",
        mime="text/csv",
    )
    uploaded = st.file_uploader("Upload CSV", type=["csv"], key="fluid_upload")
    if uploaded:
        try:
            new_df = pd.read_csv(uploaded)
            st.dataframe(new_df.head())
            mode = st.radio("Import mode", ["Replace", "Append"], key="fluid_import_mode")
            if st.button("Confirm import", key="fluid_import_confirm"):
                if mode == "Replace":
                    st.session_state.fluid_db = new_df
                else:
                    st.session_state.fluid_db = pd.concat(
                        [st.session_state.fluid_db, new_df], ignore_index=True)
                _save_fluids(st.session_state.fluid_db)
                st.success("Fluid database updated.")
        except Exception as e:
            st.error(f"Error: {e}")
