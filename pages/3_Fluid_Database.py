"""Page 3 – Fluid Database: browse solvents, explore T-dependent properties, manage custom fluids."""

import streamlit as st

import pandas as pd
import numpy as np
import pathlib

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, list_solvents, solvent_info_table,
    density, viscosity, surface_tension, diffusivity,
    is_known_solvent,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
FLUID_CSV = DATA_DIR / "fluids.csv"


# ── Custom-fluid persistence ─────────────────────────────────────────────

def _load_custom_fluids() -> pd.DataFrame:
    if FLUID_CSV.exists():
        return pd.read_csv(FLUID_CSV)
    return pd.DataFrame(columns=[
        "fluid_name", "rho_kg_m3", "mu_Pa_s", "D_mol_m2_s",
        "surface_tension_N_m", "notes",
    ])


def _save_custom_fluids(df: pd.DataFrame):
    df.to_csv(FLUID_CSV, index=False)


if "fluid_db" not in st.session_state:
    st.session_state.fluid_db = _load_custom_fluids()

st.title("💧 Fluid Database")

tab_library, tab_solvent, tab_custom, tab_blend, tab_import = st.tabs([
    "Solvent Library", "Solvent Properties (T)", "Custom Fluids", "Blend Fluids", "Import / Export",
])

# ── Solvent Library ───────────────────────────────────────────────────────
with tab_library:
    st.markdown(
        "Reference table of all built-in solvents with properties at 25 °C.  "
        "These solvents are always available in the Mixing Assessment and "
        "Reactor Comparison pages — select one and set any temperature to "
        "get properties from literature correlations."
    )

    st.dataframe(
        pd.DataFrame(solvent_info_table()),
        use_container_width=False,
        hide_index=True,
    )

    if not st.session_state.fluid_db.empty:
        st.subheader("Custom Fluids")
        st.caption("Fluids added manually (fixed properties, not temperature-dependent).")
        st.dataframe(st.session_state.fluid_db, use_container_width=False, hide_index=True)

# ── Solvent Properties at Temperature ─────────────────────────────────────
with tab_solvent:
    st.markdown(
        "Compute physical properties for common pharmaceutical solvents at any "
        "liquid-phase temperature.  Correlations use literature-fitted density, "
        "Arrhenius viscosity, linear surface tension, and Stokes-Einstein diffusivity."
    )

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
    ], vertical_spacing=0.2, horizontal_spacing=0.10)

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

# ── Custom Fluids ─────────────────────────────────────────────────────────
with tab_custom:
    st.markdown(
        "Add or edit **custom fluids** that are not in the built-in solvent library "
        "(e.g. mixtures, slurries, concentrated acids).  Custom fluids have fixed "
        "properties that do not vary with temperature."
    )

    # --- Browse & edit ---
    st.subheader("Browse & Edit")
    edited = st.data_editor(
        st.session_state.fluid_db,
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
        st.session_state.fluid_db = edited.copy()
        _save_custom_fluids(st.session_state.fluid_db)
        st.success("Custom fluid database saved.")

    # --- Add custom fluid ---
    st.subheader("Add Custom Fluid")
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
            if is_known_solvent(name):
                st.warning(f"**{name}** is already in the solvent library — no need to add it.")
            else:
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
                _save_custom_fluids(st.session_state.fluid_db)
                st.success(f"Added **{name}**.")
        elif submitted:
            st.warning("Enter a fluid name.")

# ── Blend Fluids ──────────────────────────────────────────────────────────
with tab_blend:
    st.markdown(
        "Create a **blend** from existing fluids (solvents and/or custom fluids).  "
        "Enter volumetric proportions — these are converted to mass fractions "
        "using each component's density, and mixture properties are computed as "
        "mass-weighted averages."
    )

    # Helper: get properties for any fluid name (solvent or custom) at 25 °C
    def _get_fluid_props(fname: str, T: float = 25.0) -> dict | None:
        """Return {rho, mu, D_mol, sigma} for a solvent or custom fluid."""
        if is_known_solvent(fname):
            p = get_properties(fname, T)
            return {
                "rho_kg_m3": p["rho_kg_m3"],
                "mu_Pa_s": p["mu_Pa_s"],
                "D_mol_m2_s": p["D_mol_m2_s"],
                "surface_tension_N_m": p["surface_tension_N_m"],
            }
        cust = st.session_state.fluid_db
        if not cust.empty and fname in cust["fluid_name"].values:
            row = cust[cust["fluid_name"] == fname].iloc[0]
            return {
                "rho_kg_m3": float(row["rho_kg_m3"]),
                "mu_Pa_s": float(row["mu_Pa_s"]),
                "D_mol_m2_s": float(row["D_mol_m2_s"]),
                "surface_tension_N_m": float(row["surface_tension_N_m"]),
            }
        return None

    # Available fluid names: solvents + custom
    _solvent_list = sorted(SOLVENT_DB.keys())
    _custom_list = (
        st.session_state.fluid_db["fluid_name"].tolist()
        if not st.session_state.fluid_db.empty else []
    )
    _available = _solvent_list + _custom_list

    blend_components = st.multiselect(
        "Select component fluids", _available,
        default=None, key="blend_components",
        help="Choose two or more fluids to blend.",
    )

    blend_T = st.number_input(
        "Temperature for solvent properties (°C)", value=25.0, step=1.0,
        format="%.1f", key="blend_T",
        help="Built-in solvent properties are evaluated at this temperature.  "
             "Custom fluids use their fixed values regardless.",
    )

    if len(blend_components) >= 2:
        # Collect volumetric proportions
        st.subheader("Volumetric Proportions")
        st.caption(
            "Enter the volume fraction of each component (they will be "
            "normalised automatically so they sum to 1)."
        )

        vol_fracs: dict[str, float] = {}
        cols = st.columns(min(len(blend_components), 4))
        for i, comp in enumerate(blend_components):
            with cols[i % len(cols)]:
                vol_fracs[comp] = st.number_input(
                    f"{comp}", min_value=0.0, value=1.0, step=0.1,
                    format="%.2f", key=f"blend_vf_{comp}",
                )

        total_vol = sum(vol_fracs.values())
        if total_vol <= 0:
            st.warning("Total volume must be > 0.")
        else:
            # Normalise volume fractions
            vol_norm = {k: v / total_vol for k, v in vol_fracs.items()}

            # Gather component properties
            comp_props: list[dict] = []
            _missing = []
            for comp in blend_components:
                p = _get_fluid_props(comp, blend_T)
                if p is None:
                    _missing.append(comp)
                else:
                    comp_props.append({"name": comp, "vol_frac": vol_norm[comp], **p})
            if _missing:
                st.error(f"Could not find properties for: {', '.join(_missing)}")
            else:
                # Convert volumetric to mass fractions
                # mass_i = vol_frac_i * rho_i (per unit total volume)
                masses = [cp["vol_frac"] * cp["rho_kg_m3"] for cp in comp_props]
                total_mass = sum(masses)
                for cp, m in zip(comp_props, masses):
                    cp["mass_frac"] = m / total_mass

                # Mass-weighted average properties
                blend_rho = sum(cp["mass_frac"] * cp["rho_kg_m3"] for cp in comp_props)
                blend_mu = sum(cp["mass_frac"] * cp["mu_Pa_s"] for cp in comp_props)
                blend_D = sum(cp["mass_frac"] * cp["D_mol_m2_s"] for cp in comp_props)
                blend_sig = sum(cp["mass_frac"] * cp["surface_tension_N_m"] for cp in comp_props)

                # Display composition table
                st.subheader("Blend Composition")
                comp_table = pd.DataFrame([{
                    "Component": cp["name"],
                    "Vol %": f"{cp['vol_frac'] * 100:.1f}",
                    "Mass %": f"{cp['mass_frac'] * 100:.1f}",
                    "ρ (kg/m³)": f"{cp['rho_kg_m3']:.1f}",
                    "μ (Pa·s)": f"{cp['mu_Pa_s']:.6f}",
                } for cp in comp_props])
                st.dataframe(comp_table, use_container_width=False, hide_index=True)

                # Display blended properties
                st.subheader("Blended Properties (mass-weighted)")
                bm1, bm2, bm3, bm4 = st.columns(4)
                bm1.metric("ρ (kg/m³)", f"{blend_rho:.2f}")
                bm2.metric("μ (Pa·s)", f"{blend_mu:.6f}")
                bm3.metric("σ (N/m)", f"{blend_sig:.4f}")
                bm4.metric("D_mol (m²/s)", f"{blend_D:.3e}")

                # Auto-generate name from actual normalised volume fractions
                _parts = [
                    f"{cp['vol_frac'] * 100:.0f}% {cp['name']}"
                    for cp in comp_props
                ]
                _auto_name = " / ".join(_parts) + f" ({blend_T:.0f} °C)"
                _comp_notes = ", ".join(
                    f"{cp['name']} {cp['vol_frac']*100:.0f}vol%"
                    for cp in comp_props
                )

                # Push the computed name into session state whenever it changes
                if st.session_state.get("_blend_auto_name") != _auto_name:
                    st.session_state["_blend_auto_name"] = _auto_name
                    st.session_state["blend_name"] = _auto_name

                blend_name = st.text_input(
                    "Blend name for database",
                    key="blend_name",
                )

                if st.button("➕ Add blend to custom fluids", key="blend_add_btn"):
                    if not blend_name:
                        st.warning("Enter a name for the blend.")
                    elif is_known_solvent(blend_name):
                        st.warning(f"**{blend_name}** conflicts with a built-in solvent name.")
                    else:
                        new_row = pd.DataFrame([{
                            "fluid_name": blend_name,
                            "rho_kg_m3": round(blend_rho, 2),
                            "mu_Pa_s": round(blend_mu, 8),
                            "D_mol_m2_s": blend_D,
                            "surface_tension_N_m": round(blend_sig, 5),
                            "notes": f"Blend at {blend_T:.0f} °C: {_comp_notes}",
                        }])
                        st.session_state.fluid_db = pd.concat(
                            [st.session_state.fluid_db, new_row], ignore_index=True)
                        _save_custom_fluids(st.session_state.fluid_db)
                        st.success(f"Added **{blend_name}** to the custom fluid database.")
                        st.rerun()
    elif len(blend_components) == 1:
        st.info("Select at least two fluids to create a blend.")

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.markdown("Import / export **custom fluids** (CSV).  Built-in solvents are always available.")
    st.download_button(
        "⬇️ Download custom fluids (CSV)",
        data=st.session_state.fluid_db.to_csv(index=False).encode("utf-8"),
        file_name="custom_fluids_export.csv",
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
                _save_custom_fluids(st.session_state.fluid_db)
                st.success("Custom fluid database updated.")
        except Exception as e:
            st.error(f"Error: {e}")
