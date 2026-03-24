"""Page 3 – Fluid Database: browse solvents, explore T-dependent properties, manage custom fluids."""

import streamlit as st

import pandas as pd
import numpy as np
import pathlib

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, list_solvents, solvent_info_table,
    density, viscosity, surface_tension, diffusivity,
    specific_heat, thermal_conductivity,
    is_known_solvent, vapor_pressure_mmHg, boiling_point_at_pressure,
    solvent_miscibility,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
FLUID_CSV = DATA_DIR / "fluids.csv"


# ── Custom-fluid persistence ─────────────────────────────────────────────

def _load_custom_fluids() -> pd.DataFrame:
    if FLUID_CSV.exists():
        return pd.read_csv(FLUID_CSV)
    return pd.DataFrame(columns=[
        "fluid_name", "rho_kg_m3", "mu_Pa_s", "D_mol_m2_s",
        "surface_tension_N_m", "Cp_J_per_kgK", "k_W_per_mK", "notes",
    ])


def _save_custom_fluids(df: pd.DataFrame):
    df.to_csv(FLUID_CSV, index=False)


if "fluid_db" not in st.session_state:
    st.session_state.fluid_db = _load_custom_fluids()

st.title("💧 Fluid Database")

_is_admin = st.session_state.get("admin_authenticated", False)
_ADMIN_HINT = "Log in via Admin Tools to enable editing."

tab_library, tab_solvent, tab_custom, tab_blend, tab_import = st.tabs([
    "Solvent Library", "Solvent Properties (T)", "Custom Fluids", "Blend Fluids", "Import / Export",
])

# ── Solvent Library ───────────────────────────────────────────────────────
with tab_library:
    st.markdown(
        "Reference table of all built-in solvents with **properties at 25 °C and 1 atm**.  "
        "These solvents are always available in the Mixing Assessment and "
        "Reactor Comparison pages — select one and set any temperature to "
        "get properties from literature correlations."
    )

    st.dataframe(
        pd.DataFrame(solvent_info_table()),
        width='content',
        hide_index=True,
    )

    if not st.session_state.fluid_db.empty:
        st.subheader("Custom Fluids")
        st.caption("Fluids added manually (fixed properties, not temperature-dependent).")
        st.dataframe(st.session_state.fluid_db, width='content', hide_index=True)

# ── Solvent Properties at Temperature ─────────────────────────────────────
with tab_solvent:
    st.markdown(
        "Compute physical properties for common pharmaceutical solvents at any "
        "liquid-phase temperature and pressure.  Correlations use literature-fitted density, "
        "Arrhenius viscosity, linear surface tension, and Stokes-Einstein diffusivity.  "
        "The **Antoine equation** adjusts the boiling point for non-atmospheric pressures."
    )

    # --- Solvent selector + pressure + temperature ---
    s_col1, s_col2, s_col3 = st.columns(3)
    with s_col1:
        solvent_name = st.selectbox("Solvent", list_solvents(), key="solv_sel")
    sd = SOLVENT_DB[solvent_name]
    with s_col2:
        P_atm = st.number_input(
            "Pressure (atm)", min_value=0.001, max_value=50.0,
            value=1.0, step=0.1, format="%.3f",
            key="solv_P",
            help="System pressure in atmospheres. Adjusts the boiling point "
                 "via the Antoine equation.",
        )
    bp_at_P = boiling_point_at_pressure(P_atm, sd)
    with s_col3:
        T_C = st.number_input(
            f"Temperature (°C)  [liquid range: {sd.mp_C:.1f} – {bp_at_P:.1f}]",
            min_value=-200.0, max_value=400.0,
            value=25.0, step=1.0, format="%.1f",
            key="solv_temp",
        )

    props = get_properties(solvent_name, T_C, P_atm)

    if not props["in_range"]:
        st.warning(
            f"⚠️ {T_C:.1f} °C is outside the liquid range "
            f"({sd.mp_C:.0f} – {bp_at_P:.0f} °C) for {solvent_name} "
            f"at {P_atm:.3f} atm.  Values are extrapolated and may be unreliable."
        )

    # --- Property results ---
    st.subheader(f"{solvent_name} at {T_C:.1f} °C, {P_atm:.3f} atm")
    pc1, pc2, pc3, pc4 = st.columns(4)
    pc1.metric("ρ (kg/m³)", f"{props['rho_kg_m3']:.2f}")
    pc2.metric("μ (Pa·s)", f"{props['mu_Pa_s']:.6f}")
    pc3.metric("σ (N/m)", f"{props['surface_tension_N_m']:.4f}")
    pc4.metric("D_mol (m²/s)", f"{props['D_mol_m2_s']:.3e}")

    pc5, pc6, pc7, pc8, pc9 = st.columns(5)
    pc5.metric("Cp (J/kg·K)", f"{props['Cp_J_per_kgK']:.1f}")
    pc6.metric("k (W/m·K)", f"{props['k_W_per_mK']:.4f}")
    pc7.metric("Vapor pressure (atm)", f"{props['vapor_pressure_atm']:.4f}")
    pc8.metric("b.p. at P (°C)", f"{bp_at_P:.1f}")
    pc9.metric("Normal b.p. (°C)", f"{sd.bp_C:.1f}")

    st.caption(f"MW = {props['mw']:.2f} g/mol  ·  CAS {props['cas']}  ·  "
               f"b.p.(1 atm) = {props['bp_C']:.1f} °C  ·  "
               f"b.p.({P_atm:.2f} atm) = {bp_at_P:.1f} °C  ·  "
               f"m.p. = {props['mp_C']:.1f} °C")

    # --- Property-vs-temperature curves ---
    st.subheader("Property vs Temperature")
    T_lo = sd.mp_C
    T_hi = bp_at_P
    T_arr = np.linspace(T_lo, T_hi, 200)

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=3, cols=2, subplot_titles=[
        "Density ρ (kg/m³)", "Viscosity μ (Pa·s)",
        "Surface tension σ (N/m)", "Diffusivity D (m²/s)",
        "Specific heat Cp (J/kg·K)", "Thermal conductivity k (W/m·K)",
    ], vertical_spacing=0.14, horizontal_spacing=0.10)

    rho_arr = [density(T, sd) for T in T_arr]
    mu_arr = [viscosity(T, sd) for T in T_arr]
    sig_arr = [surface_tension(T, sd) for T in T_arr]
    D_arr = [diffusivity(T, sd) for T in T_arr]
    Cp_arr = [specific_heat(T, sd) for T in T_arr]
    k_arr = [thermal_conductivity(T, sd) for T in T_arr]

    for r, c, y_arr, name in [
        (1, 1, rho_arr, "ρ"), (1, 2, mu_arr, "μ"),
        (2, 1, sig_arr, "σ"), (2, 2, D_arr, "D"),
        (3, 1, Cp_arr, "Cp"), (3, 2, k_arr, "k"),
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

    fig.update_layout(height=780, margin=dict(t=40, b=40))
    st.plotly_chart(fig, width='content')

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
        width='content',
        column_config={
            "rho_kg_m3": st.column_config.NumberColumn("ρ (kg/m³)", format="%.1f"),
            "mu_Pa_s": st.column_config.NumberColumn("μ (Pa·s)", format="%.6f"),
            "D_mol_m2_s": st.column_config.NumberColumn("D_mol (m²/s)", format="%.2e"),
            "surface_tension_N_m": st.column_config.NumberColumn("σ (N/m)", format="%.4f"),
            "Cp_J_per_kgK": st.column_config.NumberColumn("Cp (J/kg·K)", format="%.1f"),
            "k_W_per_mK": st.column_config.NumberColumn("k (W/m·K)", format="%.4f"),
            "hsp_d": st.column_config.NumberColumn("δd (MPa½)", format="%.1f"),
            "hsp_p": st.column_config.NumberColumn("δp (MPa½)", format="%.1f"),
            "hsp_h": st.column_config.NumberColumn("δh (MPa½)", format="%.1f"),
        },
        key="fluid_editor",
    )

    if st.button("💾 Save changes", key="save_fluid",
                 disabled=not _is_admin, help=None if _is_admin else _ADMIN_HINT):
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
            Cp = st.number_input("Specific heat Cp (J/kg·K)", min_value=1.0, value=4182.0, format="%.1f")
            k_val = st.number_input("Thermal conductivity k (W/m·K)", min_value=0.001, value=0.607, format="%.4f")
        st.markdown("**Hansen Solubility Parameters** _(optional — for miscibility screening)_")
        h1, h2, h3 = st.columns(3)
        hsp_d = h1.number_input("δd dispersion (MPa½)", min_value=0.0, value=0.0, format="%.1f",
                                help="Set to 0 if unknown")
        hsp_p = h2.number_input("δp polar (MPa½)", min_value=0.0, value=0.0, format="%.1f",
                                help="Set to 0 if unknown")
        hsp_h = h3.number_input("δh H-bonding (MPa½)", min_value=0.0, value=0.0, format="%.1f",
                                help="Set to 0 if unknown")
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
                    "Cp_J_per_kgK": Cp,
                    "k_W_per_mK": k_val,
                    "hsp_d": hsp_d,
                    "hsp_p": hsp_p,
                    "hsp_h": hsp_h,
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
        "Enter proportions on a **volumetric** or **mass** basis — they are converted "
        "to mass fractions using each component's density, and mixture properties "
        "are computed as mass-weighted averages."
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
                "Cp_J_per_kgK": p["Cp_J_per_kgK"],
                "k_W_per_mK": p["k_W_per_mK"],
            }
        cust = st.session_state.fluid_db
        if not cust.empty and fname in cust["fluid_name"].values:
            row = cust[cust["fluid_name"] == fname].iloc[0]
            return {
                "rho_kg_m3": float(row["rho_kg_m3"]),
                "mu_Pa_s": float(row["mu_Pa_s"]),
                "D_mol_m2_s": float(row["D_mol_m2_s"]),
                "surface_tension_N_m": float(row["surface_tension_N_m"]),
                "Cp_J_per_kgK": float(row.get("Cp_J_per_kgK", 4182.0)),
                "k_W_per_mK": float(row.get("k_W_per_mK", 0.607)),
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
        # Choose input basis
        _blend_basis = st.radio(
            "Input basis",
            ["Volume", "Mass"],
            horizontal=True,
            key="blend_basis",
            help="Enter component amounts on a volumetric or mass basis. "
                 "Properties are always averaged on a mass-fraction basis.",
        )
        _is_vol_basis = _blend_basis == "Volume"
        _basis_label = "volume" if _is_vol_basis else "mass"

        # Collect proportions
        st.subheader(f"Component Contributions ({_basis_label})")
        st.caption(
            f"Enter the {_basis_label} contribution of each component (they will be "
            "normalised automatically to sum to 1.0)."
        )

        input_fracs: dict[str, float] = {}
        cols = st.columns(min(len(blend_components), 4))
        for i, comp in enumerate(blend_components):
            with cols[i % len(cols)]:
                input_fracs[comp] = st.number_input(
                    f"{comp}", min_value=0.0, value=1.0, step=0.1,
                    format="%.2f", key=f"blend_f_{_basis_label}_{comp}",
                )

        total_input = sum(input_fracs.values())
        if total_input <= 0:
            st.warning(f"Total {_basis_label} must be > 0.")
        else:
            # Normalise input fractions
            input_norm = {k: v / total_input for k, v in input_fracs.items()}

            # Gather component properties
            comp_props: list[dict] = []
            _missing = []
            for comp in blend_components:
                p = _get_fluid_props(comp, blend_T)
                if p is None:
                    _missing.append(comp)
                else:
                    comp_props.append({"name": comp, "input_frac": input_norm[comp], **p})
            if _missing:
                st.error(f"Could not find properties for: {', '.join(_missing)}")
            else:
                if _is_vol_basis:
                    # Input is volume fractions → convert to mass fractions
                    for cp in comp_props:
                        cp["vol_frac"] = cp["input_frac"]
                    masses = [cp["vol_frac"] * cp["rho_kg_m3"] for cp in comp_props]
                    total_mass = sum(masses)
                    for cp, m in zip(comp_props, masses):
                        cp["mass_frac"] = m / total_mass
                else:
                    # Input is mass fractions → convert to volume fractions
                    for cp in comp_props:
                        cp["mass_frac"] = cp["input_frac"]
                    vols = [cp["mass_frac"] / cp["rho_kg_m3"] for cp in comp_props]
                    total_vol_calc = sum(vols)
                    for cp, v in zip(comp_props, vols):
                        cp["vol_frac"] = v / total_vol_calc

                # Mass-weighted average properties
                blend_rho = sum(cp["mass_frac"] * cp["rho_kg_m3"] for cp in comp_props)
                blend_mu = sum(cp["mass_frac"] * cp["mu_Pa_s"] for cp in comp_props)
                blend_D = sum(cp["mass_frac"] * cp["D_mol_m2_s"] for cp in comp_props)
                blend_sig = sum(cp["mass_frac"] * cp["surface_tension_N_m"] for cp in comp_props)
                blend_Cp = sum(cp["mass_frac"] * cp["Cp_J_per_kgK"] for cp in comp_props)
                blend_k = sum(cp["mass_frac"] * cp["k_W_per_mK"] for cp in comp_props)

                # Compute absolute volume (L) and mass (kg) per component
                for cp in comp_props:
                    if _is_vol_basis:
                        cp["vol_L"] = input_fracs[cp["name"]]
                        cp["mass_kg"] = cp["vol_L"] * cp["rho_kg_m3"] * 1e-3  # L→m³
                    else:
                        cp["mass_kg"] = input_fracs[cp["name"]]
                        cp["vol_L"] = cp["mass_kg"] / cp["rho_kg_m3"] * 1e3   # m³→L
                _total_vol_L = sum(cp["vol_L"] for cp in comp_props)
                _total_mass_kg = sum(cp["mass_kg"] for cp in comp_props)

                # Display total input
                st.metric(f"Total {_basis_label} entered", f"{total_input:.2f}")

                # Display composition table
                st.subheader("Blend Composition")
                blend_nu = blend_mu / blend_rho if blend_rho > 0 else 0.0
                comp_rows = [{
                    "Component": cp["name"],
                    "Added": f"{input_fracs[cp['name']]:.2f}",
                    "Volume (L)": f"{cp['vol_L']:.3g}",
                    "Mass (kg)": f"{cp['mass_kg']:.3g}",
                    "Vol %": f"{cp['vol_frac'] * 100:.1f}",
                    "Mass %": f"{cp['mass_frac'] * 100:.1f}",
                    "ρ (kg/m³)": f"{cp['rho_kg_m3']:.1f}",
                    "μ (Pa·s)": f"{cp['mu_Pa_s']:.6f}",
                    "ν (m²/s)": f"{cp['mu_Pa_s'] / cp['rho_kg_m3']:.3e}",
                    "σ (N/m)": f"{cp['surface_tension_N_m']:.4f}",
                    "D (m²/s)": f"{cp['D_mol_m2_s']:.3e}",
                    "Cp (J/kg·K)": f"{cp['Cp_J_per_kgK']:.1f}",
                    "k (W/m·K)": f"{cp['k_W_per_mK']:.4f}",
                } for cp in comp_props]
                comp_rows.append({
                    "Component": "**Blend**",
                    "Added": f"{total_input:.2f}",
                    "Volume (L)": f"{_total_vol_L:.3g}",
                    "Mass (kg)": f"{_total_mass_kg:.3g}",
                    "Vol %": "100.0",
                    "Mass %": "100.0",
                    "ρ (kg/m³)": f"{blend_rho:.1f}",
                    "μ (Pa·s)": f"{blend_mu:.6f}",
                    "ν (m²/s)": f"{blend_nu:.3e}",
                    "σ (N/m)": f"{blend_sig:.4f}",
                    "D (m²/s)": f"{blend_D:.3e}",
                    "Cp (J/kg·K)": f"{blend_Cp:.1f}",
                    "k (W/m·K)": f"{blend_k:.4f}",
                })
                st.dataframe(pd.DataFrame(comp_rows), width='content', hide_index=True)

                # ── Pairwise miscibility screening ────────────────────────
                from itertools import combinations
                _pairs = list(combinations([cp["name"] for cp in comp_props], 2))
                if _pairs:
                    st.subheader("Miscibility Screening")
                    st.caption(
                        "Pairwise assessment using experimental data (built-in solvents) "
                        "or Hansen Solubility Parameters where available."
                    )
                    _any_immiscible = False
                    _misc_rows = []
                    for _n1, _n2 in _pairs:
                        _misc = solvent_miscibility(_n1, _n2, custom_fluids=st.session_state.fluid_db)
                        _misc_rows.append({
                            "Pair": f"{_n1} / {_n2}",
                            "Assessment": _misc["assessment"],
                            "R_a (MPa½)": f"{_misc['Ra']:.1f}" if _misc["Ra"] is not None else "—",
                            "Source": _misc["source"],
                        })
                        if _misc["miscible"] is False:
                            _any_immiscible = True
                    st.dataframe(pd.DataFrame(_misc_rows), width='content', hide_index=True)

                    if _any_immiscible:
                        st.warning(
                            "⚠️ One or more component pairs are **immiscible or partially miscible**. "
                            "The blend may form multiple liquid phases. Mass-weighted property "
                            "averages shown below assume a single homogeneous phase and may "
                            "not be physically meaningful."
                        )
                    else:
                        st.success("🟢 All component pairs are miscible — single-phase blend expected.")

                # Display blended properties
                st.subheader("Blended Properties (mass-weighted)")
                bm1, bm2, bm3, bm4, bm5, bm6 = st.columns(6)
                bm1.metric("ρ (kg/m³)", f"{blend_rho:.2f}")
                bm2.metric("μ (Pa·s)", f"{blend_mu:.6f}")
                bm3.metric("σ (N/m)", f"{blend_sig:.4f}")
                bm4.metric("D_mol (m²/s)", f"{blend_D:.3e}")
                bm5.metric("Cp (J/kg·K)", f"{blend_Cp:.1f}")
                bm6.metric("k (W/m·K)", f"{blend_k:.4f}")

                # Auto-generate name from actual normalised fractions
                if _is_vol_basis:
                    _parts = [
                        f"{cp['vol_frac'] * 100:.0f}%v {cp['name']}"
                        for cp in comp_props
                    ]
                    _comp_notes = ", ".join(
                        f"{cp['name']} {cp['vol_frac']*100:.0f}vol%"
                        for cp in comp_props
                    )
                else:
                    _parts = [
                        f"{cp['mass_frac'] * 100:.0f}%w {cp['name']}"
                        for cp in comp_props
                    ]
                    _comp_notes = ", ".join(
                        f"{cp['name']} {cp['mass_frac']*100:.0f}wt%"
                        for cp in comp_props
                    )
                _auto_name = " / ".join(_parts) + f" ({blend_T:.0f} °C)"

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
                            "Cp_J_per_kgK": round(blend_Cp, 1),
                            "k_W_per_mK": round(blend_k, 4),
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
            if st.button("Confirm import", key="fluid_import_confirm",
                         disabled=not _is_admin, help=None if _is_admin else _ADMIN_HINT):
                if mode == "Replace":
                    st.session_state.fluid_db = new_df
                else:
                    st.session_state.fluid_db = pd.concat(
                        [st.session_state.fluid_db, new_df], ignore_index=True)
                _save_custom_fluids(st.session_state.fluid_db)
                st.success("Custom fluid database updated.")
        except Exception as e:
            st.error(f"Error: {e}")
