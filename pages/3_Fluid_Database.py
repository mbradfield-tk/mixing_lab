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
        "to mass fractions using each component's density.  Properties are computed "
        "with literature-recommended mixing rules (log-mixing for viscosity, "
        "volume-additive density, etc.).  Optionally add **dissolved starting "
        "material** to account for the effect of a soluble solid on the solution."
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
        help="Choose one or more fluids. Two or more will be blended.",
    )

    blend_T = st.number_input(
        "Temperature for solvent properties (°C)", value=25.0, step=1.0,
        format="%.1f", key="blend_T",
        help="Built-in solvent properties are evaluated at this temperature.  "
             "Custom fluids use their fixed values regardless.",
    )

    if len(blend_components) >= 1:
        # Choose input basis
        _blend_basis = st.radio(
            "Input basis",
            ["Volume", "Mass", "SM Ratio"],
            horizontal=True,
            key="blend_basis",
            help="**Volume / Mass** — enter amounts directly.  "
                 "**SM Ratio** — specify starting-material mass and "
                 "solvent-to-SM ratio (L solvent per kg SM); "
                 "total solvent volume is computed automatically.",
        )
        _is_sm_ratio = _blend_basis == "SM Ratio"
        _is_vol_basis = _blend_basis in ("Volume", "SM Ratio")
        _basis_label = ("relative volume" if _is_sm_ratio
                        else "volume" if _is_vol_basis else "mass")

        # ── SM-ratio sizing inputs ────────────────────────────────
        _smr_total_vol_L = 0.0
        _smr_mass_kg = 0.0
        if _is_sm_ratio:
            st.subheader("Starting-Material Sizing")
            st.caption(
                "Enter the mass of starting material and the solvent-to-SM "
                "ratio (sometimes called *relative volumes*).  The total "
                "solvent volume is SM mass × ratio."
            )
            _smr_c1, _smr_c2, _smr_c3 = st.columns(3)
            with _smr_c1:
                _smr_mass_kg = st.number_input(
                    "SM mass (kg)", min_value=0.001, value=0.010,
                    step=0.001, format="%.3f", key="smr_sm_mass",
                )
            with _smr_c2:
                _smr_ratio = st.number_input(
                    "Solvent ratio (L/kg SM)",
                    min_value=0.1, value=10.0, step=0.5,
                    format="%.1f", key="smr_ratio",
                    help="Total litres of solvent blend per kg of SM. "
                         "e.g. 10 L/kg = '10 volumes'.",
                )
            _smr_total_vol_L = _smr_mass_kg * _smr_ratio
            with _smr_c3:
                st.metric("Total solvent volume",
                          f"{_smr_total_vol_L:.3f} L")

        # Collect proportions
        st.subheader(f"Component Contributions ({_basis_label})")
        if _is_sm_ratio:
            st.caption(
                "Enter the relative volume proportion of each solvent "
                "component (they will be normalised to sum to 1.0). "
                "Absolute volumes are set by the SM ratio above."
            )
        else:
            st.caption(
                f"Enter the {_basis_label} contribution of each component "
                "(they will be normalised automatically to sum to 1.0)."
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

                # --- Mixing rules (literature-recommended) ---
                # Density: ideal mixing (volume additivity)
                #   1/ρ_mix = Σ(w_i / ρ_i)
                blend_rho = 1.0 / sum(
                    cp["mass_frac"] / cp["rho_kg_m3"]
                    for cp in comp_props)

                # Viscosity: Arrhenius log-mixing rule
                #   ln(μ_mix) = Σ(w_i · ln μ_i)
                # Ref: Irving (1977); Poling et al., Properties of Gases & Liquids
                blend_mu = np.exp(sum(
                    cp["mass_frac"] * np.log(cp["mu_Pa_s"])
                    for cp in comp_props))

                # Diffusivity: Vignes-type log mixing
                blend_D = np.exp(sum(
                    cp["mass_frac"] * np.log(cp["D_mol_m2_s"])
                    for cp in comp_props))

                # Surface tension: volume-fraction weighted (Macleod–Sugden)
                blend_sig = sum(
                    cp["vol_frac"] * cp["surface_tension_N_m"]
                    for cp in comp_props)

                # Specific heat: mass-weighted (thermodynamically exact)
                blend_Cp = sum(
                    cp["mass_frac"] * cp["Cp_J_per_kgK"]
                    for cp in comp_props)

                # Thermal conductivity: volume-fraction weighted (Li 1976)
                blend_k = sum(
                    cp["vol_frac"] * cp["k_W_per_mK"]
                    for cp in comp_props)

                # Dissolved-solid flags (set later if SM is included)
                _sm_added = False
                _sm_name_val = ""
                _sm_conc_kg_L = 0.0

                # Compute absolute volume (L) and mass (kg) per component
                for cp in comp_props:
                    if _is_sm_ratio:
                        # Distribute total volume from SM sizing
                        cp["vol_L"] = (cp["vol_frac"]
                                       * _smr_total_vol_L)
                        cp["mass_kg"] = (cp["vol_L"]
                                         * cp["rho_kg_m3"] * 1e-3)
                    elif _is_vol_basis:
                        cp["vol_L"] = input_fracs[cp["name"]]
                        cp["mass_kg"] = cp["vol_L"] * cp["rho_kg_m3"] * 1e-3  # L→m³
                    else:
                        cp["mass_kg"] = input_fracs[cp["name"]]
                        cp["vol_L"] = cp["mass_kg"] / cp["rho_kg_m3"] * 1e3   # m³→L
                _total_vol_L = sum(cp["vol_L"] for cp in comp_props)
                _total_mass_kg = sum(cp["mass_kg"] for cp in comp_props)

                # Display total input
                if _is_sm_ratio:
                    st.metric("Total solvent",
                              f"{_total_vol_L:.3f} L  "
                              f"({_total_mass_kg:.3f} kg)")
                else:
                    st.metric(f"Total {_basis_label} entered",
                              f"{total_input:.2f}")

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

                # ── Dissolved starting material (optional) ────────────
                st.subheader("Dissolved Starting Material")
                if _is_sm_ratio:
                    st.info(
                        f"SM Ratio mode: **{_smr_mass_kg:.3f} kg** of "
                        f"starting material is automatically included."
                    )
                    _include_sm = True
                else:
                    _include_sm = st.checkbox(
                        "Include dissolved starting material",
                        value=False,
                        key="blend_include_sm",
                        help="Add a soluble solid (e.g. API, reagent) that "
                             "dissolves in the solvent blend, adjusting "
                             "physical properties.",
                    )

                if _include_sm:
                    sm_c1, sm_c2, sm_c3 = st.columns(3)
                    with sm_c1:
                        _sm_name_input = st.text_input(
                            "Starting material name", value="API",
                            key="blend_sm_name",
                        )
                        if _is_sm_ratio:
                            _sm_mass_kg = _smr_mass_kg
                            st.metric("SM mass (kg)",
                                      f"{_sm_mass_kg:.3f}")
                        else:
                            _sm_mass_kg = st.number_input(
                                "Mass of dissolved solid (kg)",
                                min_value=0.0, value=0.010, step=0.001,
                                format="%.3f", key="blend_sm_mass",
                            )
                    with sm_c2:
                        _sm_rho = st.number_input(
                            "Solid density ρ_s (kg/m³)",
                            min_value=100.0, value=1300.0, step=10.0,
                            format="%.0f", key="blend_sm_rho",
                            help="True density of the solid "
                                 "(typ. 1100–1500 kg/m³ for organic APIs).",
                        )
                        _sm_mw = st.number_input(
                            "Molecular weight (g/mol, 0 = unknown)",
                            min_value=0.0, value=0.0, step=10.0,
                            format="%.1f", key="blend_sm_mw",
                        )
                    with sm_c3:
                        _sm_Cp = st.number_input(
                            "Solid Cp (J/kg·K)",
                            min_value=100.0, value=1200.0, step=50.0,
                            format="%.0f", key="blend_sm_Cp",
                            help="Specific heat of the solid. "
                                 "Typical: 1000–1500 J/kg·K for organic solids.",
                        )
                        _sm_k_visc = st.number_input(
                            "Viscosity coefficient k_μ",
                            min_value=0.0, value=2.5, step=0.1,
                            format="%.1f", key="blend_sm_kmu",
                            help="Controls viscosity increase from dissolved "
                                 "solid: μ_sol = μ_blend × exp(k_μ × φ_s). "
                                 "Default 2.5 (Einstein). Increase for "
                                 "large molecules / polymers.",
                        )

                    if _sm_mass_kg > 0:
                        _sm_vol_m3 = _sm_mass_kg / _sm_rho
                        _sm_vol_L = _sm_vol_m3 * 1000.0

                        _sol_mass_kg = _total_mass_kg + _sm_mass_kg
                        _sol_vol_L = _total_vol_L + _sm_vol_L
                        _sol_vol_m3 = _sol_vol_L / 1000.0

                        _phi_s = (_sm_vol_m3 / _sol_vol_m3
                                  if _sol_vol_m3 > 0 else 0.0)
                        _w_s = (_sm_mass_kg / _sol_mass_kg
                                if _sol_mass_kg > 0 else 0.0)

                        # Adjust properties for dissolved solid
                        blend_rho = (_sol_mass_kg / _sol_vol_m3
                                     if _sol_vol_m3 > 0 else blend_rho)
                        _mu_blend_orig = blend_mu
                        blend_mu = _mu_blend_orig * np.exp(
                            _sm_k_visc * _phi_s)
                        # Stokes–Einstein: D ∝ 1/μ
                        blend_D = blend_D * (_mu_blend_orig / blend_mu)
                        # Cp: mass-weighted with solid
                        blend_Cp = (1 - _w_s) * blend_Cp + _w_s * _sm_Cp

                        _total_mass_kg = _sol_mass_kg
                        _total_vol_L = _sol_vol_L
                        _sm_conc_kg_L = (_sm_mass_kg / _total_vol_L
                                         if _total_vol_L > 0 else 0.0)
                        _sm_added = True
                        _sm_name_val = _sm_name_input

                # Display final properties as a summary table
                _props_label = ("Solution Properties" if _sm_added
                                else "Blended Properties")
                st.subheader(_props_label)

                blend_nu = blend_mu / blend_rho if blend_rho > 0 else 0.0

                # Build rows: components → (optional SM) → blend/solution
                _summary_rows = []
                for cp in comp_props:
                    _cp_vol_pct = (cp["vol_L"] / _total_vol_L * 100
                                   if _total_vol_L > 0 else 0.0)
                    _cp_mass_pct = (cp["mass_kg"] / _total_mass_kg * 100
                                    if _total_mass_kg > 0 else 0.0)
                    _summary_rows.append({
                        "Component": cp["name"],
                        "Volume (L)": f"{cp['vol_L']:.3g}",
                        "Mass (kg)": f"{cp['mass_kg']:.3g}",
                        "Vol %": f"{_cp_vol_pct:.1f}",
                        "Mass %": f"{_cp_mass_pct:.1f}",
                        "Conc. (kg/L)": "—",
                        "Molarity (mol/L)": "—",
                        "ρ (kg/m³)": f"{cp['rho_kg_m3']:.1f}",
                        "μ (Pa·s)": f"{cp['mu_Pa_s']:.6f}",
                        "ν (m²/s)": f"{cp['mu_Pa_s'] / cp['rho_kg_m3']:.3e}",
                        "σ (N/m)": f"{cp['surface_tension_N_m']:.4f}",
                        "D (m²/s)": f"{cp['D_mol_m2_s']:.3e}",
                        "Cp (J/kg·K)": f"{cp['Cp_J_per_kgK']:.1f}",
                        "k (W/m·K)": f"{cp['k_W_per_mK']:.4f}",
                    })
                if _sm_added:
                    _sm_vol_pct = (_sm_vol_L / _total_vol_L * 100
                                   if _total_vol_L > 0 else 0.0)
                    _sm_mass_pct = (_sm_mass_kg / _total_mass_kg * 100
                                    if _total_mass_kg > 0 else 0.0)
                    _summary_rows.append({
                        "Component": f"SM: {_sm_name_val}",
                        "Volume (L)": f"{_sm_vol_L:.3g}",
                        "Mass (kg)": f"{_sm_mass_kg:.3g}",
                        "Vol %": f"{_sm_vol_pct:.2f}",
                        "Mass %": f"{_sm_mass_pct:.1f}",
                        "Conc. (kg/L)": f"{_sm_conc_kg_L:.4f}",
                        "Molarity (mol/L)": (f"{_sm_conc_kg_L * 1000.0 / _sm_mw:.3f}"
                                              if _sm_mw > 0 else "—"),
                        "ρ (kg/m³)": f"{_sm_rho:.1f}",
                        "μ (Pa·s)": "—",
                        "ν (m²/s)": "—",
                        "σ (N/m)": "—",
                        "D (m²/s)": "—",
                        "Cp (J/kg·K)": f"{_sm_Cp:.1f}",
                        "k (W/m·K)": "—",
                    })
                _sol_label = "**Solution**" if _sm_added else "**Blend**"
                _summary_rows.append({
                    "Component": _sol_label,
                    "Volume (L)": f"{_total_vol_L:.3g}",
                    "Mass (kg)": f"{_total_mass_kg:.3g}",
                    "Vol %": "100.0",
                    "Mass %": "100.0",
                    "Conc. (kg/L)": (f"{_sm_conc_kg_L:.4f}"
                                     if _sm_added else "—"),
                    "Molarity (mol/L)": (f"{_sm_conc_kg_L * 1000.0 / _sm_mw:.3f}"
                                         if _sm_added and _sm_mw > 0
                                         else "—"),
                    "ρ (kg/m³)": f"{blend_rho:.1f}",
                    "μ (Pa·s)": f"{blend_mu:.6f}",
                    "ν (m²/s)": f"{blend_nu:.3e}",
                    "σ (N/m)": f"{blend_sig:.4f}",
                    "D (m²/s)": f"{blend_D:.3e}",
                    "Cp (J/kg·K)": f"{blend_Cp:.1f}",
                    "k (W/m·K)": f"{blend_k:.4f}",
                })
                st.dataframe(
                    pd.DataFrame(_summary_rows),
                    width='content', hide_index=True)

                with st.expander("ℹ️ Mixing rules used"):
                    st.markdown(
                        "| Property | Rule | Reference |\n"
                        "|---|---|---|\n"
                        "| ρ | Ideal (volume-additive): "
                        "1/ρ = Σ(wᵢ/ρᵢ) | — |\n"
                        "| μ | Arrhenius log-mixing: "
                        "ln μ = Σ(wᵢ ln μᵢ) | Irving 1977 |\n"
                        "| σ | Volume-fraction weighted: "
                        "σ = Σ(φᵢ σᵢ) | Macleod–Sugden |\n"
                        "| D | Log-mixing: "
                        "ln D = Σ(wᵢ ln Dᵢ) | Vignes 1966 |\n"
                        "| Cp | Mass-weighted: "
                        "Cp = Σ(wᵢ Cpᵢ) | Thermodynamic |\n"
                        "| k | Volume-fraction weighted: "
                        "k = Σ(φᵢ kᵢ) | Li 1976 |"
                    )
                    if _sm_added:
                        st.markdown(
                            "**Dissolved-solid adjustments:**\n"
                            "- ρ — total mass / total volume "
                            "(additive volumes)\n"
                            "- μ — μ_blend × exp(k_μ × φ_s) "
                            "(Mooney 1951)\n"
                            "- D — Stokes–Einstein scaling (D ∝ 1/μ)\n"
                            "- Cp — mass-weighted with solid Cp"
                        )

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
                if _sm_added:
                    _auto_name += f" + {_sm_name_val} {_sm_conc_kg_L:.3f} kg/L"

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
                            "notes": (f"Blend at {blend_T:.0f} °C: {_comp_notes}"
                                      + (f" + {_sm_name_val} {_sm_conc_kg_L:.3f} kg/L"
                                         if _sm_added else "")),
                        }])
                        st.session_state.fluid_db = pd.concat(
                            [st.session_state.fluid_db, new_row], ignore_index=True)
                        _save_custom_fluids(st.session_state.fluid_db)
                        st.success(f"Added **{blend_name}** to the custom fluid database.")
                        st.rerun()
    if not blend_components:
        st.info("Select at least one fluid to get started.")

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
