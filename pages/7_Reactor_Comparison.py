"""
Page 5 – Hydrodynamic & Mixing Comparison of Selected Reactors
===============================================================
Compare key hydrodynamic parameters side-by-side across multiple reactors
for a given fluid system.  Each reactor's full operating envelope is mapped
from its RPM and volume ranges (4 corner conditions).
"""

import streamlit as st
from utils.sidebar import render_sidebar
render_sidebar()

import pandas as pd
import numpy as np
import pathlib
import plotly.graph_objects as go

from utils.calculations import (
    compute_reactor_hydro, compute_damkohler_numbers,
    settling_velocity, particle_reynolds, zwietering_njs,
    solid_liquid_mass_transfer, particle_suspension_criterion,
    reaction_rate_mol_per_s, heat_generation_rate,
    estimate_jacket_area, estimate_U, heat_removal_capacity,
    heat_balance_assessment,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── Colour palette by scale ──────────────────────────────────────────────
SCALE_COLORS = {"Lab": "#1f77b4", "Pilot": "#ff7f0e", "Manufacturing": "#2ca02c"}


def _load(key, fn):
    if key not in st.session_state:
        p = DATA_DIR / fn
        st.session_state[key] = pd.read_csv(p) if p.exists() else pd.DataFrame()
    return st.session_state[key]


def _safe_float(val, default=0.0):
    try:
        v = float(val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default


reactors = _load("reactor_db", "reactors.csv")
fluids = _load("fluid_db", "fluids.csv")
reactions = _load("reaction_db", "reactions.csv")
particles_db = _load("particle_db", "particles.csv")

st.title("📊 Reactor Comparison")

if reactors.empty:
    st.warning("Populate the Reactor Database first.")
    st.stop()

# ── Selection ────────────────────────────────────────────────────────────
st.header("Select Reactors & Conditions")

_all_reactor_names = reactors["reactor_name"].tolist()
_preferred_defaults = ["Nalas – EasyMax 102", "Cambrex – R-802", "Nalas – 20-L"]
_defaults = [n for n in _preferred_defaults if n in _all_reactor_names] or _all_reactor_names[:4]

selected_names = st.multiselect(
    "Reactors to compare",
    _all_reactor_names,
    default=_defaults,
)

if not selected_names:
    st.info("Select at least one reactor above.")
    st.stop()

col1, col2 = st.columns(2)
with col1:
    if not fluids.empty:
        fluid_name = st.selectbox("Fluid system", fluids["fluid_name"].tolist(), key="cmp_fluid")
        fluid = fluids[fluids["fluid_name"] == fluid_name].iloc[0]
        rho = float(fluid["rho_kg_m3"])
        mu = float(fluid["mu_Pa_s"])
        D_mol = float(fluid["D_mol_m2_s"])
    else:
        st.info("No fluids in database – using water defaults.")
        rho, mu, D_mol = 997.0, 0.00089, 2.3e-9
with col2:
    if not reactions.empty:
        rxn_name = st.selectbox("Reaction (for Da numbers)", reactions["reaction_name"].tolist(), key="cmp_rxn")
        rxn = reactions[reactions["reaction_name"] == rxn_name].iloc[0]
        t_rxn = float(rxn["t_rxn_s"])
        rxn_k = float(rxn["k_value"])
        rxn_C0 = float(rxn["C0_mol_L"])
        rxn_order = str(rxn.get("order", "1"))
        rxn_T_C = _safe_float(rxn.get("T_C"), 25.0)
        rxn_delta_H = _safe_float(rxn.get("delta_H_kJ_mol"), 0.0)
        if t_rxn <= 0:
            k = rxn_k
            C0 = rxn_C0
            order = rxn_order
            if order in ("1", "pseudo-1") and k > 0:
                t_rxn = np.log(2) / k
            elif order in ("2", "pseudo-2") and k * C0 > 0:
                t_rxn = 1.0 / (k * C0)
            else:
                t_rxn = 1.0
    else:
        t_rxn = 1.0
        rxn_k, rxn_C0, rxn_order, rxn_T_C, rxn_delta_H = 0.0, 0.0, "1", 25.0, 0.0

col3, col4 = st.columns(2)
with col3:
    v_s = st.number_input("Superficial gas velocity v_s (m/s)", value=0.0, min_value=0.0,
                          format="%.4f", key="cmp_vs",
                          help="Set > 0 to compute kLa. Typical: 0.001 – 0.05 m/s.")
with col4:
    coal_choice = st.selectbox("Liquid type (for kLa)",
                               ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                               key="cmp_coal")
    is_coalescing = coal_choice.startswith("Coalescing")

# ── Particle options ──────────────────────────────────────────────────────────
include_particles = st.checkbox("Include solid particles", value=False, key="cmp_include_particles")
cmp_rho_p = cmp_d50_um = cmp_phi_p = cmp_X_wt = cmp_S_zw = 0.0
if include_particles and not particles_db.empty:
    pcol1, pcol2, pcol3 = st.columns(3)
    with pcol1:
        cmp_particle_name = st.selectbox("Particle", particles_db["particle_name"].tolist(), key="cmp_particle")
        cmp_part = particles_db[particles_db["particle_name"] == cmp_particle_name].iloc[0]
        cmp_rho_p = float(cmp_part["rho_p_kg_m3"])
        cmp_d50_um = float(cmp_part["d50_um"])
        cmp_phi_p = float(cmp_part["shape_factor"])
    with pcol2:
        cmp_X_wt = st.number_input("Solids loading X (wt-%)", value=5.0, min_value=0.01,
                                    format="%.2f", key="cmp_Xwt")
    with pcol3:
        cmp_S_zw = st.number_input("Zwietering S", value=5.5, min_value=0.5,
                                    max_value=20.0, format="%.1f", key="cmp_Szw")
elif include_particles and particles_db.empty:
    st.warning("Particle database is empty.")

# ── Heat balance options ──────────────────────────────────────────────────────
include_heat = st.checkbox("Include heat balance", value=False, key="cmp_include_heat")
cmp_T_coolant = rxn_T_C - 20.0  # default 20 °C below reaction temp
if include_heat:
    hcol1, hcol2 = st.columns(2)
    with hcol1:
        cmp_T_coolant = st.number_input(
            "Coolant temperature (°C)", value=rxn_T_C - 20.0,
            format="%.1f", key="cmp_T_cool",
            help="Jacket coolant inlet temperature")
    with hcol2:
        st.markdown(f"**Reaction:** ΔH = {rxn_delta_H:.0f} kJ/mol, "
                    f"T = {rxn_T_C:.0f} °C, ΔT = {rxn_T_C - cmp_T_coolant:.1f} °C")
    if rxn_delta_H == 0:
        st.info("Selected reaction has no enthalpy value – add ΔH in the Reaction Database.")


# ── Compute 4 corners per reactor ────────────────────────────────────────

def _liquid_height(V_litres: float, D_tank: float, H_max: float) -> float:
    """Estimate liquid height for a given fill volume (cylindrical approx)."""
    A = np.pi / 4 * D_tank**2          # cross-section area  m²
    h = (V_litres / 1000.0) / A        # m
    return min(h, H_max) if h > 0 else H_max


CORNER_LABELS = [
    "min RPM / max V",
    "max RPM / max V",
    "min RPM / min V",
    "max RPM / min V",
]

envelope_rows: list[dict] = []  # one row per (reactor, corner)
skipped: list[str] = []
reactor_info: dict = {}         # stash geometry for curve precomputation

for rname in selected_names:
    r = reactors[reactors["reactor_name"] == rname].iloc[0]
    D_imp_v  = _safe_float(r.get("D_imp_m"))
    D_tank_v = _safe_float(r.get("D_tank_m"))
    H_max    = _safe_float(r.get("H_m"))
    Np_v     = _safe_float(r.get("Np"), None)
    Nq_v     = _safe_float(r.get("Nq"), None)
    scale    = r.get("scale", "")

    if D_imp_v == 0 or D_tank_v == 0 or H_max == 0:
        skipped.append(rname)
        continue

    # RPM range → rev/s
    rpm_min = _safe_float(r.get("N_rpm_min"))
    rpm_max = _safe_float(r.get("N_rpm_max"))
    n_rps_default = _safe_float(r.get("N_rps"))
    if rpm_max == 0 and n_rps_default > 0:
        rpm_max = n_rps_default * 60
        rpm_min = rpm_max          # single point
    elif rpm_max == 0:
        skipped.append(rname)
        continue
    if rpm_min == 0:
        rpm_min = rpm_max          # single point if min missing

    N_lo = rpm_min / 60.0
    N_hi = rpm_max / 60.0

    # Volume range
    V_max_L = _safe_float(r.get("V_L_max"))
    V_min_L = _safe_float(r.get("V_L_min"))
    V_default = _safe_float(r.get("V_L"))
    # Fallback: compute from geometry
    V_geo = np.pi / 4 * D_tank_v**2 * H_max * 1000  # litres
    if V_max_L == 0:
        V_max_L = V_default if V_default > 0 else V_geo
    if V_min_L == 0:
        V_min_L = V_max_L       # single volume if missing

    reactor_info[rname] = dict(
        D_imp=D_imp_v, D_tank=D_tank_v, H_max=H_max,
        Np=Np_v, Nq=Nq_v,
        N_lo=N_lo, N_hi=N_hi,
        V_max_L=V_max_L, V_min_L=V_min_L,
        rpm_max=rpm_max,
    )

    # Heat-transfer geometry for this reactor
    _r_material = str(r.get("material", ""))
    _r_bottom_dish = str(r.get("bottom_dish", ""))
    _r_U_override = _safe_float(r.get("U_W_m2K"), 0.0)
    _r_A_override = _safe_float(r.get("A_ht_m2"), 0.0)
    reactor_info[rname]["material"] = _r_material
    reactor_info[rname]["bottom_dish"] = _r_bottom_dish
    reactor_info[rname]["U_override"] = _r_U_override
    reactor_info[rname]["A_override"] = _r_A_override

    corners = [
        ("min RPM / max V", N_lo, V_max_L),
        ("max RPM / max V", N_hi, V_max_L),
        ("min RPM / min V", N_lo, V_min_L),
        ("max RPM / min V", N_hi, V_min_L),
    ]

    for label, N, V_L in corners:
        H_v = _liquid_height(V_L, D_tank_v, H_max)
        h = compute_reactor_hydro(
            N=N, D_imp=D_imp_v, D_tank=D_tank_v, H=H_v,
            rho=rho, mu=mu, Np=Np_v, Nq=Nq_v,
            v_s=v_s, coalescing=is_coalescing,
            D_mol=D_mol,
        )
        da = compute_damkohler_numbers(
            h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
            kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
        )
        # Particle parameters (if enabled)
        part_vals = {}
        if include_particles and cmp_d50_um > 0:
            _dp = cmp_d50_um * 1e-6
            _nu = mu / rho
            _drho = abs(cmp_rho_p - rho)
            _vt = settling_velocity(_dp, cmp_rho_p, rho, mu, cmp_phi_p)
            _rep = particle_reynolds(_dp, _vt, rho, mu)
            _njs = zwietering_njs(cmp_S_zw, _nu, _dp, _drho, rho, cmp_X_wt, D_imp_v)
            _ksl = solid_liquid_mass_transfer(_dp, _vt, rho, mu, D_mol)
            _njs_rpm = _njs * 60
            _n_over_njs = N / _njs if _njs > 0 else 0.0
            part_vals = {
                "N_js (RPM)": _njs_rpm,
                "N/N_js": _n_over_njs,
                "v_t (m/s)": _vt,
                "Re_p": _rep,
                "k_SL (m/s)": _ksl,
            }
        # Heat balance parameters (if enabled)
        heat_vals = {}
        if include_heat and rxn_delta_H != 0:
            _r_mol_s = reaction_rate_mol_per_s(rxn_order, rxn_k, rxn_C0, V_L)
            _Q_gen = heat_generation_rate(rxn_delta_H, _r_mol_s)
            _A_ht = _r_A_override if _r_A_override > 0 else estimate_jacket_area(D_tank_v, H_v, _r_bottom_dish)
            _U_ht = _r_U_override if _r_U_override > 0 else estimate_U(_r_material, N)
            _dT = rxn_T_C - cmp_T_coolant
            _Q_cool = heat_removal_capacity(_U_ht, _A_ht, _dT)
            heat_vals = {
                "Q_gen (W)": _Q_gen,
                "Q_cool (W)": _Q_cool,
                "U (W/m²·K)": _U_ht,
                "A_ht (m²)": _A_ht,
                "Q_gen/Q_cool": _Q_gen / _Q_cool if _Q_cool > 0 else np.inf,
            }
        envelope_rows.append({
            "Reactor": rname,
            "Scale": scale,
            "Corner": label,
            "N (rev/s)": N,
            "RPM": N * 60,
            "RPM_max": rpm_max,
            "V_L": V_L,
            **h,
            **da,
            **part_vals,
            **heat_vals,
        })

if skipped:
    st.warning(f"Skipped reactors with missing geometry/agitation data: {', '.join(skipped)}")

if not envelope_rows:
    st.info("No computable reactors in the selection.")
    st.stop()

env_df = pd.DataFrame(envelope_rows)
env_df["RPM_pct"] = env_df["RPM"] / env_df["RPM_max"] * 100

# ── Aggregate min / max per reactor for plotting ─────────────────────────
PLOT_PARAMS = [
    "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
    "Micromix time t_E (s)", "Micromix time t_E_local (s)",
    "Kolmogorov η (µm)", "Re",
    "Avg shear rate (1/s)", "Max shear rate (1/s)", "Avg shear stress (Pa)",
    "Da_macro", "Da_micro", "Da_GL", "ε_max (W/kg)",
    "kLa (1/s)", "kLa_surface (1/s)",
]
HEAT_PARAMS = ["Q_gen (W)", "Q_cool (W)", "U (W/m²·K)", "A_ht (m²)", "Q_gen/Q_cool"]
if include_heat and rxn_delta_H != 0:
    PLOT_PARAMS = PLOT_PARAMS + HEAT_PARAMS
PARTICLE_PARAMS = ["N_js (RPM)", "N/N_js", "v_t (m/s)", "Re_p", "k_SL (m/s)"]
if include_particles and cmp_d50_um > 0:
    PLOT_PARAMS = PLOT_PARAMS + PARTICLE_PARAMS

agg_dict = {p: ["min", "max"] for p in PLOT_PARAMS}
agg_dict["Scale"] = "first"
agg_dict["Volume (L)"] = ["min", "max"]
agg_df = env_df.groupby("Reactor", sort=False).agg(agg_dict)
agg_df.columns = ["_".join(c).strip("_") for c in agg_df.columns]
agg_df = agg_df.reset_index()

# ── Precompute true boundary curves at many RPM points ───────────────────
_N_INTERP = 50

curve_data: dict = {}  # rname → {pct_arr, maxV: {param: arr}, minV: {param: arr}}
for rname, info in reactor_info.items():
    N_arr = np.linspace(info["N_lo"], info["N_hi"], _N_INTERP)
    pct_arr = N_arr / info["N_hi"] * 100 if info["N_hi"] > 0 else np.zeros(_N_INTERP)

    curves: dict = {"pct_arr": pct_arr}
    for vol_key, V_L in [("maxV", info["V_max_L"]), ("minV", info["V_min_L"])]:
        H_v = _liquid_height(V_L, info["D_tank"], info["H_max"])
        param_arrs: dict = {p: np.empty(_N_INTERP) for p in PLOT_PARAMS}
        for j, N in enumerate(N_arr):
            h = compute_reactor_hydro(
                N=N, D_imp=info["D_imp"], D_tank=info["D_tank"], H=H_v,
                rho=rho, mu=mu, Np=info["Np"], Nq=info["Nq"],
                v_s=v_s, coalescing=is_coalescing, D_mol=D_mol,
            )
            da = compute_damkohler_numbers(
                h["Blend time 95% (s)"], h["Micromix time t_E (s)"], t_rxn,
                kLa=h["kLa (1/s)"], kLa_surface=h["kLa_surface (1/s)"],
            )
            vals = {**h, **da}
            # Particle parameters in curve
            if include_particles and cmp_d50_um > 0:
                _dp = cmp_d50_um * 1e-6
                _nu = mu / rho
                _drho = abs(cmp_rho_p - rho)
                _vt = settling_velocity(_dp, cmp_rho_p, rho, mu, cmp_phi_p)
                _rep = particle_reynolds(_dp, _vt, rho, mu)
                _njs = zwietering_njs(cmp_S_zw, _nu, _dp, _drho, rho, cmp_X_wt, info["D_imp"])
                _ksl = solid_liquid_mass_transfer(_dp, _vt, rho, mu, D_mol)
                vals["N_js (RPM)"] = _njs * 60
                vals["N/N_js"] = N / _njs if _njs > 0 else 0.0
                vals["v_t (m/s)"] = _vt
                vals["Re_p"] = _rep
                vals["k_SL (m/s)"] = _ksl
            # Heat balance parameters in curve
            if include_heat and rxn_delta_H != 0:
                _r_mol_s = reaction_rate_mol_per_s(rxn_order, rxn_k, rxn_C0, V_L)
                _Q_gen = heat_generation_rate(rxn_delta_H, _r_mol_s)
                _A_ht = info["A_override"] if info["A_override"] > 0 else estimate_jacket_area(info["D_tank"], H_v, info["bottom_dish"])
                _U_ht = info["U_override"] if info["U_override"] > 0 else estimate_U(info["material"], N)
                _dT = rxn_T_C - cmp_T_coolant
                _Q_cool = heat_removal_capacity(_U_ht, _A_ht, _dT)
                vals["Q_gen (W)"] = _Q_gen
                vals["Q_cool (W)"] = _Q_cool
                vals["U (W/m²·K)"] = _U_ht
                vals["A_ht (m²)"] = _A_ht
                vals["Q_gen/Q_cool"] = _Q_gen / _Q_cool if _Q_cool > 0 else np.inf
            for p in PLOT_PARAMS:
                param_arrs[p][j] = vals.get(p, np.nan)
        curves[vol_key] = param_arrs
    curve_data[rname] = curves

# ── Summary Table ─────────────────────────────────────────────────────────
st.header("Operating Envelope Summary")
st.caption("Each row shows the range across the 4 corner conditions (min/max RPM × min/max volume).")

table_rows = []
for _, a in agg_df.iterrows():
    row: dict = {"Reactor": a["Reactor"], "Scale": a["Scale_first"]}
    row["Volume (L)"] = f"{a['Volume (L)_min']:.1f} – {a['Volume (L)_max']:.1f}"
    for p in PLOT_PARAMS:
        lo, hi = a[f"{p}_min"], a[f"{p}_max"]
        if abs(lo - hi) < 1e-12:
            row[p] = f"{lo:.3g}"
        else:
            row[p] = f"{lo:.3g} – {hi:.3g}"
    table_rows.append(row)

st.dataframe(pd.DataFrame(table_rows), width="stretch", hide_index=True)

# ── Detail: 4-corner table ───────────────────────────────────────────────
with st.expander("Full 4-corner detail table", expanded=False):
    detail_cols = ["Reactor", "Corner", "N (rev/s)", "V_L", "Re",
                   "P/V (W/L)", "Tip speed (m/s)", "Blend time 95% (s)",
                   "Micromix time t_E (s)", "Micromix time t_E_local (s)",
                   "Kolmogorov η (µm)",
                   "Avg shear rate (1/s)", "Max shear rate (1/s)",
                   "Avg shear stress (Pa)", "kLa (1/s)", "kLa_surface (1/s)",
                   "Da_macro", "Da_micro", "Da_GL"]
    if include_particles and cmp_d50_um > 0:
        detail_cols += PARTICLE_PARAMS
    if include_heat and rxn_delta_H != 0:
        detail_cols += HEAT_PARAMS
    # Only include columns that exist in the dataframe
    detail_cols = [c for c in detail_cols if c in env_df.columns]
    fmt = {c: "{:.3g}" for c in detail_cols if c not in ("Reactor", "Corner")}
    st.dataframe(env_df[detail_cols].style.format(fmt), width="stretch", hide_index=True)

# ── Charts: operating envelopes ──────────────────────────────────────────
st.header("Operating Envelope Charts")

# RPM-to-percentage mapping table
st.subheader("Stir Speed Reference Table")
st.caption("Translates percentage of max RPM (chart x-axis) to actual RPM for each reactor.")

pct_steps = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
rpm_table_rows = []
for rname in agg_df["Reactor"].tolist():
    rc = env_df[env_df["Reactor"] == rname].iloc[0]
    rpm_max_val = rc["RPM_max"]
    rpm_min_val = rc["RPM"] if rc["Corner"] == "min RPM / max V" else env_df[(env_df["Reactor"] == rname) & (env_df["Corner"] == "min RPM / max V")].iloc[0]["RPM"]
    pct_min = rpm_min_val / rpm_max_val * 100
    row_data: dict = {"Reactor": rname, "RPM min": f"{rpm_min_val:.0f}", "RPM max": f"{rpm_max_val:.0f}",
                       "% min": f"{pct_min:.0f}%"}
    for pct in pct_steps:
        row_data[f"{pct}%"] = f"{rpm_max_val * pct / 100:.0f}"
    rpm_table_rows.append(row_data)

st.dataframe(pd.DataFrame(rpm_table_rows), width="stretch", hide_index=True)

st.caption(
    "Each reactor's operational region is plotted as a filled polygon. "
    "The x-axis is stir speed as a percentage of each vessel's maximum RPM. "
    "Overlapping regions indicate where reactors share comparable conditions."
)

params_to_plot = st.multiselect(
    "Parameters to plot",
    PLOT_PARAMS,
    default=["P/V (W/L)", "Blend time 95% (s)", "Da_macro", "Da_micro", "Da_GL"],
    key="env_params",
)

# Extended colour palette (one colour per reactor)
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
    "#bcbd22", "#17becf",
]

_N_INTERP = 50  # (also defined earlier for precomputation)

for param in params_to_plot:
    fig = go.Figure()

    reactor_list = agg_df["Reactor"].tolist()
    for i, rname in enumerate(reactor_list):
        color = _PALETTE[i % len(_PALETTE)]
        corners = env_df[env_df["Reactor"] == rname]

        # 4 corner values (for markers)
        c = corners.set_index("Corner")
        pct_lo = c.loc["min RPM / max V", "RPM_pct"]
        pct_hi = 100.0

        val_minN_maxV = c.loc["min RPM / max V", param]
        val_maxN_maxV = c.loc["max RPM / max V", param]
        val_maxN_minV = c.loc["max RPM / min V", param]
        val_minN_minV = c.loc["min RPM / min V", param]

        # True boundary curves from precomputed data
        curves = curve_data[rname]
        pct_arr = curves["pct_arr"]
        y_maxV = curves["maxV"][param]
        y_minV = curves["minV"][param]

        # Build polygon from actual curve points (follows nonlinear shapes)
        poly_x = np.concatenate([pct_arr, pct_arr[::-1], [pct_arr[0]]])
        poly_y = np.concatenate([y_maxV, y_minV[::-1], [y_maxV[0]]])

        # Filled polygon (shaded region)
        fig.add_trace(go.Scatter(
            x=poly_x,
            y=poly_y,
            fill="toself",
            fillcolor=color,
            opacity=0.20,
            line=dict(color=color, width=1),
            mode="lines",
            name=rname,
            showlegend=True,
            legendgroup=rname,
            hoverinfo="skip",
        ))

        # Corner markers with hover detail
        marker_x = [pct_lo, pct_hi, pct_lo, pct_hi]
        marker_y = [val_minN_maxV, val_maxN_maxV, val_minN_minV, val_maxN_minV]
        marker_labels = ["min RPM / max V", "max RPM / max V",
                         "min RPM / min V", "max RPM / min V"]
        marker_symbols = ["diamond", "circle", "square", "triangle-up"]

        fig.add_trace(go.Scatter(
            x=marker_x,
            y=marker_y,
            mode="markers",
            marker=dict(
                size=10,
                color=color,
                symbol=marker_symbols,
                line=dict(width=1, color="white"),
            ),
            text=marker_labels,
            customdata=[[rname]] * 4,
            hovertemplate=(
                "%{customdata[0]}<br>"
                "%{text}<br>"
                "RPM %%max = %{x:.0f}%%<br>"
                + param + " = %{y:.3g}<extra></extra>"
            ),
            showlegend=False,
            legendgroup=rname,
        ))

        # Max-volume boundary (solid line)
        fig.add_trace(go.Scatter(
            x=pct_arr, y=y_maxV,
            mode="lines", line=dict(color=color, width=2),
            name=f"{rname} (max V)",
            showlegend=False, legendgroup=rname,
            hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.3g}<extra>" + rname + " max V</extra>",
        ))
        # Min-volume boundary (dotted line)
        fig.add_trace(go.Scatter(
            x=pct_arr, y=y_minV,
            mode="lines", line=dict(color=color, width=2, dash="dot"),
            name=f"{rname} (min V)",
            showlegend=False, legendgroup=rname,
            hovertemplate="%% max RPM: %{x:.1f}%%<br>%{y:.3g}<extra>" + rname + " min V</extra>",
        ))

    # Reference lines for Da parameters
    if param in ("Da_macro", "Da_micro", "Da_GL"):
        import math
        for da_val, color, label in [
            (0.1, "orange", "Da=0.1 (onset of sensitivity)"),
            (1.0, "red", "Da=1 (limited)"),
        ]:
            fig.add_shape(
                type="line", x0=0, x1=1, y0=da_val, y1=da_val,
                xref="paper", yref="y",
                line=dict(color=color, width=1.5, dash="dash"),
            )
            fig.add_annotation(
                x=1.0, xref="paper", xanchor="right",
                y=math.log10(da_val), yref="y",   # log-scale: pass log10 value
                yanchor="top", yshift=-2,
                text=label, font=dict(size=11, color=color),
                showarrow=False,
            )
        fig.update_yaxes(type="log")

    # Reference line for N/N_js = 1 (just-suspended)
    if param == "N/N_js":
        fig.add_shape(
            type="line", x0=0, x1=1, y0=1.0, y1=1.0,
            xref="paper", yref="y",
            line=dict(color="red", width=1.5, dash="dash"),
        )
        fig.add_annotation(
            x=1.0, xref="paper", xanchor="right",
            y=1.0, yref="y", yanchor="bottom", yshift=2,
            text="N/N_js = 1 (just suspended)",
            font=dict(size=11, color="red"), showarrow=False,
        )

    # Reference line for Q_gen/Q_cool = 1 (insufficient cooling)
    if param == "Q_gen/Q_cool":
        fig.add_shape(
            type="line", x0=0, x1=1, y0=1.0, y1=1.0,
            xref="paper", yref="y",
            line=dict(color="red", width=1.5, dash="dash"),
        )
        fig.add_annotation(
            x=1.0, xref="paper", xanchor="right",
            y=1.0, yref="y", yanchor="bottom", yshift=2,
            text="Q_gen/Q_cool = 1 (cooling limit)",
            font=dict(size=11, color="red"), showarrow=False,
        )

    fig.update_layout(
        title=param,
        xaxis_title="Stir speed (% of vessel max RPM)",
        yaxis_title=param,
        xaxis=dict(range=[0, 105], dtick=10,
                   showspikes=True, spikemode="across",
                   spikethickness=1, spikecolor="grey", spikedash="dot"),
        yaxis=dict(showspikes=True, spikemode="across",
                   spikethickness=1, spikecolor="grey", spikedash="dot"),
        height=500,
        margin=dict(t=50, b=50),
        legend=dict(title="Reactor"),
        hovermode="closest",
    )
    st.plotly_chart(fig, width="stretch")

# ── Chart legend ─────────────────────────────────────────────────────────
with st.expander("Chart legend"):
    st.markdown("""
**Regions:** Each reactor's operational envelope is drawn as a filled polygon
spanning its RPM range (as % of max) on the x-axis.

**Boundary lines:**
- **Solid line** = max fill volume edge
- **Dotted line** = min fill volume edge

**Corner markers:**

| Symbol | Corner condition |
|--------|-----------------|
| ◆ Diamond | min RPM / max Volume |
| ● Circle | max RPM / max Volume |
| ■ Square | min RPM / min Volume |
| ▲ Triangle | max RPM / min Volume |

Overlapping shaded regions indicate where two reactors can achieve similar parameter values.
""")

# ── Heat Balance Summary ─────────────────────────────────────────────────
if include_heat and rxn_delta_H != 0:
    st.header("🌡️ Heat Balance Summary")
    st.caption(
        f"Reaction: **{rxn_name if not reactions.empty else 'N/A'}** | "
        f"ΔH = {rxn_delta_H:.0f} kJ/mol | "
        f"T_rxn = {rxn_T_C:.0f} °C | T_coolant = {cmp_T_coolant:.0f} °C | "
        f"ΔT = {rxn_T_C - cmp_T_coolant:.1f} °C"
    )

    # Build summary table: one row per reactor at max-RPM / max-V corner
    heat_summary_rows = []
    for rname in reactor_info:
        corners = env_df[(env_df["Reactor"] == rname) & (env_df["Corner"] == "max RPM / max V")]
        if corners.empty:
            continue
        c = corners.iloc[0]
        Q_gen = c.get("Q_gen (W)", 0)
        Q_cool = c.get("Q_cool (W)", 0)
        U_val = c.get("U (W/m²·K)", 0)
        A_val = c.get("A_ht (m²)", 0)
        ratio = Q_gen / Q_cool if Q_cool > 0 else np.inf
        heat_summary_rows.append({
            "Reactor": rname,
            "Volume (L)": c["V_L"],
            "U (W/m²·K)": f"{U_val:.0f}",
            "A (m²)": f"{A_val:.3f}",
            "Q_gen (W)": f"{Q_gen:.1f}",
            "Q_cool (W)": f"{Q_cool:.1f}",
            "Q_gen / Q_cool": f"{ratio:.2f}" if ratio < 100 else "∞",
            "Assessment": heat_balance_assessment(Q_gen, Q_cool),
        })

    if heat_summary_rows:
        st.dataframe(pd.DataFrame(heat_summary_rows), width="stretch", hide_index=True)

    # Bar chart: Q_gen vs Q_cool per reactor
    if heat_summary_rows:
        fig_heat = go.Figure()
        r_names = [r["Reactor"] for r in heat_summary_rows]
        q_gen_vals = [float(r["Q_gen (W)"]) for r in heat_summary_rows]
        q_cool_vals = [float(r["Q_cool (W)"]) for r in heat_summary_rows]

        fig_heat.add_trace(go.Bar(
            name="Q_gen (W)", x=r_names, y=q_gen_vals,
            marker_color="firebrick",
        ))
        fig_heat.add_trace(go.Bar(
            name="Q_cool (W)", x=r_names, y=q_cool_vals,
            marker_color="steelblue",
        ))
        fig_heat.update_layout(
            title="Heat Generation vs Cooling Capacity (max RPM / max Volume)",
            yaxis_title="Power (W)",
            barmode="group",
            height=450,
            margin=dict(t=50, b=50),
        )
        st.plotly_chart(fig_heat, width="stretch")

# ── Scale-up summary ─────────────────────────────────────────────────────
st.header("Scale-Up Impact Summary")
st.caption("Ratios use midpoint (average of 4 corners) for each parameter, relative to the first selected reactor.")

if len(agg_df) >= 2:
    # Compute midpoints per reactor
    mid_df = env_df.groupby("Reactor", sort=False)[PLOT_PARAMS + ["Volume (L)"]].mean().reset_index()

    ref = mid_df.iloc[0]
    summary_rows = []
    for _, row in mid_df.iloc[1:].iterrows():
        summary_rows.append({
            "From → To": f"{ref['Reactor']} → {row['Reactor']}",
            "Volume ratio": row["Volume (L)"] / ref["Volume (L)"] if ref["Volume (L)"] > 0 else np.nan,
            "P/V ratio": row["P/V (W/L)"] / ref["P/V (W/L)"] if ref["P/V (W/L)"] > 0 else np.nan,
            "Tip speed ratio": row["Tip speed (m/s)"] / ref["Tip speed (m/s)"] if ref["Tip speed (m/s)"] > 0 else np.nan,
            "Blend time ratio": row["Blend time 95% (s)"] / ref["Blend time 95% (s)"] if ref["Blend time 95% (s)"] > 0 else np.nan,
            "Da_macro ratio": row["Da_macro"] / ref["Da_macro"] if ref["Da_macro"] > 0 else np.nan,
        })
    summary_df = pd.DataFrame(summary_rows)
    numeric_cols = summary_df.select_dtypes(include="number").columns
    st.dataframe(summary_df.style.format({c: "{:.2f}" for c in numeric_cols}), width="stretch")
else:
    st.info("Select at least two reactors to see scale-up ratios.")
