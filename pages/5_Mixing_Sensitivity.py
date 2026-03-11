"""
Page 4 – Mixing Sensitivity Calculation Workflow
=================================================
Select a reactor, reaction, and fluid system → compute hydrodynamic parameters,
Damköhler numbers, and assess mixing sensitivity.
"""

import streamlit as st
from utils.sidebar import render_sidebar
render_sidebar()

import pandas as pd
import numpy as np
import pathlib

from utils.calculations import (
    compute_reactor_hydro,
    compute_damkohler_numbers,
    kolmogorov_length,
    batchelor_length,
    micromixing_time_engulfment,
    blend_time_turbulent,
    epsilon_max_estimate,
    impeller_power,
    reynolds_number,
    power_number_correlation,
    mixing_sensitivity_assessment,
    settling_velocity,
    particle_reynolds,
    zwietering_njs,
    solid_liquid_mass_transfer,
    particle_suspension_criterion,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

# ── helpers to load databases (reuse session state if available) ─────────

def _load_csv(key, filename, columns):
    if key not in st.session_state:
        p = DATA_DIR / filename
        if p.exists():
            st.session_state[key] = pd.read_csv(p)
        else:
            st.session_state[key] = pd.DataFrame(columns=columns)
    return st.session_state[key]


reactors = _load_csv("reactor_db", "reactors.csv", ["reactor_name"])
reactions = _load_csv("reaction_db", "reactions.csv", ["reaction_name"])
fluids = _load_csv("fluid_db", "fluids.csv", ["fluid_name"])
particles = _load_csv("particle_db", "particles.csv", ["particle_name"])

st.title("⚙️ Mixing Sensitivity Calculation Workflow")

if reactors.empty or reactions.empty or fluids.empty:
    st.warning("Please populate the Reactor, Reaction, and Fluid databases first.")
    st.stop()

# ── Step 1: Select system ────────────────────────────────────────────────
st.header("1 · Select System")
col1, col2, col3 = st.columns(3)

with col1:
    reactor_name = st.selectbox("Reactor", reactors["reactor_name"].tolist())
with col2:
    reaction_name = st.selectbox("Reaction", reactions["reaction_name"].tolist())
with col3:
    fluid_name = st.selectbox("Fluid", fluids["fluid_name"].tolist())

reactor = reactors[reactors["reactor_name"] == reactor_name].iloc[0]
reaction = reactions[reactions["reaction_name"] == reaction_name].iloc[0]
fluid = fluids[fluids["fluid_name"] == fluid_name].iloc[0]

# ── Helper: safely read a numeric field with a fallback ──────────────────
def _safe(series_val, default):
    """Return float(series_val) if non-NaN, else default."""
    try:
        v = float(series_val)
        return v if not np.isnan(v) else default
    except (ValueError, TypeError):
        return default

# Use selection names in widget keys so they reset automatically on change
_rk = reactor_name   # reactor key fragment
_fk = fluid_name     # fluid key fragment
_xk = reaction_name  # reaction key fragment

# ── Step 2: Allow overrides ──────────────────────────────────────────────
st.header("2 · Review / Override Parameters")

with st.expander("Reactor geometry & agitation", expanded=False):
    oc1, oc2, oc3, oc4 = st.columns(4)
    D_tank = oc1.number_input("D_tank (m)", value=_safe(reactor.get("D_tank_m"), 0.10), format="%.4f", key=f"ov_Dt_{_rk}")
    H = oc2.number_input("H (m)", value=_safe(reactor.get("H_m"), 0.13), format="%.4f", key=f"ov_H_{_rk}")
    D_imp = oc3.number_input("D_imp (m)", value=_safe(reactor.get("D_imp_m"), 0.05), format="%.4f", key=f"ov_Di_{_rk}")
    N_rps = oc4.number_input("N (rev/s)", value=_safe(reactor.get("N_rps"), 5.0), format="%.2f", key=f"ov_N_{_rk}")
    oc5, oc6 = st.columns(2)
    Np_in = oc5.number_input("Np", value=_safe(reactor.get("Np"), 1.27), format="%.2f", key=f"ov_Np_{_rk}")
    Nq_in = oc6.number_input("Nq", value=_safe(reactor.get("Nq"), 0.79), format="%.2f", key=f"ov_Nq_{_rk}")

with st.expander("Fluid properties", expanded=False):
    fc1, fc2, fc3 = st.columns(3)
    rho = fc1.number_input("ρ (kg/m³)", value=float(fluid["rho_kg_m3"]), format="%.1f", key=f"ov_rho_{_fk}")
    mu = fc2.number_input("μ (Pa·s)", value=float(fluid["mu_Pa_s"]), format="%.6f", key=f"ov_mu_{_fk}")
    D_mol = fc3.number_input("D_mol (m²/s)", value=float(fluid["D_mol_m2_s"]), format="%.2e", key=f"ov_Dmol_{_fk}")

# ── Particle (solid) phase ────────────────────────────────────────────────
include_particles = st.checkbox("Include solid particles", value=False, key="include_particles")
if include_particles and not particles.empty:
    with st.expander("Particle properties", expanded=True):
        def _on_particle_change():
            """Reset override widgets when a new particle is selected."""
            pname = st.session_state["sel_particle"]
            p = particles[particles["particle_name"] == pname].iloc[0]
            st.session_state["ov_rhop"] = float(p["rho_p_kg_m3"])
            st.session_state["ov_d50"] = float(p["d50_um"])
            st.session_state["ov_phi"] = float(p["shape_factor"])

        particle_name = st.selectbox(
            "Particle", particles["particle_name"].tolist(),
            key="sel_particle", on_change=_on_particle_change,
        )
        particle = particles[particles["particle_name"] == particle_name].iloc[0]
        # Initialise session state on very first render
        if "ov_rhop" not in st.session_state:
            st.session_state["ov_rhop"] = float(particle["rho_p_kg_m3"])
        if "ov_d50" not in st.session_state:
            st.session_state["ov_d50"] = float(particle["d50_um"])
        if "ov_phi" not in st.session_state:
            st.session_state["ov_phi"] = float(particle["shape_factor"])

        pc1, pc2, pc3 = st.columns(3)
        rho_p = pc1.number_input("ρ_p (kg/m³)", min_value=100.0, max_value=25000.0,
                                 step=10.0, format="%.0f", key="ov_rhop")
        d50_um = pc2.number_input("D50 (µm)", min_value=0.1, max_value=10000.0,
                                  step=1.0, format="%.1f", key="ov_d50")
        phi_p = pc3.number_input("Shape factor φ",
                                 min_value=0.01, max_value=1.0, step=0.05,
                                 format="%.2f", key="ov_phi")
        pc4, pc5 = st.columns(2)
        X_wt = pc4.number_input("Solids loading X (wt-%)", value=5.0, min_value=0.01,
                                format="%.2f", key="ov_Xwt",
                                help="Mass fraction of solids as weight-%")
        S_zw = pc5.number_input("Zwietering S constant", value=5.5, min_value=0.5,
                                max_value=20.0, format="%.1f", key="ov_Szw",
                                help="Geometry-dependent constant (typ. 1–10, PBT ≈ 5.5)")
elif include_particles and particles.empty:
    st.warning("Particle database is empty. Add particles on the Particle Database page.")

with st.expander("Gas sparging (for kLa)", expanded=False):
    gc1, gc2 = st.columns(2)
    v_s = gc1.number_input("Superficial gas velocity v_s (m/s)", value=0.0, min_value=0.0, format="%.4f",
                           key=f"ov_vs_{_rk}", help="Set > 0 to compute kLa. Typical: 0.001 – 0.05 m/s.")
    coalescing = gc2.selectbox("Liquid type", ["Coalescing (pure liquid)", "Non-coalescing (electrolyte)"],
                               key=f"ov_coal_{_fk}")
    is_coalescing = coalescing.startswith("Coalescing")

with st.expander("Reaction parameters", expanded=False):
    rc1, rc2, rc3 = st.columns(3)
    k_val = rc1.number_input("k", value=float(reaction["k_value"]), format="%.6g", key=f"ov_k_{_xk}")
    C0 = rc2.number_input("C₀ (mol/L)", value=float(reaction["C0_mol_L"]), format="%.4g", key=f"ov_C0_{_xk}")
    t_rxn_input = rc3.number_input("t_rxn (s)", value=float(reaction["t_rxn_s"]), format="%.4g", key=f"ov_trxn_{_xk}",
                                    help="Characteristic reaction time. 0 = auto-compute.")

# ── Step 3: Compute ──────────────────────────────────────────────────────
st.header("3 · Results")

# Compute reaction time if needed
t_rxn = t_rxn_input
if t_rxn <= 0 and k_val > 0:
    order = str(reaction.get("order", "1"))
    if order in ("1", "pseudo-1"):
        t_rxn = np.log(2) / k_val
    elif order in ("2", "pseudo-2") and C0 > 0:
        t_rxn = 1.0 / (k_val * C0)
    else:
        t_rxn = 1.0  # fallback

hydro = compute_reactor_hydro(
    N=N_rps, D_imp=D_imp, D_tank=D_tank, H=H,
    rho=rho, mu=mu, Np=Np_in, Nq=Nq_in,
    v_s=v_s, coalescing=is_coalescing,
    D_mol=D_mol,
)

da = compute_damkohler_numbers(
    t_blend=hydro["Blend time 95% (s)"],
    t_micro=hydro["Micromix time t_E (s)"],
    t_rxn=t_rxn,
    kLa=hydro["kLa (1/s)"],
    kLa_surface=hydro["kLa_surface (1/s)"],
)

# ── Display ──────────────────────────────────────────────────────────────

# Key metrics row
m1, m2, m3, m4 = st.columns(4)
m1.metric("Re", f"{hydro['Re']:.0f}")
m2.metric("P/V (W/L)", f"{hydro['P/V (W/L)']:.3g}")
m3.metric("Blend time (s)", f"{hydro['Blend time 95% (s)']:.2f}")
m4.metric("Micromix t_E (s)", f"{hydro['Micromix time t_E (s)']:.4g}")

m5, m6, m7, m8 = st.columns(4)
m5.metric("Tip speed (m/s)", f"{hydro['Tip speed (m/s)']:.2f}")
m6.metric("Kolmogorov η (µm)", f"{hydro['Kolmogorov η (µm)']:.1f}")
m7.metric("Da_macro", f"{da['Da_macro']:.3g}")
m8.metric("Da_micro", f"{da['Da_micro']:.3g}")

# New parameters row
m9, m10, m11, m12 = st.columns(4)
m9.metric("Avg shear rate (1/s)", f"{hydro['Avg shear rate (1/s)']:.1f}")
m10.metric("Max shear rate (1/s)", f"{hydro['Max shear rate (1/s)']:.0f}")
m11.metric("Avg shear stress (Pa)", f"{hydro['Avg shear stress (Pa)']:.3g}")
m12.metric("t_E local (s)", f"{hydro['Micromix time t_E_local (s)']:.4g}")

if v_s > 0:
    m13, m14, m15, _ = st.columns(4)
    m13.metric("kLa sparged (1/s)", f"{hydro['kLa (1/s)']:.4g}")
    m14.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
    m15.metric("Da_GL", f"{da['Da_GL']:.3g}")
else:
    m13, m14, _, _ = st.columns(4)
    m13.metric("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}")
    da_gl_val = da['Da_GL']
    if da_gl_val > 0:
        m14.metric("Da_GL (surface)", f"{da_gl_val:.3g}")
    st.caption("Sparged kLa not shown (set v_s > 0 in *Gas sparging* section). "
               "Surface kLa uses the Lamont-Scott free-surface model.")

# Assessment banner
assessment = da["Assessment"]
if "Strongly mixing-limited" in assessment or "Mixing-sensitive" in assessment or "transfer-limited" in assessment:
    st.error(f"🔴 **{assessment}**")
elif "Potentially sensitive" in assessment or "Potentially transfer" in assessment:
    st.warning(f"🟡 **{assessment}**")
else:
    st.success(f"🟢 **{assessment}**")

# Batchelor length
lam_B = batchelor_length(mu / rho, hydro["P/V (W/m³)"], D_mol)
st.info(f"Batchelor microscale λ_B = {lam_B * 1e6:.2f} µm  |  Reaction time t_rxn = {t_rxn:.4g} s")

# ── Particle (solid-phase) results ──────────────────────────────────────
particle_results = {}
particle_meta = {}
if include_particles and not particles.empty:
    st.subheader("Solid-Phase Parameters")
    d_p_m = d50_um * 1e-6  # convert µm → m
    nu = mu / rho
    delta_rho = abs(rho_p - rho)

    v_t = settling_velocity(d_p_m, rho_p, rho, mu, phi_p)
    Re_p = particle_reynolds(d_p_m, v_t, rho, mu)
    N_js = zwietering_njs(S_zw, nu, d_p_m, delta_rho, rho, X_wt, D_imp)
    k_SL = solid_liquid_mass_transfer(d_p_m, v_t, rho, mu, D_mol)
    susp = particle_suspension_criterion(N_rps, N_js)

    particle_results = {
        "d50 (µm)": d50_um,
        "ρ_p (kg/m³)": rho_p,
        "Shape factor φ": phi_p,
        "X (wt-%)": X_wt,
        "v_t (m/s)": v_t,
        "Re_p": Re_p,
        "N_js (rev/s)": N_js,
        "N_js (RPM)": N_js * 60,
        "k_SL (m/s)": k_SL,
    }
    particle_meta = {
        "Particle": particle_name,
        "Suspension": susp,
    }

    sp1, sp2, sp3, sp4 = st.columns(4)
    sp1.metric("Settling velocity (m/s)", f"{v_t:.3e}")
    sp2.metric("Re_p", f"{Re_p:.3g}")
    sp3.metric("N_js (RPM)", f"{N_js * 60:.1f}")
    sp4.metric("k_SL (m/s)", f"{k_SL:.3e}")

    # Suspension assessment
    if "Poorly" in susp:
        st.error(f"🔴 **{susp}** — current speed is below just-suspended speed")
    elif "Partially" in susp:
        st.warning(f"🟡 **{susp}** — increase speed to achieve full suspension")
    elif "Just" in susp:
        st.info(f"🟢 **{susp}** — operating near just-suspended conditions")
    else:
        st.success(f"🟢 **{susp}** — particles are well suspended")

# Full table
st.subheader("Full Hydrodynamic Parameter Table")
hydro_df = pd.DataFrame([hydro]).T
hydro_df.columns = ["Value"]
if particle_results:
    part_df = pd.DataFrame([particle_results]).T
    part_df.columns = ["Value"]
    hydro_df = pd.concat([hydro_df, part_df])
st.dataframe(hydro_df, width="stretch")
if particle_results:
    st.caption(f"Particle: **{particle_meta['Particle']}**  ·  {particle_meta['Suspension']}")

# ── Step 4: Save to recorded results ────────────────────────────────────
st.header("4 · Save Result")

if st.button("📌 Save this result to Recorded Results"):
    result_row = {
        "reactor": reactor_name,
        "reaction": reaction_name,
        "fluid": fluid_name,
        "Re": hydro["Re"],
        "P/V (W/L)": hydro["P/V (W/L)"],
        "Tip speed (m/s)": hydro["Tip speed (m/s)"],
        "Blend time (s)": hydro["Blend time 95% (s)"],
        "Micromix t_E (s)": hydro["Micromix time t_E (s)"],
        "Micromix t_E_local (s)": hydro["Micromix time t_E_local (s)"],
        "Kolmogorov η (µm)": hydro["Kolmogorov η (µm)"],
        "Batchelor λ_B (µm)": lam_B * 1e6,
        "Avg shear rate (1/s)": hydro["Avg shear rate (1/s)"],
        "Max shear rate (1/s)": hydro["Max shear rate (1/s)"],
        "Avg shear stress (Pa)": hydro["Avg shear stress (Pa)"],
        "kLa (1/s)": hydro["kLa (1/s)"],
        "kLa_surface (1/s)": hydro["kLa_surface (1/s)"],
        "t_rxn (s)": t_rxn,
        "Da_macro": da["Da_macro"],
        "Da_micro": da["Da_micro"],
        "Da_GL": da["Da_GL"],
        "Assessment": da["Assessment"],
    }
    if particle_results:
        result_row.update({
            "Particle": particle_meta.get("Particle", ""),
            "d50 (µm)": particle_results.get("d50 (µm)", ""),
            "ρ_p (kg/m³)": particle_results.get("ρ_p (kg/m³)", ""),
            "v_t (m/s)": particle_results.get("v_t (m/s)", ""),
            "Re_p": particle_results.get("Re_p", ""),
            "N_js (RPM)": particle_results.get("N_js (RPM)", ""),
            "k_SL (m/s)": particle_results.get("k_SL (m/s)", ""),
            "Suspension": particle_meta.get("Suspension", ""),
        })
    if "recorded_results" not in st.session_state:
        st.session_state.recorded_results = pd.DataFrame()
    st.session_state.recorded_results = pd.concat(
        [st.session_state.recorded_results, pd.DataFrame([result_row])],
        ignore_index=True,
    )
    # Also persist to CSV
    results_csv = DATA_DIR / "recorded_results.csv"
    st.session_state.recorded_results.to_csv(results_csv, index=False)
    st.success("Result saved! View it on the **Recorded Results** page.")
