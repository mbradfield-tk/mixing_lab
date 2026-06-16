"""
Page 8 – Bourne Protocol
====================================================
Guides the user step-by-step through the Bourne mixing-sensitivity screening
protocol as described by Bourne (2003).

The protocol systematically varies:
  Test 1 – Impeller speed  (does mixing matter at all?)
  Test 2 – Feed rate/time  (mesomixing vs micromixing)
  Test 3 – Feed location   (mesomixing vs macromixing)
  Confirmatory – # feed points / viscosity change

Based on process-response changes at each level, the decision tree identifies
which mixing scale (macro / meso / micro) is controlling.

References
----------
- Bourne, J.R. (2003). Mixing and the Selectivity of Chemical Reactions.
  Org. Process Res. Dev., 7(4), 471-508.
- Sarafinas, A. (2018). Test Process Mixing Sensitivities Using the Bourne
  Protocol. Scientific Update Webinar, 13 Nov 2018.
- Sarafinas, A. & Teich, C.I. (2016). Chapter 13 in Advances in Industrial
  Mixing (Kresta et al., Wiley).
"""

import streamlit as st

import json
import pandas as pd
import numpy as np
import pathlib

from utils.solvent_properties import (
    SOLVENT_DB, get_properties, is_known_solvent,
)
from utils.calculations import (
    compute_reactor_hydro,
    reynolds_number,
    power_number_correlation,
    impeller_power,
    power_per_volume,
    tip_speed,
    micromixing_time_engulfment,
    blend_time_turbulent,
    kolmogorov_length,
)

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"

from utils.data_helpers import load_db, safe_float, all_fluid_names, safe_iloc, reactor_search_name

reactors = load_db("reactor_db", "reactors.csv")
custom_fluids = load_db("fluid_db", "fluids.csv")

_all_fluid_names = all_fluid_names(custom_fluids)

# Response metrics relevant to reaction / crystallization screening in
# small-molecule API manufacturing. Units are captured separately (see
# RESPONSE_UNITS) so the metric name stays unit-free.
RESPONSE_METRICS = [
    "Yield",
    "Conversion",
    "Selectivity",
    "Concentration",
    "HPLC area",
    "Potency",
    "Assay",
    "Total impurities",
    "Largest single impurity",
    "Enantiomeric excess (ee)",
    "Diastereomeric ratio (dr)",
    "Residual solvent",
    "Particle size D50",
    "Particle size D90",
    "Span (D90−D10)/D50",
    "Polymorph / form purity",
    "Filtration rate (or time)",
    "Bulk density",
    "Cake density"
]

# Common units for the metrics above. Users may also enter a custom unit.
RESPONSE_UNITS = [
    "%",
    "% w/w",
    "wt%",
    "area %",
    "g/L",
    "mg/mL",
    "mol/L",
    "ppm",
    "µm",
    "kg/m²/h",
    "g/mL",
    "ratio",
    "(none)",
    "s",
    "min",
    "h",
]

_CUSTOM_UNIT = "Custom…"

# Sensible default unit for each metric (used to preselect the unit dropdown).
_METRIC_DEFAULT_UNIT = {
    "Yield": "%",
    "Conversion": "%",
    "Selectivity": "%",
    "Concentration": "g/L",
    "HPLC area": "area %",
    "Potency": "% w/w",
    "Assay": "% w/w",
    "Total impurities": "%",
    "Largest single impurity": "%",
    "Enantiomeric excess (ee)": "%",
    "Diastereomeric ratio (dr)": "ratio",
    "Residual solvent": "ppm",
    "Particle size D50": "µm",
    "Particle size D90": "µm",
    "Span (D90−D10)/D50": "ratio",
    "Polymorph / form purity": "%",
    "Filtration rate": "kg/m²/h",
    "Bulk density": "g/mL",
}


def _metric_index(name: str) -> int:
    """Return the index of a metric in RESPONSE_METRICS, or 0 if not found."""
    try:
        return RESPONSE_METRICS.index(name)
    except ValueError:
        return 0


def _unit_selector(container, key, default_unit="%", label="Unit"):
    """Render a unit dropdown (with a custom-entry option) and return the unit."""
    _opts = RESPONSE_UNITS + [_CUSTOM_UNIT]
    _idx = _opts.index(default_unit) if default_unit in _opts else 0
    _sel = container.selectbox(label, _opts, index=_idx, key=key)
    if _sel == _CUSTOM_UNIT:
        _sel = container.text_input(
            "Custom unit", value="", key=f"{key}_custom",
            placeholder="e.g. g/100 g",
        ) or ""
    return _sel


def _fmt_metric(metric: str, unit: str) -> str:
    """Combine a metric name and unit into a display label."""
    unit = (unit or "").strip()
    if unit and unit != "(none)":
        return f"{metric} ({unit})"
    return metric

# ══════════════════════════════════════════════════════════════════════════
# Constrain the page content to a centered, narrower column so widgets and
# tables don't stretch across the full (wide-layout) screen width.
st.markdown(
    """
    <style>
    section.main div.block-container {
        max-width: 900px;
        margin-left: auto;
        margin-right: auto;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🅱️ Bourne Protocol – Mixing Sensitivity Screening")
st.markdown("""
Determines **which mixing scale** (macro / meso / micro) controls your
process outcome. Applicable to any semi-batch stirred-tank process with
competitive rate processes.
""")

_BP_IMG = pathlib.Path(__file__).resolve().parent.parent / "images" / "general" / "bourne_protocol_decision_tree.png"
with st.expander("📋 Protocol overview – Decision Tree", expanded=False):
    if _BP_IMG.exists():
        st.image(str(_BP_IMG), caption="Bourne Protocol – Decision Tree")
    else:
        st.info("Decision tree image not found. Generate it from the ⚙️ Admin page.")

# ── Generation counter for reliable widget reset ─────────────────────────
if "_bp_gen" not in st.session_state:
    st.session_state["_bp_gen"] = 0

def _bk(name: str) -> str:
    """Return a widget key namespaced by the current generation counter."""
    return f"bp_{st.session_state['_bp_gen']}_{name}"

def _reset_bourne():
    """Increment generation counter and clear old Bourne Protocol state."""
    old_gen = st.session_state.get("_bp_gen", 0)
    # Remove all keys from the old generation (and non-generation legacy keys)
    for k in list(st.session_state.keys()):
        if k == "_bp_gen" or k == "bp_restart":
            continue
        if k.startswith("bp_") or k.startswith("_sel_bp_"):
            del st.session_state[k]
    st.session_state["_bp_gen"] = old_gen + 1
    # Return to the start screen so the user re-begins the protocol.
    st.session_state["_bp_started"] = False


def _json_default(o):
    """Fallback JSON encoder for numpy scalar / array values."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _collect_bourne_inputs():
    """Gather all user inputs (system + tests) into a tidy key/scope/value frame."""
    gen = st.session_state.get("_bp_gen", 0)
    prefix = f"bp_{gen}_"
    # These widgets are restored via their companion "_sel_bp_*" keys (which use
    # them as the default value/index), so we skip the raw widget keys here.
    # The button keys carry no reusable input and are skipped too.
    _skip = {"reactor", "fluid", "project", "step_no", "uop",
             "t1_assess", "t2_assess", "t3_assess", "export_pdf", "save"}
    rows = []
    for k, v in list(st.session_state.items()):
        if k.startswith(prefix):
            logical = k[len(prefix):]
            if (logical in _skip
                    or logical.startswith("gen_step")
                    or logical.startswith("dl_step")):
                continue
            try:
                vj = json.dumps(v, default=_json_default)
            except (TypeError, ValueError):
                continue
            rows.append({"key": logical, "scope": "gen", "value": vj})
        elif k.startswith("_sel_bp_"):
            try:
                vj = json.dumps(v, default=_json_default)
            except (TypeError, ValueError):
                continue
            rows.append({"key": k, "scope": "raw", "value": vj})
    df = pd.DataFrame(rows, columns=["key", "scope", "value"])
    if not df.empty:
        df = df.sort_values(["scope", "key"]).reset_index(drop=True)
    return df


def _apply_bourne_inputs(df):
    """Restore inputs from an exported frame into a fresh widget generation."""
    old_gen = st.session_state.get("_bp_gen", 0)
    for k in list(st.session_state.keys()):
        if k.startswith(f"bp_{old_gen}_") or k.startswith("_sel_bp_"):
            del st.session_state[k]
    new_gen = old_gen + 1
    for _, row in df.iterrows():
        try:
            val = json.loads(row["value"])
        except (TypeError, ValueError):
            continue
        if row["scope"] == "gen":
            st.session_state[f"bp_{new_gen}_{row['key']}"] = val
        else:
            st.session_state[str(row["key"])] = val
    st.session_state["_bp_gen"] = new_gen
    st.session_state["_bp_started"] = True


st.button("🔄 Restart protocol", key="bp_restart", on_click=_reset_bourne)

st.divider()

# ════════════════════════════════════════════════════════════════════════════
# Save / Load all protocol inputs (CSV) so a run can be reproduced later.
# ════════════════════════════════════════════════════════════════════════════
with st.expander("💾 Save / Load Protocol Inputs", expanded=False):
    st.caption(
        "Export every input you've entered (system definition plus all three "
        "tests) to a CSV, then re-import it later to reproduce and re-run the "
        "protocol from where you left off."
    )
    _exp_col, _imp_col = st.columns(2)
    with _exp_col:
        st.markdown("**Export inputs**")
        _inputs_df = _collect_bourne_inputs()
        _proj_tag = (st.session_state.get("_sel_bp_project") or "protocol").strip().replace(" ", "_") or "protocol"
        st.download_button(
            "⬇️ Download inputs CSV",
            data=_inputs_df.to_csv(index=False).encode("utf-8"),
            file_name=f"Bourne_inputs_{_proj_tag}.csv",
            mime="text/csv",
            key="bp_export_inputs",
            disabled=_inputs_df.empty,
        )
        if _inputs_df.empty:
            st.caption("Nothing entered yet — fill in the protocol first.")
    with _imp_col:
        st.markdown("**Import inputs**")
        _up = st.file_uploader("Upload inputs CSV", type=["csv"], key="bp_import_file")
        if _up is not None and st.button("📤 Load these inputs", key="bp_apply_import"):
            try:
                _imp_df = pd.read_csv(_up, dtype=str, keep_default_na=False)
                _missing = {"key", "scope", "value"} - set(_imp_df.columns)
                if _missing:
                    st.error(f"Invalid inputs file (missing columns: {', '.join(sorted(_missing))}).")
                else:
                    _apply_bourne_inputs(_imp_df)
                    st.toast("Protocol inputs loaded.", icon="✅")
                    st.rerun()
            except Exception as exc:
                st.error(f"Failed to load inputs: {exc}")

# ══════════════════════════════════════════════════════════════════════════
# SECTION 00 – Project & Step Metadata
# ══════════════════════════════════════════════════════════════════════════
st.header("Project Information")

_UNIT_OPS = [
    "Reaction",
    "Crystallization",
    "Reaction and crystallization",
    "Extraction / Wash",
    "Dissolution",
    "Distillation",
    "Other",
]

col_proj, col_step, col_uop = st.columns([2, 1, 2])
with col_proj:
    project_name = st.text_input(
        "Project name", value=st.session_state.get("_sel_bp_project", ""),
        key=_bk("project"), placeholder="e.g. Project Apollo / API-123",
    )
    st.session_state["_sel_bp_project"] = project_name
with col_step:
    step_number = st.text_input(
        "Step number", value=st.session_state.get("_sel_bp_step", ""),
        key=_bk("step_no"), placeholder="e.g. 3",
    )
    st.session_state["_sel_bp_step"] = step_number
with col_uop:
    _uop_default = st.session_state.get("_sel_bp_uop", _UNIT_OPS[0])
    unit_operation = st.selectbox(
        "Type of unit operation", _UNIT_OPS,
        index=_UNIT_OPS.index(_uop_default) if _uop_default in _UNIT_OPS else 0,
        key=_bk("uop"),
    )
    if unit_operation == "Other":
        unit_operation = st.text_input(
            "Specify unit operation", value="",
            key=_bk("uop_custom"), placeholder="e.g. Phase separation",
        ) or "Other"
    st.session_state["_sel_bp_uop"] = unit_operation

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 0 – Reactor & Fluid Selection
# ══════════════════════════════════════════════════════════════════════════
st.header("Define Your System")

# Persist selections across page navigations
_reactor_list = reactors["reactor_name"].tolist() if not reactors.empty else []

def _sel_idx(lst, key, default=0):
    val = st.session_state.get(key)
    if val in lst:
        return lst.index(val)
    return default

col_r, col_f = st.columns(2)

with col_r:
    if not reactors.empty:
        reactor_name = st.selectbox("Reactor", _reactor_list, index=_sel_idx(_reactor_list, "_sel_bp_reactor"), key=_bk("reactor"),
                                    format_func=lambda n: reactor_search_name(reactors, n))
        st.session_state["_sel_bp_reactor"] = reactor_name
        r = safe_iloc(reactors, "reactor_name", reactor_name, "Reactor")
        D_tank = safe_float(r.get("D_tank_m"), 0.10)
        H = safe_float(r.get("H_m"), 0.13)
        D_imp = safe_float(r.get("D_imp_m"), 0.05)
        # NOTE: the Bourne protocol uses a single fixed power number (the reactor's
        # tabulated Np, or 1.27 if unset) for every speed. This differs from the
        # Mixing Assessment / Reactor Comparison pages, which derive Np from a
        # Reynolds-dependent correlation or ROM. For the same reactor the power
        # (and thus P/V) reported here can therefore differ from those pages.
        Np_val = safe_float(r.get("Np"), 1.27)
        Nq_val = safe_float(r.get("Nq"), 0.79)
        V_L_min = safe_float(r.get("V_L_min"), 0.0) or safe_float(r.get("V_L"), 0.0)
        V_L_max = safe_float(r.get("V_L_max"), 0.0) or safe_float(r.get("V_L"), 0.0)
        V_L_avg = (V_L_min + V_L_max) / 2.0
        # RPM bounds from reactor database
        N_rpm_min = safe_float(r.get("N_rpm_min"), 0.0) or None
        N_rpm_max = safe_float(r.get("N_rpm_max"), 0.0) or None
        # Use average of min/max RPM as the centerpoint speed
        if N_rpm_min is not None and N_rpm_max is not None:
            N_center = (N_rpm_min + N_rpm_max) / 2.0 / 60.0  # rev/s
        else:
            N_center = safe_float(r.get("N_rps"), 5.0)
        N_rps_min = N_rpm_min / 60.0 if N_rpm_min is not None else None
        N_rps_max = N_rpm_max / 60.0 if N_rpm_max is not None else None
    else:
        st.info("No reactors in database — enter manually below.")
        D_tank = st.number_input("Tank diameter (m)", value=0.10, format="%.4f")
        H = st.number_input("Liquid height (m)", value=0.13, format="%.4f")
        D_imp = st.number_input("Impeller diameter (m)", value=0.05, format="%.4f")
        N_center = st.number_input("Centerpoint speed (rev/s)", value=5.0, format="%.2f")
        Np_val = st.number_input("Power number Np", value=1.27, format="%.2f")
        Nq_val = st.number_input("Pumping number Nq", value=0.79, format="%.2f")
        V_L_avg = np.pi / 4 * D_tank**2 * H * 1000
        V_L_min = V_L_avg
        V_L_max = V_L_avg
        N_rpm_min = None
        N_rpm_max = None
        N_rps_min = None
        N_rps_max = None

col_t, col_p = st.columns(2)

with col_f:
    fluid_name = st.selectbox("Fluid system", _all_fluid_names,
                              index=_sel_idx(_all_fluid_names, "_sel_bp_fluid"),
                              key=_bk("fluid"))
    st.session_state["_sel_bp_fluid"] = fluid_name
    _is_solvent = is_known_solvent(fluid_name)
    if _is_solvent:
        _bp_T_C = col_t.number_input("Temperature (°C)", value=25.0, step=1.0,
                                    format="%.1f", key=_bk("fluid_T"))
        _bp_P_atm = col_p.number_input("Pressure (atm)", value=1.0, min_value=0.01,
                                     step=0.1, format="%.2f", key=_bk("fluid_P"))
        _fprops = get_properties(fluid_name, _bp_T_C, _bp_P_atm)
        rho = _fprops["rho_kg_m3"]
        mu = _fprops["mu_Pa_s"]
        if not _fprops["in_range"]:
            st.warning(f"⚠️ {_bp_T_C:.1f} °C is outside the liquid range "
                       f"({_fprops['mp_C']:.0f} – {_fprops['bp_at_P_C']:.0f} °C) "
                       f"for {fluid_name}. Values are extrapolated.")
    else:
        _cust = safe_iloc(custom_fluids, "fluid_name", fluid_name, "Custom fluid")
        rho = float(_cust["rho_kg_m3"])
        mu = float(_cust["mu_Pa_s"])

nu = mu / rho if rho > 0 else 0.0

# ── Volume selection ──────────────────────────────────────────────────────
st.subheader("Working Volume")
if V_L_min != V_L_max:
    st.caption(f"Reactor volume range: {V_L_min:.3f} – {V_L_max:.3f} L  •  Average: {V_L_avg:.3f} L")

_vol_min = V_L_min if V_L_min > 0 else 0.001
_reactor_tag = reactor_name.replace(" ", "_") if not reactors.empty else "manual"
V_L = st.number_input(
    "Working volume (L)", min_value=_vol_min, value=_vol_min,
    step=0.001, format="%.3f", key=_bk(f"vol_L_{_reactor_tag}"),
    help="Defaults to the minimum volume from the reactor database.",
)
V_m3 = V_L / 1000.0  # m³

st.divider()

# ── Start gate ────────────────────────────────────────────────────────────
# The protocol walkthrough (Test 1 onward) and all results stay hidden until
# the user has defined their system and clicks Start Protocol.
if not st.session_state.get("_bp_started", False):
    st.info(
        "Define your reactor, fluid system, and working volume above, then "
        "click **Start Protocol** to begin the walkthrough."
    )
    if st.button("🚀 Start Protocol", type="primary", key="bp_start"):
        st.session_state["_bp_started"] = True
        st.rerun()
    st.stop()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 1 – TEST 1: Impeller Speed  (Does mixing matter?)
# ══════════════════════════════════════════════════════════════════════════
st.header("Test 1 - Vary Impeller Speed")
st.markdown("""
Vary **impeller speed only** (hold feed rate & location constant).\n
Target ≈**100× P/m range** across three agitation speeds.\n
If the response changes → mixing matters; proceed to Test 2.

**Speed selection:** Default centerpoint P/m ≈ **0.2 W/kg**. Set high/low at
**10×** and **0.1×** (speed ratios ≈ 2.15× and 0.46×, since P ∝ N³).
Use plant condition as centerpoint if targeting existing equipment.
""")

st.subheader("Experimental Conditions")

# Default centerpoint P/m from Sarafinas (2018); allow user to customize or use plant RPM
centerpoint_mode = st.radio(
    "Centerpoint selection method",
    ["Default (0.2 W/kg)", "Custom P/m", "Custom RPM"],
    horizontal=True, key=_bk("ctr_mode"),
)

if centerpoint_mode == "Custom P/m":
    PV_center_wkg = st.number_input("Centerpoint P/m (W/kg)", value=0.2, format="%.4g", key=_bk("pv_ctr"))
    _N_from_pv = (PV_center_wkg * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
    st.write(f"Custom centerpoint: **P/m = {PV_center_wkg:.4g} W/kg** "
             f"at N = {_N_from_pv:.2f} rev/s ({_N_from_pv*60:.0f} RPM)")
elif centerpoint_mode == "Custom RPM":
    _rpm_default = float(N_center * 60) if N_center else 300.0
    _custom_rpm = st.number_input("Centerpoint RPM", value=_rpm_default, min_value=1.0,
                                  format="%.1f", key=_bk("rpm_ctr"))
    _N_custom = _custom_rpm / 60.0
    _P_custom = impeller_power(Np_val, rho, _N_custom, D_imp)
    PV_center_wkg = power_per_volume(_P_custom, V_m3) / rho
    st.write(f"Custom centerpoint: **N = {_custom_rpm:.1f} RPM** → "
             f"**P/m = {PV_center_wkg:.4g} W/kg** "
             f"({power_per_volume(_P_custom, V_m3)/1000:.4g} W/L)")
    if N_rpm_min is not None and N_rpm_max is not None:
        if _custom_rpm < N_rpm_min or _custom_rpm > N_rpm_max:
            st.warning(f"⚠️ {_custom_rpm:.1f} RPM is outside the reactor range "
                       f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM).")
else:
    # Default centerpoint P/m = 0.2 W/kg (Sarafinas, 2018).
    # Check whether the reactor can achieve it within its RPM range.
    _DEFAULT_PV = 0.2  # W/kg
    _N_for_default = (_DEFAULT_PV * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
    _N_for_default_rpm = _N_for_default * 60

    _can_reach_default = True
    if N_rps_min is not None and _N_for_default < N_rps_min:
        _can_reach_default = False
    if N_rps_max is not None and _N_for_default > N_rps_max:
        _can_reach_default = False

    if _can_reach_default:
        PV_center_wkg = _DEFAULT_PV
        st.write(f"Centerpoint default: **P/m = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {_N_for_default_rpm:.0f} RPM")
        if N_rpm_min is not None and N_rpm_max is not None:
            st.caption(f"Default of 0.2 W/kg is within reactor RPM range "
                       f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM).")
    else:
        # Fall back to average reactor speed
        _P_ctr = impeller_power(Np_val, rho, N_center, D_imp)
        PV_center_wkg = power_per_volume(_P_ctr, V_m3) / rho
        _ctr_rpm = N_center * 60
        st.write(f"Centerpoint from reactor average speed: **P/m = {PV_center_wkg:.4g} W/kg** "
                 f"at N = {N_center:.2f} rev/s ({_ctr_rpm:.0f} RPM)")
        st.caption(f"Default of 0.2 W/kg requires {_N_for_default_rpm:.0f} RPM, "
                   f"which is outside the reactor range "
                   f"({N_rpm_min:.1f} – {N_rpm_max:.1f} RPM). "
                   f"Using average speed instead.")

# Compute the three test speeds
# P/m ∝ N³  →  N ∝ (P/m)^(1/3)
N_center_calc = (PV_center_wkg * rho * V_m3 / (Np_val * rho * D_imp**5))**(1/3)
N_high_ideal = N_center_calc * 10**(1/3)   # ≈ 2.154×
N_low_ideal = N_center_calc * 0.1**(1/3)   # ≈ 0.464×

# Clamp to reactor min/max RPM bounds
N_low = N_low_ideal
N_high = N_high_ideal
_low_clamped = False
_high_clamped = False

if N_rps_min is not None and N_low < N_rps_min:
    N_low = N_rps_min
    _low_clamped = True
if N_rps_max is not None and N_high > N_rps_max:
    N_high = N_rps_max
    _high_clamped = True

# Build the conditions table
def _hydro_row(label, N):
    Re = reynolds_number(N, D_imp, rho, mu)
    Np = Np_val
    P = impeller_power(Np, rho, N, D_imp)
    eps = power_per_volume(P, V_m3)
    eps_kg = eps / rho                          # W/kg for Kolmogorov & micromixing
    u = tip_speed(N, D_imp)
    t_blend = blend_time_turbulent(Nq_val, V_m3, D_imp, N)
    t_micro = micromixing_time_engulfment(eps_kg, nu)
    eta = kolmogorov_length(nu, eps_kg)
    pv_wkg = eps_kg
    pv_rel = pv_wkg / PV_center_wkg if PV_center_wkg > 0 else 0
    return {
        "Condition": label,
        "Volume (L)": round(V_L, 3),
        # "N (rev/s)": round(N, 3),
        "N (RPM)": round(N * 60, 1),
        "P/m (W/kg)": round(pv_wkg, 4),
        "P/m rel. to center": f"{pv_rel:.2f}×",
        "P/V (W/L)": round(eps / 1000, 4),
        "Re": round(Re, 0),
        "Tip speed (m/s)": round(u, 3),
        "Blend time (s)": round(t_blend, 2),
        "t_E micro (s)": round(t_micro, 5),
        "η Kolmogorov (µm)": round(eta * 1e6, 1),
    }

_low_label = "Low  (min RPM)" if _low_clamped else "Low  (0.1× P/m)"
_high_label = "High (max RPM)" if _high_clamped else "High (10× P/m)"

t1_rows = [
    _hydro_row(_low_label, N_low),
    _hydro_row("Center (1× P/m)", N_center_calc),
    _hydro_row(_high_label, N_high),
]
t1_df = pd.DataFrame(t1_rows)
_t1_basic_cols = ["Condition", "Volume (L)", "N (RPM)", "P/m (W/kg)", "P/m rel. to center", "P/V (W/L)", "Tip speed (m/s)"]
st.dataframe(t1_df[_t1_basic_cols], width='stretch', hide_index=True)

with st.expander("Additional hydrodynamic parameters"):
    _t1_detail_cols = ["Condition", "Re", "Blend time (s)", "t_E micro (s)", "η Kolmogorov (µm)"]
    st.dataframe(t1_df[_t1_detail_cols], width='stretch', hide_index=True)

pv_ratio = t1_rows[2]["P/m (W/kg)"] / t1_rows[0]["P/m (W/kg)"] if t1_rows[0]["P/m (W/kg)"] > 0 else 0
_speed_ratio = N_high / N_low if N_low > 0 else 0
st.caption(f"P/m ratio (high/low) = **{pv_ratio:.1f}×**  •  Speed ratio (high/low) = **{_speed_ratio:.2f}×**")
st.caption(f"Power is computed with a fixed power number Np = **{Np_val:.2f}** "
           f"(from the reactor database, or 1.27 if unset). The Mixing Assessment "
           f"and Reactor Comparison pages instead use a Reynolds-dependent / ROM "
           f"power number, so their P/V values may differ for the same reactor.")

if _low_clamped or _high_clamped:
    _msgs = []
    if _low_clamped:
        _actual_low_rel = t1_rows[0]["P/m (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"Low speed clamped to reactor minimum ({N_rpm_min:.1f} RPM) → "
                     f"actual P/m = **{_actual_low_rel:.3f}×** centerpoint "
                     f"(target was 0.1×)")
    if _high_clamped:
        _actual_high_rel = t1_rows[2]["P/m (W/kg)"] / PV_center_wkg if PV_center_wkg > 0 else 0
        _msgs.append(f"High speed clamped to reactor maximum ({N_rpm_max:.1f} RPM) → "
                     f"actual P/m = **{_actual_high_rel:.3f}×** centerpoint "
                     f"(target was 10×)")
    st.warning("⚠️ **Speed limits applied:**\n\n" + "\n\n".join(_msgs) +
              "\n\nThe achievable P/V range is narrower than the ideal 100× span. "
              "Consider whether this is sufficient to draw conclusions.")

# ── Discrete speed adjustments for fed-batch (constant P/m) ──────────────
st.subheader("Discrete Speed Adjustments (Fed-Batch)")
st.markdown(
    "For fed-batch operations, the impeller speed can be stepped at "
    "discrete volume milestones to hold **P/m constant** as the working "
    "volume grows. Enable below to plan those setpoints for each of the "
    "three Test 1 conditions."
)

_use_adj = st.checkbox(
    "Add discrete speed adjustments at higher volumes",
    key=_bk("t1_use_adj"),
)

_adj_vols: list[float] = []

if _use_adj:
    _n_adj = int(st.number_input(
        "Number of additional volume adjustment points",
        min_value=1, max_value=10, value=1, step=1,
        key=_bk("t1_n_adj"),
    ))

    _vol_upper = V_L_max if V_L_max > V_L else V_L * 2.0
    _adj_cols = st.columns(min(_n_adj, 4))
    for _i in range(_n_adj):
        _col = _adj_cols[_i % len(_adj_cols)]
        _default_vol = V_L + (_i + 1) * (_vol_upper - V_L) / _n_adj
        _v = _col.number_input(
            f"Adjustment {_i + 1} volume (L)",
            min_value=float(V_L), value=float(_default_vol),
            step=0.001, format="%.3f",
            key=_bk(f"t1_adj_vol_{_i}"),
        )
        _adj_vols.append(_v)

    _pv_targets_adj = [
        (_low_label, t1_rows[0]["P/m (W/kg)"]),
        ("Center (1× P/m)", t1_rows[1]["P/m (W/kg)"]),
        (_high_label, t1_rows[2]["P/m (W/kg)"]),
    ]

    _vol_steps = [("Initial", V_L)] + [
        (f"Adj. {_i + 1}", _v) for _i, _v in enumerate(_adj_vols)
    ]

    _adj_rows = []
    _any_clamped = False
    for _step_label, _vol_L in _vol_steps:
        _row = {"Step": _step_label, "Volume (L)": round(_vol_L, 3)}
        _vol_m3_step = _vol_L / 1000.0
        for _cond_label, _pv_wkg in _pv_targets_adj:
            # P/m = Np · N³ · D⁵ / V  →  N = (P/m · V / (Np · D⁵))^(1/3)
            _N_req = (_pv_wkg * _vol_m3_step / (Np_val * D_imp**5))**(1/3)
            _N_rpm_req = _N_req * 60
            _flag = ""
            if N_rps_min is not None and _N_req < N_rps_min:
                _N_rpm_req = N_rpm_min
                _flag = " ⚠"
                _any_clamped = True
            elif N_rps_max is not None and _N_req > N_rps_max:
                _N_rpm_req = N_rpm_max
                _flag = " ⚠"
                _any_clamped = True
            _row[f"{_cond_label} (RPM)"] = f"{_N_rpm_req:.1f}{_flag}"
        _adj_rows.append(_row)

    _adj_df = pd.DataFrame(_adj_rows)
    st.dataframe(_adj_df, width='stretch', hide_index=True)
    st.caption(
        "Speeds chosen to hold each condition's P/m constant as volume "
        "increases. Set each value as a discrete setpoint when the working "
        "volume reaches the corresponding milestone."
    )
    if _any_clamped:
        st.warning(
            "⚠ One or more required RPM values fall outside the reactor "
            "range and were clamped to the reactor min/max. The target P/m "
            "cannot be maintained at those steps."
        )

with st.expander("Practical limits to consider"):
    st.markdown("""
    **Min speed:** N_js (solids), N_jd (liquids), flooding (gas)  
    **Max speed:** vortex / surface aeration, mechanical limits, splashing
    """)

# ── P/V iso-lines plot: Speed vs Fill Volume ─────────────────────────────
if V_L_min != V_L_max:
    import plotly.graph_objects as go

    st.subheader("Speed vs Fill Volume at Constant P/m")

    _vol_points = np.linspace(V_L_min, V_L_max, 50)
    _vol_m3_points = _vol_points / 1000.0

    _pv_targets = [
        ("0.1× P/m", PV_center_wkg * 0.1, "blue"),
        ("1× P/m (center)", PV_center_wkg, "green"),
        ("10× P/m", PV_center_wkg * 10, "red"),
    ]

    fig_t1 = go.Figure()

    for _label, _pv_wkg, _color in _pv_targets:
        # P/m (W/kg) = Np * N³ * D⁵ / V  →  N = (P/m * V / (Np * D⁵))^(1/3)
        _N_vals = (_pv_wkg * _vol_m3_points / (Np_val * D_imp**5))**(1/3)
        _rpm_vals = _N_vals * 60
        fig_t1.add_trace(go.Scatter(
            x=_vol_points, y=_rpm_vals, mode='lines',
            name=_label, line=dict(color=_color, width=2),
        ))

    # Mark the user-selected centerpoint volume
    _N_ctr_mark = (PV_center_wkg * V_m3 / (Np_val * D_imp**5))**(1/3) * 60
    fig_t1.add_trace(go.Scatter(
        x=[V_L], y=[_N_ctr_mark], mode='markers',
        name=f'Centerpoint ({V_L:.3f} L)',
        marker=dict(color='black', size=12, symbol='circle'),
    ))

    # Mark min and max volume at each P/m level
    for _label, _pv_wkg, _color in _pv_targets:
        _N_min_v = (_pv_wkg * (V_L_min / 1000) / (Np_val * D_imp**5))**(1/3) * 60
        _N_max_v = (_pv_wkg * (V_L_max / 1000) / (Np_val * D_imp**5))**(1/3) * 60
        fig_t1.add_trace(go.Scatter(
            x=[V_L_min, V_L_max], y=[_N_min_v, _N_max_v], mode='markers',
            name=f'{_label} (min/max vol)',
            marker=dict(color=_color, size=9, symbol='square'),
            showlegend=False,
        ))

    # Mark discrete fed-batch adjustment points on each iso-line
    if _adj_vols:
        for _idx, (_label, _pv_wkg, _color) in enumerate(_pv_targets):
            _adj_rpm = [
                (_pv_wkg * (_v / 1000) / (Np_val * D_imp**5))**(1/3) * 60
                for _v in _adj_vols
            ]
            fig_t1.add_trace(go.Scatter(
                x=_adj_vols, y=_adj_rpm, mode='markers',
                name='Fed-batch adjustments' if _idx == 0 else None,
                marker=dict(color=_color, size=11, symbol='diamond',
                            line=dict(color='black', width=1)),
                showlegend=(_idx == 0),
                hovertemplate='%{x:.3f} L → %{y:.1f} RPM<extra></extra>',
            ))

    # RPM bounds
    if N_rpm_min is not None:
        fig_t1.add_hline(y=N_rpm_min, line_dash="dash", line_color="gray",
                         annotation_text=f"Min RPM ({N_rpm_min:.0f})",
                         annotation_position="top left")
    if N_rpm_max is not None:
        fig_t1.add_hline(y=N_rpm_max, line_dash="dash", line_color="gray",
                         annotation_text=f"Max RPM ({N_rpm_max:.0f})",
                         annotation_position="bottom left")

    fig_t1.update_layout(
        title="Impeller Speed vs Fill Volume at Constant P/m",
        xaxis_title="Fill Volume (L)",
        yaxis_title="Impeller Speed (RPM)",
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.8)"),
        height=500,
    )
    st.plotly_chart(fig_t1, width='stretch')

# ── Response entry for Test 1 ───────────────────────────────────────────
st.subheader("Record Test 1 Responses")

# Sensitivity assessment helpers (used by all three tests)
_SENSITIVITY_THRESHOLD = 5.0  # % relative change (Sarafinas, 2018)

def _assess_sensitivity(resp_values, ref_value, threshold_pct=_SENSITIVITY_THRESHOLD):
    """Return (max_pct_change, is_sensitive, detail_text)."""
    if ref_value == 0:
        # Fall back to absolute range if center is zero
        rng = max(resp_values) - min(resp_values)
        return (rng, rng > 0, f"Absolute range = {rng:.4g} (center = 0; cannot compute relative change)")
    pct_changes = [abs(v - ref_value) / abs(ref_value) * 100 for v in resp_values]
    max_pct = max(pct_changes)
    sensitive = max_pct >= threshold_pct
    return (max_pct, sensitive, pct_changes)

def _bp_centerpoint_metrics():
    """Return the centerpoint hydrodynamic metrics dict for step PDF exports."""
    _eps = power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3)
    _eps_kg = _eps / rho
    return {
        "N (RPM)": round(N_center_calc * 60, 1),
        "P/m (W/kg)": round(PV_center_wkg, 4),
        "Re": round(reynolds_number(N_center_calc, D_imp, rho, mu), 0),
        "Tip speed (m/s)": round(tip_speed(N_center_calc, D_imp), 3),
        "Blend time (s)": round(blend_time_turbulent(Nq_val, V_m3, D_imp, N_center_calc), 2),
        "Micromix t_E (s)": round(micromixing_time_engulfment(_eps_kg, nu), 5),
        "Kolmogorov eta (um)": round(kolmogorov_length(nu, _eps_kg) * 1e6, 1),
    }

def _step_export_ui(step, snap_extra):
    """Render a 'generate + download' PDF control for a single protocol step."""
    _bytes_key = f"_p6_step{step}_pdf_bytes"
    _name_key = f"_p6_step{step}_pdf_name"
    if st.button(f"📥 Create Test {step} PDF", key=_bk(f"gen_step{step}_pdf")):
        with st.spinner("Generating PDF…"):
            try:
                from utils.report_builder import build_bourne_step_pdf, report_filename
                _snap = {
                    "step": step,
                    "project_name": project_name,
                    "step_number": step_number,
                    "unit_operation": unit_operation,
                    "reactor": reactor_name if not reactors.empty else "Manual entry",
                    "fluid": fluid_name,
                    "V_L": V_L,
                    "centerpoint_metrics": _bp_centerpoint_metrics(),
                }
                _snap.update(snap_extra)
                st.session_state[_bytes_key] = build_bourne_step_pdf(_snap)
                st.session_state[_name_key] = report_filename(
                    f"Bourne_Test{step}", reactor_name if not reactors.empty else "")
            except Exception as exc:
                st.error(f"PDF generation failed: {exc}")
    if _bytes_key in st.session_state:
        st.download_button(
            f"⬇️ Download Test {step} PDF",
            data=st.session_state[_bytes_key],
            file_name=st.session_state[_name_key],
            mime="application/pdf",
            key=_bk(f"dl_step{step}_pdf"),
        )

_t1_multi = st.checkbox(
    "Track multiple KPIs",
    key=_bk("t1_multi"),
    help=("Enable to record several response metrics and combine them into "
          "an overall mixing-sensitivity assessment by majority vote."),
)

if _t1_multi:
    _t1_n_kpi = int(st.number_input(
        "Number of KPIs",
        min_value=2, max_value=10, value=2, step=1,
        key=_bk("t1_n_kpi"),
    ))
else:
    _t1_n_kpi = 1

st.caption(
    "Each KPI can be tracked **Quantitatively** (numeric values, auto-assessed "
    "at a ≥ 5% relative change) or **Qualitatively** (observations with a "
    "manual sensitivity judgment)."
)

_t1_kpi_inputs = []  # [{mode, name, low, ctr, high, [judgment]}, ...]
for _k in range(_t1_n_kpi):
    if _t1_multi:
        if _k > 0:
            st.divider()
        st.markdown(f"**KPI {_k + 1}**")
    _kpi_mode = st.radio(
        "Capture mode" if not _t1_multi else f"KPI {_k + 1} capture mode",
        ["Quantitative", "Qualitative"],
        horizontal=True, key=_bk(f"t1_mode_{_k}"),
    )
    _mcol, _ucol = st.columns([2, 1])
    _kpi_metric = _mcol.selectbox(
        "Response metric" if not _t1_multi else f"KPI {_k + 1} metric",
        RESPONSE_METRICS,
        index=0,
        key=_bk(f"t1_resp_name_{_k}"),
    )
    _kpi_unit = _unit_selector(
        _ucol, _bk(f"t1_resp_unit_{_k}"),
        default_unit=_METRIC_DEFAULT_UNIT.get(_kpi_metric, "%"),
    )
    _kpi_name = _fmt_metric(_kpi_metric, _kpi_unit)
    if _kpi_mode == "Quantitative":
        _col_a, _col_b, _col_c = st.columns(3)
        with _col_a:
            _r_low = _col_a.number_input(f"{_low_label}", value=0.0, format="%.4g",
                                         key=_bk(f"t1_resp_low_{_k}"))
        with _col_b:
            _r_ctr = _col_b.number_input("Center (1× P/m)", value=0.0, format="%.4g",
                                         key=_bk(f"t1_resp_ctr_{_k}"))
        with _col_c:
            _r_high = _col_c.number_input(f"{_high_label}", value=0.0, format="%.4g",
                                          key=_bk(f"t1_resp_high_{_k}"))
        _t1_kpi_inputs.append({"mode": "quant", "metric": _kpi_metric,
                               "unit": _kpi_unit, "name": _kpi_name,
                               "low": _r_low, "ctr": _r_ctr, "high": _r_high})
    else:
        _col_a, _col_b, _col_c = st.columns(3)
        with _col_a:
            _q_low = st.text_area(f"{_low_label}", value="", height=100,
                                  key=_bk(f"t1_qual_low_{_k}"),
                                  placeholder="e.g. hazy, slow dissolution…")
        with _col_b:
            _q_ctr = st.text_area("Center (1× P/m)", value="", height=100,
                                  key=_bk(f"t1_qual_ctr_{_k}"),
                                  placeholder="e.g. clear solution…")
        with _col_c:
            _q_high = st.text_area(f"{_high_label}", value="", height=100,
                                   key=_bk(f"t1_qual_high_{_k}"),
                                   placeholder="e.g. fully dissolved, no haze…")
        _q_judgment = st.radio(
            "Does this response indicate mixing sensitivity?",
            ["Sensitive – mixing matters", "Not sensitive – mixing not critical"],
            key=_bk(f"t1_qual_judgment_{_k}"),
            horizontal=True,
        )
        _t1_kpi_inputs.append({"mode": "qual", "metric": _kpi_metric,
                               "unit": _kpi_unit, "name": _kpi_name,
                               "low": _q_low, "ctr": _q_ctr, "high": _q_high,
                               "judgment": _q_judgment})

# back-compat: first KPI base metric used by later sections
_t1_resp_name = _t1_kpi_inputs[0]["metric"]

# Button-triggered assessment
if st.button("📊 Assess Test 1 Responses", key=_bk("t1_assess")):
    def _kpi_has_data(kp):
        if kp["mode"] == "quant":
            return not (kp["low"] == 0.0 and kp["ctr"] == 0.0 and kp["high"] == 0.0)
        return any(str(v).strip() for v in (kp["low"], kp["ctr"], kp["high"]))

    if not any(_kpi_has_data(kp) for kp in _t1_kpi_inputs):
        st.warning("Enter response values before assessing.")
    else:
        _kpi_results = []
        for kp in _t1_kpi_inputs:
            if kp["mode"] == "quant":
                _mp, _sn, _pd = _assess_sensitivity(
                    [kp["low"], kp["ctr"], kp["high"]], kp["ctr"]
                )
                _kpi_results.append({
                    "name": kp["name"],
                    "qualitative": False,
                    "max_pct": _mp,
                    "sensitive": _sn,
                    "pct_detail": _pd,
                    "resp": [kp["low"], kp["ctr"], kp["high"]],
                })
            else:
                _sn = kp["judgment"].startswith("Sensitive")
                _kpi_results.append({
                    "name": kp["name"],
                    "qualitative": True,
                    "max_pct": 0.0,
                    "sensitive": _sn,
                    "pct_detail": None,
                    "resp": [kp["low"], kp["ctr"], kp["high"]],
                    "judgment": kp["judgment"],
                })
        _n_total = len(_kpi_results)
        _n_sensitive = sum(1 for r in _kpi_results if r["sensitive"])
        if _n_sensitive == 0:
            _status = "not_sensitive"
        elif _n_sensitive > _n_total / 2:
            _status = "sensitive"
        else:
            _status = "may_be_sensitive"

        _quant_pcts = [r["max_pct"] for r in _kpi_results if not r["qualitative"]]
        st.session_state[_bk("t1_assessed")] = {
            "qualitative": all(r["qualitative"] for r in _kpi_results),
            "kpi_results": _kpi_results,
            "n_total": _n_total,
            "n_sensitive": _n_sensitive,
            "status": _status,
            "labels": [_low_label, "Center (1× P/m)", _high_label],
            # back-compat single-KPI fields (first KPI / max across KPIs)
            "max_pct": max(_quant_pcts, default=0.0),
            "sensitive": _status != "not_sensitive",
            "resp_name": _kpi_results[0]["name"],
            "resp": _kpi_results[0]["resp"],
            "pct_detail": _kpi_results[0]["pct_detail"],
        }
        # Rerun so the top-of-page input export captures the new assessment.
        st.rerun()

# Display results if assessed
if _bk("t1_assessed") not in st.session_state:
    st.info("Enter response values, then click **Assess Test 1 Responses**.")
    st.stop()

_t1a = st.session_state[_bk("t1_assessed")]
t1_max_pct = _t1a["max_pct"]
t1_sensitive = _t1a["sensitive"]
t1_status = _t1a.get("status", "sensitive" if t1_sensitive else "not_sensitive")
t1_kpi_results = _t1a.get("kpi_results", [])
_t1_qualitative = _t1a.get("qualitative", False)
_t1_labels = _t1a["labels"]

if t1_kpi_results:
    _rows = []
    for r in t1_kpi_results:
        _is_qual = r.get("qualitative", False)
        if _is_qual:
            _vals = [str(v) if str(v).strip() else "—" for v in r["resp"]]
            _delta = "—"
        else:
            _vals = [f"{v:g}" for v in r["resp"]]
            _delta = f"{r['max_pct']:.1f}%"
        _rows.append({
            "KPI": r["name"],
            "Mode": "Qual" if _is_qual else "Quant",
            _t1_labels[0]: _vals[0],
            _t1_labels[1]: _vals[1],
            _t1_labels[2]: _vals[2],
            "Max Δ from center (%)": _delta,
            "Sensitive?": "Yes" if r["sensitive"] else "No",
        })
    st.dataframe(pd.DataFrame(_rows), width='stretch', hide_index=True)
    st.caption(
        f"{_t1a['n_sensitive']} / {_t1a['n_total']} KPIs sensitive  •  "
        f"Quantitative threshold = **{_SENSITIVITY_THRESHOLD:.0f}%** relative change"
    )
elif _t1_qualitative:
    # Legacy single-KPI qualitative display (assessment from an older session)
    _t1_res_df = pd.DataFrame({
        "Condition": _t1_labels,
        _t1a["resp_name"]: [str(v) if v else "—" for v in _t1a["resp"]],
    })
    st.dataframe(_t1_res_df, width='stretch', hide_index=True)
    st.caption("Qualitative observations — sensitivity judged manually below.")
else:
    # Legacy single-KPI numeric display (assessment from an older session)
    _pd_detail = _t1a.get("pct_detail")
    if isinstance(_pd_detail, list):
        _t1_res_df = pd.DataFrame({
            "Condition": _t1_labels,
            _t1a["resp_name"]: _t1a["resp"],
            "Δ from center (%)": [f"{p:.1f}%" for p in _pd_detail],
        })
    else:
        _t1_res_df = pd.DataFrame({
            "Condition": _t1_labels,
            _t1a["resp_name"]: _t1a["resp"],
        })
    st.dataframe(_t1_res_df, width='stretch', hide_index=True)
    st.caption(
        f"Maximum relative change from center = **{t1_max_pct:.1f}%**  •  "
        f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**"
    )

_t1_single = _t1a.get("n_total", 1) == 1
_t1_single_qual = _t1_single and t1_kpi_results and t1_kpi_results[0].get("qualitative", False)

if t1_status == "not_sensitive":
    if _t1_qualitative:
        st.success(
            "✅ **Mixing not critical.** Qualitative response(s) judged "
            "insensitive to impeller speed. Protocol complete."
        )
    else:
        st.success(
            "✅ **Mixing not critical.** No KPIs indicated sensitivity. "
            "Protocol complete."
        )
elif t1_status == "sensitive":
    if not _t1_single:
        st.warning(
            f"⚠️ **Mixing matters!** Majority of KPIs "
            f"({_t1a['n_sensitive']} / {_t1a['n_total']}) indicated "
            f"sensitivity. Proceed to **Test 2**."
        )
    elif _t1_single_qual or _t1_qualitative:
        st.warning(
            "⚠️ **Mixing matters!** Qualitative response judged sensitive "
            "to impeller speed. Proceed to **Test 2**."
        )
    else:
        st.warning(
            f"⚠️ **Mixing matters!** {t1_max_pct:.1f}% change "
            f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Proceed to **Test 2**."
        )
else:  # may_be_sensitive
    st.info(
        f"🟡 **Process may be mixing-sensitive.** "
        f"{_t1a['n_sensitive']} of {_t1a['n_total']} KPIs indicated "
        f"sensitivity (fewer than half). "
        f"Proceed to **Test 2** with caution."
    )

# Allow user to override the automatic assessment
t1_override = st.checkbox("Override automatic assessment", key=_bk("t1_override"))
if t1_override:
    t1_result_manual = st.radio(
        "Manual assessment:",
        ["Sensitive – mixing matters", "Not sensitive – mixing not critical"],
        key=_bk("t1_manual"),
        horizontal=True,
    )
    t1_sensitive = t1_result_manual.startswith("Sensitive")
    t1_status = "sensitive" if t1_sensitive else "not_sensitive"

# Export this step's results to PDF
_step_export_ui(1, {
    "t1_conditions": t1_rows,
    "t1_responses": st.session_state.get(_bk("t1_assessed")),
})

if not t1_sensitive:
    st.stop()

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 2 – TEST 2: Feed Rate / Time  (Meso vs Micro)
# ══════════════════════════════════════════════════════════════════════════
st.header("Test 2 - Feed Rate / Feed Time")
st.markdown("""
Vary **feed time only** (hold speed & location at centerpoint). Test a **9× flow-rate
range** (1/3× and 3× the centerpoint feed time).

Sensitive → mesomixing-controlled. Insensitive → micromixing-controlled.

""")

st.subheader("Suggested Feed Time Conditions")

t_feed_center = st.number_input(
    "Centerpoint feed time (min)",
    min_value=0.1, value=60.0, step=5.0,
    help="Fastest safe & practical plant-scale feed time, or your current lab feed time.",
    key=_bk("tfeed_ctr"),
)

t_feed_fast = t_feed_center / 3.0
t_feed_slow = t_feed_center * 3.0

feed_vol = st.number_input("Total feed volume (mL)", min_value=0.1, value=100.0, key=_bk("feed_vol"))

t2_conditions = pd.DataFrame([
    {
        "Condition": "Fast  (1/3× feed time)",
        "Feed time (min)": round(t_feed_fast, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_fast, 2),
        "Notes": "Highest mesomixing stress — more 'bad stuff' if meso-controlled",
    },
    {
        "Condition": "Center (1× feed time)",
        "Feed time (min)": round(t_feed_center, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_center, 2),
        "Notes": "Baseline condition",
    },
    {
        "Condition": "Slow  (3× feed time)",
        "Feed time (min)": round(t_feed_slow, 2),
        "Flow rate (mL/min)": round(feed_vol / t_feed_slow, 2),
        "Notes": "Less mesomixing stress — better result if meso-controlled",
    },
])
st.dataframe(t2_conditions, width='content', hide_index=True)
st.caption(f"Flow-rate ratio (fast/slow) = **9×**  •  Hold impeller speed at N = {N_center_calc:.2f} rev/s ({N_center_calc*60:.0f} RPM)")

with st.expander("Understanding mesomixing"):
    st.markdown(r"""
    $$t_{meso} \propto \frac{Q_{feed}}{k \cdot \varepsilon_{loc}}$$

    Higher feed rate → longer plume persistence → more concentration gradients.
    If response is **insensitive** to feed rate → **micromixing-controlled**
    (plume disperses fast; molecular-scale mixing is rate-limiting).
    """)

# ── Quantitative response entry for Test 2 ──────────────────────────────
st.subheader("Record Test 2 Responses")
st.caption("Use the same response metric as Test 1.")

_t2_mcol, _t2_ucol = st.columns([2, 1])
_t2_metric = _t2_mcol.selectbox("Response metric", RESPONSE_METRICS,
                                index=_metric_index(_t1_resp_name),
                                key=_bk("t2_resp_name"))
_t2_unit = _unit_selector(_t2_ucol, _bk("t2_resp_unit"),
                          default_unit=_METRIC_DEFAULT_UNIT.get(_t2_metric, "%"))
_t2_resp_name = _fmt_metric(_t2_metric, _t2_unit)
col_t2a, col_t2b, col_t2c = st.columns(3)
with col_t2a:
    t2_resp_fast = st.number_input("Fast (1/3× feed time)", value=0.0,
                                   format="%.4g", key=_bk("t2_resp_fast"))
with col_t2b:
    t2_resp_ctr = st.number_input("Center (1× feed time)", value=0.0,
                                  format="%.4g", key=_bk("t2_resp_ctr"))
with col_t2c:
    t2_resp_slow = st.number_input("Slow (3× feed time)", value=0.0,
                                   format="%.4g", key=_bk("t2_resp_slow"))

# Button-triggered assessment
if st.button("📊 Assess Test 2 Responses", key=_bk("t2_assess")):
    if t2_resp_ctr == 0.0 and t2_resp_fast == 0.0 and t2_resp_slow == 0.0:
        st.warning("Enter all three response values before assessing.")
    else:
        t2_max_pct, t2_sensitive, t2_pct_detail = _assess_sensitivity(
            [t2_resp_fast, t2_resp_ctr, t2_resp_slow], t2_resp_ctr
        )
        st.session_state[_bk("t2_assessed")] = {
            "max_pct": t2_max_pct, "sensitive": t2_sensitive,
            "pct_detail": t2_pct_detail,
            "resp": [t2_resp_fast, t2_resp_ctr, t2_resp_slow],
            "resp_name": _t2_resp_name,
        }
        # Rerun so the top-of-page input export captures the new assessment.
        st.rerun()

if _bk("t2_assessed") not in st.session_state:
    st.info("Enter all three response values, then click **Assess Test 2 Responses**.")
    st.stop()

_t2a = st.session_state[_bk("t2_assessed")]
t2_max_pct = _t2a["max_pct"]
t2_sensitive = _t2a["sensitive"]
t2_pct_detail = _t2a["pct_detail"]

if isinstance(t2_pct_detail, list):
    _t2_res_df = pd.DataFrame({
        "Condition": ["Fast (1/3×)", "Center (1×)", "Slow (3×)"],
        _t2a["resp_name"]: _t2a["resp"],
        "Δ from center (%)": [f"{p:.1f}%" for p in t2_pct_detail],
    })
else:
    _t2_res_df = pd.DataFrame({
        "Condition": ["Fast (1/3×)", "Center (1×)", "Slow (3×)"],
        _t2a["resp_name"]: _t2a["resp"],
    })
st.dataframe(_t2_res_df, width='stretch', hide_index=True)
st.caption(f"Maximum relative change from center = **{t2_max_pct:.1f}%**  •  "
           f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**")

if not t2_sensitive:
    st.info(
        f"🔬 **Micromixing controls.** {t2_max_pct:.1f}% change "
        f"(< {_SENSITIVITY_THRESHOLD:.0f}%). Scale-up: hold **local ε** constant."
    )
    _micro_conclusion = True
else:
    st.warning(
        f"⚠️ **Mesomixing matters!** {t2_max_pct:.1f}% change "
        f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Proceed to **Test 3**."
    )
    _micro_conclusion = False

# Allow user to override
t2_override = st.checkbox("Override automatic assessment", key=_bk("t2_override"))
if t2_override:
    t2_result_manual = st.radio(
        "Manual assessment:",
        ["Sensitive – mesomixing matters", "Not sensitive – micromixing controls"],
        key=_bk("t2_manual"),
        horizontal=True,
    )
    t2_sensitive = t2_result_manual.startswith("Sensitive")
    _micro_conclusion = not t2_sensitive

# Export this step's results to PDF
_step_export_ui(2, {
    "t2_conditions": {
        "N_RPM": round(N_center_calc * 60, 1),
        "feed_vol_mL": feed_vol,
        "feed_location": "Held constant (centerpoint)",
        "rows": [
            {"Condition": "Fast (1/3x feed time)", "Feed time (min)": round(t_feed_fast, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_fast, 2)},
            {"Condition": "Center (1x feed time)", "Feed time (min)": round(t_feed_center, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_center, 2)},
            {"Condition": "Slow (3x feed time)", "Feed time (min)": round(t_feed_slow, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_slow, 2)},
        ],
    },
    "t2_responses": st.session_state.get(_bk("t2_assessed")),
})

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 3 – TEST 3: Feed Location  (Meso vs Macro)
# ══════════════════════════════════════════════════════════════════════════

# Defaults when Test 3 is not reached (micromixing controls)
t3_sensitive = False
t3_max_pct = 0.0
_macro_conclusion = False
_meso_conclusion = False

# Compute ε_avg at centerpoint (used in Test 3 and save section)
eps_avg = power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3)
eps_avg_kg = eps_avg / rho  # W/kg for micromixing calculations

if t2_sensitive:
    st.header("Test 3 - Feed Location")
    st.markdown("""
    Vary **feed location only** (hold speed & feed time at centerpoint). Move
    from surface → impeller zone. This changes **local ε** without affecting
    macromixing.

    | Feed location | ε_loc / ε_avg |
    |---|---|
    | Surface | ≈ 0.1 |
    | Sub-surface (mid-tank) | ≈ 1 |
    | Impeller zone | ≈ 2 – 5 |

    Sensitive → **mesomixing**. Insensitive → **macromixing**.
    """)

    t3_locs = pd.DataFrame([
        {
            "Feed Location": "Surface",
            "ε_loc / ε_avg": 0.1,
            "ε_loc (W/m³)": round(0.1 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(0.1 * eps_avg_kg, nu), 5),
        },
        {
            "Feed Location": "Sub-surface (mid-tank)",
            "ε_loc / ε_avg": 1.0,
            "ε_loc (W/m³)": round(1.0 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(1.0 * eps_avg_kg, nu), 5),
        },
        {
            "Feed Location": "Impeller zone",
            "ε_loc / ε_avg": 3.0,
            "ε_loc (W/m³)": round(3.0 * eps_avg, 2),
            "t_E micro (s)": round(micromixing_time_engulfment(3.0 * eps_avg_kg, nu), 5),
        },
    ])
    st.dataframe(t3_locs, width='stretch', hide_index=True)
    st.caption("Local ε estimated at centerpoint impeller speed. Actual values depend on impeller type and geometry.")

    with st.expander("Practical considerations for feed location"):
        st.markdown("""
        - Subsurface feeds → higher local ε but risk of **back-mixing** and plugging.
        - Feed velocity must exceed local bulk velocity (Jo et al., 1994).
        - *"Subsurface addition — prove it is not needed!"* — Ed Paul
        """)

    # ── Quantitative response entry for Test 3 ──────────────────────────
    st.subheader("Record Test 3 Responses")
    st.caption("Use the same response metric as Tests 1 & 2.")

    _t3_mcol, _t3_ucol = st.columns([2, 1])
    _t3_metric = _t3_mcol.selectbox("Response metric", RESPONSE_METRICS,
                                    index=_metric_index(_t2_metric),
                                    key=_bk("t3_resp_name"))
    _t3_unit = _unit_selector(_t3_ucol, _bk("t3_resp_unit"),
                              default_unit=_METRIC_DEFAULT_UNIT.get(_t3_metric, "%"))
    _t3_resp_name = _fmt_metric(_t3_metric, _t3_unit)
    col_t3a, col_t3b, col_t3c = st.columns(3)
    with col_t3a:
        t3_resp_surf = st.number_input("Surface feed", value=0.0, format="%.4g",
                                       key=_bk("t3_resp_surf"))
    with col_t3b:
        t3_resp_mid = st.number_input("Sub-surface (mid-tank)", value=0.0,
                                      format="%.4g", key=_bk("t3_resp_mid"))
    with col_t3c:
        t3_resp_imp = st.number_input("Impeller zone", value=0.0, format="%.4g",
                                      key=_bk("t3_resp_imp"))

    # Button-triggered assessment
    if st.button("📊 Assess Test 3 Responses", key=_bk("t3_assess")):
        if t3_resp_mid == 0.0 and t3_resp_surf == 0.0 and t3_resp_imp == 0.0:
            st.warning("Enter all three response values before assessing.")
        else:
            t3_max_pct, t3_sensitive, t3_pct_detail = _assess_sensitivity(
                [t3_resp_surf, t3_resp_mid, t3_resp_imp], t3_resp_mid
            )
            st.session_state[_bk("t3_assessed")] = {
                "max_pct": t3_max_pct, "sensitive": t3_sensitive,
                "pct_detail": t3_pct_detail,
                "resp": [t3_resp_surf, t3_resp_mid, t3_resp_imp],
                "resp_name": _t3_resp_name,
            }
            # Rerun so the top-of-page input export captures the new assessment.
            st.rerun()

    if _bk("t3_assessed") not in st.session_state:
        st.info("Enter all three response values, then click **Assess Test 3 Responses**.")
        st.stop()

    _t3a = st.session_state[_bk("t3_assessed")]
    t3_max_pct = _t3a["max_pct"]
    t3_sensitive = _t3a["sensitive"]
    t3_pct_detail = _t3a["pct_detail"]

    if isinstance(t3_pct_detail, list):
        _t3_res_df = pd.DataFrame({
            "Feed Location": ["Surface", "Sub-surface (mid)", "Impeller zone"],
            _t3a["resp_name"]: _t3a["resp"],
            "Δ from mid-tank (%)": [f"{p:.1f}%" for p in t3_pct_detail],
        })
    else:
        _t3_res_df = pd.DataFrame({
            "Feed Location": ["Surface", "Sub-surface (mid)", "Impeller zone"],
            _t3a["resp_name"]: _t3a["resp"],
        })
    st.dataframe(_t3_res_df, width='stretch', hide_index=True)
    st.caption(f"Maximum relative change from mid-tank = **{t3_max_pct:.1f}%**  •  "
               f"Sensitivity threshold = **{_SENSITIVITY_THRESHOLD:.0f}%**")

    if not t3_sensitive:
        st.info(
            f"🌀 **Macromixing controls.** {t3_max_pct:.1f}% change "
            f"(< {_SENSITIVITY_THRESHOLD:.0f}%). Bulk blending is rate-limiting."
        )
        st.caption("Scale-up: maintain short blend times (hydrofoils, static mixers, multiple impellers).")
        _macro_conclusion = True
        _meso_conclusion = False
    else:
        st.success(
            f"📐 **Mesomixing controls.** {t3_max_pct:.1f}% change "
            f"(≥ {_SENSITIVITY_THRESHOLD:.0f}%). Feed-plume dispersion is rate-limiting."
        )
        st.caption("Scale-up: constant impeller speed, extended feed time, or multiple feed points.")
        _macro_conclusion = False
        _meso_conclusion = True

    # Allow user to override
    t3_override = st.checkbox("Override automatic assessment", key=_bk("t3_override"))
    if t3_override:
        t3_result_manual = st.radio(
            "Manual assessment:",
            ["Sensitive – mesomixing controls", "Not sensitive – macromixing controls"],
            key=_bk("t3_manual"),
            horizontal=True,
        )
        t3_sensitive = t3_result_manual.startswith("Sensitive")
        _macro_conclusion = not t3_sensitive
        _meso_conclusion = t3_sensitive

    # Export this step's results to PDF
    _step_export_ui(3, {
        "t3_conditions": {
            "N_RPM": round(N_center_calc * 60, 1),
            "feed_time_min": round(t_feed_center, 2),
            "eps_avg_W_m3": round(eps_avg, 2),
            "rows": [
                {"Feed Location": "Surface", "eps_loc/eps_avg": 0.1, "eps_loc (W/m3)": round(0.1 * eps_avg, 2)},
                {"Feed Location": "Sub-surface (mid-tank)", "eps_loc/eps_avg": 1.0, "eps_loc (W/m3)": round(1.0 * eps_avg, 2)},
                {"Feed Location": "Impeller zone", "eps_loc/eps_avg": 3.0, "eps_loc (W/m3)": round(3.0 * eps_avg, 2)},
            ],
        },
        "t3_responses": st.session_state.get(_bk("t3_assessed")),
    })

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 4 – Confirmatory Experiments
# ══════════════════════════════════════════════════════════════════════════
st.header("4 · Confirmatory Experiments (Optional)")
st.caption("Optional tests to validate the conclusion from Tests 2–3.")

tab_fp, tab_visc = st.tabs(["A – Number of Feed Points", "B – Viscosity Change"])

with tab_fp:
    st.markdown("""
    ### A — Number of Feed Points
    Increase feed points (1 → 2–3) at same speed, feed time, and location.
    More feed points → lower local feed velocity at each point.

    | Outcome | Conclusion |
    |---|---|
    | Insensitive to # feed points | Micromixing controlled |
    | Changes with # feed points | Mesomixing controlled |
    """)

    fp_result = st.radio(
        "Was the process response sensitive to the number of feed points?",
        ["— Select —", "No – insensitive (supports micromixing control)",
         "Yes – changed (supports mesomixing control)"],
        key=_bk("fp_result"),
    )
    if fp_result.startswith("No"):
        st.info("✅ Confirms **micromixing** control.")
    elif fp_result.startswith("Yes"):
        st.info("✅ Confirms **mesomixing** control.")

with tab_visc:
    st.markdown(r"""
    ### B — Viscosity Change
    Change bulk viscosity (e.g. co-solvent ratio) while holding recipe constant.
    $t_E = 17.3 (\nu / \varepsilon)^{1/2}$

    | Outcome | Conclusion |
    |---|---|
    | Response **changes** with viscosity | Micromixing controlled |
    | Response **insensitive** | Not micromixing controlled |

    ⚠️ Large viscosity changes may shift the flow regime — keep Re turbulent.
    """)

    visc_result = st.radio(
        "Was the process response sensitive to viscosity change?",
        ["— Select —", "Yes – changed (supports micromixing control)",
         "No – insensitive (does not support micromixing control)"],
        key=_bk("visc_result"),
    )
    if visc_result.startswith("Yes"):
        st.info("✅ Supports **micromixing** control.")
    elif visc_result.startswith("No"):
        st.info("✅ Suggests micromixing is **not** the controlling mechanism.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 5 – Overall Conclusion & Decision Tree
# ══════════════════════════════════════════════════════════════════════════
st.header("5 · Protocol Summary & Conclusion")

# Build the conclusion from all recorded answers
conclusions = []
scaleup_notes = []

# Test 1 always passed if we reached here
if t1_status == "may_be_sensitive":
    conclusions.append((
        "Test 1 – Impeller Speed",
        f"**May be sensitive** ({_t1a.get('n_sensitive', '?')} / {_t1a.get('n_total', '?')} KPIs "
        f"≥ {_SENSITIVITY_THRESHOLD:.0f}%) → Mixing may matter",
        "🟡",
    ))
elif _t1a.get("qualitative", False):
    conclusions.append((
        "Test 1 – Impeller Speed",
        "**Sensitive** (qualitative judgment) → Mixing matters",
        "⚠️",
    ))
else:
    conclusions.append((
        "Test 1 – Impeller Speed",
        f"**Sensitive** (max {t1_max_pct:.1f}% change) → Mixing matters",
        "⚠️",
    ))

# Test 2
if not t2_sensitive:
    conclusions.append(("Test 2 – Feed Rate",
                        f"**Insensitive** ({t2_max_pct:.1f}% change) → Micromixing-controlled", "🔬"))
    scaleup_notes.append("Hold **local ε** constant at the feed point on scale-up.")
else:
    conclusions.append(("Test 2 – Feed Rate",
                        f"**Sensitive** ({t2_max_pct:.1f}% change) → Not purely micromixing-controlled", "⚠️"))

# Test 3
if t2_sensitive:  # only reached Test 3 if Test 2 was sensitive
    if t3_sensitive:
        conclusions.append(("Test 3 – Feed Location",
                            f"**Sensitive** ({t3_max_pct:.1f}% change) → Mesomixing-controlled", "📐"))
        scaleup_notes.append("Control feed plume dispersion: constant impeller speed, extended feed time, or multiple feed points.")
    else:
        conclusions.append(("Test 3 – Feed Location",
                            f"**Insensitive** ({t3_max_pct:.1f}% change) → Macromixing-controlled", "🌀"))
        scaleup_notes.append("Maintain short blend times on scale-up (high-efficiency impellers, static mixers).")

# Determine the dominant mixing limitation
if not t2_sensitive:
    dominant = "Micromixing"
    dominant_icon = "🔬"
    dominant_color = "blue"
elif t3_sensitive:
    dominant = "Mesomixing"
    dominant_icon = "📐"
    dominant_color = "orange"
else:
    dominant = "Macromixing"
    dominant_icon = "🌀"
    dominant_color = "red"

# Summary table
st.subheader("Test Results")
for test_name, result, icon in conclusions:
    st.markdown(f"{icon} **{test_name}:** {result}")

st.subheader("Dominant Mixing Limitation")

if dominant == "Micromixing":
    st.success(f"""
    {dominant_icon} **{dominant}** — Molecular-scale engulfment is rate-limiting.

    **Scale-up:** Maintain constant **local ε** at the feed point.
    """)
elif dominant == "Mesomixing":
    st.warning(f"""
    {dominant_icon} **{dominant}** — Feed-plume disintegration is rate-limiting.

    **Scale-up:** Constant impeller speed, extended feed time, or multiple feed points.
    """)
elif dominant == "Macromixing":
    st.error(f"""
    {dominant_icon} **{dominant}** — Bulk blending / circulation is rate-limiting.

    **Scale-up:** Reduce blend time (hydrofoil impellers, multiple impellers, static mixers).
    """)
else:
    st.info("Complete Tests 2 and 3 to determine the controlling mixing scale.")

# Scale-up notes
if scaleup_notes:
    st.subheader("Key Scale-Up Actions")
    for note in scaleup_notes:
        st.markdown(f"- {note}")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 6 – Save Protocol Record
# ══════════════════════════════════════════════════════════════════════════
st.header("6 · Save Protocol Record")

if st.button("📌 Save Bourne Protocol result to Recorded Results", key=_bk("save")):
    result_row = {
        "reactor": reactor_name if not reactors.empty else "Manual entry",
        "fluid": fluid_name,
        "reaction": "Bourne Protocol",
        "Re": round(reynolds_number(N_center_calc, D_imp, rho, mu), 0),
        "P/V (W/L)": round(power_per_volume(impeller_power(Np_val, rho, N_center_calc, D_imp), V_m3) / 1000, 4),
        "Tip speed (m/s)": round(tip_speed(N_center_calc, D_imp), 3),
        "Blend time (s)": round(blend_time_turbulent(Nq_val, V_m3, D_imp, N_center_calc), 2),
        "Micromix t_E (s)": round(micromixing_time_engulfment(eps_avg / rho, nu), 5),
        "Kolmogorov η (µm)": round(kolmogorov_length(nu, eps_avg / rho) * 1e6, 1),
        "t_rxn (s)": "N/A (protocol)",
        "Da_macro": "N/A",
        "Da_micro": "N/A",
        "Assessment": f"Bourne Protocol → {dominant}-controlled",
    }
    if "recorded_results" not in st.session_state:
        st.session_state.recorded_results = pd.DataFrame()
    st.session_state.recorded_results = pd.concat(
        [st.session_state.recorded_results, pd.DataFrame([result_row])],
        ignore_index=True,
    )
    results_csv = DATA_DIR / "recorded_results.csv"
    st.session_state.recorded_results.to_csv(results_csv, index=False)
    st.success("Protocol result saved to **Recorded Results**.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 7 – Export Report
# ══════════════════════════════════════════════════════════════════════════
st.header("7 · Export Report")

if st.button("📥 Export PDF Report", type="primary", key=_bk("export_pdf")):
    with st.spinner("Generating PDF…"):
        try:
            from utils.report_builder import build_bourne_protocol_pdf, report_filename

            _bp_snap = {
                "reactor": reactor_name if not reactors.empty else "Manual entry",
                "fluid": fluid_name,
                "V_L": V_L,
                "project_name": project_name,
                "step_number": step_number,
                "unit_operation": unit_operation,
                "dominant": dominant,
                "conclusions": conclusions,
                "scaleup_notes": scaleup_notes,
                "t1_conditions": t1_rows,
                "t1_responses": st.session_state.get(_bk("t1_assessed")),
                "t2_conditions": {
                    "N_RPM": round(N_center_calc * 60, 1),
                    "feed_vol_mL": feed_vol,
                    "feed_location": "Held constant (centerpoint)",
                    "rows": [
                        {"Condition": "Fast (1/3x feed time)", "Feed time (min)": round(t_feed_fast, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_fast, 2)},
                        {"Condition": "Center (1x feed time)", "Feed time (min)": round(t_feed_center, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_center, 2)},
                        {"Condition": "Slow (3x feed time)", "Feed time (min)": round(t_feed_slow, 2), "Flow rate (mL/min)": round(feed_vol / t_feed_slow, 2)},
                    ],
                },
                "t2_responses": st.session_state.get(_bk("t2_assessed")),
                "t3_conditions": {
                    "N_RPM": round(N_center_calc * 60, 1),
                    "feed_time_min": round(t_feed_center, 2),
                    "eps_avg_W_m3": round(eps_avg, 2),
                    "rows": [
                        {"Feed Location": "Surface", "eps_loc/eps_avg": 0.1, "eps_loc (W/m3)": round(0.1 * eps_avg, 2)},
                        {"Feed Location": "Sub-surface (mid-tank)", "eps_loc/eps_avg": 1.0, "eps_loc (W/m3)": round(1.0 * eps_avg, 2)},
                        {"Feed Location": "Impeller zone", "eps_loc/eps_avg": 3.0, "eps_loc (W/m3)": round(3.0 * eps_avg, 2)},
                    ],
                },
                "t3_responses": st.session_state.get(_bk("t3_assessed")),
                "centerpoint_metrics": {
                    "N (RPM)": round(N_center_calc * 60, 1),
                    "P/m (W/kg)": round(PV_center_wkg, 4),
                    "Re": round(reynolds_number(N_center_calc, D_imp, rho, mu), 0),
                    "Tip speed (m/s)": round(tip_speed(N_center_calc, D_imp), 3),
                    "Blend time (s)": round(blend_time_turbulent(Nq_val, V_m3, D_imp, N_center_calc), 2),
                    "Micromix t_E (s)": round(micromixing_time_engulfment(eps_avg_kg, nu), 5),
                    "Kolmogorov eta (um)": round(kolmogorov_length(nu, eps_avg_kg) * 1e6, 1),
                },
            }
            _pdf_bytes = build_bourne_protocol_pdf(_bp_snap)
            st.session_state["_p6_pdf_bytes"] = _pdf_bytes
            st.session_state["_p6_pdf_name"] = report_filename(
                "Bourne_Protocol", reactor_name if not reactors.empty else ""
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_p6_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_p6_pdf_bytes"],
        file_name=st.session_state["_p6_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")

st.divider()

# ══════════════════════════════════════════════════════════════════════════
# SECTION 7 – References
# ══════════════════════════════════════════════════════════════════════════
st.header("References")
st.markdown("""
1. Bourne, J.R. (2003). *Mixing and the Selectivity of Chemical Reactions.*
   Org. Process Res. Dev., 7(4), 471–508.
2. Paul, E.L., Atiemo-Obeng, V.A. & Kresta, S.M. (2004). *Handbook of
   Industrial Mixing.* Wiley-Interscience.
3. Baldyga, J. & Bourne, J.R. (1999). *Turbulent Mixing and Chemical
   Reactions.* Wiley.
4. Jo, M.C., Penney, W.R. & Fasano, J.B. (1994). Back-mixing into reactor
   feed pipes caused by turbulence in an agitated vessel. AIChE Symp. Ser.,
   90(299).
""")
