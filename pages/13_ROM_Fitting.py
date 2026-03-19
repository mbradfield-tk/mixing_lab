"""
Admin – ROM / Experimental Correlation Fitting
===============================================
Upload measured or CFD data, fit a reduced-order model, and register
it for use throughout Mixing Lab.

Password-protected (same credentials as Admin Tools).
"""

from __future__ import annotations

import json
import pathlib
import hmac
import textwrap

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.optimize import curve_fit

from utils.rom_registry import (
    Correlation,
    SUPPORTED_PARAMS,
    PARAM_DISPLAY,
    register,
    get_all_correlations,
    get_all_registered_reactors,
)
from utils.rom_templates import (
    TEMPLATES,
    CorrTemplate,
    compatible_templates,
    template_by_id,
    templates_for_param,
)

# ── Paths ────────────────────────────────────────────────────────────────
DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
SAVED_CORR_FILE = DATA_DIR / "fitted_correlations.json"

st.title("🔧 ROM / Experimental Correlation Fitting")

# ── Authentication gate ──────────────────────────────────────────────────
_ADMIN_USER = "admin"
_ADMIN_PASS = "admin_tak_2026"


def _check_password() -> bool:
    if st.session_state.get("admin_authenticated"):
        return True
    with st.form("rom_fit_login"):
        st.subheader("🔒 Admin login required")
        user = st.text_input("Username")
        pwd = st.text_input("Password", type="password")
        submitted = st.form_submit_button("Log in", type="primary")
    if submitted:
        if (hmac.compare_digest(user, _ADMIN_USER)
                and hmac.compare_digest(pwd, _ADMIN_PASS)):
            st.session_state["admin_authenticated"] = True
            st.rerun()
        else:
            st.error("Invalid username or password.")
    return False


if not _check_password():
    st.stop()

# ═══════════════════════════════════════════════════════════════════════════
#  Persistence helpers
# ═══════════════════════════════════════════════════════════════════════════

def _load_saved() -> list[dict]:
    """Load previously saved correlations from JSON."""
    if SAVED_CORR_FILE.exists():
        with open(SAVED_CORR_FILE, "r") as f:
            return json.load(f)
    return []


def _save_all(entries: list[dict]) -> None:
    """Write fitted correlations to JSON."""
    with open(SAVED_CORR_FILE, "w") as f:
        json.dump(entries, f, indent=2)


def _append_and_save(entry: dict) -> None:
    """Add one fitted correlation and persist."""
    entries = _load_saved()
    entries.append(entry)
    _save_all(entries)


def _register_from_entry(entry: dict) -> None:
    """Re-register a saved correlation into the live ROM registry."""
    tmpl = template_by_id(entry["template_id"])
    if tmpl is None:
        return
    coeffs = np.array(entry["coeffs"])
    func = tmpl.build_func(coeffs)
    latex = tmpl.latex_filled(coeffs)

    register(
        entry["reactor_name"],
        Correlation(
            name=entry["name"],
            param=entry["param"],
            corr_type=entry["corr_type"],
            func=func,
            latex=latex,
            source=entry.get("source", "Fitted in Mixing Lab"),
            description=entry.get("description", ""),
            input_params=tmpl.column_labels,
        ),
    )


# ── Auto-load saved correlations on first run ────────────────────────────
if "rom_fit_loaded" not in st.session_state:
    for entry in _load_saved():
        _register_from_entry(entry)
    st.session_state["rom_fit_loaded"] = True


# ═══════════════════════════════════════════════════════════════════════════
#  Fitting engine
# ═══════════════════════════════════════════════════════════════════════════

def _fit_template(
    tmpl: CorrTemplate,
    df: pd.DataFrame,
    target_col: str,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Fit *tmpl* to *df* and return (coeffs, std_errs, r2, y_pred).

    Uses ``scipy.optimize.curve_fit`` with log-transform when
    ``tmpl.log_transform`` is True (linearises power-laws and avoids
    scale-sensitivity issues).
    """
    y = df[target_col].values.astype(float)
    X = {col: df[col].values.astype(float) for col in tmpl.required_columns}

    # Build wrapper for curve_fit (expects f(x, *params) → y)
    if tmpl.log_transform:
        # Filter out non-positive values for log transform
        mask = y > 0
        for col in tmpl.required_columns:
            mask &= X[col] > 0
        y_fit = np.log(y[mask])
        X_fit = {col: vals[mask] for col, vals in X.items()}

        def _wrapper(_dummy, *params):
            c = np.array(params)
            pred = tmpl.model(c, X_fit)
            return np.log(np.clip(pred, 1e-30, None))

        x_dummy = np.arange(mask.sum())
    else:
        y_fit = y.copy()
        X_fit = X

        def _wrapper(_dummy, *params):
            c = np.array(params)
            return tmpl.model(c, X_fit)

        x_dummy = np.arange(len(y))

    p0 = tmpl.initial_guess
    lb = tmpl.bounds_lower if tmpl.bounds_lower else [-np.inf] * tmpl.n_coeffs
    ub = tmpl.bounds_upper if tmpl.bounds_upper else [np.inf] * tmpl.n_coeffs

    popt, pcov = curve_fit(
        _wrapper, x_dummy, y_fit,
        p0=p0, bounds=(lb, ub),
        maxfev=10_000,
    )
    perr = np.sqrt(np.diag(pcov))

    # Predictions on full data for R² and parity plot
    y_pred = tmpl.model(popt, X).flatten()
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return popt, perr, r2, y_pred


# ═══════════════════════════════════════════════════════════════════════════
#  Column-mapping helper
# ═══════════════════════════════════════════════════════════════════════════

# Known aliases for template columns → common data-file column names
_COLUMN_ALIASES: dict[str, list[str]] = {
    "P_V": ["P/V", "P_V", "PV", "power_per_volume", "power_per_vol", "P/V (W/m³)", "P/V (W/m3)"],
    "v_s": ["v_s", "vs", "superficial_velocity", "superficial_gas_velocity", "u_g"],
    "eps_kg": ["eps_kg", "epsilon", "eps", "dissipation", "energy_dissipation", "ε", "epsilon_avg"],
    "mu": ["mu", "viscosity", "dynamic_viscosity", "μ"],
    "Re": ["Re", "reynolds", "Reynolds", "reynolds_number"],
    "V_L": ["V_L", "volume", "vol", "V", "volume_L", "fill_volume"],
    "ND2": ["ND2", "N_D2", "ND^2"],
    "N": ["N", "speed", "rps", "impeller_speed", "stir_speed_rps"],
    "D_imp": ["D_imp", "D", "impeller_diameter", "d_imp"],
    "nu_over_eps": ["nu_over_eps", "nu/eps", "ν/ε"],
    "nu": ["nu", "kinematic_viscosity", "ν"],
}


def _auto_map_columns(
    data_cols: list[str], required: list[str]
) -> dict[str, str | None]:
    """Try to auto-map uploaded column names to template column names."""
    mapping: dict[str, str | None] = {}
    for req in required:
        found = None
        aliases = _COLUMN_ALIASES.get(req, [req])
        for alias in aliases:
            for dc in data_cols:
                if dc.strip().lower() == alias.strip().lower():
                    found = dc
                    break
            if found:
                break
        mapping[req] = found
    return mapping


# ═══════════════════════════════════════════════════════════════════════════
#  UI
# ═══════════════════════════════════════════════════════════════════════════

st.markdown("""
Upload measured or CFD-generated data, select a model form, fit
coefficients, and register the resulting correlation for use in all
Mixing Lab analysis pages.

---
""")

tab_fit, tab_manage = st.tabs(["📈 Fit New Correlation", "📋 Manage Saved Correlations"])

# ═══════════════════════════════════════════════════════════════════════════
#  Tab 1 – Fit New Correlation
# ═══════════════════════════════════════════════════════════════════════════
with tab_fit:

    # ── Step 1: Select target parameter ──────────────────────────────────
    st.header("① Select target parameter")

    param = st.selectbox(
        "Parameter to fit",
        options=SUPPORTED_PARAMS,
        format_func=lambda p: PARAM_DISPLAY.get(p, p),
    )

    # ── Step 2: Upload data ──────────────────────────────────────────────
    st.header("② Upload data")

    st.markdown(
        "Upload a **CSV** file with at least a **target column** (measured "
        "values of the parameter) and columns for the independent variables "
        "required by the chosen correlation form."
    )

    uploaded = st.file_uploader("CSV file", type=["csv"])

    _fit_ready = False  # gates downstream sections

    if uploaded is None:
        st.info("Upload a CSV file to continue.")
    else:
        df_raw = pd.read_csv(uploaded)
        st.subheader("Uploaded data preview")
        st.dataframe(df_raw.head(20), use_container_width=True)

        data_cols = list(df_raw.columns)
        numeric_cols = [c for c in data_cols if pd.api.types.is_numeric_dtype(df_raw[c])]

        # ── Step 3: Choose target column ─────────────────────────────────
        st.header("③ Map target column")

        target_col = st.selectbox(
            "Column containing measured values of the target parameter",
            options=numeric_cols,
        )

        # ── Step 4: Select correlation form ──────────────────────────────
        st.header("④ Select correlation form")

        param_templates = templates_for_param(param)

        if not param_templates:
            st.warning(f"No correlation templates defined for **{PARAM_DISPLAY.get(param, param)}**.")
        else:
            # ── Nomenclature table ───────────────────────────────────────
            _symbol_map: dict[str, str] = {}
            for _t in param_templates:
                _symbol_map.update(_t.column_labels)

            with st.expander("📖 Nomenclature – variable symbols used in the correlation forms", expanded=False):
                _nomen_rows = [{"Symbol": sym, "Description": desc} for sym, desc in _symbol_map.items()]
                st.table(pd.DataFrame(_nomen_rows))

            st.markdown("Available model forms for this parameter:")

            for _tmpl_item in param_templates:
                st.markdown(f"- **{_tmpl_item.name}** — requires: {', '.join(_tmpl_item.required_columns)}")

            selected_tmpl_id = st.selectbox(
                "Correlation form",
                options=[t.id for t in param_templates],
                format_func=lambda tid: next(t.name for t in param_templates if t.id == tid),
            )
            tmpl = next(t for t in param_templates if t.id == selected_tmpl_id)

            # ── Step 5: Map input columns ────────────────────────────────
            st.header("⑤ Map input columns")

            st.markdown(
                "Map each required input variable to a column in your uploaded data. "
                "If a column is not present, you may need to compute it before uploading."
            )

            auto_map = _auto_map_columns(data_cols, tmpl.required_columns)
            col_mapping: dict[str, str] = {}

            for req_col in tmpl.required_columns:
                label = tmpl.column_labels.get(req_col, req_col)
                default_idx = 0
                if auto_map[req_col] and auto_map[req_col] in numeric_cols:
                    default_idx = numeric_cols.index(auto_map[req_col])
                chosen = st.selectbox(
                    f"**{req_col}** — {label}",
                    options=numeric_cols,
                    index=default_idx,
                    key=f"map_{req_col}",
                )
                col_mapping[req_col] = chosen

            # Rename columns in a working copy
            df_work = df_raw.copy()
            reverse_map = {v: k for k, v in col_mapping.items()}
            df_work = df_work.rename(columns=reverse_map)

            missing = [c for c in tmpl.required_columns if c not in df_work.columns]
            if missing:
                st.error(f"Missing mapped columns: {missing}")
            else:
                if target_col not in df_work.columns:
                    df_work[target_col] = df_raw[target_col].values

                relevant = tmpl.required_columns + [target_col]
                df_clean = df_work.dropna(subset=relevant)
                n_dropped = len(df_work) - len(df_clean)
                if n_dropped > 0:
                    st.info(f"Dropped {n_dropped} rows with missing values.")

                if len(df_clean) < tmpl.n_coeffs + 1:
                    st.error(
                        f"Need at least {tmpl.n_coeffs + 1} data points for "
                        f"{tmpl.n_coeffs} coefficients. Only {len(df_clean)} usable rows."
                    )
                else:
                    _fit_ready = True

    # ── Step 6: Fit (only if data is ready) ──────────────────────────────
    if _fit_ready:
        st.header("⑥ Fit model")

        if st.button("Run fit", type="primary"):
            try:
                popt, perr, r2, y_pred = _fit_template(tmpl, df_clean, target_col)
            except Exception as exc:
                st.error(f"Fitting failed: {exc}")
                _fit_ready = False
            else:
                st.session_state["rom_fit_result"] = {
                    "popt": popt.tolist(),
                    "perr": perr.tolist(),
                    "r2": r2,
                    "y_pred": y_pred.tolist(),
                    "y_obs": df_clean[target_col].values.tolist(),
                    "template_id": tmpl.id,
                    "param": param,
                    "target_col": target_col,
                    "n_points": len(df_clean),
                }
                st.rerun()

        # ── Show results if available ────────────────────────────────────
        result = st.session_state.get("rom_fit_result")

        if result and result.get("template_id") == tmpl.id and result.get("param") == param:
            popt = np.array(result["popt"])
            perr = np.array(result["perr"])
            r2 = result["r2"]
            y_pred = np.array(result["y_pred"])
            y_obs = np.array(result["y_obs"])

            st.subheader("Fitted coefficients")
            coeff_df = pd.DataFrame({
                "Coefficient": tmpl.coeff_names,
                "Value": [f"{c:.5g}" for c in popt],
                "Std Error": [f"{e:.3g}" for e in perr],
            })
            st.table(coeff_df)

            st.subheader("Fitted equation")
            st.latex(tmpl.latex_filled(popt))
            st.metric("R²", f"{r2:.4f}")
            st.caption(f"Fitted on {result['n_points']} data points.")

            # ── Parity plot ──────────────────────────────────────────────
            st.subheader("Parity plot")

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=y_obs, y=y_pred,
                mode="markers",
                marker=dict(size=8, color="#1f77b4"),
                name="Data",
            ))
            lo = min(y_obs.min(), y_pred.min()) * 0.9
            hi = max(y_obs.max(), y_pred.max()) * 1.1
            fig.add_trace(go.Scatter(
                x=[lo, hi], y=[lo, hi],
                mode="lines",
                line=dict(dash="dash", color="grey"),
                name="1:1",
            ))
            fig.add_trace(go.Scatter(
                x=[lo, hi], y=[lo * 1.2, hi * 1.2],
                mode="lines",
                line=dict(dash="dot", color="lightgrey"),
                name="+20%",
            ))
            fig.add_trace(go.Scatter(
                x=[lo, hi], y=[lo * 0.8, hi * 0.8],
                mode="lines",
                line=dict(dash="dot", color="lightgrey"),
                name="−20%",
            ))
            fig.update_layout(
                xaxis_title="Observed",
                yaxis_title="Predicted",
                width=600, height=500,
                showlegend=True,
            )
            st.plotly_chart(fig, use_container_width=False)

            # ── Residual plot ────────────────────────────────────────────
            st.subheader("Residuals")
            residuals = y_obs - y_pred
            fig_res = go.Figure()
            fig_res.add_trace(go.Scatter(
                x=y_pred, y=residuals,
                mode="markers",
                marker=dict(size=8, color="#ff7f0e"),
            ))
            fig_res.add_hline(y=0, line_dash="dash", line_color="grey")
            fig_res.update_layout(
                xaxis_title="Predicted",
                yaxis_title="Residual (Obs − Pred)",
                width=600, height=400,
            )
            st.plotly_chart(fig_res, use_container_width=False)

            # ── Step 7: Register the correlation ─────────────────────────
            st.header("⑦ Register correlation")

            st.markdown(
                "Save the fitted correlation so it becomes available in the "
                "**Mixing Sensitivity**, **Reactor Comparison**, and other "
                "analysis pages."
            )

            with st.form("register_corr"):
                _reactor_csv = DATA_DIR / "reactors.csv"
                if _reactor_csv.exists():
                    _reactor_df = pd.read_csv(_reactor_csv)
                    _reactor_names = _reactor_df["reactor_name"].dropna().tolist()
                else:
                    _reactor_names = []
                if _reactor_names:
                    reactor_name = st.selectbox(
                        "Reactor",
                        options=_reactor_names,
                        key="reg_reactor_name",
                    )
                else:
                    reactor_name = st.text_input(
                        "Reactor name",
                        placeholder="e.g. Nalas – EasyMax 102",
                    )
                corr_type = st.radio(
                    "Correlation type",
                    options=["Experimental", "ROM"],
                    horizontal=True,
                )
                corr_name = st.text_input(
                    "Correlation name",
                    value=f"{PARAM_DISPLAY.get(param, param)} – fitted {corr_type}",
                )
                corr_source = st.text_input(
                    "Source / reference",
                    value="Fitted in Mixing Lab",
                )
                corr_desc = st.text_area(
                    "Description (optional)",
                    value=f"Fitted from uploaded data using template '{tmpl.name}'. R² = {r2:.4f}.",
                )
                register_btn = st.form_submit_button("Register & save", type="primary")

            if register_btn:
                if not reactor_name.strip():
                    st.error("Reactor name is required.")
                else:
                    func = tmpl.build_func(popt)
                    latex = tmpl.latex_filled(popt)
                    register(
                        reactor_name.strip(),
                        Correlation(
                            name=corr_name,
                            param=param,
                            corr_type=corr_type,
                            func=func,
                            latex=latex,
                            source=corr_source,
                            description=corr_desc,
                            input_params=tmpl.column_labels,
                        ),
                    )

                    entry = {
                        "reactor_name": reactor_name.strip(),
                        "name": corr_name,
                        "param": param,
                        "corr_type": corr_type,
                        "template_id": tmpl.id,
                        "coeffs": popt.tolist(),
                        "r2": r2,
                        "source": corr_source,
                        "description": corr_desc,
                    }
                    _append_and_save(entry)

                    st.success(
                        f"✅ Registered **{corr_name}** for **{reactor_name}** "
                        f"({corr_type}).  It is now available on analysis pages."
                    )
                    st.balloons()

# ═══════════════════════════════════════════════════════════════════════════
#  Tab 2 – Manage Saved Correlations
# ═══════════════════════════════════════════════════════════════════════════
with tab_manage:
    saved = _load_saved()
    if not saved:
        st.info("No fitted correlations saved yet.")
    else:
        st.caption(f"{len(saved)} saved correlation(s).")
        for i, entry in enumerate(saved):
            tmpl_ref = template_by_id(entry["template_id"])
            tmpl_name = tmpl_ref.name if tmpl_ref else entry["template_id"]
            coeffs = np.array(entry["coeffs"])
            label = (
                f"{entry['corr_type']} · {entry['reactor_name']} · "
                f"{PARAM_DISPLAY.get(entry['param'], entry['param'])}"
            )
            with st.expander(label):
                st.markdown(f"**{entry['name']}**")
                if tmpl_ref:
                    st.latex(tmpl_ref.latex_filled(coeffs))
                st.markdown(f"R² = {entry.get('r2', 'N/A')}")
                st.markdown(f"Template: {tmpl_name}")
                st.markdown(f"Source: {entry.get('source', '—')}")
                if entry.get("description"):
                    st.markdown(entry["description"])
                st.caption(f"Coefficients: {coeffs.tolist()}")

                if st.button("🗑️ Delete", key=f"del_{i}"):
                    saved.pop(i)
                    _save_all(saved)
                    st.success("Deleted. Restart the app to fully unregister.")
                    st.rerun()
