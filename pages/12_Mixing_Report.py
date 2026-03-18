"""
Page 12 – Mixing Assessment Report
====================================
Generate a PDF report summarising the latest Mixing Assessment (Page 5) results.
Includes key metrics, operating envelope plots, and sensitivity assessment.
"""

import streamlit as st
import pathlib
import io
import datetime
import numpy as np
import plotly.graph_objects as go

from fpdf import FPDF

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_LOGO = _ROOT / "images" / "general" / "logo.png"

# DejaVu Sans – Unicode-capable TTF bundled with matplotlib
def _find_dejavu() -> pathlib.Path | None:
    """Locate DejaVuSans.ttf from the matplotlib package."""
    try:
        import matplotlib
        d = pathlib.Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"
        if (d / "DejaVuSans.ttf").exists():
            return d
    except ImportError:
        pass
    return None

_DEJAVU_DIR = _find_dejavu()

# Unicode → ASCII substitution map (used when falling back to Helvetica)
_UNICODE_MAP = str.maketrans({
    "\u2014": "--",   # em dash
    "\u2013": "-",    # en dash
    "\u00b7": ".",    # middle dot
    "\u00b0": "deg",  # degree
    "\u00b2": "2",    # superscript 2
    "\u00b3": "3",    # superscript 3
    "\u00b5": "u",    # micro sign
    "\u03b7": "eta",  # eta
    "\u03b5": "eps",  # epsilon
    "\u03bb": "lambda",
    "\u03c1": "rho",
    "\u03c6": "phi",
    "\u0394": "D",    # Delta
})

def _safe_text(text: str, is_unicode_font: bool) -> str:
    """If using a non-Unicode font, replace special characters with ASCII."""
    if is_unicode_font:
        return text
    return text.translate(_UNICODE_MAP).encode("latin-1", errors="replace").decode("latin-1")

st.title("📄 Mixing Assessment Report")

# ── Check for results ────────────────────────────────────────────────────
_snap = st.session_state.get("_ms_report_snapshot")
if not _snap:
    st.warning(
        "No mixing assessment results available.\n\n"
        "Run a **Mixing Assessment** (Page 5) first — results are automatically "
        "captured when the computation completes."
    )
    st.stop()

# ── Unpack snapshot ──────────────────────────────────────────────────────
reactor_name = _snap["reactor"]
reaction_name = _snap["reaction"]
fluid_name = _snap["fluid"]
fluid_T_C = _snap["fluid_T_C"]
N_rpm = _snap["N_rpm"]
V_L = _snap["V_L"]
hydro = _snap["hydro"]
da = _snap["da"]
t_rxn = _snap["t_rxn"]
heat_results = _snap.get("heat_results", {})
particle_results = _snap.get("particle_results", {})
particle_meta = _snap.get("particle_meta", {})
batchelor_um = _snap.get("batchelor_um", 0.0)
envelope = _snap.get("envelope")  # may be None

# ── Helpers ──────────────────────────────────────────────────────────────
_DISPLAY_NAMES = {
    "Da_macro": "Macromixing (Da_macro)",
    "Da_micro": "Micromixing (Da_micro)",
    "Da_GL": "Gas-Liquid Mass Transfer (Da_GL)",
    "Q_gen/Q_cool (%)": "Heat Capacity (Q_gen/Q_cool %)",
}

def _da_text(Da: float) -> str:
    if Da < 0.01:
        return "Not sensitive"
    if Da < 0.1:
        return "Likely not sensitive"
    if Da < 1:
        return "Potentially sensitive"
    if Da < 10:
        return "Likely sensitive"
    return "Highly sensitive"

def _da_symbol(Da: float) -> str:
    if Da < 0.1:
        return "GREEN"
    if Da < 1:
        return "AMBER"
    return "RED"


# ── Build envelope chart images (PNG bytes via Plotly) ───────────────────
_REPORT_PARAMS = ["Da_micro", "Da_macro", "Da_GL", "Blend time 95% (s)", "P/V (W/L)"]
if heat_results:
    _REPORT_PARAMS.append("Q_gen/Q_cool (%)")

_MODE_COLORS = {"Literature": "#3366CC", "ROM": "#33AA66", "Experimental": "#FF8800"}

def _build_envelope_fig(param: str) -> go.Figure | None:
    """Build a single operating-envelope Plotly figure for *param*."""
    if envelope is None:
        return None
    curve_data = envelope["curve_data"]
    pct_arr = np.array(envelope["pct_arr"])
    active_modes = envelope["active_modes"]
    priority_mode_label = envelope["priority_mode_label"]
    current_pct = envelope["current_pct"]
    env_V_max = envelope["env_V_max"]
    env_V_min = envelope["env_V_min"]
    rpm_max = envelope["rpm_max"]

    # Check param exists in curve data
    first_mode = list(curve_data.keys())[0]
    if param not in curve_data[first_mode]["maxV"]:
        return None

    fig = go.Figure()
    for mode_label in active_modes:
        if mode_label not in curve_data:
            continue
        color = _MODE_COLORS.get(mode_label, "#999999")
        mc = curve_data[mode_label]
        y_max = np.array(mc["maxV"][param])
        y_min = np.array(mc["minV"][param])
        # Envelope fill
        poly_x = np.concatenate([pct_arr, pct_arr[::-1], [pct_arr[0]]])
        poly_y = np.concatenate([y_max, y_min[::-1], [y_max[0]]])
        fig.add_trace(go.Scatter(
            x=poly_x, y=poly_y, fill="toself", fillcolor=color, opacity=0.15,
            line=dict(color=color, width=1), mode="lines",
            name=f"{mode_label} envelope", legendgroup=mode_label,
            hoverinfo="skip",
        ))
        fig.add_trace(go.Scatter(
            x=pct_arr, y=y_max, mode="lines",
            line=dict(color=color, width=2),
            name=f"{mode_label} max V ({env_V_max:.1f} L)",
            legendgroup=mode_label,
        ))
        fig.add_trace(go.Scatter(
            x=pct_arr, y=y_min, mode="lines",
            line=dict(color=color, width=2, dash="dot"),
            name=f"{mode_label} min V ({env_V_min:.1f} L)",
            legendgroup=mode_label,
        ))

    # Operating point
    if priority_mode_label in curve_data:
        pc = curve_data[priority_mode_label]
        y_maxp = np.array(pc["maxV"][param])
        y_minp = np.array(pc["minV"][param])
        if abs(env_V_max - env_V_min) > 1e-6:
            frac = max(0.0, min(1.0, (V_L - env_V_min) / (env_V_max - env_V_min)))
            y_interp = (np.interp(current_pct, pct_arr, y_minp) * (1 - frac)
                        + np.interp(current_pct, pct_arr, y_maxp) * frac)
        else:
            y_interp = np.interp(current_pct, pct_arr, y_maxp)
        fig.add_trace(go.Scatter(
            x=[current_pct], y=[y_interp],
            mode="markers", marker=dict(size=12, color="red", symbol="star",
                                         line=dict(width=1, color="white")),
            name="Current",
        ))

    # Reference lines
    if param in ("Da_macro", "Da_micro", "Da_GL"):
        import math
        for da_val, da_color, label in [
            (0.1, "orange", "Da=0.1"), (1.0, "red", "Da=1"),
        ]:
            fig.add_shape(type="line", x0=0, x1=1, y0=da_val, y1=da_val,
                          xref="paper", yref="y",
                          line=dict(color=da_color, width=1.5, dash="dash"))
        fig.update_yaxes(type="log")
    if param == "Q_gen/Q_cool (%)":
        fig.add_shape(type="line", x0=0, x1=1, y0=100, y1=100,
                      xref="paper", yref="y",
                      line=dict(color="red", width=1.5, dash="dash"))

    display = _DISPLAY_NAMES.get(param, param)
    fig.update_layout(
        title=display, xaxis_title=f"Stir speed (% of max RPM = {rpm_max:.0f})",
        yaxis_title=display,
        xaxis=dict(range=[0, 105], dtick=10),
        height=400, width=700, margin=dict(t=50, b=50),
    )
    return fig


def _fig_to_png_bytes(fig: go.Figure) -> bytes:
    return fig.to_image(format="png", scale=2)


# ── PDF builder ──────────────────────────────────────────────────────────
class MixingReport(FPDF):
    """Custom FPDF subclass with header/footer branding."""

    _FONT = "DejaVu"   # registered family name
    _unicode_font = True

    def __init__(self, logo_path: str | None, report_title: str, **kw):
        super().__init__(**kw)
        self._logo_path = logo_path
        self._report_title = report_title
        self.set_auto_page_break(auto=True, margin=25)
        # Register DejaVu Sans (Unicode-capable), fall back to Helvetica
        if _DEJAVU_DIR:
            try:
                self.add_font("DejaVu", "",  str(_DEJAVU_DIR / "DejaVuSans.ttf"))
                self.add_font("DejaVu", "B", str(_DEJAVU_DIR / "DejaVuSans-Bold.ttf"))
                self.add_font("DejaVu", "I", str(_DEJAVU_DIR / "DejaVuSans-Oblique.ttf"))
                self.add_font("DejaVu", "BI", str(_DEJAVU_DIR / "DejaVuSans-BoldOblique.ttf"))
            except Exception:
                self._FONT = "Helvetica"
                self._unicode_font = False
        else:
            self._FONT = "Helvetica"
            self._unicode_font = False

    def header(self):
        if self._logo_path and pathlib.Path(self._logo_path).exists():
            self.image(self._logo_path, x=10, y=8, h=12)
        self.set_font(self._FONT, "B", 10)
        self.set_text_color(100, 100, 100)
        self.cell(0, 10, self._s(self._report_title), align="R")
        self.ln(14)
        # Thin line under header
        self.set_draw_color(200, 200, 200)
        self.line(10, self.get_y(), self.w - 10, self.get_y())
        self.ln(4)

    def footer(self):
        self.set_y(-20)
        self.set_draw_color(200, 200, 200)
        self.line(10, self.get_y(), self.w - 10, self.get_y())
        self.set_y(-15)
        self.set_font(self._FONT, "I", 8)
        self.set_text_color(140, 140, 140)
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        self.cell(0, 10, f"Generated {ts}", align="L")
        self.cell(0, 10, f"Page {self.page_no()}/{{nb}}", align="R")

    # ── convenience methods ──────────────────────────────────────────────
    def _s(self, text: str) -> str:
        """Sanitise text for the active font."""
        return _safe_text(text, self._unicode_font)

    def section_title(self, text: str):
        self.set_font(self._FONT, "B", 13)
        self.set_text_color(30, 30, 80)
        self.cell(0, 9, self._s(text))
        self.ln(10)

    def sub_title(self, text: str):
        self.set_font(self._FONT, "B", 11)
        self.set_text_color(50, 50, 50)
        self.cell(0, 8, self._s(text))
        self.ln(8)

    def kv(self, key: str, value: str, bold_val: bool = False):
        self.set_font(self._FONT, "", 10)
        self.set_text_color(60, 60, 60)
        self.cell(70, 6, self._s(key))
        style = "B" if bold_val else ""
        self.set_font(self._FONT, style, 10)
        self.set_text_color(20, 20, 20)
        self.cell(0, 6, self._s(value))
        self.ln(6)

    def metric_table(self, rows: list[tuple[str, str]], cols: int = 2):
        """Render key-value pairs in a multi-column table."""
        col_w = (self.w - 20) / cols / 2
        self.set_font(self._FONT, "", 9)
        for i, (k, v) in enumerate(rows):
            if i > 0 and i % cols == 0:
                self.ln(6)
            self.set_text_color(80, 80, 80)
            self.cell(col_w, 6, self._s(k))
            self.set_text_color(20, 20, 20)
            self.set_font(self._FONT, "B", 9)
            self.cell(col_w, 6, self._s(v))
            self.set_font(self._FONT, "", 9)
        self.ln(8)

    def assessment_box(self, text: str, colour: str = "GREEN"):
        """Draw a coloured assessment banner."""
        if colour == "GREEN":
            r, g, b = 34, 139, 34
            br, bg, bb = 220, 245, 220
        elif colour == "AMBER":
            r, g, b = 180, 130, 0
            br, bg, bb = 255, 245, 210
        else:
            r, g, b = 200, 30, 30
            br, bg, bb = 255, 220, 220
        self.set_fill_color(br, bg, bb)
        self.set_text_color(r, g, b)
        self.set_font(self._FONT, "B", 10)
        self.cell(0, 8, f"  {self._s(text)}", fill=True)
        self.ln(9)
        self.set_text_color(0, 0, 0)


def _build_pdf() -> bytes:
    logo = str(_LOGO) if _LOGO.exists() else None
    title = f"Mixing Assessment \u2014 {reactor_name}"
    pdf = MixingReport(logo, title, orientation="P", unit="mm", format="A4")
    pdf.alias_nb_pages()

    # ── Page 1: Title & System Info ──────────────────────────────────────
    pdf.add_page()
    pdf.set_font(pdf._FONT, "B", 20)
    pdf.set_text_color(30, 30, 80)
    pdf.cell(0, 14, pdf._s("Mixing Assessment Report"), align="C")
    pdf.ln(18)

    pdf.section_title("System Configuration")
    pdf.kv("Reactor", reactor_name, bold_val=True)
    pdf.kv("Reaction", reaction_name, bold_val=True)
    pdf.kv("Fluid", f"{fluid_name}  ({fluid_T_C:.1f} °C)", bold_val=True)
    pdf.kv("Stir speed", f"{N_rpm:.0f} RPM")
    pdf.kv("Liquid volume", f"{V_L:.2f} L")
    pdf.kv("Reaction time (t_rxn)", f"{t_rxn:.4g} s")
    pdf.ln(4)

    # ── Hydrodynamic metrics ─────────────────────────────────────────────
    pdf.section_title("Hydrodynamic & Mixing Metrics")
    metrics = [
        ("Re", f"{hydro['Re']:.0f}"),
        ("P/V (W/L)", f"{hydro['P/V (W/L)']:.3g}"),
        ("Blend time 95% (s)", f"{hydro['Blend time 95% (s)']:.2f}"),
        ("Micromix t_E (s)", f"{hydro['Micromix time t_E (s)']:.4g}"),
        ("Tip speed (m/s)", f"{hydro['Tip speed (m/s)']:.2f}"),
        ("Kolmogorov eta (um)", f"{hydro['Kolmogorov η (µm)']:.1f}"),
        ("Batchelor lambda_B (um)", f"{batchelor_um:.2f}"),
        ("Circulation time (s)", f"{hydro['Circulation time (s)']:.2f}"),
        ("Avg shear rate (1/s)", f"{hydro['Avg shear rate (1/s)']:.1f}"),
        ("Max shear rate (1/s)", f"{hydro['Max shear rate (1/s)']:.0f}"),
        ("Avg shear stress (Pa)", f"{hydro['Avg shear stress (Pa)']:.3g}"),
        ("EDCF (W/kg/s)", f"{hydro['EDCF (W/kg/s)']:.3g}"),
        ("Torque (N.m)", f"{hydro['Torque (N·m)']:.3g}"),
        ("Froude number", f"{hydro['Froude number']:.4g}"),
    ]
    pdf.metric_table(metrics, cols=2)
    pdf.ln(2)

    # Mass transfer
    pdf.sub_title("Mass Transfer")
    mt_rows = [
        ("kLa surface (1/s)", f"{hydro['kLa_surface (1/s)']:.4g}"),
    ]
    if hydro.get("kLa (1/s)", 0) > 0:
        mt_rows.insert(0, ("kLa sparged (1/s)", f"{hydro['kLa (1/s)']:.4g}"))
    pdf.metric_table(mt_rows, cols=2)

    # Heat transfer
    if heat_results:
        pdf.sub_title("Heat Transfer")
        ht_rows = [
            ("Q_gen (W)", f"{heat_results['Q_gen (W)']:.1f}"),
            ("Q_cool (W)", f"{heat_results['Q_cool (W)']:.1f}"),
            ("U (W/m2.K)", f"{heat_results['U (W/m²·K)']:.0f}"),
            ("A_ht (m2)", f"{heat_results['A_ht (m²)']:.3f}"),
            ("Q_gen/Q_cool (%)", f"{heat_results['Q_gen/Q_cool (%)']:.1f}%"),
        ]
        pdf.metric_table(ht_rows, cols=2)

    # Particles
    if particle_results:
        pdf.sub_title("Solid Particles")
        sp_rows = [
            ("Particle", particle_meta.get("Particle", "")),
            ("d50 (um)", f"{particle_results['d50 (µm)']:.1f}"),
            ("rho_p (kg/m3)", f"{particle_results['ρ_p (kg/m³)']:.0f}"),
            ("N_js design (RPM)", f"{particle_results['N_js (RPM)']:.1f}"),
            ("v_t (m/s)", f"{particle_results['v_t (m/s)']:.3e}"),
            ("k_SL (m/s)", f"{particle_results['k_SL (m/s)']:.3e}"),
        ]
        pdf.metric_table(sp_rows, cols=2)
        susp = particle_meta.get("Suspension", "")
        _susp_col = "GREEN" if "Well" in susp or "Just" in susp else ("AMBER" if "Partial" in susp else "RED")
        pdf.assessment_box(susp, _susp_col)

    # ── Page 2: Sensitivity Assessment ───────────────────────────────────
    pdf.add_page()
    pdf.section_title("Mixing Sensitivity Assessment")

    da_macro = da["Da_macro"]
    da_micro = da["Da_micro"]
    da_gl = da["Da_GL"]

    pdf.kv("Da_macro", f"{da_macro:.3g}")
    pdf.assessment_box(
        f"Macromixing: {_da_text(da_macro)} (Da = {da_macro:.3g})",
        _da_symbol(da_macro),
    )
    pdf.kv("Da_micro", f"{da_micro:.3g}")
    pdf.assessment_box(
        f"Micromixing: {_da_text(da_micro)} (Da = {da_micro:.3g})",
        _da_symbol(da_micro),
    )
    if da_gl > 0:
        pdf.kv("Da_GL", f"{da_gl:.3g}")
        pdf.assessment_box(
            f"Gas-liquid: {_da_text(da_gl)} (Da = {da_gl:.3g})",
            _da_symbol(da_gl),
        )
    if heat_results:
        ratio = heat_results["Q_gen/Q_cool (%)"]
        heat_col = "GREEN" if ratio < 100 else "RED"
        pdf.assessment_box(
            f"Heat balance: Q_gen/Q_cool = {ratio:.1f}%",
            heat_col,
        )

    overall = da.get("Assessment", "")
    if overall:
        pdf.ln(4)
        pdf.sub_title("Overall Assessment")
        if "not" in overall.lower():
            _oc = "GREEN"
        elif "potentially" in overall.lower():
            _oc = "AMBER"
        else:
            _oc = "RED"
        pdf.assessment_box(overall, _oc)

    # ── Page 3+: Operating Envelope Charts ───────────────────────────────
    if envelope is not None:
        pdf.add_page()
        pdf.section_title("Operating Envelopes")
        pdf.set_font(pdf._FONT, "", 9)
        pdf.set_text_color(80, 80, 80)
        env_info = envelope
        pdf.cell(0, 6, pdf._s(
                 f"RPM range: {env_info['rpm_min']:.0f} - {env_info['rpm_max']:.0f}  |  "
                 f"Volume range: {env_info['env_V_min']:.1f} - {env_info['env_V_max']:.1f} L"))
        pdf.ln(10)

        _chart_count_on_page = 0
        for param in _REPORT_PARAMS:
            fig = _build_envelope_fig(param)
            if fig is None:
                continue
            try:
                png = _fig_to_png_bytes(fig)
            except Exception:
                continue
            if not png or png[:4] != b"\x89PNG":
                continue
            # Each chart: ~90mm tall. Fit 2 per page
            if _chart_count_on_page >= 2:
                pdf.add_page()
                _chart_count_on_page = 0
            pdf.image(io.BytesIO(png), x=15, w=180)
            pdf.ln(5)
            _chart_count_on_page += 1

    return bytes(pdf.output())


# ── Preview section ──────────────────────────────────────────────────────
st.header("Report Preview")

st.subheader("System")
c1, c2, c3 = st.columns(3)
c1.metric("Reactor", reactor_name)
c2.metric("Reaction", reaction_name)
c3.metric("Fluid", f"{fluid_name} ({fluid_T_C:.1f} °C)")

c4, c5, c6 = st.columns(3)
c4.metric("RPM", f"{N_rpm:.0f}")
c5.metric("Volume (L)", f"{V_L:.2f}")
c6.metric("t_rxn (s)", f"{t_rxn:.4g}")

st.divider()
st.subheader("Key Metrics")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Re", f"{hydro['Re']:.0f}")
m2.metric("P/V (W/L)", f"{hydro['P/V (W/L)']:.3g}")
m3.metric("Blend time (s)", f"{hydro['Blend time 95% (s)']:.2f}")
m4.metric("Micromix t_E (s)", f"{hydro['Micromix time t_E (s)']:.4g}")

m5, m6, m7, m8 = st.columns(4)
m5.metric("Tip speed (m/s)", f"{hydro['Tip speed (m/s)']:.2f}")
m6.metric("Kolmogorov (µm)", f"{hydro['Kolmogorov η (µm)']:.1f}")
m7.metric("Da_macro", f"{da['Da_macro']:.3g}")
m8.metric("Da_micro", f"{da['Da_micro']:.3g}")

st.divider()
st.subheader("Sensitivity Assessment")

def _banner(label, Da):
    if Da < 0.1:
        st.success(f"🟢 **{label}:** {_da_text(Da)} (Da = {Da:.3g})")
    elif Da < 1:
        st.warning(f"🟡 **{label}:** {_da_text(Da)} (Da = {Da:.3g})")
    else:
        st.error(f"🔴 **{label}:** {_da_text(Da)} (Da = {Da:.3g})")

_banner("Macromixing", da["Da_macro"])
_banner("Micromixing", da["Da_micro"])
if da["Da_GL"] > 0:
    _banner("Gas-liquid", da["Da_GL"])
if heat_results:
    ratio = heat_results["Q_gen/Q_cool (%)"]
    if ratio < 100:
        st.success(f"🟢 **Heat balance:** Q_gen/Q_cool = {ratio:.1f}%")
    else:
        st.error(f"🔴 **Heat balance:** Q_gen/Q_cool = {ratio:.1f}%")

overall = da.get("Assessment", "")
if overall:
    st.info(f"**Overall:** {overall}")

# Envelope preview charts
if envelope is not None:
    st.divider()
    st.subheader("Operating Envelopes (included in report)")
    for param in _REPORT_PARAMS:
        fig = _build_envelope_fig(param)
        if fig is not None:
            st.plotly_chart(fig, use_container_width=True)

# ── Export button ────────────────────────────────────────────────────────
st.divider()
st.header("Export")

if st.button("📥 Export PDF Report", type="primary"):
    with st.spinner("Generating PDF…"):
        try:
            pdf_bytes = _build_pdf()
            st.session_state["_ms_pdf_bytes"] = pdf_bytes
            st.session_state["_ms_pdf_name"] = (
                f"Mixing_Report_{reactor_name.replace(' ', '_')}_"
                f"{datetime.datetime.now():%Y%m%d_%H%M}.pdf"
            )
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")

if "_ms_pdf_bytes" in st.session_state:
    st.download_button(
        "⬇️ Download PDF",
        data=st.session_state["_ms_pdf_bytes"],
        file_name=st.session_state["_ms_pdf_name"],
        mime="application/pdf",
    )
    st.success("PDF ready for download.")
