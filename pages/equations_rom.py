"""Equations – Reactor-Specific Correlations (ROM & Experimental).

Dynamically displays all registered reduced-order-model (ROM) and
experimental correlations from the ROM registry.
"""

import streamlit as st
from utils.rom_registry import (
    get_all_correlations,
    get_all_registered_reactors,
    PARAM_DISPLAY,
    SUPPORTED_PARAMS,
)

st.title("📐 Reactor-Specific Correlations (ROM & Experimental)")

st.markdown("""
In addition to the standard **literature correlations** used throughout
Mixing Lab, reactor-specific correlations can be registered from:

- **ROM (Reduced Order Models)** – fitted from CFD parametric studies or
  surrogate models.
- **Experimental** – fitted directly to measured data (e.g. kLa from
  gassing-out experiments, blend time from PLIF/decolourisation).

These override the literature correlation for a specific parameter when
the corresponding mode is selected on the **Mixing Sensitivity** or
**Reactor Comparison** pages.

---
""")

# ── Which parameters can be overridden? ──────────────────────────────────
st.header("Supported Parameters")
st.markdown(
    "The following parameters can be overridden by ROM or Experimental "
    "correlations.  All other parameters continue to use literature "
    "correlations.\n\n"
    + "\n".join(f"- {PARAM_DISPLAY[p]}" for p in SUPPORTED_PARAMS)
)

# ── Registered correlations ──────────────────────────────────────────────
all_corrs = get_all_correlations()
registered_reactors = get_all_registered_reactors()

if not registered_reactors:
    st.info("No reactor-specific correlations have been registered yet.")
    st.stop()

st.header("Registered Correlations")

for reactor_name in registered_reactors:
    corrs = all_corrs[reactor_name]
    st.subheader(f"🔧 {reactor_name}")

    for corr in corrs:
        label = f"{corr.corr_type} · {PARAM_DISPLAY.get(corr.param, corr.param)}"
        with st.expander(label, expanded=True):
            st.markdown(f"**{corr.name}**")
            st.latex(corr.latex)

            if corr.description:
                st.markdown(corr.description)

            if corr.input_params:
                rows = "\n".join(
                    f"| `{k}` | {v} |" for k, v in corr.input_params.items()
                )
                st.markdown(
                    "| Input | Description |\n"
                    "|-------|-------------|\n"
                    + rows
                )

            st.caption(f"Source: {corr.source}")

st.markdown("---")
st.markdown(
    "💡 **Adding new correlations:** Register them in "
    "`utils/rom_registry.py` using the `register()` function.  "
    "They will automatically appear here and become selectable on "
    "the calculation pages."
)
