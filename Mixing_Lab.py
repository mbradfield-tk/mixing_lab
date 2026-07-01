"""
Mixing Lab – Hydrodynamics & Mixing-Sensitivity Assessment Tool
================================================================
Main entry point.  Run with:  streamlit run Mixing_Lab.py
"""

import streamlit as st
from utils.sidebar import render_sidebar
from utils.usage import log_access

st.set_page_config(
    page_title="Mixing Lab",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Record one access event per browser session (client IP, timestamp, etc.).
log_access()

# Render branding first so the logo sits at the top of the sidebar,
# above the navigation menu.
render_sidebar()

# ── Page navigation ───────────────────────────────────────────────────────
home = st.Page("pages/home.py", title="App Overview", default=True)

databases = [
    st.Page("pages/1_Reactor_Database.py",  title="Reactor Database",  icon="⚗️"),
    st.Page("pages/2_Reaction_Database.py", title="Reaction Database", icon="🧪"),
    st.Page("pages/3_Fluid_Database.py",    title="Fluid Database",    icon="💧"),
    st.Page("pages/4_Particle_Database.py", title="Particle Database", icon="❇️"),
]

analysis = [
    st.Page("pages/10_Mixing_Sensitivity_Protocol.py", title="Reaction Sensitivity Protocol", icon="🧭"),
    st.Page("pages/11_Crystallization_Sensitivity_Protocol.py", title="Crystallization Sensitivity Protocol (WIP)", icon="🚧"),
    st.Page("pages/5_Mixing_Sensitivity.py", title="Mixing Assessment", icon="🌀"),
    st.Page("pages/6_Bourne_Protocol.py",    title="Bourne Protocol",    icon="🅱️"),
    st.Page("pages/7_Reactor_Comparison.py", title="Reactor Comparison", icon="📈"),
]

tools = [
    st.Page("pages/9_Unit_Converter.py", title="Unit Converter", icon="🔄"),
    st.Page("pages/14_Heat_Transfer.py", title="Heat Transfer", icon="🔥"),
]

results = [
    st.Page("pages/8_Recorded_Results.py", title="Recorded Results", icon="📋"),
]

reference = [
    st.Page("pages/equations_reference.py", title="Equations Reference", icon="📐"),
]

admin = [
    st.Page("pages/0_Admin_Import.py", title="Admin Tools", icon="🛠️"),
    st.Page("pages/13_ROM_Fitting.py", title="ROM / Correlation Fitting", icon="🔧"),
]

nav = st.navigation({
    "": [home],
    "Databases": databases,
    "Analysis": analysis,
    "Tools": tools,
    "Results": results,
    "Reference": reference,
    "Admin": admin,
},
position="sidebar",
expanded=False)

nav.run()
