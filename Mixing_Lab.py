"""
Mixing Lab – Hydrodynamics & Mixing-Sensitivity Assessment Tool
================================================================
Main entry point.  Run with:  streamlit run Mixing_Lab.py
"""

import streamlit as st
from utils.sidebar import render_sidebar

st.set_page_config(
    page_title="Mixing Lab",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

render_sidebar()

# ── Page navigation ───────────────────────────────────────────────────────
home = st.Page("pages/home.py", title="Mixing Lab", icon="🔬", default=True)

databases = [
    st.Page("pages/1_Reactor_Database.py",  title="Reactor Database",  icon="⚗️"),
    st.Page("pages/2_Reaction_Database.py", title="Reaction Database", icon="🧪"),
    st.Page("pages/3_Fluid_Database.py",    title="Fluid Database",    icon="💧"),
    st.Page("pages/4_Particle_Database.py", title="Particle Database", icon="❇️"),
]

analysis = [
    st.Page("pages/10_Mixing_Sensitivity_Protocol.py", title="Sensitivity Protocol", icon="🧭"),
    st.Page("pages/5_Mixing_Sensitivity.py", title="Mixing Assessment", icon="🌀"),
    st.Page("pages/6_Bourne_Protocol.py",    title="Bourne Protocol",    icon="🧐"),
    st.Page("pages/7_Reactor_Comparison.py", title="Reactor Comparison", icon="📈"),
]

tools = [
    st.Page("pages/9_Unit_Converter.py", title="Unit Converter", icon="🔄"),
]

results = [
    st.Page("pages/8_Recorded_Results.py", title="Recorded Results", icon="📋"),
]

reference = [
    st.Page("pages/equations_hydrodynamics.py",  title="Hydrodynamics & Shear",  icon="📐"),
    st.Page("pages/equations_mixing.py",         title="Mixing & Damköhler",     icon="📐"),
    st.Page("pages/equations_mass_transfer.py",  title="Mass Transfer & kLa",    icon="📐"),
    st.Page("pages/equations_heat.py",           title="Heat Balance",            icon="📐"),
    st.Page("pages/equations_particles.py",      title="Solid Particles",         icon="📐"),
    st.Page("pages/equations_scaleup.py",        title="Scale-Up & Properties",   icon="📐"),
    st.Page("pages/equations_rom.py",            title="ROM & Experimental",      icon="📐"),
]

admin = [
    st.Page("pages/0_Admin_Import.py", title="Admin Tools", icon="🛠️"),
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
