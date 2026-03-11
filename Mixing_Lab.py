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
home = st.Page("pages/home.py", title="🔬 Mixing Lab", default=True)

databases = [
    st.Page("pages/1_Reactor_Database.py",  title="⚗️ Reactor Database"),
    st.Page("pages/2_Reaction_Database.py", title="🧪 Reaction Database"),
    st.Page("pages/3_Fluid_Database.py",    title="💧 Fluid Database"),
    st.Page("pages/4_Particle_Database.py", title="❉ Particle Database"),
]

analysis = [
    st.Page("pages/5_Mixing_Sensitivity.py", title="⚙️ Mixing Sensitivity"),
    st.Page("pages/6_Bourne_Protocol.py",    title="🧫 Bourne Protocol"),
    st.Page("pages/7_Reactor_Comparison.py", title="📊 Reactor Comparison"),
]

results = [
    st.Page("pages/8_Recorded_Results.py", title="📋 Recorded Results"),
]

reference = [
    st.Page("pages/equations_hydrodynamics.py",  title="📐 Hydrodynamics & Shear"),
    st.Page("pages/equations_mixing.py",         title="📐 Mixing & Damköhler"),
    st.Page("pages/equations_mass_transfer.py",  title="📐 Mass Transfer & kLa"),
    st.Page("pages/equations_heat.py",           title="📐 Heat Balance"),
    st.Page("pages/equations_particles.py",      title="📐 Solid Particles"),
    st.Page("pages/equations_scaleup.py",        title="📐 Scale-Up & Properties"),
]

admin = [
    st.Page("pages/0_Admin_Import.py", title="🛠️ Admin Import"),
]

nav = st.navigation({
    "": [home],
    "Databases": databases,
    "Analysis": analysis,
    "Results": results,
    "Reference": reference,
    "Admin": admin,
})

nav.run()
