"""Page 2 – Reaction Database: browse, add, import/export reaction kinetics."""

import streamlit as st

import pandas as pd
import pathlib

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
REACTION_CSV = DATA_DIR / "reactions.csv"


def _load_reactions() -> pd.DataFrame:
    if REACTION_CSV.exists():
        return pd.read_csv(REACTION_CSV)
    return pd.DataFrame(columns=[
        "reaction_name", "type", "order", "k_value", "k_units",
        "C0_mol_L", "t_rxn_s", "T_C", "solvent", "delta_H_kJ_mol", "notes",
    ])


def _save_reactions(df: pd.DataFrame):
    df.to_csv(REACTION_CSV, index=False)


if "reaction_db" not in st.session_state:
    st.session_state.reaction_db = _load_reactions()

st.title("🧪 Reaction Database")

tab_browse, tab_add, tab_import = st.tabs(["Browse & Edit", "Add Reaction", "Import / Export"])

# ── Browse & Edit ─────────────────────────────────────────────────────────
with tab_browse:
    st.markdown("Edit kinetic data inline.  Click **Save changes** when done.")

    c1, c2 = st.columns(2)
    with c1:
        filt_type = st.multiselect("Reaction type", st.session_state.reaction_db["type"].dropna().unique().tolist())
    with c2:
        filt_solvent = st.multiselect("Solvent", st.session_state.reaction_db["solvent"].dropna().unique().tolist())

    df = st.session_state.reaction_db.copy()
    if filt_type:
        df = df[df["type"].isin(filt_type)]
    if filt_solvent:
        df = df[df["solvent"].isin(filt_solvent)]

    edited = st.data_editor(df, num_rows="dynamic", use_container_width=False, key="rxn_editor")

    btn_col1, btn_col2, _ = st.columns([1, 1, 4])
    with btn_col1:
        if st.button("💾 Save changes", key="save_rxn"):
            st.session_state.reaction_db = edited.copy()
            _save_reactions(st.session_state.reaction_db)
            st.success("Reaction database saved.")
    with btn_col2:
        if st.button("🔄 Refresh", key="refresh_rxn"):
            st.session_state.reaction_db = _load_reactions()
            st.rerun()

# ── Add Reaction ──────────────────────────────────────────────────────────
with tab_add:
    # Optional: pre-fill from an existing reaction
    rxn_db = st.session_state.reaction_db
    template_options = ["— blank —"] + rxn_db["reaction_name"].dropna().tolist()
    template_choice = st.selectbox(
        "Pre-fill from existing reaction",
        template_options,
        key="rxn_template_choice",
        help="Select an existing reaction to populate the fields below, then adjust as needed.",
    )
    tpl: dict = {}
    if template_choice != "— blank —":
        tpl_row = rxn_db[rxn_db["reaction_name"] == template_choice]
        if not tpl_row.empty:
            tpl = tpl_row.iloc[0].to_dict()

    def _tpl(col, default):
        """Return template value for *col*, falling back to *default*."""
        v = tpl.get(col, default)
        if isinstance(v, float) and pd.isna(v):
            return default
        return v

    with st.form("add_rxn"):
        c1, c2, c3 = st.columns(3)
        with c1:
            name = st.text_input("Reaction name *",
                                 value=f"{_tpl('reaction_name', '')} (copy)" if tpl else "")
            rxn_type = st.text_input("Type (e.g. Cross-coupling)",
                                     value=str(_tpl("type", "")))
            order_options = ["1", "2", "pseudo-1", "pseudo-2", "n/a"]
            tpl_order = str(_tpl("order", "1"))
            order_idx = order_options.index(tpl_order) if tpl_order in order_options else 0
            order = st.selectbox("Kinetic order", order_options, index=order_idx)
        with c2:
            k_val = st.number_input("Rate constant k", min_value=0.0,
                                    value=float(_tpl("k_value", 0.01)), format="%.6g")
            k_units = st.text_input("k units",
                                    value=str(_tpl("k_units", "1/s")))
            C0 = st.number_input("C₀ (mol/L)", min_value=0.0,
                                 value=float(_tpl("C0_mol_L", 0.1)), format="%.4g")
        with c3:
            t_rxn = st.number_input("Characteristic reaction time (s)", min_value=0.0,
                                    value=float(_tpl("t_rxn_s", 0.0)), format="%.4g",
                                    help="Leave 0 to auto-compute from k and C₀")
            T = st.number_input("Temperature (°C)",
                                value=float(_tpl("T_C", 25.0)))
            solvent = st.text_input("Solvent",
                                    value=str(_tpl("solvent", "THF")))
        c4, c5 = st.columns(2)
        with c4:
            delta_H = st.number_input("ΔH_rxn (kJ/mol)",
                                      value=float(_tpl("delta_H_kJ_mol", 0.0)), format="%.1f",
                                      help="Heat of reaction (negative = exothermic)")
        with c5:
            notes = st.text_area("Notes",
                                 value=str(_tpl("notes", "")))
        submitted = st.form_submit_button("Add reaction")

        if submitted and name:
            # Auto-compute t_rxn if not provided
            if t_rxn == 0 and k_val > 0:
                if order in ("1", "pseudo-1"):
                    import numpy as np
                    t_rxn = float(np.log(2) / k_val)
                elif order in ("2", "pseudo-2") and C0 > 0:
                    t_rxn = 1.0 / (k_val * C0)

            new = pd.DataFrame([{
                "reaction_name": name, "type": rxn_type, "order": order,
                "k_value": k_val, "k_units": k_units, "C0_mol_L": C0,
                "t_rxn_s": t_rxn, "T_C": T, "solvent": solvent,
                "delta_H_kJ_mol": delta_H, "notes": notes,
            }])
            st.session_state.reaction_db = pd.concat(
                [st.session_state.reaction_db, new], ignore_index=True)
            _save_reactions(st.session_state.reaction_db)
            st.success(f"Added **{name}**.")
        elif submitted:
            st.warning("Enter a reaction name.")

# ── Import / Export ───────────────────────────────────────────────────────
with tab_import:
    st.download_button(
        "⬇️ Download reaction database (CSV)",
        data=st.session_state.reaction_db.to_csv(index=False).encode("utf-8"),
        file_name="reactions_export.csv",
        mime="text/csv",
    )
    uploaded = st.file_uploader("Upload CSV", type=["csv"], key="rxn_upload")
    if uploaded:
        try:
            new_df = pd.read_csv(uploaded)
            st.dataframe(new_df.head())
            mode = st.radio("Import mode", ["Replace", "Append"], key="rxn_import_mode")
            if st.button("Confirm import", key="rxn_import_confirm"):
                if mode == "Replace":
                    st.session_state.reaction_db = new_df
                else:
                    st.session_state.reaction_db = pd.concat(
                        [st.session_state.reaction_db, new_df], ignore_index=True)
                _save_reactions(st.session_state.reaction_db)
                st.success("Reaction database updated.")
        except Exception as e:
            st.error(f"Error: {e}")
