"""
Admin – Import / Sync Reactor Master File
==========================================
Converts the transposed *reactors_alt.csv* master file (properties as rows,
reactors as columns) into the row-per-reactor format used by *reactors.csv*.

Workflow
--------
1. Upload or reload ``reactors_alt.csv`` from ``data/``.
2. Preview the transposed source and the converted output.
3. Edit the field mapping if columns change.
4. Choose **merge** (update existing + add new) or **replace** (overwrite).
5. Click **Apply** to write ``data/reactors.csv``.
"""

import streamlit as st
import hmac

import pandas as pd
import numpy as np
import pathlib
import graphviz

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
ALT_CSV = DATA_DIR / "reactors_alt.csv"
REACTOR_CSV = DATA_DIR / "reactors.csv"

st.title("🛠️ Admin Tools")

# ── Authentication gate ──────────────────────────────────────────────────
_ADMIN_USER = "admin"
_ADMIN_PASS = "admin_tak_2026"


def _check_password() -> bool:
    """Return True if the user has entered valid credentials."""
    if st.session_state.get("admin_authenticated"):
        return True

    with st.form("admin_login"):
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

st.caption(
    "Convert the transposed master file **reactors_alt.csv** into the "
    "row-per-reactor format used by the Reactor Database."
)

with st.expander("Reactor Import Tool", expanded=False):
    # ─────────────────────────────────────────────────────────────────────────
    # 1. Load the transposed source
    # ─────────────────────────────────────────────────────────────────────────
    st.header("1 · Source File")

    upload = st.file_uploader(
        "Upload a replacement reactors_alt.csv (or leave blank to use the one in data/)",
        type=["csv"],
        key="alt_upload",
    )

    if upload is not None:
        raw = pd.read_csv(upload, header=None)
        st.success(f"Uploaded file: {upload.name}  ({raw.shape[0]} rows × {raw.shape[1]-2} reactors)")
    elif ALT_CSV.exists():
        raw = pd.read_csv(ALT_CSV, header=None)
        st.info(f"Using **{ALT_CSV.name}** from data/  ({raw.shape[0]} rows × {raw.shape[1]-2} reactors)")
    else:
        st.warning("No reactors_alt.csv found. Upload one above.")
        st.stop()

    # ── Parse ────────────────────────────────────────────────────────────────
    # Col 0 = property label, Col 1 = units, Cols 2+ = reactor values
    prop_labels = raw.iloc[:, 0].astype(str).str.strip()
    units = raw.iloc[:, 1].astype(str).str.strip()
    n_reactors = raw.shape[1] - 2

    # Build reactor name from Owner + Tag rows
    owners = raw.iloc[0, 2:].astype(str).str.strip()
    tags = raw.iloc[1, 2:].astype(str).str.strip()
    reactor_names = [
        f"{o} – {t}" if o != "nan" and t != "nan" else None
        for o, t in zip(owners, tags)
    ]

    # Show raw preview
    with st.expander("Raw source preview (first 20 rows)", expanded=False):
        preview = raw.iloc[:20].copy()
        preview.columns = ["Property", "Units"] + [n or f"Col {i}" for i, n in enumerate(reactor_names)]
        st.dataframe(preview, width='content', hide_index=True)

    # ─────────────────────────────────────────────────────────────────────────
    # 2. Transpose into row-per-reactor dataframe
    # ─────────────────────────────────────────────────────────────────────────

    def _val(series_idx, col_idx):
        """Extract a scalar value from the raw dataframe, cleaning strings."""
        v = raw.iloc[series_idx, col_idx]
        if pd.isna(v):
            return np.nan
        s = str(v).strip()
        if s in ("", "nan", "NaN", "-", "#REF!", "#DIV/0!", "#N/A"):
            return np.nan
        return s


    def _num(series_idx, col_idx):
        """Extract a numeric value, returning NaN for non-numeric."""
        v = _val(series_idx, col_idx)
        if v is np.nan or v is None:
            return np.nan
        try:
            return float(v)
        except (ValueError, TypeError):
            return np.nan


    # Row index lookup (property label → row number)
    _row_idx = {lbl: i for i, lbl in enumerate(prop_labels)}


    def _r(label: str) -> int:
        """Get row index for a property label (fuzzy-ish)."""
        if label in _row_idx:
            return _row_idx[label]
        # Try stripped trailing spaces
        for k, v in _row_idx.items():
            if k.strip() == label.strip():
                return v
        return -1


    # ── Scale mapping ────────────────────────────────────────────────────────
    SCALE_MAP = {
        "Commercial": "Manufacturing",
        "Pilot": "Pilot",
        "Lab": "Lab",
        "Kilo": "Kilo",
    }

    # ── Build converted rows ────────────────────────────────────────────────
    converted_rows = []
    extra_fields_all = []  # For the "all fields" view

    for c in range(2, 2 + n_reactors):
        name = reactor_names[c - 2]
        if name is None:
            continue

        # Core fields for reactors.csv
        owner_val = _val(0, c) or ""
        scale_raw = _val(_r("Scale"), c) or ""
        scale = SCALE_MAP.get(scale_raw, scale_raw)

        row = {
            "reactor_name": name,
            "owner": owner_val.replace(" CC", "").replace(" HP", ""),
            "type": "Batch",
            "scale": scale,
            "D_tank_m": _num(_r("Internal diameter (ID)"), c) / 1000 if not np.isnan(_num(_r("Internal diameter (ID)"), c)) else np.nan,
            "H_m": _num(_r("Straight wall height (tan-tan)"), c) / 1000 if not np.isnan(_num(_r("Straight wall height (tan-tan)"), c)) else np.nan,
            "H_max_m": _num(_r("Max height"), c) / 1000 if not np.isnan(_num(_r("Max height"), c)) else np.nan,
            "D_imp_m": _num(_r("Impeller 1 diameter"), c) / 1000 if not np.isnan(_num(_r("Impeller 1 diameter"), c)) else np.nan,
            "impeller_type": _val(_r("Impeller 1 type"), c),
            "Np": _num(_r("Impeller 1 power number"), c) if not np.isnan(_num(_r("Impeller 1 power number"), c)) else np.nan,  # Not in alt file – must be set manually
            "Nq": np.nan,
            "N_rpm_min": _num(_r("Min agitation speed"), c),
            "N_rpm_max": _num(_r("Max agitation speed"), c),
            "N_rps": _num(_r("Max agitation speed"), c) / 60 if not np.isnan(_num(_r("Max agitation speed"), c)) else np.nan,
            "V_L_min": _num(_r("Volume min"), c),
            "V_L_max": _num(_r("Volume max"), c),
            "V_L": _num(_r("Volume (nominal)"), c),
            "material": _val(_r("Shell material"), c),
            "baffles": _num(_r("Baffle count"), c),
            "bottom_dish": _val(_r("Bottom dish type"), c),
            "top_dish": _val(_r("Top dish type"), c),
            "impeller_count": _num(_r("Impeller count"), c),
            "imp1_clearance_m": _num(_r("Impeller 1 clearance"), c) / 1000 if not np.isnan(_num(_r("Impeller 1 clearance"), c)) else np.nan,
            "imp1_height_m": _num(_r("Impeller 1 height"), c) / 1000 if not np.isnan(_num(_r("Impeller 1 height"), c)) else np.nan,
            "D_imp2_m": _num(_r("Impeller 2 diameter"), c) / 1000 if not np.isnan(_num(_r("Impeller 2 diameter"), c)) else np.nan,
            "Np2": _num(_r("Impeller 2 power number"), c) if not np.isnan(_num(_r("Impeller 2 power number"), c)) else np.nan,
            "imp2_clearance_m": _num(_r("Impeller 2 clearance"), c) / 1000 if not np.isnan(_num(_r("Impeller 2 clearance"), c)) else np.nan,
            "imp2_height_m": _num(_r("Impeller 2 height"), c) / 1000 if not np.isnan(_num(_r("Impeller 2 height"), c)) else np.nan,
            "D_imp3_m": _num(_r("Impeller 3 diameter"), c) / 1000 if not np.isnan(_num(_r("Impeller 3 diameter"), c)) else np.nan,
            "Np3": _num(_r("Impeller 3 power number"), c) if not np.isnan(_num(_r("Impeller 3 power number"), c)) else np.nan,
            "imp3_clearance_m": _num(_r("Impeller 3 clearance"), c) / 1000 if not np.isnan(_num(_r("Impeller 3 clearance"), c)) else np.nan,
            "imp3_height_m": _num(_r("Impeller 3 height"), c) / 1000 if not np.isnan(_num(_r("Impeller 3 height"), c)) else np.nan,
            "Zwietering_S": np.nan,
            "GMB_z": np.nan,
            "wall_thickness_mm": _num(_r("Shell thickness"), c),
            "OD_m": np.nan,
            "knuckle_radius_m": np.nan,
            "notes": "",
        }
        converted_rows.append(row)

        # Extra fields (everything from alt that doesn't map directly)
        extras = {"reactor_name": name}
        for ri in range(len(prop_labels)):
            lbl = prop_labels[ri]
            if lbl == "nan":
                continue
            extras[lbl] = _val(ri, c)
        extra_fields_all.append(extras)

    converted_df = pd.DataFrame(converted_rows)
    extras_df = pd.DataFrame(extra_fields_all)

    # ─────────────────────────────────────────────────────────────────────────
    # 3. Preview converted data
    # ─────────────────────────────────────────────────────────────────────────
    st.header("2 · Converted Preview")

    st.subheader("Core fields (reactors.csv format)")
    st.dataframe(converted_df, width='content', hide_index=True)

    with st.expander("All fields from source (including unmapped)", expanded=False):
        st.dataframe(extras_df, width='content', hide_index=True)

    # ─────────────────────────────────────────────────────────────────────────
    # 4. Load current reactor DB for comparison
    # ─────────────────────────────────────────────────────────────────────────
    st.header("3 · Current Reactor Database")

    if REACTOR_CSV.exists():
        current_df = pd.read_csv(REACTOR_CSV)
        st.dataframe(current_df, width='content', hide_index=True)
        st.caption(f"{len(current_df)} reactors in current database.")
    else:
        current_df = pd.DataFrame()
        st.info("No reactors.csv found – a new one will be created.")

    # ─────────────────────────────────────────────────────────────────────────
    # 5. Merge or Replace
    # ─────────────────────────────────────────────────────────────────────────
    st.header("4 · Apply Changes")

    mode = st.radio(
        "Import mode",
        ["Merge (update existing + add new)", "Replace (overwrite entirely)"],
        help=(
            "**Merge** matches reactors by name: existing rows are updated with "
            "non-NaN values from the import, new reactors are appended, and "
            "reactors only in the current DB are kept.  "
            "**Replace** overwrites reactors.csv completely with the imported data."
        ),
    )

    preserve_np = st.checkbox(
        "Preserve manually-set Np / Nq / Zwietering_S / GMB_z values from current DB",
        value=True,
        help="The alt file doesn't contain these; keep the values you've already entered.",
    )

    # Manual fields to preserve
    _PRESERVE_COLS = ["Np", "Nq", "Np2", "Np3", "Zwietering_S", "GMB_z",
                      "V_L_min", "V_L_max"]

    if st.button("🚀 Apply to reactors.csv", type="primary"):
        if "Replace" in mode:
            result = converted_df.copy()
            # Preserve manual fields from current DB if requested
            if preserve_np and not current_df.empty:
                for _, cur_row in current_df.iterrows():
                    mask = result["reactor_name"] == cur_row["reactor_name"]
                    if mask.any():
                        for col in _PRESERVE_COLS:
                            if col in cur_row and pd.notna(cur_row[col]):
                                result.loc[mask, col] = cur_row[col]
        else:
            # Merge mode
            if current_df.empty:
                result = converted_df.copy()
            else:
                result = current_df.copy()
                for _, new_row in converted_df.iterrows():
                    name = new_row["reactor_name"]
                    mask = result["reactor_name"] == name
                    if mask.any():
                        # Update existing: overwrite only non-NaN imported values
                        idx = result.index[mask][0]
                        for col in converted_df.columns:
                            if col == "reactor_name":
                                continue
                            new_val = new_row[col]
                            if pd.notna(new_val) and str(new_val).strip() != "":
                                if preserve_np and col in _PRESERVE_COLS:
                                    # Keep existing if it has a value
                                    if pd.notna(result.at[idx, col]):
                                        continue
                                result.at[idx, col] = new_val
                    else:
                        # Append new reactor
                        result = pd.concat([result, new_row.to_frame().T], ignore_index=True)

        # Ensure all expected columns exist
        for col in converted_df.columns:
            if col not in result.columns:
                result[col] = np.nan

        # Write
        result.to_csv(REACTOR_CSV, index=False)
        # Update session state if cached
        if "reactor_db" in st.session_state:
            st.session_state.reactor_db = result

        st.success(f"✅ Wrote {len(result)} reactors to **{REACTOR_CSV.name}**.")
        st.dataframe(result, width='content', hide_index=True)

# ─────────────────────────────────────────────────────────────────────────
# 6. Field mapping reference
# ─────────────────────────────────────────────────────────────────────────
with st.expander("📋 Field mapping reference", expanded=False):
    st.markdown("""
| reactors.csv column | reactors_alt.csv row | Transform |
|---|---|---|
| `reactor_name` | Owner + Tag | `"{Owner} – {Tag}"` |
| `owner` | Owner | Strip site suffix (CC/HP) |
| `scale` | Scale | Commercial → Manufacturing |
| `D_tank_m` | Internal diameter (ID) | mm → m ÷ 1000 |
| `H_m` | Straight wall height (tan-tan) | mm → m ÷ 1000 |
| `H_max_m` | Max height | mm → m ÷ 1000 |
| `D_imp_m` | Impeller 1 diameter | mm → m ÷ 1000 |
| `impeller_type` | Impeller 1 type | direct |
| `N_rpm_min` | Min agitation speed | direct (rpm) |
| `N_rpm_max` | Max agitation speed | direct (rpm) |
| `V_L` | Volume (nominal) | direct (L) |
| `impeller_count` | Impeller count | direct |
| `imp1_clearance_m` | Impeller 1 clearance | mm → m ÷ 1000 |
| `imp1_height_m` | Impeller 1 height | mm → m ÷ 1000 |
| `D_imp2_m` | Impeller 2 diameter | mm → m ÷ 1000 |
| `imp2_clearance_m` | Impeller 2 clearance | mm → m ÷ 1000 |
| `imp2_height_m` | Impeller 2 height | mm → m ÷ 1000 |
| `D_imp3_m` | Impeller 3 diameter | mm → m ÷ 1000 |
| `imp3_clearance_m` | Impeller 3 clearance | mm → m ÷ 1000 |
| `top_dish` | Top dish type | direct |
| `bottom_dish` | Bottom dish type | direct |
| `wall_thickness_mm` | Shell thickness | direct (mm) |
| `baffles` | Baffle count | direct |
| `material` | Shell material | direct |

**Not in alt file** (must be set manually in Reactor Database):
`Np`, `Nq`, `Np2`, `Np3`, `Zwietering_S`, `GMB_z`, `V_L_min`, `V_L_max`,
`OD_m`, `knuckle_radius_m`
""")

# ─────────────────────────────────────────────────────────────────────────
# 7. Bourne Protocol Decision Tree
# ─────────────────────────────────────────────────────────────────────────
st.divider()
st.header("5 · Bourne Protocol – Decision Tree")
st.caption(
    "Generate an interactive decision-tree diagram of the Bourne mixing-sensitivity "
    "screening protocol (Bourne 2003 / Sarafinas 2018)."
)

_BOURNE_DOT = """
    digraph bourne_protocol {
        rankdir=TB
        fontname="Arial"
        node [fontname="Arial" fontsize=10 style=filled shape=box
              fillcolor="#F5F5F5" color="#BDBDBD" margin="0.15,0.08"]
        edge [fontname="Arial" fontsize=9 color="#616161"]

        /* Start */
        START [label="🧫  Bourne Protocol\nMixing Sensitivity Screening" shape=Mrecord
               fillcolor="#4CAF50" fontcolor=white fontsize=11]

        /* Test 1 */
        T1 [label="Test 1\nVary Impeller Speed\n(~100× change in P/V)" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9" fontsize=10]
        T1_SETUP [label="3 speeds:\nLow (0.1× P/V)\nCenter (1× P/V)\nHigh (10× P/V)"
                  fillcolor="#FFF3E0" color="#FFB74D"]
        T1_NO [label="✅ Mixing is NOT critical\nover tested P/V range\n— Protocol Complete"
               fillcolor="#E8F5E9" color="#81C784" fontsize=10]

        /* Test 2 */
        T2 [label="Test 2\nVary Feed Rate / Feed Time\n(9× range on flow rate)" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9" fontsize=10]
        T2_SETUP [label="3 feed times:\nFast (1/3× t_feed)\nCenter (1× t_feed)\nSlow (3× t_feed)"
                  fillcolor="#FFF3E0" color="#FFB74D"]
        T2_NO [label="🔬 Micromixing\nControls Process\n\n→ Hold local ε constant\nat feed point on scale-up"
               fillcolor="#E8EAF6" color="#7986CB" fontsize=10]

        /* Test 3 */
        T3 [label="Test 3\nVary Feed Location\n(surface → impeller zone)" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9" fontsize=10]
        T3_SETUP [label="3 locations:\nSurface (ε_loc ≈ 0.1 ε_avg)\nSub-surface (ε_loc ≈ 1 ε_avg)\nImpeller zone (ε_loc ≈ 3 ε_avg)"
                  fillcolor="#FFF3E0" color="#FFB74D"]
        T3_MACRO [label="🌀 Macromixing\nControls Process\n\n→ Reduce blend time\n→ High-efficiency impellers\n→ Multiple impellers / static mixers"
                  fillcolor="#FFEBEE" color="#EF9A9A" fontsize=10]
        T3_MESO [label="📐 Mesomixing\nControls Process\n\n→ Keep impeller speed constant\n→ Extend feed time\n→ Multiple feed points"
                 fillcolor="#FFF8E1" color="#FFD54F" fontsize=10]

        /* Confirmatory */
        CONFIRM [label="Confirmatory Experiments\n(Optional)" shape=Mrecord
                 fillcolor="#F3E5F5" color="#CE93D8" fontsize=10]
        CONF_A [label="A – Number of Feed Points\nInsensitive → Micromixing\nSensitive → Mesomixing"
                fillcolor="#F3E5F5" color="#CE93D8"]
        CONF_B [label="B – Viscosity Change\nSensitive → Micromixing\nInsensitive → Not micromixing"
                fillcolor="#F3E5F5" color="#CE93D8"]

        /* Edges */
        START -> T1
        T1 -> T1_SETUP [label="Experimental\nconditions" style=dashed]

        T1 -> T1_NO     [label="No change in\nprocess response"]
        T1 -> T2         [label="Process response\nchanged"]

        T2 -> T2_SETUP [label="Experimental\nconditions" style=dashed]

        T2 -> T2_NO     [label="No change\n(feed-rate insensitive)"]
        T2 -> T3         [label="Process response\nchanged"]

        T3 -> T3_SETUP [label="Experimental\nconditions" style=dashed]

        T3 -> T3_MACRO  [label="No change\n(location insensitive)"]
        T3 -> T3_MESO   [label="Process response\nchanged"]

        /* Confirmatory links */
        T2_NO -> CONFIRM [label="Optional\nvalidation" style=dashed]
        T3_MESO -> CONFIRM [label="Optional\nvalidation" style=dashed]
        CONFIRM -> CONF_A
        CONFIRM -> CONF_B
    }
"""

if st.button("📊 Generate Bourne Protocol Decision Tree", key="gen_bourne_tree"):
    st.session_state["_show_bourne_tree"] = True

if st.session_state.get("_show_bourne_tree"):
    st.graphviz_chart(_BOURNE_DOT, width='content')

    # Export buttons
    _bourne_graph = graphviz.Source(_BOURNE_DOT)
    _col_png, _col_svg, _ = st.columns([1, 1, 4])
    with _col_png:
        st.download_button(
            "⬇️ Download PNG",
            data=_bourne_graph.pipe(format="png"),
            file_name="bourne_protocol_decision_tree.png",
            mime="image/png",
        )
    with _col_svg:
        st.download_button(
            "⬇️ Download SVG",
            data=_bourne_graph.pipe(format="svg"),
            file_name="bourne_protocol_decision_tree.svg",
            mime="image/svg+xml",
        )

# ─────────────────────────────────────────────────────────────────────────
# 8. Mixing Sensitivity Protocol Decision Tree
# ─────────────────────────────────────────────────────────────────────────
st.divider()
st.header("6 · Mixing Sensitivity Protocol – Decision Tree")
st.caption(
    "Generate an interactive decision-tree diagram of the Mixing Sensitivity "
    "Protocol (pre-screening → kinetics → phases → competing reactions → "
    "heat transfer → mixing-time comparison)."
)

_MSP_DOT = """
    digraph protocol {
        rankdir=TB
        fontname="Arial"
        node [fontname="Arial" fontsize=10 style=filled shape=box
              fillcolor="#F5F5F5" color="#BDBDBD" margin="0.15,0.08"]
        edge [fontname="Arial" fontsize=9 color="#616161"]

        /* start / end nodes */
        START [label="🧭  Start Protocol" shape=Mrecord
               fillcolor="#4CAF50" fontcolor=white fontsize=11]
        SUM   [label="6 · Summary & \\nRecommendations" shape=Mrecord
               fillcolor="#2196F3" fontcolor=white fontsize=11]

        /* Step 0 – Bourne pre-screening */
        BOURNE [label="0 · Bourne Protocol\\nPart 1 suggested pre-screen" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        BOURNE_DO [label="Run Bourne Protocol\\nPart 1 (quick screen)\\n→ vary P/V (i.e stir speed)"
                   fillcolor="#FFF3E0" color="#FFB74D"]
        BOURNE_OK [label="✅ Pre-screen\\ncomplete"
                   fillcolor="#E8F5E9" color="#81C784"]

        /* Step 1 – Kinetics */
        K     [label="1 · Kinetics\\navailable?" shape=diamond
               fillcolor="#E3F2FD" color="#90CAF9"]
        KACT  [label="Measure kinetics\\n& calorimetry\\n→ add to Reaction DB"
               fillcolor="#FFF3E0" color="#FFB74D"]
        KAPPROX [label="Select approximate\\nkinetics from DB\\n(similar reaction)"
                 fillcolor="#FFF8E1" color="#FFD54F"]
        KSEL  [label="Select reaction\\nfrom database"
               fillcolor="#E8F5E9" color="#81C784"]

        /* Step 2 – Phases */
        PH    [label="2 · How many\\nphases?" shape=diamond
               fillcolor="#E3F2FD" color="#90CAF9"]
        PH_OK [label="✅ Mass transfer\\nnot a factor"
               fillcolor="#E8F5E9" color="#81C784"]
        PH_MT [label="⚠️ Interphase mass\\ntransfer may limit\\n→ characterise kLa, k_SL"
               fillcolor="#FFF8E1" color="#FFD54F"]

        /* Step 3 – Competing reactions */
        CR      [label="3 · Competing\\nreactions?" shape=diamond
                 fillcolor="#E3F2FD" color="#90CAF9"]
        MESO    [label="⚠️ Micro/meso-mixing limitation possible\\n→ Bourne Protocol"
                 fillcolor="#FFF8E1" color="#EF9A9A"]
        MESO_OK [label="✅ Micro/meso-mixing\\nnot a factor"
                 fillcolor="#E8F5E9" color="#81C784"]

        /* Step 4 – Heat */
        HT     [label="4 · ΔH data\\navailable?" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        HT_MEASURE [label="Perform calorimetry\\n→ measure ΔH"
                    fillcolor="#FFF3E0" color="#FFB74D"]
        HT_EST [label="Estimate ΔH from\\nsimilar reaction"
                fillcolor="#FFF8E1" color="#FFD54F"]
        HT_CHK [label="|ΔH| ≥ 50\\nkJ/mol?" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        HT_HOT [label="🔴 Heat-sensitive\\n→ run heat balance"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        HT_OK  [label="🟢 Modest\\nexothermicity"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Step 5 – Mixing time */
        MIX    [label="5 · Reaction\\ntime t_rxn" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        MICRO  [label="🔴 Micromixing\\nlikely sensitive"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        MACRO  [label="🟡 Macromixing\\nmay matter at scale"
                fillcolor="#FFF8E1" color="#FFD54F"]
        MIX_OK [label="🟢 Mixing unlikely\\nto be limiting"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Edges — Step 0 */
        START -> BOURNE

        BOURNE -> BOURNE_DO [label="Not done\\nyet"]
        BOURNE -> BOURNE_OK [label="Done /\\nSkip"]
        BOURNE_DO -> BOURNE_OK [label="Complete"]
        BOURNE_OK -> K

        /* Edges — Step 1 */
        K -> KACT    [label="No kinetics"]
        K -> KAPPROX [label="Approximate"]
        K -> KSEL    [label="Yes"]
        KACT -> K    [label="Data obtained\\n→ repeat Step 1"]
        KAPPROX -> PH
        KSEL -> PH

        /* Edges — Step 2 */
        PH -> PH_OK [label="Single\\nliquid"]
        PH -> PH_MT [label="Multi-phase\\n(G / L / S)"]

        PH_OK -> CR
        PH_MT -> CR

        /* Edges — Step 3 */
        CR -> MESO    [label="Yes /\\nNot sure"]
        CR -> MESO_OK [label="No"]

        MESO    -> HT
        MESO_OK -> HT

        /* Edges — Step 4 */
        HT -> HT_MEASURE [label="No ΔH"]
        HT -> HT_EST     [label="Estimate"]
        HT -> HT_CHK     [label="Yes"]
        HT_MEASURE -> HT [label="Data obtained\\n→ repeat Step 4"]
        HT_EST -> HT_CHK
        HT_CHK -> HT_HOT [label="Yes"]
        HT_CHK -> HT_OK  [label="No"]

        /* Edges — Step 5 */
        HT_HOT -> MIX
        HT_OK  -> MIX

        MIX -> MICRO  [label="< 1 s"]
        MIX -> MACRO  [label="1 – 10 s"]
        MIX -> MIX_OK [label="> 10 s"]

        MICRO  -> SUM
        MACRO  -> SUM
        MIX_OK -> SUM
    }
"""

if st.button("🧭 Generate Mixing Sensitivity Protocol Decision Tree", key="gen_msp_tree"):
    st.session_state["_show_msp_tree"] = True

if st.session_state.get("_show_msp_tree"):
    st.graphviz_chart(_MSP_DOT, width='content')

    # Export buttons
    _msp_graph = graphviz.Source(_MSP_DOT)
    _col_png2, _col_svg2, _ = st.columns([1, 1, 4])
    with _col_png2:
        st.download_button(
            "⬇️ Download PNG",
            data=_msp_graph.pipe(format="png"),
            file_name="mixing_sensitivity_protocol.png",
            mime="image/png",
        )
    with _col_svg2:
        st.download_button(
            "⬇️ Download SVG",
            data=_msp_graph.pipe(format="svg"),
            file_name="mixing_sensitivity_protocol.svg",
            mime="image/svg+xml",
        )

# ─────────────────────────────────────────────────────────────────────────
# 9. Crystallization Sensitivity Protocol Decision Tree
# ─────────────────────────────────────────────────────────────────────────
st.divider()
st.header("7 · Crystallization Sensitivity Protocol – Decision Tree")
st.caption(
    "Generate an interactive decision-tree diagram of the Crystallization "
    "Sensitivity Protocol (pre-screening → crystallization parameters → "
    "type & supersaturation → nucleation/growth competition → heat transfer → "
    "Damköhler analysis → summary)."
)

_CSP_DOT = """
    digraph cryst_protocol {
        rankdir=TB
        fontname="Arial"
        node [fontname="Arial" fontsize=10 style=filled shape=box
              fillcolor="#F5F5F5" color="#BDBDBD" margin="0.15,0.08"]
        edge [fontname="Arial" fontsize=9 color="#616161"]

        /* Start / End */
        START [label="💎  Start Crystallization\\nSensitivity Protocol" shape=Mrecord
               fillcolor="#4CAF50" fontcolor=white fontsize=11]
        SUM   [label="6 · Summary &\\nRecommendations" shape=Mrecord
               fillcolor="#2196F3" fontcolor=white fontsize=11]

        /* Step 0 – Bourne pre-screening */
        BOURNE [label="0 · Bourne Pre-Screen\\n(vary impeller speed/feed rate)" shape=diamond
                fillcolor="#E3F2FD" color="#90CAF9"]
        B_SENS [label="🔴 Mixing sensitivity\\nconfirmed experimentally"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        B_OK   [label="🟢 No sensitivity\\nobserved at lab scale"
                fillcolor="#E8F5E9" color="#81C784"]
        B_SKIP [label="⚪ Skipped\\n(proceed with theory)"
                fillcolor="#F5F5F5" color="#BDBDBD"]

        /* Step 1 – Crystallization parameters */
        CP [label="1 · Crystallization\\nParameters" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9"]
        CP_DB  [label="Select from\\nCrystallization DB"
                fillcolor="#E8F5E9" color="#81C784"]
        CP_MAN [label="Enter parameters\\nmanually"
                fillcolor="#FFF3E0" color="#FFB74D"]
        CP_OK  [label="Parameters loaded:\\nt_ind, k_g, σ, MSZW,\\nΔH_cryst, polymorph info"
                fillcolor="#E8F5E9" color="#81C784"]
        CP_DER [label="Derived times:\\nt_G = L_target / G\\nG = k_g · σ^g"
                fillcolor="#FFF8E1" color="#FFD54F"]

        /* Step 2 – Type & supersaturation */
        TY [label="2 · Crystallization Type\\n& Supersaturation" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9"]
        TY_COOL [label="🟢 Cooling / Evap\\nLow inherent risk\\n→ macromixing concern"
                 fillcolor="#E8F5E9" color="#81C784"]
        TY_FEED [label="🔴 Anti-solvent /\\nReactive / pH-shift\\nHigh inherent risk\\n→ micro/mesomixing"
                 fillcolor="#FFEBEE" color="#EF9A9A"]
        SIG [label="σ / σ_max\\nratio?" shape=diamond
             fillcolor="#E3F2FD" color="#90CAF9"]
        SIG_HI [label="⚠️ Near metastable\\nlimit — σ spike risk"
                fillcolor="#FFF8E1" color="#FFD54F"]
        SIG_OK [label="🟢 Well within\\nmetastable zone"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Step 3 – Nucleation / growth competition */
        NG [label="3 · Nucleation / Growth\\nCompetition" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9"]
        NG_INFO [label="B ∝ σ^n  (n ~ 2–10)\\nG ∝ σ^g  (g ~ 1–2)\\nn >> g → nucleation\\nmore σ-sensitive"
                 fillcolor="#FFF8E1" color="#FFD54F"]
        NG_MICRO [label="🔴 Feed-sensitive:\\nmicro/mesomixing\\ncontrols local σ"
                  fillcolor="#FFEBEE" color="#EF9A9A"]
        NG_FAST  [label="🟡 Fast t_ind\\n(< 10 s) — sensitive\\neven without feed"
                  fillcolor="#FFF8E1" color="#FFD54F"]
        NG_OK    [label="🟢 Slow nucleation\\n→ low mixing risk"
                  fillcolor="#E8F5E9" color="#81C784"]
        POLY [label="Polymorphism?" shape=diamond
              fillcolor="#E3F2FD" color="#90CAF9"]
        POLY_Y [label="⚠️ Ostwald's rule:\\nmetastable form may\\nnucleate at high σ"
                fillcolor="#FFF8E1" color="#FFD54F"]
        SEED [label="Seeded?" shape=diamond
              fillcolor="#E3F2FD" color="#90CAF9"]
        SEED_Y [label="🟢 Seed surface area\\nreduces nucleation\\ndependence (2–5× less\\nmixing-sensitive)"
                fillcolor="#E8F5E9" color="#81C784"]
        SEED_N [label="🟡 Unseeded:\\nall nucleation is\\nprimary/secondary"
                fillcolor="#FFF8E1" color="#FFD54F"]

        /* Step 4 – Heat transfer */
        HT [label="4 · Heat Transfer\\nScreening" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9"]
        HT_HI  [label="🔴 |ΔH| ≥ 40 kJ/mol\\nHeat removal may limit\\ncooling/feed rates"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        HT_MED [label="🟡 |ΔH| 20–40 kJ/mol\\nCheck at scale"
                fillcolor="#FFF8E1" color="#FFD54F"]
        HT_LO  [label="🟢 |ΔH| < 20 kJ/mol\\nHeat not limiting"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Step 5 – Damköhler analysis */
        DA [label="5 · Damköhler Analysis\\nt_mix vs t_cryst" shape=diamond
            fillcolor="#E3F2FD" color="#90CAF9"]
        DA_INFO [label="Da_micro = t_E / t_ind\\nDa_macro = θ₉₅ / t_G"
                 fillcolor="#FFF8E1" color="#FFD54F"]
        DA_HI  [label="🔴 Da > 1\\nMixing-sensitive"
                fillcolor="#FFEBEE" color="#EF9A9A"]
        DA_MED [label="🟡 Da 0.1 – 1\\nPotentially sensitive"
                fillcolor="#FFF8E1" color="#FFD54F"]
        DA_LO  [label="🟢 Da < 0.1\\nNot sensitive"
                fillcolor="#E8F5E9" color="#81C784"]

        /* Edges — Step 0 */
        START -> BOURNE

        BOURNE -> B_SENS [label="Sensitivity\\nconfirmed"]
        BOURNE -> B_OK   [label="No change\\nobserved"]
        BOURNE -> B_SKIP [label="Skip"]
        B_SENS -> CP
        B_OK   -> CP
        B_SKIP -> CP

        /* Edges — Step 1 */
        CP -> CP_DB  [label="From DB"]
        CP -> CP_MAN [label="Manual"]
        CP_DB  -> CP_OK
        CP_MAN -> CP_OK
        CP_OK  -> CP_DER
        CP_DER -> TY

        /* Edges — Step 2 */
        TY -> TY_COOL [label="Cooling /\\nEvaporative"]
        TY -> TY_FEED [label="Anti-solvent /\\nReactive / pH-shift"]
        TY_COOL -> SIG
        TY_FEED -> SIG
        SIG -> SIG_HI [label="σ/σ_max > 0.4"]
        SIG -> SIG_OK [label="σ/σ_max ≤ 0.4"]
        SIG_HI -> NG
        SIG_OK -> NG

        /* Edges — Step 3 */
        NG -> NG_INFO [style=dashed label="Background"]
        NG -> NG_MICRO [label="Feed-sensitive\\ntype"]
        NG -> NG_FAST  [label="t_ind < 10 s"]
        NG -> NG_OK    [label="Slow nucleation\\n& no feed"]
        NG_MICRO -> POLY
        NG_FAST  -> POLY
        NG_OK    -> POLY
        POLY -> POLY_Y [label="Yes"]
        POLY -> SEED   [label="No"]
        POLY_Y -> SEED
        SEED -> SEED_Y [label="Yes"]
        SEED -> SEED_N [label="No"]
        SEED_Y -> HT
        SEED_N -> HT

        /* Edges — Step 4 */
        HT -> HT_HI  [label="|ΔH| ≥ 40"]
        HT -> HT_MED [label="20–40"]
        HT -> HT_LO  [label="< 20 or\\nunknown"]
        HT_HI  -> DA
        HT_MED -> DA
        HT_LO  -> DA

        /* Edges — Step 5 */
        DA -> DA_INFO [style=dashed label="Definitions"]
        DA -> DA_HI  [label="Da > 1"]
        DA -> DA_MED [label="0.1 – 1"]
        DA -> DA_LO  [label="Da < 0.1"]

        DA_HI  -> SUM
        DA_MED -> SUM
        DA_LO  -> SUM
    }
"""

if st.button("💎 Generate Crystallization Sensitivity Protocol Decision Tree", key="gen_csp_tree"):
    st.session_state["_show_csp_tree"] = True

if st.session_state.get("_show_csp_tree"):
    st.graphviz_chart(_CSP_DOT, width='content')

    # Export buttons
    _csp_graph = graphviz.Source(_CSP_DOT)
    _col_png3, _col_svg3, _ = st.columns([1, 1, 4])
    with _col_png3:
        st.download_button(
            "⬇️ Download PNG",
            data=_csp_graph.pipe(format="png"),
            file_name="crystallization_sensitivity_protocol.png",
            mime="image/png",
        )
    with _col_svg3:
        st.download_button(
            "⬇️ Download SVG",
            data=_csp_graph.pipe(format="svg"),
            file_name="crystallization_sensitivity_protocol.svg",
            mime="image/svg+xml",
        )
