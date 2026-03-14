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

import pandas as pd
import numpy as np
import pathlib

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
ALT_CSV = DATA_DIR / "reactors_alt.csv"
REACTOR_CSV = DATA_DIR / "reactors.csv"

st.title("⚙️ Admin – Reactor Import")
st.caption(
    "Convert the transposed master file **reactors_alt.csv** into the "
    "row-per-reactor format used by the Reactor Database."
)

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
    st.dataframe(preview, use_container_width=False, hide_index=True)

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
st.dataframe(converted_df, use_container_width=False, hide_index=True)

with st.expander("All fields from source (including unmapped)", expanded=False):
    st.dataframe(extras_df, use_container_width=False, hide_index=True)

# ─────────────────────────────────────────────────────────────────────────
# 4. Load current reactor DB for comparison
# ─────────────────────────────────────────────────────────────────────────
st.header("3 · Current Reactor Database")

if REACTOR_CSV.exists():
    current_df = pd.read_csv(REACTOR_CSV)
    st.dataframe(current_df, use_container_width=False, hide_index=True)
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
    st.dataframe(result, use_container_width=False, hide_index=True)

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
