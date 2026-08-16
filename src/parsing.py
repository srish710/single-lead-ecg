from __future__ import annotations

import pandas as pd

REQUIRED_COLUMNS = [
    "Well",
    "Sample_Name",
    "Group",
    "Target_Name",
    "Task",
    "Technical_Replicate",
    "Ct",
]


def load_ct_data(filepath_or_buffer) -> pd.DataFrame:
    name = getattr(filepath_or_buffer, "name", str(filepath_or_buffer))
    if str(name).lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(filepath_or_buffer)
    else:
        df = pd.read_csv(filepath_or_buffer)

    df.columns = [c.strip() for c in df.columns]
    validate_columns(df)

    df["Ct"] = pd.to_numeric(df["Ct"], errors="coerce")

    return df[REQUIRED_COLUMNS]


def validate_columns(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Uploaded file is missing required column(s): {missing}. "
            f"Expected columns: {REQUIRED_COLUMNS}"
        )