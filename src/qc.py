from __future__ import annotations

import numpy as np
import pandas as pd

REPLICATE_STD_THRESHOLD = 0.5
NTC_CONTAMINATION_CT = 35.0


def flag_undetermined(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["flag_undetermined"] = df["Ct"].isna()
    return df


def flag_replicate_variance(df: pd.DataFrame, threshold: float = REPLICATE_STD_THRESHOLD) -> pd.DataFrame:
    df = df.copy()
    group_std = (
        df.groupby(["Sample_Name", "Target_Name"])["Ct"]
        .transform(lambda s: s.std(ddof=1))
    )
    df["replicate_std"] = group_std
    df["flag_high_variance"] = group_std > threshold
    return df


def flag_ntc_contamination(df: pd.DataFrame, ct_threshold: float = NTC_CONTAMINATION_CT) -> pd.DataFrame:
    df = df.copy()
    is_ntc = df["Task"].str.upper().eq("NTC")
    df["flag_ntc_contamination"] = is_ntc & df["Ct"].notna() & (df["Ct"] < ct_threshold)
    return df


def run_qc(df: pd.DataFrame) -> pd.DataFrame:
    df = flag_undetermined(df)
    df = flag_replicate_variance(df)
    df = flag_ntc_contamination(df)
    df["any_flag"] = (
        df["flag_undetermined"] | df["flag_high_variance"] | df["flag_ntc_contamination"]
    )
    return df


def qc_summary(df: pd.DataFrame) -> pd.DataFrame:
    flagged = df[df["any_flag"]].copy()
    cols = [
        "Well", "Sample_Name", "Target_Name", "Ct",
        "flag_undetermined", "flag_high_variance", "flag_ntc_contamination",
    ]
    return flagged[cols].reset_index(drop=True)