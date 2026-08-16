from __future__ import annotations

import numpy as np
import pandas as pd


def compute_mean_ct(qc_df: pd.DataFrame, exclude_flagged: bool = True) -> pd.DataFrame:
    df = qc_df.copy()
    if exclude_flagged:
        df = df[~df["any_flag"]]

    mean_ct = (
        df.groupby(["Sample_Name", "Group", "Target_Name"], as_index=False)["Ct"]
        .mean()
        .rename(columns={"Ct": "Mean_Ct"})
    )
    return mean_ct


def compute_dct(mean_ct_df: pd.DataFrame, reference_gene: str) -> pd.DataFrame:
    pivot = mean_ct_df.pivot_table(
        index=["Sample_Name", "Group"], columns="Target_Name", values="Mean_Ct"
    ).reset_index()

    if reference_gene not in pivot.columns:
        raise ValueError(f"Reference gene '{reference_gene}' not found in data.")

    target_genes = [c for c in pivot.columns if c not in ("Sample_Name", "Group", reference_gene)]

    records = []
    for _, row in pivot.iterrows():
        for gene in target_genes:
            records.append(
                {
                    "Sample_Name": row["Sample_Name"],
                    "Group": row["Group"],
                    "Target_Name": gene,
                    "dCt": row[gene] - row[reference_gene],
                }
            )
    return pd.DataFrame(records)


def compute_ddct(dct_df: pd.DataFrame, control_group: str) -> pd.DataFrame:
    df = dct_df.copy()
    control_means = (
        df[df["Group"] == control_group]
        .groupby("Target_Name")["dCt"]
        .mean()
        .rename("control_mean_dCt")
    )
    df = df.merge(control_means, on="Target_Name", how="left")
    df["ddCt"] = df["dCt"] - df["control_mean_dCt"]
    df["Fold_Change"] = 2 ** (-df["ddCt"])
    return df.drop(columns=["control_mean_dCt"])


def run_ddct_pipeline(qc_df: pd.DataFrame, reference_gene: str, control_group: str) -> pd.DataFrame:
    mean_ct = compute_mean_ct(qc_df)
    dct = compute_dct(mean_ct, reference_gene)
    ddct = compute_ddct(dct, control_group)
    return ddct