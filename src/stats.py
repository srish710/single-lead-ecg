from __future__ import annotations

import pandas as pd
from scipy import stats


def run_ttest(dct_df: pd.DataFrame, group_a: str, group_b: str) -> pd.DataFrame:
    results = []
    for gene, sub in dct_df.groupby("Target_Name"):
        a = sub[sub["Group"] == group_a]["dCt"].dropna()
        b = sub[sub["Group"] == group_b]["dCt"].dropna()
        if len(a) < 2 or len(b) < 2:
            results.append({"Target_Name": gene, "p_value": None, "significant": None})
            continue
        t_stat, p_value = stats.ttest_ind(a, b, equal_var=False)
        results.append(
            {
                "Target_Name": gene,
                "t_statistic": round(t_stat, 3),
                "p_value": round(p_value, 4),
                "significant": bool(p_value < 0.05),
            }
        )
    return pd.DataFrame(results)


def summarize(ddct_df: pd.DataFrame, ttest_df: pd.DataFrame) -> pd.DataFrame:
    fc_summary = (
        ddct_df.groupby(["Group", "Target_Name"], as_index=False)["Fold_Change"]
        .mean()
    )
    return fc_summary.merge(ttest_df, on="Target_Name", how="left")