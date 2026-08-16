import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.ddct import compute_mean_ct, compute_dct, compute_ddct
from src.qc import run_qc


def make_toy_df():
    rows = [
        {"Well": "A1", "Sample_Name": "Control_A", "Group": "Control", "Target_Name": "GAPDH", "Task": "UNKNOWN", "Technical_Replicate": 1, "Ct": 20.0},
        {"Well": "A2", "Sample_Name": "Control_A", "Group": "Control", "Target_Name": "GAPDH", "Task": "UNKNOWN", "Technical_Replicate": 2, "Ct": 20.0},
        {"Well": "A3", "Sample_Name": "Control_A", "Group": "Control", "Target_Name": "TargetX", "Task": "UNKNOWN", "Technical_Replicate": 1, "Ct": 25.0},
        {"Well": "A4", "Sample_Name": "Control_A", "Group": "Control", "Target_Name": "TargetX", "Task": "UNKNOWN", "Technical_Replicate": 2, "Ct": 25.0},
        {"Well": "B1", "Sample_Name": "Treated_A", "Group": "Treated", "Target_Name": "GAPDH", "Task": "UNKNOWN", "Technical_Replicate": 1, "Ct": 20.0},
        {"Well": "B2", "Sample_Name": "Treated_A", "Group": "Treated", "Target_Name": "GAPDH", "Task": "UNKNOWN", "Technical_Replicate": 2, "Ct": 20.0},
        {"Well": "B3", "Sample_Name": "Treated_A", "Group": "Treated", "Target_Name": "TargetX", "Task": "UNKNOWN", "Technical_Replicate": 1, "Ct": 27.0},
        {"Well": "B4", "Sample_Name": "Treated_A", "Group": "Treated", "Target_Name": "TargetX", "Task": "UNKNOWN", "Technical_Replicate": 2, "Ct": 27.0},
    ]
    return pd.DataFrame(rows)


def test_mean_ct():
    df = run_qc(make_toy_df())
    mean_ct = compute_mean_ct(df)
    row = mean_ct[(mean_ct.Sample_Name == "Control_A") & (mean_ct.Target_Name == "GAPDH")]
    assert row["Mean_Ct"].iloc[0] == 20.0


def test_dct():
    df = run_qc(make_toy_df())
    mean_ct = compute_mean_ct(df)
    dct = compute_dct(mean_ct, reference_gene="GAPDH")
    row = dct[(dct.Sample_Name == "Control_A") & (dct.Target_Name == "TargetX")]
    assert row["dCt"].iloc[0] == 5.0


def test_ddct_and_fold_change():
    df = run_qc(make_toy_df())
    mean_ct = compute_mean_ct(df)
    dct = compute_dct(mean_ct, reference_gene="GAPDH")
    ddct = compute_ddct(dct, control_group="Control")

    control_row = ddct[ddct.Sample_Name == "Control_A"].iloc[0]
    treated_row = ddct[ddct.Sample_Name == "Treated_A"].iloc[0]

    assert abs(control_row["ddCt"] - 0.0) < 1e-9
    assert abs(control_row["Fold_Change"] - 1.0) < 1e-9
    assert abs(treated_row["ddCt"] - 2.0) < 1e-9
    assert abs(treated_row["Fold_Change"] - 0.25) < 1e-9


def test_qc_flags_undetermined():
    toy = make_toy_df()
    toy.loc[0, "Ct"] = None
    df = run_qc(toy)
    assert df.loc[0, "flag_undetermined"] == True  # noqa: E712