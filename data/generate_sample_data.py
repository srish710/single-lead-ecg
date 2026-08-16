"""
Generates a SYNTHETIC qPCR Ct-value dataset shaped like a real instrument
export (e.g. Applied Biosystems QuantStudio / Bio-Rad CFX Maestro).
"""

import numpy as np
import pandas as pd

rng = np.random.default_rng(seed=42)

REFERENCE_GENE = "GAPDH"
TARGET_GENE = "BRCA1"
GROUPS = {
    "Control": ["Control_1", "Control_2", "Control_3"],
    "Olaparib": ["Olaparib_1", "Olaparib_2", "Olaparib_3"],
}
N_TECH_REPLICATES = 3

CT_PARAMS = {
    (REFERENCE_GENE, "Control"): (19.0, 0.25),
    (REFERENCE_GENE, "Olaparib"): (19.1, 0.25),
    (TARGET_GENE, "Control"): (25.0, 0.25),
    (TARGET_GENE, "Olaparib"): (27.5, 0.25),
}

rows = []
well_letters = "ABCDEFGH"
well_counter = 0


def next_well():
    global well_counter
    row = well_letters[well_counter // 12]
    col = (well_counter % 12) + 1
    well_counter += 1
    return f"{row}{col}"


for group, samples in GROUPS.items():
    for sample in samples:
        for gene in (REFERENCE_GENE, TARGET_GENE):
            mean_ct, sd = CT_PARAMS[(gene, group)]
            for rep in range(1, N_TECH_REPLICATES + 1):
                ct = rng.normal(mean_ct, sd)
                rows.append(
                    {
                        "Well": next_well(),
                        "Sample_Name": sample,
                        "Group": group,
                        "Target_Name": gene,
                        "Task": "UNKNOWN",
                        "Technical_Replicate": rep,
                        "Ct": round(ct, 2),
                    }
                )

for r in rows:
    if r["Sample_Name"] == "Olaparib_2" and r["Target_Name"] == TARGET_GENE and r["Technical_Replicate"] == 1:
        r["Ct"] = np.nan

for r in rows:
    if r["Sample_Name"] == "Control_3" and r["Target_Name"] == TARGET_GENE and r["Technical_Replicate"] == 2:
        r["Ct"] = round(r["Ct"] + 1.8, 2)

ntc_rows = [
    {"Well": next_well(), "Sample_Name": "NTC", "Group": "NTC", "Target_Name": REFERENCE_GENE,
     "Task": "NTC", "Technical_Replicate": 1, "Ct": np.nan},
    {"Well": next_well(), "Sample_Name": "NTC", "Group": "NTC", "Target_Name": TARGET_GENE,
     "Task": "NTC", "Technical_Replicate": 1, "Ct": 29.8},
]
rows.extend(ntc_rows)

df = pd.DataFrame(rows)
out_path = "data/sample_qpcr_run.csv"
df.to_csv(out_path, index=False)
print(f"Wrote {len(df)} rows to {out_path}")