import os
import sys

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.parsing import load_ct_data
from src.qc import run_qc, qc_summary
from src.ddct import compute_mean_ct, compute_dct, compute_ddct
from src.stats import run_ttest, summarize

st.set_page_config(page_title="qPCR Analysis Dashboard", layout="wide")

st.title("qPCR / RT-qPCR Analysis Dashboard")
st.caption(
    "Upload a raw Ct export, flag QC issues, and compute relative gene "
    "expression via the ddCt method."
)

with st.sidebar:
    st.header("1. Data")
    uploaded = st.file_uploader("Upload Ct export (CSV or Excel)", type=["csv", "xlsx", "xls"])
    use_sample = st.checkbox("Use bundled sample dataset instead", value=uploaded is None)

if use_sample or uploaded is None:
    sample_path = os.path.join(os.path.dirname(__file__), "..", "data", "sample_qpcr_run.csv")
    if not os.path.exists(sample_path):
        st.error(
            "Sample dataset not found. Run `python data/generate_sample_data.py` first, "
            "or upload your own file."
        )
        st.stop()
    raw_df = load_ct_data(sample_path)
    st.sidebar.info("Using bundled synthetic sample dataset.")
else:
    raw_df = load_ct_data(uploaded)

qc_df = run_qc(raw_df)
flagged = qc_summary(qc_df)

st.subheader("Raw data")
st.dataframe(raw_df, use_container_width=True, height=200)

st.subheader("QC flags")
if flagged.empty:
    st.success("No QC issues detected.")
else:
    st.warning(f"{len(flagged)} well(s) flagged. These are excluded from the mean-Ct calculation below.")
    st.dataframe(flagged, use_container_width=True)

st.subheader("Plate layout (colored by Ct, flagged wells outlined in red)")


def render_plate(df: pd.DataFrame):
    plate_rows = list("ABCDEFGH")
    plate_cols = list(range(1, 13))
    z = [[None for _ in plate_cols] for _ in plate_rows]
    text = [["" for _ in plate_cols] for _ in plate_rows]

    for _, r in df.iterrows():
        well = r["Well"]
        if len(well) < 2:
            continue
        row_i = plate_rows.index(well[0]) if well[0] in plate_rows else None
        try:
            col_i = int(well[1:]) - 1
        except ValueError:
            col_i = None
        if row_i is None or col_i is None or col_i >= 12:
            continue
        z[row_i][col_i] = r["Ct"]
        text[row_i][col_i] = f"{r['Sample_Name']}<br>{r['Target_Name']}<br>Ct={r['Ct']}"

    fig = go.Figure(
        data=go.Heatmap(
            z=z, x=plate_cols, y=plate_rows, colorscale="Viridis",
            hoverongaps=False, text=text, hoverinfo="text",
        )
    )
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(height=350, margin=dict(l=20, r=20, t=20, b=20))
    return fig


st.plotly_chart(render_plate(qc_df), use_container_width=True)

st.subheader("Relative expression (ddCt)")

genes = sorted(qc_df["Target_Name"].unique())
groups = sorted(qc_df["Group"].unique())
groups = [g for g in groups if g.upper() != "NTC"]

col1, col2 = st.columns(2)
with col1:
    reference_gene = st.selectbox("Reference (housekeeping) gene", genes, index=genes.index("GAPDH") if "GAPDH" in genes else 0)
with col2:
    control_group = st.selectbox("Control group", groups, index=0)

try:
    mean_ct = compute_mean_ct(qc_df)
    dct = compute_dct(mean_ct, reference_gene=reference_gene)
    ddct = compute_ddct(dct, control_group=control_group)

    other_groups = [g for g in groups if g != control_group]
    ttest_results = pd.DataFrame()
    if other_groups:
        ttest_results = run_ttest(dct, group_a=control_group, group_b=other_groups[0])

    summary = summarize(ddct, ttest_results) if not ttest_results.empty else ddct

    left, right = st.columns([2, 1])
    with left:
        fig = px.bar(
            ddct, x="Sample_Name", y="Fold_Change", color="Group",
            facet_col="Target_Name", title="Fold change by sample",
        )
        st.plotly_chart(fig, use_container_width=True)
    with right:
        st.markdown("**Summary (mean fold change + significance)**")
        st.dataframe(summary, use_container_width=True)

    st.download_button(
        "Download full results (CSV)",
        data=ddct.to_csv(index=False),
        file_name="ddct_results.csv",
        mime="text/csv",
    )

except Exception as e:
    st.error(f"Could not compute ddCt with the current settings: {e}")

st.divider()
st.caption(
    "Bundled sample data is synthetic and generated for demonstration purposes only. "
    "Methodology: Livak & Schmittgen (2001) ddCt relative quantification."
)