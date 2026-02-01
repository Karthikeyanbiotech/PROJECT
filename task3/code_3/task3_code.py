import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import os

# --- DYNAMIC PATH SETUP ---
# This looks for the 'input' folder relative to where the script is
# It handles the "../" automatically so you don't get FileNotfound errors
base_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(base_dir, "..", "input")
OUTPUT_DIR = os.path.join(base_dir, "..", "output_3")

# Create output_3 directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# STEP 1: TCGA TRANSCRIPTOMICS
# =========================

# Loading from task3/input
expr_path = os.path.join(INPUT_DIR, "TCGA.LUAD.sampleMap_HiSeqV2")
expr = pd.read_csv(expr_path, sep="\t")
expr = expr.rename(columns={"sample": "Gene"})

egfr_row = expr[expr["Gene"] == "EGFR"]
egfr_expr = egfr_row.set_index("Gene").T
egfr_expr.columns = ["EGFR_expression"]

egfr_expr["Sample_Type"] = egfr_expr.index.str[-2:]
egfr_expr["Group"] = egfr_expr["Sample_Type"].map({"01": "Tumor", "11": "Normal"})
egfr_expr = egfr_expr.dropna(subset=["Group"])

egfr_final = egfr_expr.reset_index().rename(columns={"index": "Sample_ID"})
egfr_final = egfr_final[["Sample_ID", "EGFR_expression", "Group"]]

# Saving to task3/output_3
egfr_final.to_csv(os.path.join(OUTPUT_DIR, "EGFR_LUAD_Tumor_vs_Normal.csv"), index=False)

tumor = egfr_expr[egfr_expr["Group"] == "Tumor"]["EGFR_expression"]
normal = egfr_expr[egfr_expr["Group"] == "Normal"]["EGFR_expression"]

stat, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")

print("TCGA-LUAD Transcriptomics")
print(f"P-value: {p_value}")

plt.figure(figsize=(6,6))
plt.boxplot([tumor, normal], labels=["Tumor", "Normal"])
plt.ylabel("EGFR expression (log2)")
plt.title("EGFR expression in TCGA-LUAD")
plt.savefig(os.path.join(OUTPUT_DIR, "EGFR_expression_TCGA.png"), dpi=300)
plt.show()

# =========================
# STEP 2: DEPMAP FUNCTIONAL GENOMICS
# =========================

# Loading from task3/input
gene_effect = pd.read_csv(os.path.join(INPUT_DIR, "CRISPRGeneEffect.csv"), index_col=0, low_memory=False)
model_info = pd.read_csv(os.path.join(INPUT_DIR, "Model.csv"))

egfr_col = [c for c in gene_effect.columns if c.startswith("EGFR")][0]
egfr_scores = gene_effect[egfr_col].reset_index()
egfr_scores.columns = ["DepMap_ID", "EGFR_gene_effect"]

model_info = model_info.rename(columns={"ModelID": "DepMap_ID"})
egfr_df = egfr_scores.merge(model_info, on="DepMap_ID", how="inner")
luad_df = egfr_df[egfr_df["OncotreeCode"] == "LUAD"]

plt.figure(figsize=(8,5))
plt.hist(luad_df["EGFR_gene_effect"], bins=25, color='seagreen')
plt.axvline(-0.5, color='red', linestyle="--")
plt.xlabel("EGFR gene-effect score")
plt.title("EGFR CRISPR dependency in LUAD (DepMap)")
plt.savefig(os.path.join(OUTPUT_DIR, "EGFR_CRISPR_dependency.png"), dpi=300)
plt.show()

print(f"\nTask 3 Complete. Results saved in: {OUTPUT_DIR}")