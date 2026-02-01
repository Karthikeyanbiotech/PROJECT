The objective of Task 3 is to validate the prioritized target gene (EGFR) using two independent biological evidence layers:

Transcriptomics (TCGA-LUAD) — to examine tumor vs normal expression

Functional genomics (DepMap CRISPR) — to evaluate gene dependency in LUAD cell lines

This task helps determine whether EGFR is:

expression-driven, or

mutation-driven with functional dependency in specific tumor subsets.

📂 Input Files
1. TCGA Expression Dataset

File name:
TCGA.LUAD.sampleMap_HiSeqV2

Source:
UCSC Xena (TCGA project)

Description:
RNA-seq expression matrix for LUAD samples.

Structure:

Rows → genes

Columns → patient samples

Values → log2-normalized expression

Sample identification:

01 → Tumor sample

11 → Normal tissue

2. DepMap Functional Genomics Dataset

File name:
CRISPRGeneEffect.csv

Source:
DepMap (Cancer Dependency Map)

Description:
Genome-wide CRISPR knockout gene-effect scores.

Interpretation:

More negative score → higher gene dependency

3. DepMap Cell Line Metadata

File name:
Model.csv

Source:
DepMap Portal

Description:
Contains cancer-type annotations for each cell line, including Oncotree codes.

🧪 Workflow Overview

Load TCGA LUAD transcriptomic data

Extract EGFR expression

Classify tumor and normal samples

Perform statistical comparison

Visualize expression differences

Load DepMap CRISPR dependency data

Filter LUAD-specific cell lines

Quantify EGFR dependency

Interpret integrated results

🧩 Step-by-Step Methodology
Step 1: Import Required Libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


These libraries are used for data handling, statistics, and visualization.

Step 2: Load TCGA LUAD Expression Data
expr = pd.read_csv("TCGA.LUAD.sampleMap_HiSeqV2", sep="\t")
expr = expr.rename(columns={"sample": "Gene"})


Loads the RNA-seq expression matrix and prepares the gene column for filtering.

Step 3: Extract EGFR Expression
egfr_row = expr[expr["Gene"] == "EGFR"]
egfr_expr = egfr_row.set_index("Gene").T
egfr_expr.columns = ["EGFR_expression"]


Selects EGFR and reshapes the data so each row corresponds to one patient sample.

Step 4: Identify Tumor and Normal Samples
egfr_expr["Sample_Type"] = egfr_expr.index.str[-2:]
egfr_expr["Group"] = egfr_expr["Sample_Type"].map({
    "01": "Tumor",
    "11": "Normal"
})
egfr_expr = egfr_expr.dropna(subset=["Group"])


Uses TCGA barcode information to classify samples.

Step 5: Save Processed Expression File
egfr_expr_reset = egfr_expr.reset_index()
egfr_expr_reset = egfr_expr_reset.rename(columns={"index": "Sample_ID"})
egfr_final = egfr_expr_reset[["Sample_ID", "EGFR_expression", "Group"]]
egfr_final.to_csv("EGFR_LUAD_Tumor_vs_Normal.csv", index=False)


Creates a clean output file containing EGFR expression with sample labels.

Step 6: Statistical Comparison
tumor = egfr_expr[egfr_expr["Group"] == "Tumor"]["EGFR_expression"]
normal = egfr_expr[egfr_expr["Group"] == "Normal"]["EGFR_expression"]

stat, p_value = mannwhitneyu(tumor, normal, alternative="two-sided")


Performs a Mann–Whitney U test to compare tumor and normal expression values.

Step 7: Visualization
plt.boxplot([tumor, normal], labels=["Tumor", "Normal"])
plt.ylabel("EGFR expression (log2)")
plt.title("EGFR expression in TCGA-LUAD")
plt.show()


Visualizes expression distribution between tumor and normal samples.

🔬 Functional Genomics Analysis (DepMap)
Step 8: Load CRISPR Dependency Data
gene_effect = pd.read_csv("CRISPRGeneEffect.csv", index_col=0, low_memory=False)


Loads genome-wide CRISPR knockout dependency scores.

Step 9: Extract EGFR Dependency Scores
egfr_col = [c for c in gene_effect.columns if c.startswith("EGFR")][0]
egfr_scores = gene_effect[egfr_col].reset_index()
egfr_scores.columns = ["DepMap_ID", "EGFR_gene_effect"]


Selects EGFR-specific knockout scores.

Step 10: Load and Merge Cell Line Metadata
model_info = pd.read_csv("Model.csv")
model_info = model_info.rename(columns={"ModelID": "DepMap_ID"})

egfr_df = egfr_scores.merge(model_info, on="DepMap_ID", how="inner")


Adds cancer-type information to dependency data.

Step 11: Filter LUAD Cell Lines
luad_df = egfr_df[egfr_df["OncotreeCode"] == "LUAD"]


Restricts analysis to lung adenocarcinoma models only.

Step 12: Dependency Classification
strong_dep = (luad_df["EGFR_gene_effect"] <= -0.5).sum()
moderate_dep = ((luad_df["EGFR_gene_effect"] > -0.5) &
                (luad_df["EGFR_gene_effect"] <= -0.3)).sum()
weak_dep = (luad_df["EGFR_gene_effect"] > -0.3).sum()


Classifies dependency strength using standard DepMap thresholds.

Step 13: Visualization
plt.hist(luad_df["EGFR_gene_effect"], bins=25)
plt.axvline(-0.5, linestyle="--")
plt.axvline(-0.3, linestyle="--")
plt.xlabel("EGFR gene-effect score")
plt.ylabel("Number of LUAD cell lines")
plt.title("EGFR CRISPR dependency in LUAD")
plt.show()


Displays distribution of EGFR dependency across LUAD cell lines.

✅ Final Outcome

EGFR does not show strong transcriptional upregulation in LUAD tumors.

A subset of LUAD cell lines shows strong dependency on EGFR.

This pattern is consistent with mutation-driven oncogene addiction.

Final decision:
GO — precision oncology context
