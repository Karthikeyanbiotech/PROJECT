import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- 1. SETUP PATHS ---
# This ensures the code knows where to look regardless of who runs it
INPUT_FILE = "../task1_input/openTarget_dataset.tsv"
OUTPUT_DIR = "../output"

# Create output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# --- 2. DATA LOADING & CLEANING ---
df = pd.read_csv(INPUT_FILE, na_values=["No data"], sep='\t')

# Selecting specific columns for analysis
cols_to_keep = [
    "symbol", "globalScore", "cancerGeneCensus", "intogen", 
    "evaSomatic", "cancerBiomarkers", "chembl", "reactome", "europepmc"
]
df = df[cols_to_keep]
df = df.fillna(0)

# Filter strong disease-associated targets and sort
df = df[df["globalScore"] >= 0.2]
df = df.sort_values(by="globalScore", ascending=False)

# Save intermediate filtered data
df.to_csv(f"{OUTPUT_DIR}/SelectedTargets.csv", index=False)

# --- 3. PRIORITY SCORING ---
weights = {
    "globalScore": 0.30,
    "cancerGeneCensus": 0.15,
    "intogen": 0.15,
    "evaSomatic": 0.10,
    "cancerBiomarkers": 0.10,
    "chembl": 0.10,
    "reactome": 0.05,
    "europepmc": 0.05
}

df["priorityScore"] = (
    df["globalScore"] * weights["globalScore"] +
    df["cancerGeneCensus"] * weights["cancerGeneCensus"] +
    df["intogen"] * weights["intogen"] +
    df["evaSomatic"] * weights["evaSomatic"] +
    df["cancerBiomarkers"] * weights["cancerBiomarkers"] +
    df["chembl"] * weights["chembl"] +
    df["reactome"] * weights["reactome"] +
    df["europepmc"] * weights["europepmc"]
)

# Rank and save final CSV
df = df.sort_values("priorityScore", ascending=False)
df.to_csv(f"{OUTPUT_DIR}/RankedTargets.csv", index=False)

# --- 4. VISUALIZATION: TOP 10 BAR CHART ---
top10 = df.head(10)
plt.figure(figsize=(10, 6))
plt.bar(top10["symbol"], top10["priorityScore"], color='skyblue', edgecolor='navy')
plt.xlabel("Target Genes", fontweight='bold')
plt.ylabel("Priority Score", fontweight='bold')
plt.title("Top 10 Prioritized Cancer Targets", fontsize=14)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/top_10_cancer_targets.png", dpi=300)
plt.show()

# --- 5. VISUALIZATION: FEATURE HEATMAP ---
features = [
    "globalScore", "cancerGeneCensus", "intogen", "evaSomatic", 
    "cancerBiomarkers", "chembl", "reactome", "europepmc"
]
top5 = df.head(5)
heatmap_data = top5.set_index("symbol")[features]

plt.figure(figsize=(10, 5))
sns.heatmap(heatmap_data, annot=True, cmap="viridis", linewidths=0.5)
plt.title("Feature Contribution Heatmap (Top 5 Targets)", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/feature_contribution_heatmap.png", dpi=300)
plt.show()

print(f"Analysis complete. All files saved to: {OUTPUT_DIR}")