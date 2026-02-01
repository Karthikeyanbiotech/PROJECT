# lung adenocarcinoma data

# OPEN TARGET

import pandas as pd
import numpy as np
df = pd.read_csv("task1_2/task1_input/openTarget_dataset.tsv", na_values=["No data"], sep='\t')
df = df[
    [
        "symbol",
        "globalScore",
        "cancerGeneCensus",
        "intogen",
        "evaSomatic",
        "cancerBiomarkers",
        "chembl",
        "reactome",
        "europepmc"
    ]
]
df = df.fillna(0);

# filter strong disease-associated targets
df = df[df["globalScore"] >= 0.2]

# sort by association strength
df = df.sort_values(by="globalScore", ascending=False)

# save step  output
df.to_csv("task1_2/output/SelectedTargets.csv", index=False)
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

# calculate priority score
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

# rank targets
df = df.sort_values("priorityScore", ascending=False)

# save ranked targets
df.to_csv("task1_2/output/RankedTargets.csv", index=False)
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

df = pd.read_csv("task1_2/output/RankedTargets.csv")

top10 = df.head(10)

plt.figure(figsize=(8,5))
plt.bar(top10["symbol"], top10["priorityScore"])
plt.xlabel("Target genes")
plt.ylabel("Priority score")
plt.title("Top 10 prioritized cancer targets")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("task1_2/output/top_10_cancer_targets.png", dpi=300)
plt.show()


features = [
    "globalScore",
    "cancerGeneCensus",
    "intogen",
    "evaSomatic",
    "cancerBiomarkers",
    "chembl",
    "reactome",
    "europepmc"
]

top5 = df.head(5)
heatmap_data = top5.set_index("symbol")[features]

plt.figure(figsize=(8,4))
sns.heatmap(
    heatmap_data,
    annot=True,
    cmap="viridis",
    linewidths=0.5
)
plt.title("Feature contribution heatmap (Top 5 targets)")
plt.tight_layout()
plt.savefig("task1_2/output/feature_contribution_heatmap.png", dpi=300) 
plt.show()