import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Specialized Chemoinformatics & ML libraries
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

# --- DYNAMIC PATH SETUP ---
base_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(base_dir, "..", "task6_input")
OUTPUT_DIR = os.path.join(base_dir, "..", "output_6")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# =====================================================
# 1. LOAD RAW ChEMBL DATA
# =====================================================
input_file = os.path.join(INPUT_DIR, "chembl_dataset.tsv")

df_raw = pd.read_csv(
    input_file,
    sep="\t",
    on_bad_lines="skip",
    engine="python",
    encoding="latin1" 
)

print(f"Raw dataset loaded. Shape: {df_raw.shape}")

# =====================================================
# 2. FILTER & STANDARDIZE ACTIVITY (Including Names)
# =====================================================
# Keeping Molecule Name and ID columns
df_activity = df_raw[
    (df_raw["Standard Type"].isin(["IC50", "Ki", "Kd"])) &
    (df_raw["Smiles"].notna()) &
    (df_raw["Standard Value"].notna())
].copy()

# Select SMILES, Activity, and Identity columns
cols = ["Smiles", "Standard Value", "Standard Units", "Molecule Name", "Molecule ChEMBL ID"]
qsar_df = df_activity[cols].copy()
qsar_df.rename(columns={"Smiles": "SMILES"}, inplace=True)
qsar_df["Activity_nM"] = pd.to_numeric(qsar_df["Standard Value"], errors="coerce")
qsar_df = qsar_df.dropna(subset=["Activity_nM"])

# Convert nM to pActivity
qsar_df["pActivity"] = -np.log10(qsar_df["Activity_nM"] * 1e-9)

# Group by SMILES but keep the Names/IDs
# We use 'first' to keep the name associated with that structure
qsar_df = qsar_df.groupby("SMILES", as_index=False).agg({
    "pActivity": "mean",
    "Molecule Name": "first",
    "Molecule ChEMBL ID": "first"
})

# Fill empty names with ChEMBL IDs so the column isn't blank
qsar_df["Molecule Name"] = qsar_df["Molecule Name"].fillna(qsar_df["Molecule ChEMBL ID"])

# =====================================================
# 3. RDKit DESCRIPTOR GENERATION
# =====================================================
def calc_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    return [
        Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)
    ]

descriptor_names = ["MolWt", "LogP", "HBD", "HBA", "TPSA", "RotB"]
desc_data, valid_idx = [], []

print("Generating molecular descriptors...")
for i, smi in enumerate(qsar_df["SMILES"]):
    d = calc_descriptors(smi)
    if d is not None:
        desc_data.append(d)
        valid_idx.append(i)

desc_df = pd.DataFrame(desc_data, columns=descriptor_names)
qsar_df = qsar_df.iloc[valid_idx].reset_index(drop=True)
qsar_df["Class"] = (qsar_df["pActivity"] >= 6).astype(int)

# =====================================================
# 4. ML MODELING
# =====================================================
X_train, X_test, y_train, y_test = train_test_split(
    desc_df, qsar_df["Class"], test_size=0.2, random_state=42, stratify=qsar_df["Class"]
)

models = {
    "RandomForest": RandomForestClassifier(n_estimators=300, random_state=42),
    "GradientBoosting": GradientBoostingClassifier(random_state=42)
}

# Training the best model (Random Forest is usually best for QSAR)
best_model = models["RandomForest"]
best_model.fit(desc_df, qsar_df["Class"])

# =====================================================
# 5. OUTPUTS & VISUALIZATION (With Names)
# =====================================================
qsar_df["Hit_Probability"] = best_model.predict_proba(desc_df)[:, 1]

# Sort by probability and keep the names in the final display
top_hits = qsar_df.sort_values("Hit_Probability", ascending=False).head(20)

# Reordering columns for better CSV readability
final_cols = ["Molecule Name", "Molecule ChEMBL ID", "Hit_Probability", "pActivity", "SMILES"]
top_hits = top_hits[final_cols]

top_hits.to_csv(os.path.join(OUTPUT_DIR, "Top_AI_Prioritized_Hits.csv"), index=False)

print("\nTOP 5 PRIORITIZED HITS:")
print(top_hits[["Molecule Name", "Hit_Probability"]].head())

# Save Feature Importance Plot
if hasattr(best_model, "feature_importances_"):
    importance = pd.Series(best_model.feature_importances_, index=descriptor_names).sort_values(ascending=False)
    plt.figure(figsize=(8,5))
    importance.plot(kind="bar", color='teal')
    plt.title("Feature Importance (QSAR Model)")
    plt.savefig(os.path.join(OUTPUT_DIR, "QSAR_Feature_Importance.png"), dpi=300)

print(f"\nAnalysis Complete. Results saved to: {OUTPUT_DIR}")