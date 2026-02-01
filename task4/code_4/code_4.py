import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import requests
import os

# --- DYNAMIC PATH SETUP ---
base_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(base_dir, "..", "task4_input")
OUTPUT_DIR = os.path.join(base_dir, "..", "output_4")

# Create output_4 directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------------------------------
# STEP 1: Load Reactome pathway molecules (nodes)
# -----------------------------------------------------
# File name matches your folder exactly
input_file = os.path.join(INPUT_DIR, "Participating Molecules [R-HSA-1643713].tsv")
pathway = pd.read_csv(input_file, sep="\t")

# Extract gene symbols from 'MoleculeName'
genes = pathway["MoleculeName"].dropna().apply(lambda x: x.split(' ')[-1]).unique().tolist()

print(f"Number of pathway genes identified: {len(genes)}")

# -----------------------------------------------------
# STEP 2: Retrieve high-confidence interactions (STRING)
# -----------------------------------------------------
string_api = "https://string-db.org/api"
edges = []

print("Fetching interactions from STRING API... (this may take a moment)")
for gene in genes:
    params = {
        "identifiers": gene,
        "species": 9606,
        "required_score": 700
    }
    url = f"{string_api}/tsv-no-header/network"
    try:
        response = requests.post(url, data=params)
        if not response.text.strip() or response.text.strip().startswith("Error"):
            continue

        for line in response.text.strip().split("\n"):
            l = line.split("\t")
            if len(l) > 5:
                protein1, protein2, score = l[2], l[3], float(l[5])
                if protein1 in genes and protein2 in genes:
                    edges.append((protein1, protein2, score))
    except Exception as e:
        print(f"Connection error for {gene}: {e}")

edge_df = pd.DataFrame(edges, columns=["Protein1", "Protein2", "Score"]).drop_duplicates()
print(f"Number of high-confidence interactions: {edge_df.shape[0]}")

# -----------------------------------------------------
# STEP 3: Build protein–protein interaction network
# -----------------------------------------------------
G = nx.Graph()
for _, row in edge_df.iterrows():
    G.add_edge(row["Protein1"], row["Protein2"], weight=row["Score"])

# -----------------------------------------------------
# STEP 4: Network visualization
# -----------------------------------------------------
plt.figure(figsize=(12,10))
pos = nx.spring_layout(G, k=0.15, seed=42) # k adjusts spacing

nx.draw(
    G, pos,
    with_labels=True,
    node_size=1000,
    node_color="skyblue",
    edge_color="silver",
    font_size=8,
    font_weight='bold'
)

plt.title("EGFR Pathway Interaction Network (Reactome + STRING)", fontsize=15)
# Save to output_4
plt.savefig(os.path.join(OUTPUT_DIR, "EGFR_PPI_Network.png"), dpi=300, bbox_inches='tight')
plt.show()

# -----------------------------------------------------
# STEP 5: Centrality analysis
# -----------------------------------------------------
degree_centrality = nx.degree_centrality(G)
centrality_df = (
    pd.DataFrame.from_dict(degree_centrality, orient="index", columns=["DegreeCentrality"])
    .sort_values("DegreeCentrality", ascending=False)
)

# Save the centrality analysis to CSV
centrality_df.to_csv(os.path.join(OUTPUT_DIR, "Network_Centrality_Analysis.csv"))

print("\nTop interacting proteins (Hubs):")
print(centrality_df.head(10))
print(f"\nTask 4 analysis complete. Results saved in: {OUTPUT_DIR}")
