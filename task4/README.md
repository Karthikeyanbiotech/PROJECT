STEP 1: Load Reactome pathway file
pathway = pd.read_csv(
    "Participating Molecules [R-HSA-1643713].tsv", sep="\t"
)

What you did

You loaded a Reactome pathway TSV file into Python.

What this file contains

Each row = one molecule participating in a pathway

Column MoleculeName = protein or complex name

This file defines which genes belong to the EGFR pathway.

Why this step is needed

You must first know which genes are part of the pathway before building any interaction network.

This step defines the nodes of the network.

STEP 2: Extract gene symbols
genes = pathway["MoleculeName"].dropna() \
        .apply(lambda x: x.split(" ")[-1]) \
        .unique().tolist()

What you did

You extracted gene symbols from the MoleculeName column.

Example:

"EGFR [plasma membrane]" → EGFR

Why you used this logic

Reactome stores molecules as descriptive strings.

STRING database accepts gene symbols only, not full descriptions.

So you cleaned the names to obtain:

EGFR, GRB2, SOS1, KRAS, PIK3CA ...


This creates the final pathway gene list.

STEP 3: Prepare STRING API connection
string_api = "https://string-db.org/api"
edges = []

What you did

You prepared to query the STRING database programmatically.

STRING provides protein–protein interaction (PPI) data.

Why STRING is used

Reactome gives:

pathway membership

STRING gives:

physical and functional interactions between proteins

You are combining:

curated pathway genes + interaction evidence

STEP 4: Query STRING for each gene
for gene in genes:

What you did

You looped through each pathway gene one by one.

For each gene, you requested:

its interacting proteins

only in humans (species = 9606)

STEP 5: Set interaction confidence threshold
"required_score": 700

What this means

STRING interaction scores range from 0–1000.

Score ≥ 700 = high-confidence interactions only.

Why this is important

You intentionally removed:

weak

predicted

noisy interactions

Only strong interactions are kept.

STEP 6: Send request to STRING
response = requests.post(url, data=params)

What you did

You sent a POST request to STRING and received interaction data.

Each response contains multiple interaction pairs.

STEP 7: Validate response
if not response.text.strip():

Why you added this

STRING does not always return data for every gene.

This prevents:

crashes

empty parsing errors

This is good defensive coding.

STEP 8: Parse interaction lines
protein1, protein2, score = l[2], l[3], float(l[5])

What you extracted

From STRING output:

interacting protein A

interacting protein B

interaction confidence score

Each row becomes one interaction (edge).

STEP 9: Filter interactions
if protein1 in genes and protein2 in genes:

What this does

You kept only interactions within your pathway genes.

You removed:

outside proteins

unrelated signaling noise

So the final network is:

Reactome pathway–restricted PPI network

This is very important.

STEP 10: Store interactions
edges.append((protein1, protein2, score))


Now you have:

node list → pathway genes

edge list → STRING interactions

STEP 11: Create interaction table
edge_df = pd.DataFrame(edges)


This table contains:

Protein1

Protein2

Interaction score

This is the network edge list.

STEP 12: Build graph using NetworkX
G = nx.Graph()


You created an undirected protein–protein interaction network.

STEP 13: Add edges to network
G.add_edge(row["Protein1"], row["Protein2"], weight=row["Score"])


Each interaction becomes:

a connection between two proteins

weighted by confidence score

Now the network structure exists.

STEP 14: Count nodes and edges
G.number_of_nodes()
G.number_of_edges()


This tells:

how many pathway proteins are connected

how dense the signaling network is

STEP 15: Network visualization
pos = nx.spring_layout(G)
nx.draw(...)

What this does

Creates a force-directed layout

Highly connected proteins move toward center

Weakly connected proteins move outward

This gives a visual signaling map.

STEP 16: Centrality analysis
degree_centrality = nx.degree_centrality(G)

What degree centrality means

For each protein:

how many direct interactions it has with others

Higher value = more connected = more important in signaling flow.

STEP 17: Rank proteins
centrality_df.sort_values(ascending=False)


This identifies:

hub proteins

key signaling mediators

potential secondary targets

FINAL OUTPUT OF TASK 4

From this workflow you obtained:

✅ Pathway gene list (Reactome)
✅ High-confidence interaction network (STRING)
✅ Visual network map
✅ Hub proteins using centralityThis folder contains files related to Task 4.
