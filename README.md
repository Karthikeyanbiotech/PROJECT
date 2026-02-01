# PROJECT
Project-Task wise organization
# Multi-Omic Drug Target Discovery & Hit Prioritization
## Project-Task Wise Organization

This project demonstrates an end-to-end computational drug discovery pipeline focused on **EGFR** as a therapeutic target in **Lung Adenocarcinoma (LUAD)**.

---

## 📂 Project Structure
* **task1_2/**: Target identification via Open Targets and Python-based prioritization scoring.
* **task3/**: Multi-omics validation using TCGA (transcriptomics) and DepMap (CRISPR dependency).
* **task4/**: PPI Network analysis using Reactome and STRING DB.
* **task6/**: QSAR Machine Learning model for bioactivity prediction using ChEMBL data.

---

## 🔬 Task Summaries

### Task 1 & 2: Disease Landscape & AI Prioritization
* **Objective**: Identify and rank high-potential targets for Lung Adenocarcinoma.
* **Methodology**: Built a Python scoring framework using features like genetic association and druggability.
* **Key Result**: EGFR was identified as the top-ranked target for further validation.

### Task 3: Multi-Omics Target Validation
* **Data Sources**: TCGA-LUAD (Transcriptomics) and DepMap Public 24Q2 (Functional Genomics).
* **Findings**: 
    * **Expression**: EGFR expression is high in tumors, though the difference vs. normal was not statistically significant ($p > 0.05$).
    * **Dependency**: CRISPR knockout data confirmed a "Go" signal, with ~17% of cell lines showing strong dependency.

### Task 4: Pathway & Network Analysis
* **Tools**: NetworkX, STRING API, Reactome.
* **Insight**: Constructed a Protein-Protein Interaction (PPI) network identifying key signaling partners (e.g., GRB2, SHC1) involved in oncogenic signaling.

### Task 6: AI-Based Hit Scoring (QSAR)
* **Objective**: Predict bioactivity of small molecules against EGFR.
* **Model**: Random Forest Classifier trained on ChEMBL bioactivity data ($IC_{50}$ / $K_i$).
* **Features**: RDKit descriptors (MolWt, LogP, HBD, HBA, TPSA).
* **Performance**: Successfully prioritized compounds with high predicted hit probability.

---

## 🛠️ Installation & Usage
### Prerequisites
* Python 3.8+
* RDKit (`pip install rdkit`)

### Installation
```bash
pip install pandas numpy matplotlib seaborn scipy networkx scikit-learn requests
Running the Pipeline
Navigate to any task folder and run the python script:

Bash
cd task6/code_6
python3 code_6.py
Task,Output,Outcome
1 & 2,Priority Table,EGFR ranked as Top Target
3,Boxplot & Histogram,Confirmed functional dependency
4,PPI Network Graph,Identified key pathway hubs
6,Hit Probabilities,Prioritized top 20 lead compounds
