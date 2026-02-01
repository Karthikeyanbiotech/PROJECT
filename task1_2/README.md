# Disease & Target Prioritization Analysis – Lung Adenocarcinoma (LUAD)

## Project Overview

This project focuses on systematic target identification and prioritization for **Lung Adenocarcinoma (LUAD)** using public biological databases and a simple AI/ML-inspired scoring framework.

The objective is to simulate an early-stage drug discovery workflow by integrating:

- Disease biology understanding  
- Multi-source target evidence  
- Feature-based computational ranking  

The workflow is divided into two parts:

1. **Task 1 – Disease & Target Landscape Analysis**
2. **Task 2 – Target Prioritization Using AI/ML Signals**

---

## Task 1: Disease & Target Landscape Analysis

### 1. Disease Rationale

Lung Adenocarcinoma (LUAD) is the most common subtype of non-small cell lung cancer (NSCLC).  
Despite advances in targeted therapies, long-term disease control remains limited due to:

- High molecular heterogeneity  
- Development of acquired drug resistance  
- Activation of bypass signaling pathways  
- Presence of non-actionable tumors  

These challenges highlight the need for **robust target prioritization strategies** that integrate genetic, expression, and pathway-level evidence.

---

### 2. Database Selection and Justification

To ensure biological relevance and reduce bias, multiple public databases were used, each contributing a distinct evidence layer.

| Database | Purpose |
|--------|---------|
| Open Targets | Integrated disease–target association and evidence scoring |
| DisGeNET | Curated disease–gene associations |
| GWAS Catalog | Human genetic susceptibility evidence |
| UniProt | Protein function, annotation, and druggability context |

This multi-database approach reflects industry-standard target discovery workflows.

---

### 3. Shortlisting of Candidate Targets

Using LUAD-specific evidence from Open Targets and supporting databases, biologically relevant genes involved in tumor growth and survival were shortlisted.

Selected targets include:

- **EGFR** – classical oncogenic driver with therapeutic relevance  
- **MET** – bypass and resistance-associated receptor tyrosine kinase  
- **KRAS** – mutation-driven oncogene with emerging druggability  
- **HER3 (ERBB3)** – signaling amplifier involved in resistance mechanisms  

These targets were selected based on:

- Genetic and somatic mutation evidence  
- Pathway involvement (MAPK, PI3K–AKT)  
- Clinical relevance in LUAD biology  

---

### 4. Evidence-Based Evaluation

Each shortlisted target was qualitatively evaluated using:

- Genetic association strength  
- Oncogenic driver status  
- Expression and pathway involvement  
- Known or potential druggability  

This step established a **biological justification** before moving to computational prioritization.

---

## Task 2: Target Prioritization Using AI/ML Signals

### 1. Objective

The goal of Task 2 is to design a **simple scoring framework** that ranks candidate targets using biologically meaningful features derived from public datasets.

Rather than building a predictive ML model, this task focuses on:

- Feature engineering  
- Weighted evidence integration  
- Transparent and interpretable scoring  

This approach mirrors early-stage AI-assisted decision-making used in pharmaceutical R&D.

---

### 2. Input Dataset

**File:** `openTarget_dataset.tsv`  
**Source:** Open Targets Platform  
**Disease:** Lung Adenocarcinoma (LUAD)

Each row represents a candidate gene, and each column represents a specific evidence modality.

---

### 3. Feature Definition

The following features were selected for target scoring:

| Feature | Description |
|------|-------------|
| globalScore | Overall disease–target association score |
| cancerGeneCensus | Curated cancer driver annotation |
| intogen | Statistically significant driver mutations |
| evaSomatic | Somatic mutation evidence |
| cancerBiomarkers | Clinical biomarker relevance |
| chembl | Drug-target interaction and tractability |
| reactome | Pathway-level biological involvement |
| europepmc | Literature support |

These features collectively capture **causality, function, clinical relevance, and feasibility**.

---

### 4. Data Preprocessing

The following steps were applied:

1. Load LUAD target dataset  
2. Retain biologically relevant evidence columns  
3. Replace missing values with zero (absence of evidence)  
4. Filter weak associations using `globalScore ≥ 0.2`  

This ensures consistent scoring across targets.

---

### 5. Scoring Framework Design

A weighted scoring model was implemented using Python.

Weights were assigned based on biological importance:

| Evidence Category | Weight |
|------------------|--------|
| Genetic association | 30% |
| Cancer driver evidence | 30% |
| Clinical / somatic evidence | 20% |
| Druggability & pathway | 15% |
| Literature support | 5% |

This weighting prioritizes **causal biology over publication frequency**.

---

### 6. Priority Score Calculation

A composite priority score was computed as a weighted sum of all features:

priorityScore = Σ (feature × assigned weight)

Targets were then ranked in descending order of priority score.

---

### 7. Visualization

The ranked targets were visualized using:

- Bar plots to show relative priority scores  
- Heatmaps to display evidence contribution patterns  

This enables intuitive interpretation of why certain targets rank higher.

---

## Model Assumptions

- Public database scores are treated as normalized biological proxies  
- Multiple independent evidence sources increase confidence  
- Higher composite scores indicate stronger translational potential  

---

## Model Limitations

- Bias toward well-annotated genes  
- Limited representation of novel or poorly studied targets  
- Does not model tumor microenvironment or patient heterogeneity  

These limitations are inherent to public-data–driven prioritization frameworks.

---

## Final Conclusion

This project demonstrates a structured approach to early-stage target prioritization by integrating:

- Disease biology understanding  
- Multi-source evidence aggregation  
- Interpretable AI/ML-inspired scoring  

The framework enables rational decision-making and reduces early discovery risk by highlighting targets with strong biological validity and therapeutic feasibility.

---

## Tools Used

- Python (pandas, numpy, matplotlib, seaborn)
- Open Targets Platform
- Public genetic and pathway databases

---

## Author

Karthikeyan  
M.Tech Bioinformatics  

