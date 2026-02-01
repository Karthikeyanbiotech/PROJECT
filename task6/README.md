AI-Based Hit Scoring / QSAR Prototype (EGFR)
Overview

This task implements a machine learning–based QSAR (Quantitative Structure–Activity Relationship) prototype to support early-stage hit prioritization for the validated oncology target EGFR.

The objective is to learn relationships between molecular structure–derived descriptors and experimental bioactivity data, enabling ranking of compounds by predicted likelihood of activity. The model is designed to assist decision-making during early drug discovery by reducing experimental screening burden.

Objective

Predict bioactivity trends using molecular descriptors

Prioritize candidate compounds based on hit probability

Evaluate regression vs classification modeling strategies

Select the most reliable model for early-stage screening

This task focuses on hit prioritization, not exact potency prediction.

Input Data
Primary Input File

File name:

EGFR_ChEMBL_Bioactivity.tsv

Data Source

Database: ChEMBL

Target: EGFR (Human)

Data type: experimentally measured bioactivity values

Bioactivity Types Included

IC50

Ki

Kd

Input File Description

Each row in the input file represents a compound–activity measurement for EGFR.

Key columns used in the analysis include:

Column	Description
Molecule ChEMBL ID	Unique compound identifier
Smiles	Canonical SMILES string
Standard Type	Activity measurement type (IC50, Ki, Kd)
Standard Value	Experimental activity value
Standard Units	Measurement unit (nM)
Target ChEMBL ID	EGFR target annotation
Assay Information	Experimental context metadata

Due to GitHub file size limitations, the full dataset is hosted externally.
Download links are provided in the Data Availability section.

Data Preprocessing

The following preprocessing steps were applied:

Retained only quantitative activity measurements (IC50, Ki, Kd)

Removed entries with missing SMILES or activity values

Standardized all activity values to nanomolar (nM)

Converted activity values to pActivity scale (−log10 molar)

This transformation stabilizes variance, improves numerical learning behavior, and allows comparability across heterogeneous assays.

Descriptor Generation

Molecular descriptors were generated using RDKit, converting SMILES strings into numerical representations of chemical structure.

Descriptors included:

Molecular weight (MolWt)

Lipophilicity (LogP)

Topological polar surface area (TPSA)

Hydrogen bond donors (HBD)

Hydrogen bond acceptors (HBA)

Number of rotatable bonds

These descriptors capture physicochemical properties known to influence kinase inhibitor binding.

Modeling Strategy

Two modeling paradigms were evaluated:

Regression

Predict continuous pActivity values

Classification

Predict active vs inactive compounds using an activity threshold

Due to experimental noise, assay heterogeneity, and strong class imbalance in public datasets, classification was selected as the primary hit-scoring approach, as it better reflects early-stage drug discovery decisions.

Machine Learning Models Evaluated

The following models were benchmarked:

Logistic Regression

Random Forest Classifier

Gradient Boosting Classifier

All models were trained using stratified train–test splits to preserve class balance.

Model Performance
Regression Benchmarking
Model	RMSE	R²
Linear Regression	1.30	0.16
Random Forest Regressor	1.08	0.42
Gradient Boosting Regressor	1.21	0.28

Regression performance was limited by experimental variability and assay noise.

Classification Performance
Model	ROC-AUC
Logistic Regression	0.76
Random Forest Classifier	0.87
Gradient Boosting Classifier	0.81

The Random Forest classifier demonstrated the strongest discriminatory power.

Final Model Selection

Selected Model: Random Forest Classifier

Rationale:

Highest ROC-AUC performance

Robust to noisy experimental data

Handles nonlinear chemical relationships

Suitable for imbalanced datasets

Provides interpretable feature importance

Hit Scoring and Prioritization

The final model was retrained on the complete dataset and used to compute hit probability scores for all compounds.

Each compound was assigned:

Experimental pActivity

Binary activity label

Predicted probability of being active (Hit Probability)

Compounds were ranked based on predicted hit probability to identify high-confidence candidates for experimental validation.

Model Reliability and Data Bias
Applicability Domain

The model is reliable primarily for compounds chemically similar to those present in the training dataset. Predictions outside this chemical space may be unreliable.

Activity Cliffs

Small chemical changes can lead to large bioactivity differences, limiting the ability of machine learning models to perfectly capture structure–activity relationships.

Decoy Bias

If inactive compounds differ strongly in physicochemical properties from actives, models may learn superficial patterns rather than true binding interactions.

Therefore, predictions should be interpreted as prioritization guidance, not definitive potency estimates.

Output

Ranked compound list with hit probabilities

Performance metrics (ROC-AUC, RMSE, R²)

Feature importance analysis

Prioritized candidates for downstream validation
