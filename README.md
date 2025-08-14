# CYCLOPEpDB

CYCLOPEpDB is a database containing data from over **10,000 Molecular Dynamics (MD) simulations** of cyclic peptide nanotubes interacting with a bacterial membrane model.  Currently, the database includes data for **3,855 different cyclic peptide (CP) sequences** and more than **1,000 MD descriptors** per simulation.

The generated data has been:
- **Clustered** into three groups based on the ability of CPs to interact with the membrane.  
- Used to compute a **pseudo-permeability score**, derived from principal components obtained through PCA of selected MD descriptors.  
- Leveraged to train a **Message Passing Neural Network (MPNN)** for prediction tasks.

---

## Getting Started

To use the code provided in this repository, you need to create a Python environment and install the required libraries. A `.yml` file is provided for this:

```bash
conda env create -f cyclopepdb_env.yml
```

---

## Run Predictions

This repository contains code for predicting MD results of **D,L-α-cyclic octapeptides**. The main script, **`MD_Predictor.py`**, takes a sequence of amino acids as input and outputs both a classification and a pseudo-permeability score:

```bash
python MD_Predictor.py -s RRKWLWLW

CP with sequence dRRdKWdLWdLW belongs to the Lipophilic cluster 
with a predicted pseudo-permeability of 9.04
--- TOTAL TIME: 17.26 seconds ---
```

The output classifies the CP sequence into **Lipophilic**, **Intermediate**, or **Lipophobic**, and provides a pseudo-permeability score ranging approximately from **0 to 10**:  
- **0** → weaker CP–membrane interaction  
- **10** → strong CP–membrane interaction

---

## Other Applications

This repository also includes **`CP_PDB_Generator.py`** and **`GraphGenerator.py`**, which are used internally by `MD_Predictor.py` but can also be run independently.

- **`CP_PDB_Generator.py`** — Adapted from code in [CYCLOPEpBuilder](https://cyclopep.com/builder). Generates `.pdb` and `.itp` files for the MA(R/S)TINI3 **forcefield** from a CP sequence:  

  ```bash
  python CP_PDB_Generator.py -s RRKWLWLW
  ```

- **`GraphGenerator.py`** — Creates a *coarse-grain*-based graph using the Deep Graph Library (DGL).  
  - **Node features**: mass, charge, 3D coordinates, Lennard–Jones parameters (per bead)  
  - **Edge features**: bead–bead distances between connected nodes  

  Example usage in Python:
  ```python
  from GraphGenerator import CP_Graph
  CP = CP_Graph("RRKWLWLW", "marstini3")

---
## Tools Overview

| Script | Inputs (CLI / Python) | Outputs | Purpose | Example |
|---|---|---|---|---|
| `MD_Predictor.py` | `-s <SEQ>` (one-letter CP sequence) | Cluster label (**Lipophilic / Intermediate / Lipophobic**), pseudo-permeability score (≈0–10) | Predict CP–membrane interaction class and pseudo-permeability using an MPNN trained on MD-derived descriptors | `python MD_Predictor.py -s RRKWLWLW` |
| `CP_PDB_Generator.py` | `-s <SEQ>` (one-letter CP sequence) | `.pdb` and `.itp` files for the MA(R/S)TINI3 **forcefield** | Generate CP structure/topology from a sequence (adapted from code in CYCLOPEpBuilder) | `python CP_PDB_Generator.py -s RRKWLWLW` |
| `GraphGenerator.py` | Python API: `CP_Graph(sequence, "marstini3")` | DGL graph | Build a *coarse-grain*-based graph for model training: node features = mass, charge, 3D coordinates, Lennard–Jones params; edge features = bead–bead distances | `from GraphGenerator import CP_Graph; CP = CP_Graph("RRKWLWLW", "marstini3")` |


> **Notes**
> - Sequence input uses the standard **one-letter code** (e.g., `RRKWLWLW`).  
> - Graphs are built with **Deep Graph Library (DGL)**.  
> - Throughout the repo, we use **forcefield** as a single word for consistency.
