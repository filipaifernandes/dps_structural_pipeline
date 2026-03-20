# dps_structural_pipeline

An automated and reproducible structural bioinformatics pipeline for retrieving protein structures and performing structural alignment and dendrogram generation.

---

## Overview

This project implements a fully automated workflow using **Snakemake** to:

* Retrieve protein structures (PDB files)
* Preprocess and organize structural data
* Perform structural alignment using **MODELLER**
* Generate a structural dendrogram and alignment file

The pipeline combines **containerized steps** for reproducibility with a **local execution step** for MODELLER due to licensing constraints.

---

## Pipeline Architecture

| Step                 | Description                      | Execution                      |
| -------------------- | -------------------------------- | ------------------------------ |
| Data Retrieval       | Fetch PDB structures             | Container (Docker/Singularity) |
| Preprocessing        | Organize and validate structures | Container                      |
| Structural Alignment | SALIGN via MODELLER              | Local environment              |
| Output Generation    | Dendrogram + alignment           | Local                          |

---

## How It Works

The workflow is running using **Snakemake**, which automatically manages dependencies between steps.

Most steps run inside a container for reproducibility. However, the structural alignment step:

* Calls a Python script (`scripts/salign.py`)
* Uses **MODELLER**
* Runs **locally** (outside the container)

This hybrid design ensures both reproducibility and compatibility with licensed software.

---

## Important Note on MODELLER

**MODELLER** is not included in the container because:

* It requires a license key
* It must be installed locally

### Setup

1. Install MODELLER from the official website
2. Obtain and activate your license key
3. Create a dedicated environment:

```bash
conda create -n modeller python=3.10
conda activate modeller
pip install modeller
```

---

## Running the Pipeline

From the project root directory:

```bash
snakemake --use-singularity --cores 4
```

This will:

* Execute all containerized steps automatically
* Run the MODELLER step locally via the defined Snakemake rule

---

##  Output

Results are generated in:

```
data/tree/
├── structural.tree   # Structural dendrogram (Newick format)
├── structural.ali    # Structural alignment (PIR format)
```

---

##  Methodology

Structural alignment is performed using the SALIGN algorithm implemented in **MODELLER**, which:

* Aligns protein structures based on 3D coordinates
* Computes similarity scores
* Generates a dendrogram representing structural relationships

---

##  Reproducibility

This pipeline ensures reproducibility through:

* Workflow management with **Snakemake**
* Containerization of non-licensed steps
* Fixed directory structure and deterministic execution

Note: Full reproducibility requires a valid local installation of MODELLER.

---

##  Project Structure

```
project/
├── Snakefile
├── scripts/
│   └── salign.py
├── data/
│   ├── raw/
│   └── tree/
├── envs/
└── README.md
```

![Pipeline DAG](dag.png)

---

##  Final Notes

This pipeline reflects real-world bioinformatics practices by:

* Integrating multiple tools into a single workflow
* Handling software licensing constraints
* Balancing reproducibility with practical limitations

---
