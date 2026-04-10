# dps_structural_pipeline

An automated and reproducible structural bioinformatics pipeline for retrieving protein structures and performing structural alignment.

---

## Overview

This project implements a fully automated workflow using **Snakemake** to:

* Retrieve protein structures (PDB files)
* Filter structures based on organism and protein family
* Preprocess and organize structural data
* Perform structural alignment using **MODELLER (SALIGN)**

The pipeline combines **containerized steps** for reproducibility with a **local execution step** for MODELLER due to licensing constraints.

---

## Pipeline Architecture


| Step                 | Description                      | Execution                      |
| -------------------- | -------------------------------- | ------------------------------ |
| Data Retrieval       | Query and fetch PDB structures   | Container (Docker/Singularity) |
| Preprocessing        | Organize and validate structures | Container                      |
| Structural Alignment | SALIGN via MODELLER              | Local environment              |

---

## How It Works

The workflow is running using **Snakemake**, which  manages dependencies between steps.

1. A structured query retrieves protein structures based on:
   * Organism (e.g. *Deinococcus*)
   * Protein family keywords (e.g. Dps, ferritin)
2. Structures are downloaded and standardized
3. A structural alignment is performed using **MODELLER SALIGN**

Most steps run inside a container for reproducibility. The alignment step runs locally due to MODELLER licensing.

---

## Query Configuration

The dataset is fully defined in `config.yaml`:

```yaml
threads: 4

query:
  organism: "Deinococcus"
  keywords:
    - "dps"
    - "dna protection protein"
    - "ferritin"
```

This ensures:
* Reproducibility
* Transparency of dataset selection
* Easy modification for other organisms or protein families

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

### What happens:

* Containerized steps:
  * Query RCSB
  * Download PDB files

* Local step:
  * Runs `scripts/salign.py` using MODELLER

---

##  Output

Results are generated in:

```
data/alignment/
└── structural.ali    # Structural alignment (PIR format)
```

---

##  Methodology

Structural alignment is performed using the **SALIGN algorithm** in MODELLER, which:
* Aligns protein structures based on 3D coordinates
* Optimizes residue correspondences
* Produces a multiple structural alignment

---

##  Reproducibility

This pipeline ensures reproducibility through:
* Workflow management with **Snakemake**
* Containerization of non-licensed steps
* Config-driven dataset definition (`config.yaml`)

Note: Full reproducibility requires a valid local installation of MODELLER.

---

##  Project Structure

```
project/
├── Snakefile
├── config.yaml
├── scripts/
│   ├── query_rcsb.py
│   └── salign.py
├── data/
│   ├── raw/
│   └── alignment/
├── envs/
└── README.md
```

![Pipeline DAG](dag.png)

---

## Notes

* The number of retrieved structures depends on availability in the PDB
* Dps proteins in *Deinococcus* are relatively limited
* The pipeline may return a small but biologically relevant dataset

---
