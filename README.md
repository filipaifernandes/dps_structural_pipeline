# DPS Structural Pipeline

**Automated and reproducible workflow for retrieving Dps protein structures, performing structural alignment, constructing phylogenetic trees, and visualizing structural divergence via RMSD heatmaps.**

Built with Snakemake and executed inside Docker containers — no manual steps, fully portable.

---

1. [Overview](#overview)
2. [Features](#features)
3. [Installation](#installation)
4. [MODELLER Setup](#modeller-setup)
5. [Configuration](#configuration)
6. [Pipeline Steps](#pipeline-steps)
7. [Output Structure](#output-structure)
8. [DAG](#dag)
9. [Reproducibility](#reproducibility)
10. [Troubleshooting](#troubleshooting)
11. [References](#references)
12. [Contact](#contact)

---

## Overview

The DPS Structural Pipeline provides a **reproducible and automated workflow for structural phylogenetics of Dps proteins across all species in the PDB**. It queries RCSB, selects the best structure per species, aligns them in 3D using MODELLER's SALIGN algorithm, builds a phylogenetic tree, generates an RMSD heatmap to visualize structural divergence, and produces a ready-to-use iTOL label file that replaces PDB codes with species names in the tree.

Three distinct analyses are combined:
- **Structural alignment** — SALIGN optimizes residue correspondences using 3D atomic coordinates, more robust than sequence-only methods for divergent proteins
- **RMSD heatmap** — pairwise Cα RMSD matrix across all structures, visualizing structural distance at a glance
- **iTOL annotation** — species names are fetched from RCSB for each PDB ID and written to a label file ready to upload to [iTOL](https://itol.embl.de/)

> !! Most steps run inside containers. The structural alignment step runs **locally** because MODELLER requires a license key.

---

## Features

- **Automated RCSB querying** — retrieves all Dps structures from the PDB
- **Smart structure selection** — keeps only the best structure per species (highest resolution)
- **3D structural alignment** — SALIGN algorithm via MODELLER
- **Phylogenetic tree** — maximum-likelihood inference with FastTree (Newick output)
- **RMSD heatmap** — pairwise Cα RMSD matrix visualized with seaborn
- **iTOL label file** — automatically maps PDB IDs to species names for publication-ready tree visualization
- **Containerized** — retrieval, conversion, tree-building, and label generation run in identical Docker environments
- **Config-driven** — swap keywords to analyze any protein family

---

## Installation

**Requirements:** Snakemake, Apptainer, and MODELLER (local)

```bash
# 1. Install Snakemake
conda install -c conda-forge -c bioconda snakemake

# 2. Install Apptainer (Ubuntu)
# Download apptainer_1.4.5_amd64.deb from https://github.com/apptainer/apptainer/releases/tag/v1.4.5
sudo apt install ./apptainer_1.4.5_amd64.deb

# 3. Clone the repo
git clone https://github.com/yourname/dps_structural_pipeline.git
cd dps_structural_pipeline

# 4. Run (with modeller environment active — see below)
conda activate modeller
snakemake --use-singularity --cores 4
```

The container image (`docker://filipafernandes/dps_structural_pipeline:010`) is pulled automatically.

---

## MODELLER Setup

MODELLER must be installed **locally** (not in the container) because it requires a personal license key.

```bash
# 1. Register and get a free academic license at https://salilab.org/modeller/registration.html

# 2. Create a dedicated environment
conda create -n modeller python=3.10
conda activate modeller
pip install modeller

# 3. Activate your license key (follow MODELLER's documentation)

# 4. Verify
python -c "from modeller import *; print('MODELLER ready!')"
```

Always activate this environment before running the pipeline:

```bash
conda activate modeller
snakemake --use-singularity --cores 4
```

---

## Configuration

All pipeline behaviour is controlled via `config.yaml`:

```yaml
threads: 4

query:
  keywords:
    - "Dps"
```

| Parameter | Description |
|---|---|
| `threads` | CPU cores to use |
| `keywords` | Protein names to search in RCSB PDB |

To repurpose the pipeline for another protein, just swap the keyword:

```yaml
query:
  keywords:
    - "catalase"
```

---

## Pipeline Steps

| Step | Tool | Execution | Output |
|---|---|---|---|
| RCSB query + best-per-species selection | `query_rcsb.py` | Container | `data/pdb_ids.txt` |
| Structure download | `download_pdbs.py` (Biopython) | Container | `data/raw/*.pdb` |
| Structural alignment | `salign.py` (MODELLER SALIGN) | **Local** | `data/alignment/structural.ali` |
| Format conversion PIR → FASTA | `ali_to_fasta.py` | Container | `data/alignment/structural.fasta` |
| Phylogenetic tree | FastTree | Container | `data/tree/tree.nwk` |
| RMSD heatmap | `rmsd_heatmap.py` (Biopython + seaborn) | Container | `data/heatmap/` |
| iTOL label file | `itol_labels.py` (RCSB API) | Container | `data/itol/labels.txt` |

### iTOL Label Generation

The `itol_labels` rule reads the structural alignment file (`structural.ali`), extracts each PDB ID, and queries the RCSB REST API to retrieve the scientific name of the source organism. The result is written as an iTOL-compatible label file with the format:

```
LABELS
SEPARATOR TAB
DATA
1e9y_A	Deinococcus_radiodurans
8ci9_C	Escherichia_coli
...
```

To visualize the tree with species names:
1. Go to [iTOL](https://itol.embl.de/) and upload `data/tree/tree.nwk`
2. Drag and drop `data/itol/labels.txt` onto the tree
3. Export as SVG or PDF for publication

> If any PDB ID cannot be resolved (e.g. due to a network timeout), the label falls back to `Unknown`. Re-running the rule will retry those entries.

---

## Output Structure

```
data/
├── pdb_ids.txt                  # Selected PDB IDs (one per species, best resolution)
├── raw/
│   └── *.pdb                    # Downloaded structure files
├── alignment/
│   ├── structural.ali           # Structural alignment (PIR format)
│   └── structural.fasta         # Alignment in FASTA format
├── tree/
│   └── tree.nwk                 # Phylogenetic tree (Newick format)
├── heatmap/
│   ├── rmsd_matrix.csv          # Pairwise Cα RMSD matrix
│   └── rmsd_heatmap.png         # Heatmap visualization
└── itol/
    └── labels.txt               # iTOL annotation file with species names
```

---

## DAG

![Pipeline DAG](dag.png)

Generate your own:
```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## Reproducibility

- Snakemake tracks dependencies and only reruns changed rules
- Docker containers pin all tool versions for retrieval and tree-building steps
- `config.yaml` makes the pipeline reusable for any protein/keyword
- MODELLER is explicitly managed in a dedicated local environment

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'modeller'`**
→ Activate the modeller environment: `conda activate modeller`

**`Apptainer not found`**
→ Verify with `apptainer --version` and reinstall if needed (see Installation)

**No structures retrieved / empty `pdb_ids.txt`**
→ Try broader keywords in `config.yaml` or check RCSB PDB directly for your protein

**License key error from MODELLER**
→ Re-register at [salilab.org/modeller](https://salilab.org/modeller/registration.html) and re-activate in the modeller environment

**Alignment fails with structural errors**
→ Inspect `data/raw/` to verify PDB files are valid; try `head -20 data/raw/*.pdb`. Very distant structures may need manual review.

**`labels.txt` contains many `Unknown` entries**
→ The RCSB API request timed out for some IDs. Re-run just the iTOL rule with `snakemake data/itol/labels.txt --use-singularity --cores 1` to retry.

For verbose output: `snakemake --use-singularity --cores 4 -v`

---

## References

- **Snakemake** — Köster & Rahmann, *Bioinformatics* 2012
- **MODELLER / SALIGN** — Šali & Blundell, *J Mol Biol* 1993; Madhusudhan et al., *Bioinformatics* 2006
- **FastTree** — Price et al., *PLoS ONE* 2010
- **Biopython** — Cock et al., *Bioinformatics* 2009
- **RCSB PDB** — Burley et al., *Nucleic Acids Res* 2021

---

## Contact

**Filipa Fernandes** — Bioinformatics Student
📧 [filipaifernandes.2005@gmail.com](mailto:filipaifernandes.2005@gmail.com)
