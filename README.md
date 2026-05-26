# DPS Structural Pipeline

**Automated and reproducible workflow for retrieving Dps protein structures, performing structural alignment, constructing phylogenetic trees, and visualizing structural divergence via RMSD heatmaps.**

Built with Snakemake and executed inside Docker/Apptainer containers — no manual steps, fully portable.

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

The DPS Structural Pipeline provides a **reproducible and automated workflow for structural phylogenetics of Dps proteins across all species in the PDB**. It queries InterPro and EBI databases, selects the best structure per species (with paralog-aware detection), aligns them in 3D using MODELLER's SALIGN algorithm, builds a phylogenetic tree, generates an RMSD heatmap to visualize structural divergence, and produces a ready-to-use iTOL label file with species names.

Three distinct analyses are combined:

- **Structural alignment** — SALIGN optimizes residue correspondences using 3D atomic coordinates, more robust than sequence-only methods for divergent proteins
- **RMSD heatmap** — pairwise Cα RMSD matrix across all structures, visualizing structural distance at a glance
- **iTOL annotation** — species names are mapped to each PDB ID and written to a label file ready to upload to [iTOL](https://itol.embl.de/)

> !! Most steps run inside containers. The structural alignment step runs **locally** because MODELLER requires a license key.

---

## Features

- **InterPro-driven querying** — retrieves all Dps family structures using HMM-validated family membership (IPR002177, IPR014490), not keyword search — biologically rigorous
- **SIFTS supplement** — cross-references UniProt accessions against the EBI SIFTS flat file to catch structures that InterPro's structure endpoint missed; SIFTS cache is automatically refreshed every 30 days
- **Paralog-aware selection** — groups structures by UniProt accession per species; species with multiple distinct DPS proteins (e.g. *Deinococcus radiodurans* DrDps1/DrDps2, *Lactococcus lactis* DpsA/DpsB) keep one representative per paralog
- **Contamination filter** — validates every structure's UniProt accession against the InterPro-derived DPS family set, removing any mis-annotated non-DPS entries
- **Organism resolution** — unknown organisms are resolved from PDB file headers after download, avoiding any reliance on RCSB API calls
- **Manual include support** — structures not yet annotated in InterPro can be added via `config.yaml` with a documented reason; fetched automatically from PDBe
- **Selection report** — full audit trail in `data/selection_report.tsv` listing every selected structure with PDB ID, organism, resolution, UniProt accession, source, and paralog flag
- **3D structural alignment** — SALIGN algorithm via MODELLER
- **Phylogenetic tree** — maximum-likelihood inference with FastTree (Newick output)
- **RMSD heatmap** — pairwise Cα RMSD matrix visualized with seaborn
- **iTOL label file** — automatically maps PDB IDs to species names for publication-ready tree visualization
- **Containerized** — all steps except SALIGN run in identical Docker/Apptainer environments
- **Config-driven** — swap InterPro IDs to analyze any protein family with zero code changes

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

The container image (`docker://filipafernandes/dps_structural_pipeline:010`) is pulled automatically on first run.

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
  interpro_ids:
    - IPR002177   # DPS family — DNA-binding protein from starved cells
    - IPR014490   # DPS-like family — archaeal/alternative DPS
  manual_include:
    - pdb_id: 1zuj
      reason: "Lactococcus lactis DpsA - genuine DPS paralog, not annotated in InterPro"
    - pdb_id: 1zs3
      reason: "Lactococcus lactis DpsB - genuine DPS paralog, not annotated in InterPro"

sifts_max_age_days: 30
batch_size: 50
```

| Parameter | Description |
|---|---|
| `threads` | CPU cores to use |
| `interpro_ids` | InterPro family accessions to query — drives all structure retrieval |
| `manual_include` | Structures to force-include regardless of InterPro annotation, with mandatory reason |
| `sifts_max_age_days` | How many days to cache the SIFTS flat file before re-downloading (default: 30) |
| `batch_size` | Batch size for API requests |

### Repurposing for another protein family

To analyze a different protein family, replace the InterPro IDs:

```yaml
query:
  interpro_ids:
    - IPR001133   # ferritin
```

Zero code changes needed.

---

## Pipeline Steps

| Step | Rule | Tool | Execution | Output |
|---|---|---|---|---|
| Structure retrieval | `query_interpro` | InterPro API + SIFTS + PDBe | Container | `data/pdb_ids.txt`, `data/species_map.json`, `data/selection_report.tsv` |
| Structure download | `download_pdbs` | Biopython PDBList | Container | `data/raw/*.pdb` |
| Structural alignment | `salign_alignment` | MODELLER SALIGN | **Local** | `data/alignment/structural.ali` |
| Format conversion | `ali_to_fasta` | custom script | Container | `data/alignment/structural.fasta` |
| Phylogenetic tree | `build_tree` | FastTree | Container | `data/tree/tree.nwk` |
| RMSD heatmap | `rmsd_heatmap` | Biopython + seaborn | Container | `data/heatmap/` |
| iTOL labels | `itol_labels` | custom script | Container | `data/itol/labels.txt` |

### Structure Retrieval Strategy

Structure retrieval is the most scientifically critical step. It uses a three-source approach:

**Source 1 — InterPro structure endpoint**
Queries `/api/structure/pdb/entry/interpro/{IPR_ID}/` for each configured InterPro ID. Returns PDB entries whose sequences have been matched to the family HMM profile — the most biologically rigorous definition of family membership. Fast (~130 structures for IPR002177 + IPR014490).

**Source 2 — SIFTS flat file**
Downloads `pdb_chain_uniprot.csv.gz` from EBI FTP (cached locally, refreshed every `sifts_max_age_days` days). Cross-references the UniProt accessions collected in Source 1 against every PDB chain in SIFTS — catches structures that InterPro's structure endpoint missed due to incomplete annotation. Also patches missing UniProt accessions for structures whose organism couldn't be resolved.

**Source 3 — Manual include**
Structures explicitly listed in `config.yaml` under `manual_include`. Used for historically unannotated but scientifically validated structures (e.g. `1zuj`, `1zs3` — *Lactococcus lactis* DpsA/DpsB from 2005, predating comprehensive InterPro annotation). Each entry requires a documented reason. Resolution and organism are fetched automatically from PDBe — nothing is truly manual.

### Paralog-Aware Selection

Most species have one DPS protein — the pipeline keeps its best-resolution structure. Some species have genuine DPS paralogs with distinct UniProt accessions (different proteins, not crystal forms of the same protein). These are kept separately:

- *Deinococcus radiodurans* — DrDps1 (`Q9RS64`) + DrDps2 (`Q9RZN1`)
- *Pseudomonas aeruginosa* — canonical DPS (`Q9I4Z7`, IPR002177) + DPS-like (`Q9HUT3`, IPR014490)
- *Lactococcus lactis* — DpsA (`1zuj`) + DpsB (`1zs3`)

Structures with the same UniProt accession (same protein, different crystal forms or conditions) always collapse to the best resolution.

### Organism Resolution

Organism names are resolved in order of priority:
1. InterPro structure response (`source_organism.scientificName`)
2. InterPro protein endpoint (fallback for structures where Source 1 has no organism)
3. PDB file header `SOURCE` record (parsed after download — catches structures with no API organism)

Species names are normalised to binomial (*Genus species*), dropping strain and subspecies information.

### iTOL Label Generation

The `itol_labels` rule reads `data/alignment/structural.ali`, extracts each PDB ID, and maps it to a species name using `data/species_map.json` (built during retrieval). The result is an iTOL-compatible label file:

```
LABELS
SEPARATOR TAB
DATA
2c2u	Deinococcus radiodurans
8ouc	Escherichia coli
...
```

To visualize:
1. Go to [iTOL](https://itol.embl.de/) and upload `data/tree/tree.nwk`
2. Drag and drop `data/itol/labels.txt` onto the tree
3. Export as SVG or PDF

---

## Output Structure

```
data/
├── pdb_ids.txt                  # Selected PDB IDs (one per paralog per species)
├── species_map.json             # PDB ID -> species name mapping
├── selection_report.tsv         # Full audit trail of structure selection
├── pdb_chain_uniprot.csv.gz     # Cached SIFTS flat file (auto-refreshed)
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

### Selection Report (`selection_report.tsv`)

Full audit trail of every structure selected, with columns:

| Column | Description |
|---|---|
| `pdb_id` | PDB accession (uppercase) |
| `organism` | Source organism (binomial) |
| `resolution_A` | Crystal resolution in Ångströms |
| `uniprot` | UniProt accession |
| `source` | Where it was found: `IPR002177`, `IPR014490`, `sifts`, or `manual` |
| `is_paralog` | `yes` if the species has multiple DPS paralogs retained |

---

## DAG

![Pipeline DAG](dag.png)

Generate your own:
```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## Reproducibility

- Snakemake tracks all file dependencies — only reruns steps whose inputs changed
- Docker/Apptainer containers pin all tool versions for every containerized step
- `config.yaml` is the single source of truth for all scientific decisions — InterPro IDs, manual inclusions, SIFTS refresh frequency
- `selection_report.tsv` provides a complete audit trail of every structure in the analysis
- SIFTS cache age is configurable and logged on every run
- Manual inclusions are documented with reasons directly in `config.yaml` — no silent curation

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'modeller'`**
→ Activate the modeller environment: `conda activate modeller`

**`Apptainer not found`**
→ Verify with `apptainer --version` and reinstall if needed (see Installation)

**Folder path with spaces or special characters causes Singularity error**
→ Move the pipeline to a path without parentheses or spaces, e.g. `~/dps_pipeline/`

**Empty `pdb_ids.txt`**
→ Check InterPro is reachable: `curl https://www.ebi.ac.uk/interpro/api/entry/interpro/IPR002177/`
→ Verify InterPro IDs in `config.yaml` are valid

**SIFTS download fails**
→ Check EBI FTP is reachable: `curl -I https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz`
→ Delete `data/pdb_chain_uniprot.csv.gz` and rerun to force fresh download

**Structures missing from expected set**
→ Check `data/selection_report.tsv` — the `source` column shows where each structure came from
→ If a known DPS structure is absent, add it to `manual_include` in `config.yaml` with a reason

**Alignment fails with structural errors**
→ Inspect `data/raw/` to verify PDB files are valid: `head -5 data/raw/*.pdb`
→ Very distant structures may require manual review

**`labels.txt` has PDB IDs instead of species names**
→ The organism could not be resolved from any source
→ Check `data/species_map.json` for that PDB ID
→ The PDB file header parser in `download_pdbs.py` should have resolved it — check the download log

For verbose Snakemake output: `snakemake --use-singularity --cores 4 -v`

---

## References

- **Snakemake** — Köster & Rahmann, *Bioinformatics* 2012
- **MODELLER / SALIGN** — Šali & Blundell, *J Mol Biol* 1993; Madhusudhan et al., *Bioinformatics* 2006
- **FastTree** — Price et al., *PLoS ONE* 2010
- **Biopython** — Cock et al., *Bioinformatics* 2009
- **InterPro** — Paysan-Lafosse et al., *Nucleic Acids Res* 2023
- **SIFTS** — Dana et al., *Nucleic Acids Res* 2019
- **PDBe** — Armstrong et al., *Nucleic Acids Res* 2020
- **RCSB PDB** — Burley et al., *Nucleic Acids Res* 2021

---

## Contact

**Filipa Fernandes** — Bioinformatics Student
📧 [filipaifernandes.2005@gmail.com](mailto:filipaifernandes.2005@gmail.com)
