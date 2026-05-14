configfile: "config.yaml"

rule all:
    input:
        "data/pdb_ids.txt",
        "data/raw/.done",
        "data/alignment/structural.ali",
        "data/tree/tree.nwk",
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png",
        "data/itol/labels.txt"

# -----------------------------------------------------------------------
# Step 1: Query InterPro (IPR002177) -> SIFTS -> best PDB per species
# -----------------------------------------------------------------------
rule query_interpro:
    output:
        "data/pdb_ids.txt"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python scripts/query_rcsb.py"

# -----------------------------------------------------------------------
# Step 2: Download PDB files
# -----------------------------------------------------------------------
rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python scripts/download_pdbs.py && touch {output}"

# -----------------------------------------------------------------------
# Step 3: Structural alignment with MODELLER SALIGN
# -----------------------------------------------------------------------
rule salign_alignment:
    input:
        "data/raw/.done"
    output:
        "data/alignment/structural.ali"
    singularity: None
    shell:
        """
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate modeller
        python scripts/salign.py
        """

# -----------------------------------------------------------------------
# Step 4: Convert PIR alignment to FASTA
# -----------------------------------------------------------------------
rule ali_to_fasta:
    input:
        "data/alignment/structural.ali"
    output:
        "data/alignment/structural.fasta"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python3 scripts/ali_to_fasta.py {input} {output}"

# -----------------------------------------------------------------------
# Step 5: Phylogenetic tree
# -----------------------------------------------------------------------
rule build_tree:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/tree/tree.nwk"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "fasttree {input} > {output}"

# -----------------------------------------------------------------------
# Step 6: RMSD heatmap
# -----------------------------------------------------------------------
rule rmsd_heatmap:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python scripts/rmsd_heatmap.py data/raw/ {output[0]} {output[1]}"

# -----------------------------------------------------------------------
# Step 7: iTOL species labels
# -----------------------------------------------------------------------
rule itol_labels:
    input:
        "data/alignment/structural.ali"
    output:
        "data/itol/labels.txt"
    script:
        "scripts/itol_labels.py"
