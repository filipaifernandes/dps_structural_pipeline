configfile: "config.yaml"

rule all:
    input:
        "data/pdb_ids.txt",
        "data/raw/.done",
        "data/alignment/structural.ali",
        "data/alignment/structural.fasta",
        "data/tree/tree.nwk",
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png",
        "data/itol/labels.txt"


# 1. Get PDB list
rule query_interpro:
    output:
        "data/pdb_ids.txt"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python3 scripts/query_rcsb.py"


# 2. Download PDBs
rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python3 scripts/download_pdbs.py && touch {output}"


# 3. Alignment
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

# 4. PIR → FASTA
rule ali_to_fasta:
    input:
        "data/alignment/structural.ali"
    output:
        "data/alignment/structural.fasta"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    script:
    	"scripts/itol_labels.py"


# 5. Tree
rule build_tree:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/tree/tree.nwk"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "fasttree {input} > {output}"


# 6. RMSD
rule rmsd_heatmap:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python3 scripts/rmsd_heatmap.py data/raw/ {output[0]} {output[1]}"


# 7. iTOL
rule itol_labels:
    input:
        "data/alignment/structural.ali"
    output:
        "data/itol/labels.txt"
    container:
        "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python3 scripts/itol_labels.py"
