configfile: "config.yaml"

rule all:
    input:
        "data/tree/structural_tree.nwk"

rule query_rcsb:
    output:
        "data/pdb_ids.txt"
    container: "docker://filipafernandes/dps_structural_pipeline:006"
    shell:
        "python scripts/query_rcsb.py config.yaml"

rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    container: "docker://filipafernandes/dps_structural_pipeline:006"
    shell:
        "python scripts/download_pdbs.py && touch {output}"

rule salign_tree:
    input:
        "data/raw/.done"
    output:
        "data/tree/structural.tree"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate modeller
        python scripts/salign.py
        """
        
rule compute_distances:
    input:
        "data/tree/structural.ali"
    output:
        "data/tree/distance_matrix.txt"
    container:
        "docker://filipafernandes/dps_structural_pipeline:006"
    shell:
        "scripts/compute_distances.py"

rule build_tree:
    input:
        "data/tree/distance_matrix.txt"
    output:
        "data/tree/structural_tree.nwk"
    container:
        "docker://filipafernandes/dps_structural_pipeline:006"
    shell:
        "scripts/build_tree.py"
