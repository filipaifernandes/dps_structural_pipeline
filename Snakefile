configfile: "config.yaml"
container: "docker://filipafernandes/dps_structural_pipeline:006"

rule all:
    input:
        "data/tree/structural.tree"

rule query_rcsb:
    output:
        "data/pdb_ids.txt"
    shell:
        "python scripts/query_rcsb.py config.yaml"

rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    shell:
        "python scripts/download_pdbs.py && touch {output}"
l
rule salign_tree:
    input:
        "data/raw/.done"
    output:
        "data/tree/structural.tree"
    shell:
        "python scripts/salign.py"
