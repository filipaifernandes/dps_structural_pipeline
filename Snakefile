configfile: "config.yaml"
container: "docker://filipafernandes/dps_structural_pipeline:002"

rule all:
    input:
        "data/tree/dps_struct_tree.nwk"

rule query_rcsb:
    output:
        "data/pdb_ids.txt"
    shell:
        "python scripts/query_rcsb.py config.yaml > {output}"

rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    shell:
        "python scripts/download_pdbs.py && touch {output}"

rule modeller_struct_alignment:
    input:
        "data/raw/.done"
    output:
        "data/tree/dps_struct_tree.nwk"
    shell:
        "python scripts/modeller_struct_alignment.py {input} {output}"
