configfile: "config.yaml"
container: "docker://filipafernandes/dps_structural_pipeline:006"

rule all:
    input:
        "data/tree/dps_struct_tree.nwk"

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

# 👇 STOP HERE FOR MANUAL MODELLER STEP
rule prepare_for_modeller:
    input:
        "data/raw/.done"
    output:
        "data/aligned/READY_FOR_MODELLER.txt"
    shell:
        """
        mkdir -p data/aligned
        echo "Run MODELLER manually and place structures_aligned.fasta here" > {output}
        """

# 👇 resumes AFTER you run MODELLER
rule build_tree:
    input:
        "data/aligned/structures_aligned.fasta"
    output:
        "data/tree/dps_struct_tree.nwk"
    shell:
        """
        python scripts/build_tree.py {input} {output}
        """
