configfile: "config.yaml"

rule all:
    input:
        "data/pdb_ids.txt",
        "data/raw/.done",
        "data/alignment/structural.ali"

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
