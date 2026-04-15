configfile: "config.yaml"

rule all:
    input:
        "data/pdb_ids.txt",
        "data/raw/.done",
        "data/alignment/structural.ali",
        "data/tree/tree.nwk",
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png"

rule query_rcsb:
    output:
        "data/pdb_ids.txt"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        "python scripts/query_rcsb.py config.yaml"

rule download_pdbs:
    input:
        "data/pdb_ids.txt"
    output:
        "data/raw/.done"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
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
    
rule ali_to_fasta:
    input:
        "data/alignment/structural.ali"
    output:
        "data/alignment/structural.fasta"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        """
        python3 scripts/ali_to_fasta.py {input} {output}
        """

rule build_tree:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/tree/tree.nwk"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        """
        fasttree {input} > {output}
        """

rule rmsd_heatmap:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png"
    container: "docker://filipafernandes/dps_structural_pipeline:010"
    shell:
        """
        python scripts/rmsd_heatmap.py data/raw/ {output[0]} {output[1]}
        """
