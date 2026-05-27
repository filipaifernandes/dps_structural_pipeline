configfile: "config.yaml"

rule all:
    input:
        "data/pdb_ids.txt",
        "data/selection_report.tsv",
        "data/raw/.done",
        "data/alignment/structural.ali",
        "data/alignment/structural.fasta",
        "data/tree/tree.nwk",
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png",
        "data/itol/labels.txt"


# 1. Query InterPro + SIFTS -> pdb_ids.txt + species_map.json + report.tsv
rule query_interpro:
    output:
        "data/pdb_ids.txt",
        "data/species_map.json",
        "data/selection_report.tsv",
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        "python3 scripts/query_rcsb.py"


# 2. Download PDBs + patch unknown organisms from PDB headers
rule download_pdbs:
    input:
        "data/pdb_ids.txt",
        "data/species_map.json",
    output:
        "data/raw/.done"
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        "python3 scripts/download_pdbs.py && touch {output}"


# 3. Structural alignment with MODELLER SALIGN
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

# 4. PIR -> FASTA
rule ali_to_fasta:
    input:
        "data/alignment/structural.ali"
    output:
        "data/alignment/structural.fasta"
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        "python3 scripts/ali_to_fasta.py {input} {output}"


# 5. Phylogenetic tree
rule build_tree:
    input:
        "data/alignment/structural.fasta"
    output:
        "data/tree/tree.nwk"
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        """
        iqtree2 \
            -s {input} \
            -m LG+G4 \
            -blmin 1e-6 \
            -bnni \
            -nt AUTO \
            -redo \
            --prefix data/tree/tree

        cp data/tree/tree.treefile {output}
        """

# 6. RMSD heatmap
rule rmsd_heatmap:
    input:
        alignment="data/alignment/structural.fasta",
        done="data/raw/.done"
    output:
        "data/heatmap/rmsd_matrix.csv",
        "data/heatmap/rmsd_heatmap.png"
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        """
    	 python3 scripts/rmsd_heatmap.py \
         data/raw/ \
         {input.alignment} \
         {output[0]} \
         {output[1]}
    	"""

# 7. iTOL labels
rule itol_labels:
    input:
        "data/alignment/structural.ali",
        "data/species_map.json",
    output:
        "data/itol/labels.txt"
    container:
        "docker://filipafernandes/dps_structural_pipeline:012"
    shell:
        "python3 scripts/itol_labels.py {input[0]} {output}"
