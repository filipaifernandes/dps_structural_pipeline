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

rule modeller_struct_alignment:
    input:
        "data/raw/.done"
    output:
        "data/tree/dps_struct_tree.nwk"
    container:
        "docker://filipafernandes/dps_structural_pipeline:006"
    shell:
        r"""
        python - <<'EOF'
import glob
import os

key = "MODELIRANJE"

patterns = [
    "/opt/conda/envs/snakemake/lib/modeller-*/modlib/modeller/config.py",
    "/opt/conda/envs/snakemake/lib/python*/site-packages/modeller/config.py",
]

patched = []

for pattern in patterns:
    for path in glob.glob(pattern):
        with open(path, "r") as f:
            txt = f.read()

        txt2 = txt.replace("license = 'XXXX'", f"license = '{key}'")
        txt2 = txt2.replace('license = "XXXX"', f'license = "{key}"')

        with open(path, "w") as f:
            f.write(txt2)

        patched.append(path)

print("Patched:", patched)
EOF

        python scripts/modeller_struct_alignment.py {output}
        """
