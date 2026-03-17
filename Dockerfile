FROM snakemake/snakemake:latest

ARG MODELLER_KEY
ENV MODELLER_KEY=${MODELLER_KEY}

RUN conda install -y -c bioconda -c conda-forge \
    biopython requests

RUN conda install -y -c salilab -c conda-forge modeller

# Patch every likely MODELLER config.py
RUN python - <<EOF
import os, glob

key = os.environ.get("MODELLER_KEY", "").strip()
if not key:
    raise RuntimeError("MODELLER_KEY not provided")

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

if not patched:
    raise RuntimeError("No MODELLER config.py found")

print("Patched MODELLER config files:")
for p in patched:
    print(" -", p)
EOF

WORKDIR /workflow
