from modeller import *
import os
from Bio.PDB import PDBParser

log.verbose()
env = environ()
env.io.atom_files_directory = [os.path.abspath("data/raw")]

aln = alignment(env)
parser = PDBParser(QUIET=True)

def get_first_chain(pdb_path):
    structure = parser.get_structure("struct", pdb_path)
    for model_obj in structure:
        for chain in model_obj:
            return chain.id
    return None

# Read PDB IDs
with open("data/pdb_ids.txt") as f:
    pdbs = [line.strip().lower() for line in f if line.strip()]

if len(pdbs) < 2:
    raise ValueError("Need at least 2 PDBs")

# -------------------------------
# 🔥 Reference structure
# -------------------------------
ref_code = pdbs[0]
ref_file = f"data/raw/{ref_code}.pdb"

if not os.path.exists(ref_file):
    raise FileNotFoundError(f"Missing reference PDB: {ref_file}")

ref_chain = get_first_chain(ref_file)

print(f"REFERENCE: {ref_code} (chain {ref_chain})")

ref_model = model(
    env,
    file=ref_file,
    model_segment=(f'FIRST:{ref_chain}', f'LAST:{ref_chain}')
)

aln.append_model(
    ref_model,
    atom_files=ref_file,
    align_codes=f"{ref_code}_{ref_chain}"
)

# -------------------------------
# 🔥 Add remaining structures
# -------------------------------
loaded = 1

for code in pdbs[1:]:
    pdb_file = f"data/raw/{code}.pdb"

    if not os.path.exists(pdb_file):
        print(f"Missing file: {pdb_file}")
        continue

    try:
        chain = get_first_chain(pdb_file)
        if chain is None:
            print(f"No chain found in {code}")
            continue

        print(f"{code} -> using chain {chain}")

        mdl = model(
            env,
            file=pdb_file,
            model_segment=(f'FIRST:{chain}', f'LAST:{chain}')
        )

        aln.append_model(
            mdl,
            atom_files=pdb_file,
            align_codes=f"{code}_{chain}"
        )

        loaded += 1

    except Exception as e:
        print(f"Skipping {code}: {e}")

if loaded < 2:
    raise ValueError("Not enough valid structures")

# -------------------------------
# 🔥 SALIGN (STRUCTURAL PASS)
# -------------------------------
aln.salign(
    alignment_type='tree',
    feature_weights=(1., 0., 0., 0., 1., 0.),
    gap_penalties_1d=(-600, -100),
    gap_penalties_3d=(0, 3),
    overhang=30,
    improve_alignment=True,
    fit=True
)

# -------------------------------
# 🔥 SALIGN (FORCE GAP ALIGNMENT)
# -------------------------------
aln.salign(
    alignment_type='progressive',
    gap_penalties_1d=(-600, -100)
)

# -------------------------------
# OUTPUT
# -------------------------------
os.makedirs("data/alignment", exist_ok=True)

aln.write(
    file='data/alignment/structural.ali',
    alignment_format='PIR'
)

print("✅ SALIGN DONE")
