from modeller import *
import os
from Bio.PDB import PDBParser

log.verbose()
env = environ()
env.io.atom_files_directory = [os.path.abspath("data/raw")]

aln = alignment(env)
parser = PDBParser(QUIET=True)

# Detect first chain
def get_first_chain(pdb_path):
    structure = parser.get_structure("struct", pdb_path)
    for model_obj in structure:
        for chain in model_obj:
            return chain.id
    return None

# Read PDB IDs
with open("data/pdb_ids.txt") as f:
    pdbs = [line.strip().lower() for line in f if line.strip()]

loaded = 0

for code in pdbs:
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
    raise ValueError("Need at least 2 valid structures for SALIGN")

os.makedirs("data/alignment", exist_ok=True)

# Structural alignment
aln.salign(
    rms_cutoff=3.5,
    normalize_pp_scores=False,
    rr_file='$(LIB)/as1.sim.mat',
    overhang=30,
    gap_penalties_1d=(-450, -50),
    gap_penalties_3d=(0, 3),
    gap_gap_score=0,
    gap_residue_score=0,
    alignment_type='progressive',
    feature_weights=(1., 0., 0., 0., 1., 0.),
    improve_alignment=True,
    fit=True,
    write_fit=True,
    write_whole_pdb=False,
    output='ALIGNMENT QUALITY'
)

aln.write(file='data/alignment/structural.ali', alignment_format='PIR')

print("SALIGN DONE")
