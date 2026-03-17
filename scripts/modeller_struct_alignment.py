import os
import sys
import subprocess
from modeller import *

# -----------------------------
# Arguments
# -----------------------------
output_tree = sys.argv[1]
threads = sys.argv[2]

# -----------------------------
# Paths
# -----------------------------
raw_dir = "data/raw"
align_dir = "data/aligned"
tree_dir = os.path.dirname(output_tree)

os.makedirs(align_dir, exist_ok=True)
os.makedirs(tree_dir, exist_ok=True)

pir_file = os.path.join(align_dir, "structures.ali")
fasta_file = os.path.join(align_dir, "structures_aligned.fasta")
dendrogram_file = os.path.join(tree_dir, "modeller_dendrogram.tree")

# -----------------------------
# Find structure files
# -----------------------------
pdb_files = sorted([
    f for f in os.listdir(raw_dir)
    if f.endswith(".ent") or f.endswith(".pdb")
])

if len(pdb_files) < 2:
    raise ValueError("Need at least 2 structure files to build a structural tree.")

# -----------------------------
# MODELLER environment
# -----------------------------
env = environ()
env.io.atom_files_directory = [raw_dir]

aln = alignment(env)

# -----------------------------
# Load all structures
# -----------------------------
loaded = 0

for pdb_file in pdb_files:
    code = os.path.splitext(pdb_file)[0]

    try:
        mdl = model(env, file=code, model_segment=('FIRST:@', 'LAST:'))
        aln.append_model(mdl, atom_files=code, align_codes=code)
        loaded += 1
        print(f"Loaded structure: {code}")
    except Exception as e:
        print(f"Warning: failed to load {pdb_file}: {e}")

if loaded < 2:
    raise ValueError("Fewer than 2 structures could be loaded into MODELLER.")

# -----------------------------
# Structural alignment with SALIGN
# -----------------------------
aln.salign(
    rms_cutoff=3.5,
    normalize_pp_scores=False,
    rr_file='$(LIB)/as1.sim.mat',
    overhang=30,
    gap_penalties_1d=(-450, -50),
    gap_penalties_3d=(0, 3),
    gap_gap_score=0,
    gap_residue_score=0,
    dendrogram_file=dendrogram_file,
    alignment_type='tree',
    feature_weights=(1., 0., 0., 0., 1., 0.),
    improve_alignment=True,
    fit=True,
    write_fit=True,
    write_whole_pdb=False,
    output='ALIGNMENT QUALITY'
)

# -----------------------------
# Save aligned sequences
# -----------------------------
aln.write(file=pir_file, alignment_format='PIR')
aln.write(file=fasta_file, alignment_format='FASTA')

print(f"Structural alignment written to {pir_file}")
print(f"Aligned FASTA written to {fasta_file}")
print(f"MODELLER dendrogram written to {dendrogram_file}")

# -----------------------------
# Build tree from structural alignment
# -----------------------------
cmd = f"FastTree {fasta_file} > {output_tree}"
subprocess.run(cmd, shell=True, check=True)

print(f"Structural tree written to {output_tree}")
