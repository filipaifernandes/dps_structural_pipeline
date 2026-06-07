import sys
import os
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import AlignIO
from Bio.PDB import PDBParser, Superimposer

plt.switch_backend("Agg")

pdb_dir = sys.argv[1]
alignment_file = sys.argv[2]
output_matrix = sys.argv[3]
output_plot = sys.argv[4]
# Optional 5th arg: path to species_map.json for human-readable axis labels
species_map = {}
if len(sys.argv) > 5 and os.path.exists(sys.argv[5]):
    with open(sys.argv[5]) as _f:
        species_map = json.load(_f)
    print(f"Loaded species map: {len(species_map)} entries")

parser = PDBParser(QUIET=True)

# -----------------------------
# Load alignment
# -----------------------------

alignment = AlignIO.read(alignment_file, "fasta")

aligned_sequences = {}

for record in alignment:

    pdb_id = record.id.split("_")[0].lower()

    aligned_sequences[pdb_id] = str(record.seq)

# -----------------------------
# Load structures
# -----------------------------

alignment_ids = sorted(aligned_sequences.keys())

pdb_files = [
    f"{pdb_id}.pdb"
    for pdb_id in alignment_ids
    if os.path.exists(os.path.join(pdb_dir, f"{pdb_id}.pdb"))
]

names = alignment_ids

all_names = [
    os.path.splitext(f)[0].lower()
    for f in pdb_files
]

names = [
    n for n in all_names
    if n in aligned_sequences
]

missing = sorted(set(all_names) - set(names))

if missing:
    print("\nWARNING: missing from alignment:")
    for m in missing:
        print(m)

structures = {
    name: parser.get_structure(name, os.path.join(pdb_dir, f))
    for name, f in zip(names, pdb_files)
}

# -----------------------------
# Extract ordered CA atoms
# -----------------------------

def get_ca_atoms(structure):
    atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:

                if residue.id[0] != " ":
                    continue

                if "CA" in residue:
                    atoms.append(residue["CA"])

    return atoms

# -----------------------------
# Alignment-aware atom mapping
# -----------------------------

def aligned_atom_pairs(name1, name2):

    seq1 = aligned_sequences[name1]
    seq2 = aligned_sequences[name2]

    atoms1_all = get_ca_atoms(structures[name1])
    atoms2_all = get_ca_atoms(structures[name2])

    idx1 = 0
    idx2 = 0

    atoms1 = []
    atoms2 = []

    for a, b in zip(seq1, seq2):

        atom1 = None
        atom2 = None

        if a != "-":
            if idx1 < len(atoms1_all):
                atom1 = atoms1_all[idx1]
            idx1 += 1

        if b != "-":
            if idx2 < len(atoms2_all):
                atom2 = atoms2_all[idx2]
            idx2 += 1

        if atom1 is not None and atom2 is not None:
            atoms1.append(atom1)
            atoms2.append(atom2)

    return atoms1, atoms2

# -----------------------------
# RMSD matrix
# -----------------------------

n = len(names)

matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):

        atoms1, atoms2 = aligned_atom_pairs(
            names[i],
            names[j]
        )

        if len(atoms1) < 3:
            matrix[i, j] = np.nan
            continue

        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)

        matrix[i, j] = sup.rms

df = pd.DataFrame(
    matrix,
    index=names,
    columns=names
)

df.to_csv(output_matrix)

# -----------------------------
# Heatmap
# Use species names if available, otherwise fall back to PDB IDs
# -----------------------------

def make_label(pdb_id):
    species = species_map.get(pdb_id)
    if species and species.lower() != "unknown":
        return f"{species} ({pdb_id.upper()})"
    return pdb_id.upper()

display_labels = [make_label(n) for n in names]

df_display = df.copy()
df_display.index   = display_labels
df_display.columns = display_labels

figsize = max(10, len(names) * 0.4)
plt.figure(figsize=(figsize, figsize * 0.8))

sns.heatmap(
    df_display,
    cmap="magma",
    square=True,
    xticklabels=True,
    yticklabels=True
)

plt.xticks(fontsize=7, rotation=45, ha="right")
plt.yticks(fontsize=7, rotation=0)
plt.title("Alignment-aware Cα RMSD")

plt.tight_layout()

plt.savefig(output_plot)
