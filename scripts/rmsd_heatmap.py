import sys
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Superimposer

# Required for headless Docker
plt.switch_backend("Agg")

pdb_dir = sys.argv[1]
output_matrix = sys.argv[2]
output_plot = sys.argv[3]

parser = PDBParser(QUIET=True)

pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith(".pdb")])
names = [os.path.splitext(f)[0] for f in pdb_files]

structures = [parser.get_structure(name, os.path.join(pdb_dir, f))
              for name, f in zip(names, pdb_files)]

def get_ca_atoms(structure):
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    atoms.append(residue["CA"])
    return atoms

n = len(structures)
matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        atoms1 = get_ca_atoms(structures[i])
        atoms2 = get_ca_atoms(structures[j])

        min_len = min(len(atoms1), len(atoms2))
        atoms1 = atoms1[:min_len]
        atoms2 = atoms2[:min_len]

        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)
        matrix[i, j] = sup.rms

df = pd.DataFrame(matrix, index=names, columns=names)

# Save matrix
df.to_csv(output_matrix)

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(df, cmap="magma")
plt.title("Dps Structural RMSD Heatmap")
plt.tight_layout()
plt.savefig(output_plot)
