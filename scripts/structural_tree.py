import sys
import os
import numpy as np

from Bio import AlignIO
from Bio.PDB import PDBParser, Superimposer
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import to_tree

pdb_dir = sys.argv[1]
alignment_file = sys.argv[2]
output_tree = sys.argv[3]

parser = PDBParser(QUIET=True)

alignment = AlignIO.read(alignment_file, "fasta")

aligned_sequences = {}

for record in alignment:
    pdb_id = record.id.split("_")[0].lower()
    aligned_sequences[pdb_id] = str(record.seq)

names = sorted(aligned_sequences.keys())

structures = {
    name: parser.get_structure(
        name,
        os.path.join(pdb_dir, f"{name}.pdb")
    )
    for name in names
}


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
            atom1 = atoms1_all[idx1]
            idx1 += 1

        if b != "-":
            atom2 = atoms2_all[idx2]
            idx2 += 1

        if atom1 is not None and atom2 is not None:
            atoms1.append(atom1)
            atoms2.append(atom2)

    return atoms1, atoms2


n = len(names)

matrix = np.zeros((n, n))

for i in range(n):
    for j in range(i + 1, n):

        atoms1, atoms2 = aligned_atom_pairs(
            names[i],
            names[j]
        )

        if len(atoms1) < 3:
            rmsd = 100.0
        else:
            sup = Superimposer()
            sup.set_atoms(atoms1, atoms2)
            rmsd = sup.rms

        matrix[i, j] = rmsd
        matrix[j, i] = rmsd


condensed = squareform(matrix)

Z = linkage(
    condensed,
    method="average"
)


def build_newick(node, parent_dist, leaf_names):

    if node.is_leaf():
        return "%s:%.4f" % (
            leaf_names[node.id],
            parent_dist - node.dist
        )

    left = build_newick(
        node.left,
        node.dist,
        leaf_names
    )

    right = build_newick(
        node.right,
        node.dist,
        leaf_names
    )

    return "(%s,%s):%.4f" % (
        left,
        right,
        parent_dist - node.dist
    )


tree = to_tree(Z)

newick = build_newick(
    tree,
    tree.dist,
    names
) + ";"

with open(output_tree, "w") as f:
    f.write(newick)
