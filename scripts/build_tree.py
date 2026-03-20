from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo

# Read matrix
with open("data/tree/distance_matrix.txt") as f:
    lines = f.readlines()

names = lines[0].strip().split()[1:]
matrix = []

for line in lines[1:]:
    parts = line.strip().split()
    matrix.append([float(x) for x in parts[1:]])

# Convert to lower triangle
lower_triangle = []
for i in range(len(matrix)):
    lower_triangle.append(matrix[i][:i+1])

dm = DistanceMatrix(names, lower_triangle)

# Build Neighbor-Joining tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Save Newick
Phylo.write(tree, "data/tree/structural_tree.nwk", "newick")

print("Tree written to data/tree/structural_tree.nwk")
