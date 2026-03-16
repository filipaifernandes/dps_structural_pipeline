import os

os.makedirs("data/tree", exist_ok=True)

# List all PDBs
pdb_files = [f for f in os.listdir("data/raw") if f.endswith(".pdb")]
print("PDB files for alignment:", pdb_files)

# Dummy tree (replace with real Modeller alignment)
with open("data/tree/dps_struct_tree.nwk", "w") as f:
    f.write("(dummy_tree);")

print("Created structural tree at data/tree/dps_struct_tree.nwk")
