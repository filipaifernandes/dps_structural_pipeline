import os
from Bio.PDB import PDBList

os.makedirs("data/raw", exist_ok=True)

with open("data/pdb_ids.txt") as f:
    pdbs = [line.strip() for line in f]

pdbl = PDBList()
for pdb in pdbs:
    pdbl.retrieve_pdb_file(pdb, pdir="data/raw", file_format="pdb", file_name=f"{pdb}.pdb")
