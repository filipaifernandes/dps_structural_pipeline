import os
import re
from Bio.PDB import PDBList

os.makedirs("data/raw", exist_ok=True)

with open("data/pdb_ids.txt") as f:
    pdbs = [line.strip() for line in f if line.strip()]

valid_pdbs = [p.lower() for p in pdbs if re.fullmatch(r"[A-Za-z0-9]{4}", p)]

if not valid_pdbs:
    raise ValueError("No valid PDB IDs found in data/pdb_ids.txt")

pdbl = PDBList()

for pdb in valid_pdbs:
    print(f"Downloading {pdb}...")
    pdbl.retrieve_pdb_file(pdb, pdir="data/raw", file_format="pdb", overwrite=True)
