# scripts/download_pdbs.py
import os
import re
from Bio.PDB import PDBList

os.makedirs("data/raw", exist_ok=True)

# Read PDB IDs
with open("data/pdb_ids.txt") as f:
    pdbs = [line.strip() for line in f if line.strip()]

# Keep only valid 4-character PDB IDs
valid_pdbs = sorted([p.lower() for p in pdbs if re.fullmatch(r"[A-Za-z0-9]{4}", p)])

pdbl = PDBList()

for pdb in valid_pdbs:
    print(f"Downloading {pdb}...")
    file_path = pdbl.retrieve_pdb_file(
        pdb,
        pdir="data/raw",
        file_format="pdb",
        overwrite=True
    )

    # Rename to clean lowercase
    new_name = os.path.join("data/raw", f"{pdb}.pdb")
    if os.path.exists(file_path):
        os.rename(file_path, new_name)

print(f"Downloaded {len(valid_pdbs)} PDB files")
