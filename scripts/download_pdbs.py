import os
from Bio.PDB import PDBList

# Create output directory
out_dir = "data/raw"
os.makedirs(out_dir, exist_ok=True)

# Read PDB IDs
with open("data/pdb_ids.txt") as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

pdbl = PDBList()

successful = []
failed = []

for pdb_id in pdb_ids:
    try:
        print(f"Downloading {pdb_id}...")

        file_path = pdbl.retrieve_pdb_file(
            pdb_id,
            pdir=out_dir,
            file_format="pdb"
        )

        # Handle failed download (returns None or missing file)
        if not file_path or not os.path.exists(file_path):
            print(f"Failed: {pdb_id}")
            failed.append(pdb_id)
            continue

        # Rename to clean filename
        new_path = os.path.join(out_dir, f"{pdb_id}.pdb")
        os.rename(file_path, new_path)

        successful.append(pdb_id)

    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        failed.append(pdb_id)

print("\n--- SUMMARY ---")
print(f"Successful: {len(successful)}")
print(f"Failed: {len(failed)}")

# Save failed IDs for debugging
with open("data/failed_pdb_ids.txt", "w") as f:
    for pid in failed:
        f.write(pid + "\n")
