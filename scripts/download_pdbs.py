import os
import json
import re
from Bio.PDB import PDBList

out_dir = "data/raw"
os.makedirs(out_dir, exist_ok=True)

with open("data/pdb_ids.txt") as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

# Load species map to patch unknowns after download
species_map_path = "data/species_map.json"
with open(species_map_path) as f:
    species_map = json.load(f)

pdbl = PDBList()
successful = []
failed     = []

def extract_organism_from_pdb_header(pdb_path):
    """
    Parse ORGANISM_SCIENTIFIC from PDB file SOURCE record.
    PDB files always have this — it's our last resort for unknown organisms.
    """
    try:
        with open(pdb_path) as f:
            in_source = False
            for line in f:
                record = line[:6].strip()
                if record == "SOURCE":
                    in_source = True
                elif in_source and record not in ("SOURCE", ""):
                    in_source = False

                if in_source:
                    # Look for ORGANISM_SCIENTIFIC: Name;
                    match = re.search(
                        r"ORGANISM_SCIENTIFIC:\s*([^;]+);",
                        line, re.IGNORECASE
                    )
                    if match:
                        name = match.group(1).strip().title()
                        return name
    except Exception:
        pass
    return None


for pdb_id in pdb_ids:
    try:
        print(f"Downloading {pdb_id}...")

        file_path = pdbl.retrieve_pdb_file(
            pdb_id,
            pdir=out_dir,
            file_format="pdb"
        )

        if not file_path or not os.path.exists(file_path):
            print(f"Failed: {pdb_id}")
            failed.append(pdb_id)
            continue

        new_path = os.path.join(out_dir, f"{pdb_id}.pdb")
        os.rename(file_path, new_path)
        successful.append(pdb_id)

        # Patch unknown organisms from PDB header
        current_label = species_map.get(pdb_id, "")
        if not current_label or current_label == pdb_id.upper():
            organism = extract_organism_from_pdb_header(new_path)
            if organism:
                print(f"  Resolved organism from PDB header: {pdb_id} -> {organism}")
                species_map[pdb_id] = organism

    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        failed.append(pdb_id)

# Save patched species map
with open(species_map_path, "w") as f:
    json.dump(species_map, f, indent=2)

print("\n--- SUMMARY ---")
print(f"Successful: {len(successful)}")
print(f"Failed:     {len(failed)}")

with open("data/failed_pdb_ids.txt", "w") as f:
    for pid in failed:
        f.write(pid + "\n")
