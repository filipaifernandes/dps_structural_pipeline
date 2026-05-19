import sys
import json
import os

input_file  = sys.argv[1]   # data/alignment/structural.ali
output_file = sys.argv[2]   # data/itol/labels.txt

# =========================================================
# Load species map saved by query_rcsb.py
# This avoids any RCSB calls (blocked inside container)
# =========================================================

try:
    with open("data/species_map.json") as f:
        species_map = json.load(f)
    print(f"Loaded species map: {len(species_map)} entries", flush=True)
except FileNotFoundError:
    print("WARNING: data/species_map.json not found — labels will be PDB IDs only")
    species_map = {}

# =========================================================
# Parse alignment file for PDB IDs
# =========================================================

labels = []

with open(input_file) as f:
    for line in f:
        if line.startswith(">P1;"):
            try:
                # Extract PDB ID from header e.g. >P1;1dpsA
                raw = line.strip().split(";")[1]
                pdb_id = raw[:4].lower()

                species = species_map.get(pdb_id, "Unknown")

                print(f"  {pdb_id} -> {species}", flush=True)
                labels.append(f"{pdb_id}\t{species}")

            except Exception as e:
                print(f"  Error parsing line '{line.strip()}': {e}", flush=True)

# =========================================================
# Write iTOL labels file
# =========================================================

os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out:
    out.write("LABELS\n")
    out.write("SEPARATOR TAB\n")
    out.write("DATA\n")
    for line in labels:
        out.write(line + "\n")

print(f"\nSaved {len(labels)} labels -> {output_file}\n")
