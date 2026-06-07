import sys
import json
import os

input_file = sys.argv[1]
output_file = sys.argv[2]

# =========================================================
# Load species map
# =========================================================

try:
    with open("data/species_map.json") as f:
        species_map = json.load(f)

    print(f"Loaded species map: {len(species_map)} entries", flush=True)

except FileNotFoundError:
    print(
        "WARNING: data/species_map.json not found — labels will use PDB IDs only",
        flush=True
    )
    species_map = {}

# =========================================================
# Parse alignment file
# =========================================================

labels = []

with open(input_file) as f:
    for line in f:

        if not line.startswith(">P1;"):
            continue

        try:
            raw = line.strip().split(";")[1]

            # Full node name exactly as it appears in the tree
            node_id = raw

            # Species lookup uses only the PDB accession
            pdb_id = raw[:4].lower()

            species = species_map.get(pdb_id, pdb_id)

            print(f"{node_id} -> {species}", flush=True)

            labels.append(f"{node_id}\t{species}")

        except Exception as e:
            print(
                f"Error parsing line '{line.strip()}': {e}",
                flush=True
            )

# =========================================================
# Write iTOL labels dataset
# =========================================================

os.makedirs(os.path.dirname(output_file), exist_ok=True)

with open(output_file, "w") as out:

    out.write("LABELS\n")
    out.write("SEPARATOR TAB\n")
    out.write("DATA\n")

    for label in labels:
        out.write(label + "\n")

print(f"\nSaved {len(labels)} labels -> {output_file}\n")
