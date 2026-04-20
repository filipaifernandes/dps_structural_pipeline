import requests
import time

input_file = snakemake.input[0]
output_file = snakemake.output[0]

def get_species(pdb_id):
    try:
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            return "Unknown"

        data = r.json()
        orgs = data.get("rcsb_entity_source_organism", [])

        if orgs:
            return orgs[0].get("scientific_name", "Unknown")

        return "Unknown"

    except:
        return "Unknown"


labels = []

# extract IDs from .ali
with open(input_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            seq_id = line.split(";")[1]
            pdb_id = seq_id.split("_")[0].lower()
            labels.append((seq_id, pdb_id))

mapping = []

for seq_id, pdb_id in labels:
    species = get_species(pdb_id)
    species_clean = species.replace(" ", "_")

    print(f"{seq_id} -> {species}")
    mapping.append((seq_id, species_clean))

    time.sleep(0.2)

# write iTOL file
with open(output_file, "w") as out:
    out.write("LABELS\nSEPARATOR TAB\nDATA\n")
    for seq_id, species in mapping:
        out.write(f"{seq_id}\t{species}\n")
