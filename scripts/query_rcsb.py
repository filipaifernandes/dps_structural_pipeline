# scripts/query_rcsb.py
import requests
import os
import yaml
import sys

# -------------------------------
# Load config
# -------------------------------
config_file = sys.argv[1]
with open(config_file) as f:
    config = yaml.safe_load(f)

keywords = [k.lower() for k in config["query"]["keywords"]]

os.makedirs("data", exist_ok=True)

# -------------------------------
# RCSB search query
# -------------------------------
url = "https://search.rcsb.org/rcsbsearch/v2/query"

query = {
    "query": {
        "type": "group",
        "logical_operator": "or",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "value": "DNA-binding protein from starved cells"
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "value": "Dps"
                }
            }
        ]
    },
    "return_type": "polymer_entity"
}


# -------------------------------
# Query RCSB
# -------------------------------
response = requests.post(url, json=query)
if response.status_code != 200:
    raise RuntimeError(f"RCSB query failed: {response.status_code}")

data = response.json()
pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
print(f"Initial candidate hits: {len(pdb_ids)}")

# -------------------------------
# Filter by keyword & select best per species
# -------------------------------
best_per_species = {}

for pdb_id in pdb_ids:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    r = requests.get(entry_url)
    if r.status_code != 200:
        continue

    entry_data = r.json()
    text_blob = str(entry_data).lower()
    if not any(k in text_blob for k in keywords):
        continue  # skip if no keyword match

    # Get species
    sources = entry_data.get("rcsb_entity_source_organism", [])
    if not sources:
        continue

    # Use the first scientific name available
    sci_name = sources[0].get("scientific_name", "")
    if not sci_name:
        continue

    # Normalize species to "Genus species" only
    species_words = sci_name.split()
    if len(species_words) < 2:
        continue
    species = " ".join(species_words[:2])

    # Get best resolution
    res = entry_data.get("rcsb_entry_info", {}).get("resolution_combined")
    if res is None:
        continue

    # Keep best resolution per species
    if species not in best_per_species or res < best_per_species[species][1]:
        best_per_species[species] = (pdb_id, res)

# -------------------------------
# Output
# -------------------------------
if not best_per_species:
    raise ValueError("No matching structures found")

print(f"Selected one best-resolution structure per species: {len(best_per_species)}")
for species in sorted(best_per_species):
    pdb_id, res = best_per_species[species]
    print(f"{species}: {pdb_id} ({res} Å)")

with open("data/pdb_ids.txt", "w") as f:
    for species in sorted(best_per_species):
        pdb_id, _ = best_per_species[species]
        f.write(pdb_id + "\n")

print(f"Saved {len(best_per_species)} PDB IDs")
