import requests
import os
import yaml
import sys

# Load config
config_file = sys.argv[1]

with open(config_file) as f:
    config = yaml.safe_load(f)

organism = config["query"]["organism"]
keywords = [k.lower() for k in config["query"]["keywords"]]

url = "https://search.rcsb.org/rcsbsearch/v2/query"

# Structured query: organism + protein
query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                    "operator": "contains_words",
                    "value": organism
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": "Protein"
                }
            }
        ]
    },
    "return_type": "entry",
    "request_options": {
        "return_all_hits": True
    }
}

response = requests.post(url, json=query)

if response.status_code != 200:
    raise RuntimeError(f"RCSB query failed: {response.status_code}")

data = response.json()
pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]

print(f"Initial {organism} protein hits: {len(pdb_ids)}")

filtered_ids = set()

# Filter by keywords
for pdb_id in pdb_ids:
    entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    r = requests.get(entry_url)

    if r.status_code != 200:
        continue

    entry_data = r.json()
    text_blob = str(entry_data).lower()

    if any(keyword in text_blob for keyword in keywords):
        filtered_ids.add(pdb_id)

print(f"Filtered structures: {len(filtered_ids)}")

if not filtered_ids:
    raise ValueError("No matching structures found")

os.makedirs("data", exist_ok=True)

with open("data/pdb_ids.txt", "w") as f:
    for pdb_id in sorted(filtered_ids):
        f.write(pdb_id + "\n")

print(f"Saved {len(filtered_ids)} PDB IDs")
