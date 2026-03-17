import requests
import os
import sys
import yaml

config_file = sys.argv[1]

with open(config_file) as f:
    config = yaml.safe_load(f)

queries = config.get("queries", [])

if not queries:
    raise ValueError("No 'queries' found in config.yaml")

url = "https://search.rcsb.org/rcsbsearch/v2/query"
all_pdb_ids = set()

for q in queries:
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": q
            }
        },
        "return_type": "entry"
    }

    response = requests.post(url, json=payload)

    if response.status_code != 200:
        print(f"Warning: query '{q}' failed with status {response.status_code}")
        continue

    data = response.json()

    for entry in data.get("result_set", []):
        all_pdb_ids.add(entry["identifier"])

if not all_pdb_ids:
    raise ValueError("No PDB structures found for the configured queries.")

os.makedirs("data", exist_ok=True)

with open("data/pdb_ids.txt", "w") as f:
    for pdb_id in sorted(all_pdb_ids):
        f.write(pdb_id + "\n")

print(f"Saved {len(all_pdb_ids)} PDB IDs to data/pdb_ids.txt")
