import requests
import sys
import yaml

# Load config
config_file = sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"
with open(config_file) as f:
    config = yaml.safe_load(f)

queries = config.get("queries", [])
if not queries:
    sys.exit("No queries found in config.yaml")

pdb_ids = set()

# Query each term separately, strip whitespace
for q in queries:
    q = q.strip()
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "operator": "exact_match",
                "value": q
            }
        },
        "return_type": "entry"
    }

    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=payload)
    if r.status_code != 200:
        print(f"Warning: query '{q}' failed with status {r.status_code}")
        continue

    results = r.json().get("result_set", [])
    for entry in results:
        pdb_ids.add(entry["identifier"])

# Output sorted PDB IDs
for pdb in sorted(pdb_ids):
    print(pdb)
