import requests
import sys
import itertools
import yaml

config_file = sys.argv[1]

with open(config_file) as f:
    config = yaml.safe_load(f)

queries = config.get("queries", [])

if not queries:
    sys.exit("No queries found in config.yaml")

# Build a search string by joining query terms
search_string = " AND ".join(queries)  

# RCSB search API
url = "https://search.rcsb.org/rcsbsearch/v2/query"

payload = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "operator": "contains_words",
            "value": search_string
        }
    },
    "return_type": "entry"
}

r = requests.post(url, json=payload)
r.raise_for_status()
results = r.json()["result_set"]

pdb_ids = [entry["identifier"] for entry in results]

# Print PDB IDs to stdout (one per line)
for pdb in sorted(pdb_ids):
    print(pdb)
