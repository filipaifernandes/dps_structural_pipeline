import os
import requests
import yaml

os.makedirs("data", exist_ok=True)

with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

queries = config["queries"]

nodes = []
for q in queries:
    nodes.append({
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "struct.title",
            "operator": "contains_words",
            "value": q
        }
    })

query_json = {
    "query": {
        "type": "group",
        "logical_operator": "or",
        "nodes": nodes
    },
    "return_type": "entry"
}

url = "https://search.rcsb.org/rcsbsearch/v2/query"
response = requests.post(url, json=query_json)
data = response.json()

pdb_ids = [r["identifier"] for r in data["result_set"]]

with open("data/pdb_ids.txt", "w") as f:
    for pdb in pdb_ids:
        f.write(pdb + "\n")

print(f"Saved {len(pdb_ids)} PDB IDs to data/pdb_ids.txt")
