# scripts/query_rcsb.py
import requests
import os
import sys
import time

os.makedirs("data", exist_ok=True)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

# -------------------------------
# Step 1 — Search for DPS (like website)
# -------------------------------
query = {
    "query": {
        "type": "terminal",
        "service": "full_text",
        "parameters": {
            "value": "dps"   # <- fixed, as you want
        }
    },
    "return_type": "entry",
    "request_options": {
        "return_all_hits": True
    }
}

response = requests.post(SEARCH_URL, json=query)
if response.status_code != 200:
    raise RuntimeError(f"Search failed: {response.status_code}")

data = response.json()
entry_ids = [r["identifier"] for r in data.get("result_set", [])]

print(f"Initial entry hits: {len(entry_ids)}")

if not entry_ids:
    raise ValueError("No results from search")

# -------------------------------
# Step 2 — Batch fetch via GraphQL
# -------------------------------
def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i+size]

def fetch_batch(entry_batch):
    query = """
    query ($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        rcsb_entry_info {
          resolution_combined
        }
        polymer_entities {
          rcsb_entity_source_organism {
            scientific_name
          }
        }
      }
    }
    """

    response = requests.post(
        GRAPHQL_URL,
        json={"query": query, "variables": {"ids": entry_batch}}
    )

    if response.status_code != 200:
        return []

    return response.json().get("data", {}).get("entries", [])

# -------------------------------
# Step 3 — Process data
# -------------------------------
best_per_species = {}

batch_size = 50

for batch in chunk_list(entry_ids, batch_size):
    entries = fetch_batch(batch)

    for entry in entries:
        try:
            pdb_id = entry.get("rcsb_id")

            # --- resolution ---
            res_list = entry.get("rcsb_entry_info", {}).get("resolution_combined", [])
            if not res_list:
                continue

            resolution = min(res_list)

            # --- species extraction ---
            for entity in entry.get("polymer_entities", []):
                orgs = entity.get("rcsb_entity_source_organism", [])
                if not orgs:
                    continue

                sci_name = orgs[0].get("scientific_name", "")
                if not sci_name:
                    continue

                parts = sci_name.split()
                if len(parts) < 2:
                    continue

                species = " ".join(parts[:2])

                # --- keep best ---
                if species not in best_per_species or resolution < best_per_species[species][1]:
                    best_per_species[species] = (pdb_id, resolution)

        except Exception:
            continue

    time.sleep(0.1)

# -------------------------------
# Step 4 — Output
# -------------------------------
if not best_per_species:
    raise ValueError("No matching structures found")

print(f"\nSelected structures: {len(best_per_species)}\n")

for species in sorted(best_per_species):
    pdb_id, res = best_per_species[species]
    print(f"{species}: {pdb_id} ({res} Å)")

# save file
with open("data/pdb_ids.txt", "w") as f:
    for species in sorted(best_per_species):
        pdb_id, _ = best_per_species[species]
        f.write(pdb_id + "\n")

print(f"\nSaved {len(best_per_species)} PDB IDs → data/pdb_ids.txt")
