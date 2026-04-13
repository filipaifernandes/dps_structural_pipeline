# scripts/query_rcsb.py

import requests
import os
import time
import yaml
import sys

print("Script started", flush=True)

# -------------------------------
# Setup
# -------------------------------
os.makedirs("data", exist_ok=True)

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

# -------------------------------
# Load config
# -------------------------------
try:
    CONFIG_PATH = os.path.join(os.path.dirname(__file__), "..", "config.yaml")
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    print("Config loaded", flush=True)
except Exception as e:
    print(f"Failed to load config.yaml: {e}", flush=True)
    sys.exit(1)

keywords = config.get("query", {}).get("keywords", ["dps"])
batch_size = config.get("batch_size", 50)

if not keywords:
    print("No keywords found in config.yaml", flush=True)
    sys.exit(1)

search_value = " ".join(keywords)
print(f"🔎 Using search query: {search_value}", flush=True)

# -------------------------------
# Step 1 — Search RCSB
# -------------------------------
print("📡 Sending search request...", flush=True)

query = {
    "query": {
        "type": "terminal",
        "service": "full_text",
        "parameters": {
            "value": search_value
        }
    },
    "return_type": "entry",
    "request_options": {
        "return_all_hits": True
    }
}

try:
    response = requests.post(SEARCH_URL, json=query)
except Exception as e:
    print(f"Request failed: {e}", flush=True)
    sys.exit(1)

print(f"STATUS: {response.status_code}", flush=True)

# --- DEBUG RAW RESPONSE ---
print("RAW RESPONSE (first 500 chars):", flush=True)
print(response.text[:500], flush=True)

# --- Parse JSON safely ---
try:
    data = response.json()
except Exception as e:
    print("Failed to parse JSON:", e, flush=True)
    sys.exit(1)

print("🔑 Keys in response:", list(data.keys()), flush=True)

# --- Validate structure ---
result_set = data.get("result_set")

if result_set is None:
    print("No 'result_set' in response!", flush=True)
    print("FULL RESPONSE:", data, flush=True)
    sys.exit(1)

entry_ids = [r.get("identifier") for r in result_set if "identifier" in r]

print(f"Initial entry hits: {len(entry_ids)}", flush=True)

if not entry_ids:
    print("No results from search", flush=True)
    sys.exit(1)

# -------------------------------
# Step 2 — Batch fetch via GraphQL
# -------------------------------
def chunk_list(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i+size]

def fetch_batch(entry_batch, retries=3):
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

    for attempt in range(retries):
        try:
            response = requests.post(
                GRAPHQL_URL,
                json={"query": query, "variables": {"ids": entry_batch}}
            )

            if response.status_code == 200:
                return response.json().get("data", {}).get("entries", [])

        except Exception as e:
            print(f"Batch request error: {e}", flush=True)

        time.sleep(1)

    print(f"Failed batch after retries: {entry_batch[:3]}...", flush=True)
    return []

# -------------------------------
# Step 3 — Process data
# -------------------------------
print("Processing entries...", flush=True)

best_per_species = {}

for i, batch in enumerate(chunk_list(entry_ids, batch_size), start=1):
    print(f"Processing batch {i}", flush=True)

    entries = fetch_batch(batch)

    for entry in entries:
        try:
            pdb_id = entry.get("rcsb_id")

            res_list = entry.get("rcsb_entry_info", {}).get("resolution_combined", [])
            if not res_list:
                continue

            resolution = min(res_list)

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

                if (
                    species not in best_per_species
                    or resolution < best_per_species[species][1]
                ):
                    best_per_species[species] = (pdb_id, resolution)

        except Exception as e:
            print(f"Error processing entry: {e}", flush=True)
            continue

    time.sleep(0.1)

# -------------------------------
# Step 4 — Output
# -------------------------------
if not best_per_species:
    print("No matching structures found", flush=True)
    sys.exit(1)

print(f"\nSelected structures: {len(best_per_species)}\n", flush=True)

for species in sorted(best_per_species):
    pdb_id, res = best_per_species[species]
    print(f"{species}: {pdb_id} ({res} Å)", flush=True)

output_file = "data/pdb_ids.txt"

with open(output_file, "w") as f:
    for species in sorted(best_per_species):
        pdb_id, _ = best_per_species[species]
        f.write(pdb_id + "\n")

print(f"\nSaved {len(best_per_species)} PDB IDs → {output_file}", flush=True)
