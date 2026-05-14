# scripts/query_rcsb.py

import os
import requests
import yaml
import time

# =========================================================
# LOAD CONFIG
# =========================================================

with open("config.yaml") as f:
    config = yaml.safe_load(f)

INTERPRO_ID = config["query"]["interpro_id"]

# =========================================================
# OUTPUT DIR
# =========================================================

os.makedirs("data", exist_ok=True)

# =========================================================
# API URL
# =========================================================

URL = (
    f"https://www.ebi.ac.uk/interpro/api/"
    f"structure/pdb/entry/interpro/{INTERPRO_ID}/"
)

print(f"\nFetching InterPro structures for {INTERPRO_ID}\n")

# =========================================================
# FETCH ALL PAGES
# =========================================================

pdb_ids = set()

next_url = URL

page = 1

while next_url:

    print(f"Fetching page {page} ...")

    r = requests.get(next_url, timeout=60)

    r.raise_for_status()

    data = r.json()

    results = data.get("results", [])

    print(f"  Results on page: {len(results)}")

    for item in results:

        metadata = item.get("metadata", {})

        # THIS is the important fix
        accession = metadata.get("accession")

        if accession:

            pdb_ids.add(accession.lower())

    next_url = data.get("next")

    page += 1

    time.sleep(0.2)

# =========================================================
# SAVE
# =========================================================

with open("data/pdb_ids.txt", "w") as f:

    for pdb in sorted(pdb_ids):

        f.write(pdb + "\n")

print(f"\nUnique PDB IDs collected: {len(pdb_ids)}")

print("Saved -> data/pdb_ids.txt")

print("\nDone.\n")
