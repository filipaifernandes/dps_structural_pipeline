import os
import requests
import yaml
import time
import json
from collections import Counter

# =========================================================
# LOAD CONFIG
# =========================================================

with open("config.yaml") as f:
    config = yaml.safe_load(f)

INTERPRO_IDS = config["query"]["interpro_ids"]
print(f"Loaded {len(INTERPRO_IDS)} InterPro IDs: {INTERPRO_IDS}", flush=True)

os.makedirs("data", exist_ok=True)

best_per_species = {}

# =========================================================
# HELPER: fetch one page with retry
# =========================================================

def fetch_page(url, retries=5):
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 408:
                print(f"  Server timeout (408), waiting 61s ...")
                time.sleep(61)
                continue
            if r.status_code == 204:
                return {"results": [], "next": None, "count": 0}
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f"  Request error (attempt {attempt+1}): {e}")
            time.sleep(5)
    return None

# =========================================================
# GET ORGANISM VIA PROTEIN ENDPOINT
# Confirmed response shape:
#   results[0].metadata.source_organism.scientificName
# =========================================================

def get_organism_for_pdb(pdb_id, ipr_id):
    url = (
        f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/"
        f"entry/interpro/{ipr_id}/structure/pdb/{pdb_id.lower()}/"
        f"?format=json&page_size=1"
    )
    data = fetch_page(url)
    if not data:
        return "Unknown"
    results = data.get("results", [])
    if not results:
        return "Unknown"

    source_org = results[0].get("metadata", {}).get("source_organism", {})

    # Confirmed field name: scientificName (camelCase)
    name = source_org.get("scientificName", "")
    if name:
        return name

    # Fallback spellings just in case
    name = source_org.get("scientific_name", "") or source_org.get("fullName", "")
    return name or "Unknown"


def extract_organism_from_structure_item(item):
    """
    Try to get organism from the structure endpoint response.
    These items have a different shape to the protein endpoint —
    organism is often absent here entirely, hence the fallback.
    """
    metadata = item.get("metadata", {})

    for org_block in [
        metadata.get("source_organism", {}),
        metadata.get("taxonomy", {}),
    ]:
        if not isinstance(org_block, dict):
            continue
        for key in ("scientificName", "scientific_name", "fullName"):
            name = org_block.get(key, "")
            if name and name.lower() != "unknown":
                return name

    proteins = item.get("proteins", [])
    if proteins:
        org_block = proteins[0].get("organism", {})
        if isinstance(org_block, dict):
            for key in ("scientificName", "scientific_name", "fullName"):
                name = org_block.get(key, "")
                if name and name.lower() != "unknown":
                    return name

    return None  # not found — use protein endpoint fallback


# =========================================================
# MAIN LOOP
# =========================================================

for INTERPRO_ID in INTERPRO_IDS:

    print(f"\n{'='*60}", flush=True)
    print(f"Processing: {INTERPRO_ID}", flush=True)
    print(f"{'='*60}\n", flush=True)

    next_url = (
        f"https://www.ebi.ac.uk/interpro/api/"
        f"structure/pdb/entry/interpro/{INTERPRO_ID}/"
        f"?page_size=200"
    )

    page = 1
    page_hits = 0

    while next_url:

        print(f"  Page {page} ...", flush=True)
        data = fetch_page(next_url)

        if data is None:
            print(f"  Failed to fetch page {page}, skipping {INTERPRO_ID}", flush=True)
            break

        results = data.get("results", [])

        if page == 1:
            print(f"  Total structures reported: {data.get('count', '?')}\n", flush=True)

        seen_on_page = set()

        for item in results:
            metadata = item.get("metadata", {})
            pdb_id = metadata.get("accession", "").upper()

            if not pdb_id or pdb_id in seen_on_page:
                continue
            seen_on_page.add(pdb_id)

            # Resolution
            resolution = metadata.get("resolution")
            try:
                resolution = float(resolution) if resolution is not None else 999.0
            except (ValueError, TypeError):
                resolution = 999.0

            # Organism — try structure item first, then protein endpoint
            organism = extract_organism_from_structure_item(item)
            if organism is None:
                organism = get_organism_for_pdb(pdb_id, INTERPRO_ID)
                time.sleep(0.2)

            # Normalise to binomial (Genus species)
            parts = organism.strip().split()
            species = " ".join(parts[:2]) if len(parts) >= 2 else organism.strip()

            if not species or species.lower() == "unknown":
                # Don't skip — use PDB ID as key so structure is retained
                # itol_labels will show PDB ID instead of species name
                species = pdb_id.upper()
                print(f"    {pdb_id} | organism unknown — keeping with PDB ID as label", flush=True)
            else:
                print(f"    {pdb_id} | {species} | {resolution} Å", flush=True)
            page_hits += 1

            if (
                species not in best_per_species
                or resolution < best_per_species[species]["resolution"]
            ):
                best_per_species[species] = {
                    "pdb":        pdb_id.lower(),
                    "resolution": resolution,
                    "interpro":   INTERPRO_ID,
                    "organism":   species,
                }

        next_url = data.get("next")
        page += 1
        time.sleep(0.3)

    print(f"\n  {INTERPRO_ID} done — {page_hits} unique structures processed", flush=True)

# =========================================================
# SAVE pdb_ids.txt
# =========================================================

selected = sorted([v["pdb"] for v in best_per_species.values()])

with open("data/pdb_ids.txt", "w") as f:
    for pdb in selected:
        f.write(pdb + "\n")

# =========================================================
# SAVE species_map.json — reused by itol_labels.py
# =========================================================

species_map = {v["pdb"]: v["organism"] for v in best_per_species.values()}

with open("data/species_map.json", "w") as f:
    json.dump(species_map, f, indent=2)

print(f"\n{'='*60}")
print(f"Total species retained : {len(selected)}")
print(f"Saved                  -> data/pdb_ids.txt")
print(f"Saved                  -> data/species_map.json")
print(f"{'='*60}\n")

ipr_counts = Counter(v["interpro"] for v in best_per_species.values())
print("Species per IPR entry:")
for ipr, count in sorted(ipr_counts.items()):
    print(f"  {ipr}: {count} species")

print("\nDone.\n")
