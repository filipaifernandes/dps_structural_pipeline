import os
import requests
import yaml
import time

# =========================================================
# LOAD CONFIG
# =========================================================

with open("config.yaml") as f:
    config = yaml.safe_load(f)

INTERPRO_IDS = config["query"]["interpro_ids"]

# =========================================================
# OUTPUT DIR
# =========================================================

os.makedirs("data", exist_ok=True)

# =========================================================
# KEEP BEST PDB PER SPECIES (GLOBAL ACROSS ALL IPRs)
# =========================================================

best_per_species = {}

# =========================================================
# MAIN LOOP OVER INTERPRO IDS
# =========================================================

for INTERPRO_ID in INTERPRO_IDS:

    print(f"\nFetching InterPro structures for {INTERPRO_ID}\n")

    URL = (
        f"https://www.ebi.ac.uk/interpro/api/"
        f"structure/pdb/entry/interpro/{INTERPRO_ID}/"
    )

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
            pdb_id = metadata.get("accession")

            if not pdb_id:
                continue

            pdb_id = pdb_id.lower()

            try:

                # =====================================================
                # FETCH ENTRY DATA (RCSB)
                # =====================================================

                rcsb_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                rcsb = requests.get(rcsb_url, timeout=30)

                if rcsb.status_code != 200:
                    continue

                rcsb_data = rcsb.json()

                # =====================================================
                # GET RESOLUTION
                # =====================================================

                resolution = 999.0

                if "rcsb_entry_info" in rcsb_data:

                    resolutions = rcsb_data["rcsb_entry_info"].get(
                        "resolution_combined",
                        []
                    )

                    if resolutions:
                        resolution = min(resolutions)

                # =====================================================
                # GET SPECIES
                # =====================================================

                species = "Unknown"

                polymer_entities = (
                    rcsb_data.get(
                        "rcsb_entry_container_identifiers",
                        {}
                    ).get(
                        "polymer_entity_ids",
                        []
                    )
                )

                if polymer_entities:

                    entity_id = polymer_entities[0]

                    entity_url = (
                        f"https://data.rcsb.org/rest/v1/core/polymer_entity/"
                        f"{pdb_id}/{entity_id}"
                    )

                    entity_response = requests.get(entity_url, timeout=30)

                    if entity_response.status_code == 200:

                        entity_data = entity_response.json()

                        organisms = entity_data.get(
                            "rcsb_entity_source_organism",
                            []
                        )

                        if organisms:
                            species = organisms[0].get(
                                "scientific_name",
                                "Unknown"
                            )

                print(f"{pdb_id} | {species} | {resolution} Å")

                # =====================================================
                # KEEP BEST STRUCTURE PER SPECIES
                # =====================================================

                if species not in best_per_species:

                    best_per_species[species] = {
                        "pdb": pdb_id,
                        "resolution": resolution
                    }

                else:

                    current_best = best_per_species[species]["resolution"]

                    if resolution < current_best:

                        best_per_species[species] = {
                            "pdb": pdb_id,
                            "resolution": resolution
                        }

                time.sleep(0.1)

            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")

        next_url = data.get("next")
        page += 1

# =========================================================
# SAVE OUTPUT
# =========================================================

selected = sorted([
    v["pdb"]
    for v in best_per_species.values()
])

with open("data/pdb_ids.txt", "w") as f:
    for pdb in selected:
        f.write(pdb + "\n")

print(f"\nSpecies retained: {len(selected)}")
print("Saved -> data/pdb_ids.txt")
print("\nDone.\n")
