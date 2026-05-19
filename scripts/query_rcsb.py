import os
import requests
import yaml
import time
from collections import Counter

# =========================================================
# LOAD CONFIG
# =========================================================

with open("config.yaml") as f:
    config = yaml.safe_load(f)

INTERPRO_IDS = config["query"]["interpro_ids"]

os.makedirs("data", exist_ok=True)

best_per_species = {}

# =========================================================
# HELPER: fetch one page with retry on 408 timeout
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
                # No content — valid but empty
                return {"results": [], "next": None, "count": 0}
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f"  Request error (attempt {attempt+1}): {e}")
            time.sleep(5)
    return None

# =========================================================
# GET ENTRY TYPE FROM INTERPRO
# =========================================================

def get_entry_type(ipr_id):
    """
    Returns the entry type string e.g. 'family', 'domain',
    'homologous_superfamily', or None if unreachable.
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{ipr_id}/"
    data = fetch_page(url)
    if not data:
        return None, "unknown", "unknown"

    metadata = data.get("metadata", {})
    entry_type = metadata.get("type", "unknown")
    name = metadata.get("name", {})
    name = name.get("name", "unknown") if isinstance(name, dict) else str(name)
    counters = metadata.get("counters", {})
    structures = counters.get("structures", 0)
    return entry_type, name, structures

# =========================================================
# FETCH STRUCTURES — tries multiple endpoint patterns
# =========================================================

def fetch_structures_for_ipr(ipr_id, entry_type):
    """
    InterPro API uses different URL patterns depending on entry type.
    We try all known patterns and use whichever returns results.

    Known working patterns:
      /structure/pdb/entry/interpro/{ID}/          <- families, domains
      /structure/pdb/protein/UniProt/entry/InterPro/{ID}/  <- superfamilies (indirect)
    
    Returns list of (pdb_id, resolution, organism) tuples.
    """

    candidate_urls = [
        # Standard — works for most families and domains
        f"https://www.ebi.ac.uk/interpro/api/structure/pdb/entry/interpro/{ipr_id}/?page_size=200",
        # Alternate capitalisation InterPro vs interpro
        f"https://www.ebi.ac.uk/interpro/api/structure/pdb/entry/InterPro/{ipr_id}/?page_size=200",
    ]

    for start_url in candidate_urls:

        print(f"  Trying: {start_url}", flush=True)

        next_url = start_url
        page = 1
        hits = []

        while next_url:
            print(f"    Page {page} ...", flush=True)
            data = fetch_page(next_url)

            if data is None:
                print(f"    Failed — moving to next URL pattern")
                break

            results = data.get("results", [])
            count   = data.get("count", 0)

            if page == 1:
                print(f"    Total structures reported: {count}")

            if not results and page == 1:
                print(f"    No results from this endpoint pattern")
                break

            for item in results:
                metadata = item.get("metadata", {})
                pdb_id = metadata.get("accession", "").upper()
                if not pdb_id:
                    continue

                # Resolution
                resolution = metadata.get("resolution")
                try:
                    resolution = float(resolution) if resolution is not None else 999.0
                except (ValueError, TypeError):
                    resolution = 999.0

                # Organism — check multiple locations in response
                organism = "Unknown"

                source_org = metadata.get("source_organism", {})
                if isinstance(source_org, dict):
                    organism = source_org.get("scientific_name", "Unknown") or "Unknown"

                if organism == "Unknown":
                    proteins = item.get("proteins", [])
                    if proteins:
                        org_block = proteins[0].get("organism", {})
                        if isinstance(org_block, dict):
                            organism = org_block.get("scientific_name", "Unknown") or "Unknown"

                if organism == "Unknown":
                    # last resort: check 'taxonomy' block
                    tax = metadata.get("taxonomy", {})
                    if isinstance(tax, dict):
                        organism = tax.get("scientific_name", "Unknown") or "Unknown"

                hits.append((pdb_id, resolution, organism))

            next_url = data.get("next")
            page += 1
            time.sleep(0.3)

        if hits:
            print(f"  -> Got {len(hits)} structures from this endpoint", flush=True)
            return hits

        print(f"  -> 0 structures, trying next pattern...", flush=True)

    # If standard endpoints both failed, try going via proteins
    # /protein/UniProt/entry/interpro/{ID}/ -> collect accessions -> SIFTS per accession
    # (only triggered as last resort, limited to reviewed proteins to keep it fast)
    print(f"  Standard endpoints returned nothing — trying protein->structure fallback ...", flush=True)
    return fetch_via_proteins(ipr_id)


def fetch_via_proteins(ipr_id):
    """
    Fallback: get reviewed UniProt accessions for this IPR,
    then query SIFTS best_structures for each.
    Limited to reviewed (SwissProt) to keep it manageable.
    """
    protein_url = (
        f"https://www.ebi.ac.uk/interpro/api/protein/reviewed/"
        f"entry/interpro/{ipr_id}/?page_size=200"
    )

    accessions = []
    next_url = protein_url
    page = 1

    print(f"  Collecting reviewed UniProt accessions ...", flush=True)

    while next_url:
        data = fetch_page(next_url)
        if not data:
            break
        for item in data.get("results", []):
            acc = item.get("metadata", {}).get("accession")
            if acc:
                accessions.append(acc)
        next_url = data.get("next")
        page += 1
        time.sleep(0.3)

    print(f"  Found {len(accessions)} reviewed accessions", flush=True)

    if not accessions:
        return []

    hits = []
    SIFTS = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"

    for i, acc in enumerate(accessions, 1):
        if i % 20 == 0:
            print(f"  SIFTS {i}/{len(accessions)} ...", flush=True)
        try:
            r = requests.get(f"{SIFTS}/{acc}", timeout=15)
            if r.status_code == 404:
                continue
            r.raise_for_status()
            data = r.json()
            for entry in data.get(acc, []):
                pdb_id    = entry.get("pdb_id", "").upper()
                resolution = entry.get("resolution")
                organism  = entry.get("organism_scientific_name", "Unknown")
                if pdb_id and resolution is not None:
                    hits.append((pdb_id, float(resolution), organism))
        except Exception as e:
            print(f"  SIFTS error for {acc}: {e}", flush=True)
        time.sleep(0.05)

    print(f"  Fallback found {len(hits)} structure hits", flush=True)
    return hits

# =========================================================
# MAIN LOOP
# =========================================================

for INTERPRO_ID in INTERPRO_IDS:

    print(f"\n{'='*60}")
    print(f"Processing: {INTERPRO_ID}")
    print(f"{'='*60}")

    entry_type, entry_name, structure_count = get_entry_type(INTERPRO_ID)

    if entry_type is None:
        print(f"  ERROR: Could not reach InterPro for {INTERPRO_ID} — skipping")
        continue

    print(f"  Name       : {entry_name}")
    print(f"  Type       : {entry_type}")
    print(f"  Structures : {structure_count} (as reported by InterPro)")

    hits = fetch_structures_for_ipr(INTERPRO_ID, entry_type)

    for (pdb_id, resolution, organism) in hits:
        parts = organism.strip().split()
        species = " ".join(parts[:2]) if len(parts) >= 2 else organism.strip() or "Unknown"

        if species == "Unknown":
            continue

        if (
            species not in best_per_species
            or resolution < best_per_species[species]["resolution"]
        ):
            best_per_species[species] = {
                "pdb": pdb_id.lower(),
                "resolution": resolution,
                "interpro": INTERPRO_ID,
            }

    print(f"  {INTERPRO_ID} done — {len(hits)} hits processed", flush=True)

# =========================================================
# SAVE OUTPUT
# =========================================================

selected = sorted([v["pdb"] for v in best_per_species.values()])

with open("data/pdb_ids.txt", "w") as f:
    for pdb in selected:
        f.write(pdb + "\n")

print(f"\n{'='*60}")
print(f"Total species retained : {len(selected)}")
print(f"Saved                  -> data/pdb_ids.txt")
print(f"{'='*60}\n")

ipr_counts = Counter(v["interpro"] for v in best_per_species.values())
print("Species per IPR entry:")
for ipr, count in sorted(ipr_counts.items()):
    print(f"  {ipr}: {count} species")

print("\nDone.\n")
