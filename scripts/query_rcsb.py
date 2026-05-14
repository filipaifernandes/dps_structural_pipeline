# scripts/query_rcsb.py
#
# DPS family structure retrieval — optimised pipeline
# Strategy: InterPro IPR002177 (?has_structure=true) → PDBe batch API → best resolution per species
#
#   1. Query InterPro for UniProt accessions in IPR002177 that HAVE PDB structures
#      (skips the ~22k accessions that have no structure — they're useless downstream)
#   2. Batch-POST all accessions to PDBe graph-api in chunks of 1000
#      (one request per chunk instead of one per accession)
#   3. Select best-resolution structure per species
#   4. Write data/pdb_ids.txt

import requests
import time
import os
import sys
import yaml

print("Script started", flush=True)
os.makedirs("data", exist_ok=True)

# -----------------------------------------------------------------------
# Load config
# -----------------------------------------------------------------------
try:
    CONFIG_PATH = os.path.join(os.path.dirname(__file__), "..", "config.yaml")
    with open(CONFIG_PATH) as f:
        config = yaml.safe_load(f)
    print("Config loaded", flush=True)
except Exception as e:
    print(f"Failed to load config.yaml: {e}", flush=True)
    sys.exit(1)

INTERPRO_ID = config.get("query", {}).get("interpro_id", "IPR002177")
print(f"Using InterPro family: {INTERPRO_ID}", flush=True)

# -----------------------------------------------------------------------
# API endpoints
# -----------------------------------------------------------------------
INTERPRO_API   = "https://www.ebi.ac.uk/interpro/api"
PDBE_BATCH_API = "https://www.ebi.ac.uk/pdbe/graph-api/uniprot/uniprot_mappings"

# -----------------------------------------------------------------------
# Step 1 — InterPro: only accessions that have at least one PDB structure
# -----------------------------------------------------------------------
def fetch_interpro_members_with_structures(interpro_id):
    accessions = []
    url = (
        f"{INTERPRO_API}/protein/UniProt/entry/InterPro/{interpro_id}/"
        f"?format=json&page_size=200&has_structure=true"
    )

    page = 0
    while url:
        page += 1
        print(f"  Fetching InterPro page {page} ...", flush=True)
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"  InterPro request failed: {e}", flush=True)
            break

        for entry in data.get("results", []):
            acc = entry.get("metadata", {}).get("accession")
            if acc:
                accessions.append(acc)

        url = data.get("next")
        if url:
            time.sleep(0.3)

    return accessions


print("\n[Step 1] Querying InterPro for DPS members WITH structures ...", flush=True)
uniprot_accessions = fetch_interpro_members_with_structures(INTERPRO_ID)
print(f"  UniProt accessions with PDB structures: {len(uniprot_accessions)}", flush=True)

if not uniprot_accessions:
    print("No accessions returned. Exiting.", flush=True)
    sys.exit(1)

# -----------------------------------------------------------------------
# Step 2 — PDBe batch API: get all PDB mappings in chunks of 1000
# -----------------------------------------------------------------------
def fetch_pdbe_batch(accession_chunk, retries=3):
    """
    POST a list of UniProt accessions (up to 1000) to PDBe graph-api.
    Returns dict: { uniprot_acc: [ {pdb_id, resolution, organism, ...} ] }
    """
    for attempt in range(retries):
        try:
            r = requests.post(
                PDBE_BATCH_API,
                json={"accessions": accession_chunk},
                timeout=60,
            )
            if r.status_code == 404:
                return {}
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f"  PDBe batch error (attempt {attempt+1}): {e}", flush=True)
            time.sleep(2)
    return {}


def parse_pdbe_response(data):
    """
    Flatten PDBe graph-api response into a list of structure dicts.
    Response shape:
      { uniprot_acc: { mappings: [ { pdb_id, resolution, organism, ... } ] } }
    """
    results = []
    for uniprot_acc, content in data.items():
        # graph-api wraps in a list for each accession
        entries = content if isinstance(content, list) else [content]
        for entry in entries:
            mappings = entry.get("mappings", [])
            for m in mappings:
                pdb_id     = m.get("pdb_id", "").upper()
                resolution = m.get("resolution")
                organism   = (
                    m.get("organism_scientific_name")
                    or m.get("organism", {}).get("scientific_name", "")
                    if isinstance(m.get("organism"), dict)
                    else m.get("organism", "")
                )
                if pdb_id and resolution is not None:
                    try:
                        results.append({
                            "pdb_id":     pdb_id,
                            "resolution": float(resolution),
                            "organism":   str(organism).strip(),
                            "uniprot":    uniprot_acc,
                        })
                    except (ValueError, TypeError):
                        continue
    return results


CHUNK_SIZE = 1000
all_structures = []
chunks = [
    uniprot_accessions[i : i + CHUNK_SIZE]
    for i in range(0, len(uniprot_accessions), CHUNK_SIZE)
]

print(f"\n[Step 2] Fetching PDB mappings in {len(chunks)} batch request(s) ...", flush=True)

for i, chunk in enumerate(chunks, start=1):
    print(f"  Batch {i}/{len(chunks)} ({len(chunk)} accessions) ...", flush=True)
    raw = fetch_pdbe_batch(chunk)
    hits = parse_pdbe_response(raw)
    all_structures.extend(hits)
    time.sleep(0.5)     # one small pause between batches — not per accession

print(f"  Total PDB structure-hits: {len(all_structures)}", flush=True)

if not all_structures:
    print("No PDB structures found. Exiting.", flush=True)
    sys.exit(1)

# -----------------------------------------------------------------------
# Step 3 — Best-resolution structure per species
# -----------------------------------------------------------------------
print("\n[Step 3] Selecting best-resolution structure per species ...", flush=True)

best_per_species = {}

for hit in all_structures:
    organism = hit["organism"]
    if not organism:
        continue

    # Normalise to binomial — drop strain/subspecies noise
    parts = organism.split()
    if len(parts) < 2:
        continue
    species = " ".join(parts[:2])

    if (
        species not in best_per_species
        or hit["resolution"] < best_per_species[species]["resolution"]
    ):
        best_per_species[species] = hit

print(f"  Unique species with structures: {len(best_per_species)}", flush=True)

if not best_per_species:
    print("No species could be assigned. Exiting.", flush=True)
    sys.exit(1)

# -----------------------------------------------------------------------
# Step 4 — Report & write output
# -----------------------------------------------------------------------
print("\n[Step 4] Selected structures:\n", flush=True)
for species in sorted(best_per_species):
    info = best_per_species[species]
    print(
        f"  {species:<40}  {info['pdb_id']}  "
        f"{info['resolution']:.2f} A  [{info['uniprot']}]",
        flush=True,
    )

output_file = "data/pdb_ids.txt"
with open(output_file, "w") as f:
    for species in sorted(best_per_species):
        f.write(best_per_species[species]["pdb_id"] + "\n")

print(f"\nSaved {len(best_per_species)} PDB IDs -> {output_file}", flush=True)
