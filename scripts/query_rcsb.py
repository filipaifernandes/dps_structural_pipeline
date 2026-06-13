# scripts/query_rcsb.py
#
# DPS family structure retrieval — InterPro + SIFTS + paralog-aware deduplication
#
# Strategy:
# 1. InterPro structure endpoint -> known PDB IDs + resolution + organism + UniProt
# 2. SIFTS flat file (cached, 30-day expiry) -> two jobs:
#    a. Fill missing UniProt accessions for structures that InterPro didn't annotate
#    b. Find extra PDB entries for our UniProt set that InterPro missed
# 3. PDBe summary API -> resolution + organism for SIFTS extras
# 4. Paralog-aware deduplication:
#    - group by species -> UniProt accession -> keep best resolution
#    - species with multiple UniProt accessions = paralogs, keep one per paralog
#    - unknowns grouped by UniProt (not PDB ID) to correctly deduplicate same-protein variants
# 5. Write pdb_ids.txt + species_map.json + selection_report.tsv

import os
import csv
import gzip
import json
import time
import datetime
import requests
import yaml
from collections import defaultdict

# =========================================================
# LOAD CONFIG
# =========================================================

with open("config.yaml") as f:
    config = yaml.safe_load(f)

INTERPRO_IDS       = config["query"]["interpro_ids"]
SIFTS_MAX_AGE_DAYS = config.get("sifts_max_age_days", 30)

print(f"InterPro IDs     : {INTERPRO_IDS}", flush=True)
print(f"SIFTS max age    : {SIFTS_MAX_AGE_DAYS} days\n", flush=True)

os.makedirs("data", exist_ok=True)

# =========================================================
# HELPERS
# =========================================================

def fetch_page(url, retries=5):
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 408:
                print(f"  408 timeout — waiting 61s ...")
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


def get_organism_and_uniprot(pdb_id, ipr_id):
    """
    InterPro protein endpoint for a specific PDB+IPR combination.
    Returns (organism_name, uniprot_accession).
    Confirmed field: results[0].metadata.source_organism.scientificName
    """
    url = (
        f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/"
        f"entry/interpro/{ipr_id}/structure/pdb/{pdb_id.lower()}/"
        f"?format=json&page_size=1"
    )
    data = fetch_page(url)
    if not data or not data.get("results"):
        return "Unknown", None
    metadata   = data["results"][0].get("metadata", {})
    source_org = metadata.get("source_organism", {})
    uniprot    = metadata.get("accession", None)
    for key in ("scientificName", "scientific_name", "fullName"):
        name = source_org.get(key, "")
        if name and name.lower() != "unknown":
            return name, uniprot
    return "Unknown", uniprot


def extract_organism_from_structure_item(item):
    """Try every known location in InterPro /structure/ response item."""
    metadata = item.get("metadata", {})
    for block in [metadata.get("source_organism", {}), metadata.get("taxonomy", {})]:
        if not isinstance(block, dict):
            continue
        for key in ("scientificName", "scientific_name", "fullName"):
            name = block.get(key, "")
            if name and name.lower() != "unknown":
                return name
    for protein in item.get("proteins", []):
        org = protein.get("organism", {})
        if isinstance(org, dict):
            for key in ("scientificName", "scientific_name", "fullName"):
                name = org.get(key, "")
                if name and name.lower() != "unknown":
                    return name
    return None


def normalise_species(organism):
    """Reduce to binomial (Genus species), dropping strain info."""
    if not organism:
        return "Unknown"
    parts = organism.strip().split()
    return " ".join(parts[:2]) if len(parts) >= 2 else organism.strip()


def sifts_needs_refresh(path, max_age_days):
    if not os.path.exists(path):
        return True
    age = datetime.datetime.now() - datetime.datetime.fromtimestamp(
        os.path.getmtime(path)
    )
    return age.days >= max_age_days


# =========================================================
# STEP 1 — InterPro structure endpoint
#          Fast: only ~130 structures across both IPR IDs
# =========================================================

print("=" * 60, flush=True)
print("STEP 1  InterPro structure endpoint", flush=True)
print("=" * 60, flush=True)

all_structures     = []   # list of dicts
uniprot_accessions = set() # all UniProt accs seen — used for SIFTS lookup
seen_pdbs          = set() # dedup across IPR IDs

for ipr_id in INTERPRO_IDS:
    print(f"\n  {ipr_id}", flush=True)
    next_url = (
        f"https://www.ebi.ac.uk/interpro/api/structure/pdb/"
        f"entry/interpro/{ipr_id}/?page_size=200"
    )
    page = 1
    hits = 0

    while next_url:
        data = fetch_page(next_url)
        if not data:
            break
        if page == 1:
            print(f"  Total structures reported: {data.get('count', '?')}", flush=True)

        for item in data.get("results", []):
            metadata = item.get("metadata", {})
            pdb_id   = metadata.get("accession", "").upper()
            if not pdb_id or pdb_id in seen_pdbs:
                continue
            seen_pdbs.add(pdb_id)

            # Resolution
            try:
                resolution = float(metadata["resolution"])
            except (KeyError, TypeError, ValueError):
                resolution = 999.0

            # Organism — try structure response first, then protein endpoint
            organism = extract_organism_from_structure_item(item)
            uniprot  = None

            # Collect UniProt from entries[] block regardless
            for entry in item.get("entries", []):
                protein_acc = entry.get("protein", "")
                if protein_acc:
                    acc = protein_acc.upper()
                    uniprot_accessions.add(acc)
                    if uniprot is None:
                        uniprot = acc

            if organism is None:
                # Fall back to protein endpoint — also gets us the UniProt
                organism, uni2 = get_organism_and_uniprot(pdb_id, ipr_id)
                if uni2 and uniprot is None:
                    uniprot = uni2
                    uniprot_accessions.add(uni2)
                time.sleep(0.2)

            species = normalise_species(organism)

            print(f"    {pdb_id:<6}  {species:<40}  {resolution:.2f} Å  {uniprot}", flush=True)
            hits += 1

            all_structures.append({
                "pdb_id"    : pdb_id.lower(),
                "resolution": resolution,
                "organism"  : species,
                "uniprot"   : uniprot,
                "source"    : ipr_id,
            })

        next_url = data.get("next")
        page += 1
        time.sleep(0.3)

    print(f"  {ipr_id} done — {hits} structures", flush=True)

print(f"\nStep 1 total: {len(all_structures)} structures, "
      f"{len(uniprot_accessions)} UniProt accessions\n", flush=True)

# =========================================================
# STEP 2 — SIFTS flat file
#          Two jobs in one pass:
#          (a) patch missing UniProt for structures already in all_structures
#          (b) find additional PDB entries for known UniProt accessions
# =========================================================

print("=" * 60, flush=True)
print("STEP 2  SIFTS flat file", flush=True)
print("=" * 60, flush=True)

SIFTS_URL   = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz"
SIFTS_CACHE = "data/pdb_chain_uniprot.csv.gz"

if sifts_needs_refresh(SIFTS_CACHE, SIFTS_MAX_AGE_DAYS):
    reason = "cache expired" if os.path.exists(SIFTS_CACHE) else "no cache"
    print(f"  Downloading SIFTS flat file ({reason}) ...", flush=True)
    r = requests.get(SIFTS_URL, stream=True, timeout=300)
    r.raise_for_status()
    with open(SIFTS_CACHE, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    print(f"  Saved -> {SIFTS_CACHE}", flush=True)
else:
    age = (datetime.datetime.now() -
           datetime.datetime.fromtimestamp(os.path.getmtime(SIFTS_CACHE))).days
    print(f"  Using cached SIFTS file (age: {age} days)", flush=True)

# Build lookup: pdb_id -> UniProt from SIFTS
# (used to patch missing UniProt and find extras)
# Load blacklists from config.yaml — nothing hardcoded here
UNIPROT_BLACKLIST = {
    entry["accession"]
    for entry in config.get("query", {}).get("uniprot_blacklist", [])
}
PDB_ID_BLACKLIST = {
    entry["pdb_id"].lower()
    for entry in config.get("query", {}).get("pdb_blacklist", [])
}

if UNIPROT_BLACKLIST:
    print(f"  UniProt blacklist loaded: {UNIPROT_BLACKLIST}", flush=True)
if PDB_ID_BLACKLIST:
    print(f"  PDB ID blacklist loaded : {PDB_ID_BLACKLIST}", flush=True)
removed_from_set = uniprot_accessions & UNIPROT_BLACKLIST
if removed_from_set:
    print(f"  Purging blacklisted accessions before SIFTS: {removed_from_set}", flush=True)
    uniprot_accessions -= UNIPROT_BLACKLIST

# Also remove any structures already collected whose UniProt is blacklisted
# or whose PDB ID is explicitly blacklisted
all_structures = [
    s for s in all_structures
    if s.get("uniprot") not in UNIPROT_BLACKLIST
    and s["pdb_id"] not in PDB_ID_BLACKLIST
]
print(f"  Structures after blacklist filter: {len(all_structures)}", flush=True)

# =========================================================
# pdb_force — for a given UniProt, always use a specific PDB
#
# Removes all other PDB IDs for that UniProt from the pool,
# and auto-fetches the forced PDB from PDBe/RCSB if it isn't
# already present. No need to also list it in manual_include.
# =========================================================
PDB_FORCE = {
    entry["uniprot"].strip(): entry["pdb_id"].lower().strip()
    for entry in config.get("query", {}).get("pdb_force", [])
}

if PDB_FORCE:
    print(f"\n  pdb_force entries: {PDB_FORCE}", flush=True)
    PDBE_SUMMARY = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/"

    for uniprot_acc, forced_pdb in PDB_FORCE.items():
        # remove all other PDBs for this UniProt
        before = len(all_structures)
        all_structures = [
            s for s in all_structures
            if not (s.get("uniprot") == uniprot_acc and s["pdb_id"] != forced_pdb)
        ]
        removed = before - len(all_structures)
        if removed:
            print(f"  pdb_force: removed {removed} alternative(s) for {uniprot_acc}, keeping {forced_pdb.upper()}", flush=True)

        # auto-fetch forced PDB if not already in pool
        if not any(s["pdb_id"] == forced_pdb for s in all_structures):
            print(f"  pdb_force: {forced_pdb.upper()} not in pool — fetching automatically...", flush=True)
            try:
                resolution = None
                organism   = "Unknown"

                # try PDBe first
                r = requests.get(PDBE_SUMMARY + forced_pdb, timeout=15)
                if r.status_code == 200:
                    data = r.json()
                    pdbe_entry = data.get(forced_pdb, [None])[0]
                    if pdbe_entry:
                        for key in ("resolution", "resolution_high"):
                            val = pdbe_entry.get(key)
                            if val is not None:
                                try:
                                    resolution = float(val)
                                    break
                                except (ValueError, TypeError):
                                    pass
                        for src in pdbe_entry.get("source", []):
                            name = src.get("organism_scientific_name", "")
                            if name and name.lower() != "unknown":
                                organism = name
                                break

                # RCSB fallback if PDBe failed
                if resolution is None:
                    print(f"  pdb_force: PDBe had no data for {forced_pdb.upper()} — trying RCSB...", flush=True)
                    rcsb_r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{forced_pdb.upper()}", timeout=15)
                    if rcsb_r.status_code == 200:
                        refine = rcsb_r.json().get("refine", [{}])
                        if refine:
                            val = refine[0].get("ls_d_res_high")
                            if val is not None:
                                try:
                                    resolution = float(val)
                                except (ValueError, TypeError):
                                    pass
                    poly_r = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{forced_pdb.upper()}/1", timeout=15)
                    if poly_r.status_code == 200:
                        src_orgs = poly_r.json().get("rcsb_entity_source_organism", [])
                        if src_orgs:
                            organism = src_orgs[0].get("ncbi_scientific_name") or src_orgs[0].get("scientific_name") or "Unknown"

                if resolution is None:
                    print(f"  pdb_force: WARNING — could not get resolution for {forced_pdb.upper()}, skipping", flush=True)
                else:
                    species = normalise_species(organism)
                    print(f"  pdb_force: fetched {forced_pdb.upper()} | {species} | {resolution} Å", flush=True)
                    all_structures.append({
                        "pdb_id"    : forced_pdb,
                        "resolution": resolution,
                        "organism"  : species,
                        "uniprot"   : uniprot_acc,
                        "source"    : "pdb_force",
                    })
            except Exception as e:
                print(f"  pdb_force: ERROR fetching {forced_pdb.upper()}: {e}", flush=True)
        else:
            print(f"  pdb_force: {forced_pdb.upper()} already in pool ✓", flush=True)

known_pdbs        = {s["pdb_id"] for s in all_structures}

# Two SIFTS lookups built in one pass:
# sifts_pdb_uniprot : pdb_id -> set of DPS-family UniProt accs (for extras + patching known)
# sifts_all_pdb_uni : pdb_id -> set of ALL UniProt accs (for patching unknowns like 3ak8)
sifts_pdb_uniprot = defaultdict(set)
sifts_all_pdb_uni  = defaultdict(set)

with gzip.open(SIFTS_CACHE, "rt") as f:
    # SIFTS file has a comment line starting with '#' before the real header
    # We must skip it to get correct column names
    raw_lines = (line for line in f if not line.startswith("#"))
    reader = csv.DictReader(raw_lines)
    for i, row in enumerate(reader):
        if i == 0:
            print(f"  SIFTS columns: {list(row.keys())}", flush=True)

        uniprot_acc = (row.get("SP_PRIMARY", "") or row.get("UNIPROT", "")).strip().upper()
        pdb_id      = (row.get("PDB", "")        or row.get("pdb_id", "")).strip().lower()

        if not uniprot_acc or not pdb_id:
            continue

        # Track ALL mappings for known PDB IDs (needed to patch unknowns)
        if pdb_id in known_pdbs:
            sifts_all_pdb_uni[pdb_id].add(uniprot_acc)

        # Track DPS-family mappings (for finding extras)
        if uniprot_acc in uniprot_accessions:
            sifts_pdb_uniprot[pdb_id].add(uniprot_acc)

# Snapshot of DPS-family UniProt accessions collected from InterPro ONLY
# Used later for contamination filter — must be taken BEFORE any patching
original_dps_uniprots = set(uniprot_accessions)

# (a) Patch missing UniProt for existing structures
# First try: match against known DPS UniProt accessions
# Second try: use any UniProt SIFTS maps to this PDB (catches 3ak8/3ak9)
patched = 0
for s in all_structures:
    if s["uniprot"] is not None:
        continue
    pdb_id = s["pdb_id"]

    # Try DPS-family match first
    if pdb_id in sifts_pdb_uniprot:
        acc = next(iter(sifts_pdb_uniprot[pdb_id]))
        s["uniprot"] = acc
        uniprot_accessions.add(acc)
        patched += 1
        print(f"  Patched UniProt (DPS match): {pdb_id.upper()} -> {acc}", flush=True)

    # Fallback: any UniProt from SIFTS for this PDB
    elif pdb_id in sifts_all_pdb_uni:
        # Pick the reviewed (shortest, most likely SwissProt) accession
        accs = sorted(sifts_all_pdb_uni[pdb_id])
        acc  = accs[0]
        s["uniprot"] = acc
        uniprot_accessions.add(acc)
        patched += 1
        print(f"  Patched UniProt (SIFTS all): {pdb_id.upper()} -> {acc}", flush=True)

print(f"  UniProt patched for {patched} structures", flush=True)

# Validate all structures — remove any whose UniProt is not a real DPS protein
# This catches contaminating entries like 8FA2 (P0DTC2 = SARS-CoV-2 spike)
# A structure is valid if:
#   (a) its UniProt was in our original InterPro-derived set, OR
#   (b) it has no UniProt (unknown, keep for now, PDB header will resolve)

valid_structures = []
for s in all_structures:
    uni = s["uniprot"]
    if uni is None or uni in original_dps_uniprots:
        valid_structures.append(s)
    else:
        print(f"  REMOVED non-DPS entry: {s['pdb_id'].upper()} "
              f"(UniProt {uni} not in DPS family)", flush=True)

removed = len(all_structures) - len(valid_structures)
print(f"  Removed {removed} non-DPS entries", flush=True)
all_structures = valid_structures

# (b) Extra PDB entries not found by InterPro
extra_pdb_to_uniprot = {
    pdb_id: next(iter(unis))
    for pdb_id, unis in sifts_pdb_uniprot.items()
    if pdb_id not in known_pdbs
}

print(f"  Extra PDB entries via SIFTS: {len(extra_pdb_to_uniprot)}", flush=True)
if extra_pdb_to_uniprot:
    print(f"  {sorted(extra_pdb_to_uniprot)}", flush=True)

# =========================================================
# STEP 3 — PDBe summary for SIFTS extras
# =========================================================

if extra_pdb_to_uniprot:
    print("\n" + "=" * 60, flush=True)
    print("STEP 3  PDBe summary for SIFTS extras", flush=True)
    print("=" * 60, flush=True)

    PDBE = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/"
    extra_list = sorted(extra_pdb_to_uniprot)

    for i in range(0, len(extra_list), 100):
        batch = extra_list[i:i+100]
        try:
            r = requests.get(PDBE + ",".join(batch), timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"  PDBe error: {e}", flush=True)
            continue

        for pdb_id, entries in data.items():
            if not entries:
                continue
            entry = entries[0]

            resolution = None
            for key in ("resolution", "resolution_high"):
                val = entry.get(key)
                if val is not None:
                    try:
                        resolution = float(val)
                        break
                    except (ValueError, TypeError):
                        pass
            if resolution is None:
                continue   # skip NMR / no resolution

            organism = "Unknown"
            for src in entry.get("source", []):
                name = src.get("organism_scientific_name", "")
                if name and name.lower() != "unknown":
                    organism = name
                    break

            species = normalise_species(organism)
            uniprot = extra_pdb_to_uniprot[pdb_id.lower()]
            uniprot_accessions.add(uniprot)

            print(f"  {pdb_id.upper():<6}  {species:<40}  {resolution:.2f} Å  {uniprot}", flush=True)

            all_structures.append({
                "pdb_id"    : pdb_id.lower(),
                "resolution": resolution,
                "organism"  : species,
                "uniprot"   : uniprot,
                "source"    : "sifts",
            })
        time.sleep(0.3)


# =========================================================
# STEP 3b -- Manual includes (from config.yaml)
#            Structures not annotated in InterPro/SIFTS but
#            scientifically validated. Fetched automatically
#            from PDBe so nothing is truly manual.
# =========================================================

manual_entries = config.get("query", {}).get("manual_include", [])

if manual_entries:
    print("\n" + "=" * 60, flush=True)
    print("STEP 3b  Manual includes from config.yaml", flush=True)
    print("=" * 60, flush=True)

    PDBE_SUMMARY = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/"
    already_have = {s["pdb_id"] for s in all_structures}

    for entry in manual_entries:
        pdb_id = entry.get("pdb_id", "").lower().strip()
        reason = entry.get("reason", "no reason given")

        if not pdb_id:
            continue
        if pdb_id in already_have:
            print(f"  {pdb_id.upper()} already found automatically -- skipping manual add", flush=True)
            continue

        print(f"  Adding {pdb_id.upper()} -- {reason}", flush=True)

        try:
            r = requests.get(PDBE_SUMMARY + pdb_id, timeout=15)
            r.raise_for_status()
            try:
                data = r.json()
                entries_list = data.get(pdb_id, [])
            except Exception:
                print(f"  WARNING: PDBe returned unparseable response for {pdb_id.upper()} — trying RCSB fallback", flush=True)
                entries_list = []
            if not entries_list:
                print(f"  WARNING: PDBe returned no data for {pdb_id.upper()} — trying RCSB fallback", flush=True)
                pdbe_entry = None
            else:
                pdbe_entry = entries_list[0]

            resolution = None

            if pdbe_entry is not None:
                for key in ("resolution", "resolution_high", "r_factor"):
                    val = pdbe_entry.get(key)
                    if val is not None:
                        try:
                            resolution = float(val)
                            break
                        except (ValueError, TypeError):
                            pass

                # Fallback: try PDBe experiment endpoint which has resolution for older entries
                if resolution is None:
                    try:
                        exp_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/{pdb_id}"
                        exp_r = requests.get(exp_url, timeout=15)
                        if exp_r.status_code == 200:
                            exp_data = exp_r.json().get(pdb_id, [])
                            for exp in exp_data:
                                val = exp.get("resolution") or exp.get("resolution_high")
                                if val is not None:
                                    try:
                                        resolution = float(val)
                                        break
                                    except (ValueError, TypeError):
                                        pass
                    except Exception:
                        pass

            # RCSB fallback — used when PDBe has no record (e.g. older entries like 1ZUJ)
            organism_rcsb = None
            if resolution is None or pdbe_entry is None:
                try:
                    print(f"  Trying RCSB for {pdb_id.upper()}...", flush=True)
                    rcsb_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
                    rcsb_r = requests.get(rcsb_url, timeout=15)
                    if rcsb_r.status_code == 200:
                        rcsb_data = rcsb_r.json()
                        if resolution is None:
                            refine = rcsb_data.get("refine", [{}])
                            if refine:
                                val = refine[0].get("ls_d_res_high")
                                if val is not None:
                                    try:
                                        resolution = float(val)
                                    except (ValueError, TypeError):
                                        pass
                        # get organism from RCSB polymer entity
                        poly_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/1"
                        poly_r = requests.get(poly_url, timeout=15)
                        if poly_r.status_code == 200:
                            poly_data = poly_r.json()
                            src_orgs = poly_data.get("rcsb_entity_source_organism", [])
                            if src_orgs:
                                organism_rcsb = src_orgs[0].get("ncbi_scientific_name") or src_orgs[0].get("scientific_name")
                except Exception as rcsb_e:
                    print(f"  RCSB fallback also failed for {pdb_id.upper()}: {rcsb_e}", flush=True)

            if resolution is None:
                print(f"  WARNING: no resolution for {pdb_id.upper()} from PDBe or RCSB -- skipping", flush=True)
                continue

            organism = "Unknown"
            if organism_rcsb:
                organism = organism_rcsb
            elif pdbe_entry is not None:
                for src in pdbe_entry.get("source", []):
                    name = src.get("organism_scientific_name", "")
                    if name and name.lower() != "unknown":
                        organism = name
                        break

            species = normalise_species(organism)
            print(f"  {pdb_id.upper()} | {species} | {resolution} A | manual", flush=True)

            all_structures.append({
                "pdb_id"    : pdb_id,
                "resolution": resolution,
                "organism"  : species,
                "uniprot"   : None,
                "source"    : "manual",
            })

        except Exception as e:
            print(f"  ERROR fetching {pdb_id.upper()} from PDBe: {e}", flush=True)

        time.sleep(0.2)

# =========================================================
# STEP 4 — Paralog-aware deduplication
#
# Group by: species_key -> uniprot_key -> best resolution
#
# species_key:
#   - known organism  -> "Genus species"
#   - unknown org but known UniProt -> "UNKNOWN::UNIPROTACC"
#     (two crystal forms of same protein collapse correctly)
#   - truly unknown   -> "UNKNOWN::pdb_id" (last resort)
#
# uniprot_key:
#   - known UniProt -> UniProt accession
#     (same protein = same UniProt = same key -> keep best resolution)
#   - unknown UniProt -> pdb_id (can't deduplicate further)
#
# Paralog detection:
#   A species_key with >1 distinct uniprot_key = has paralogs
#   Keep best resolution per uniprot_key
# =========================================================

print("\n" + "=" * 60, flush=True)
print("STEP 4  Paralog-aware deduplication", flush=True)
print("=" * 60, flush=True)

grouped = defaultdict(dict)  # species_key -> { uniprot_key -> best struct }

for s in all_structures:
    organism = s["organism"]
    uniprot  = s["uniprot"]
    pdb_id   = s["pdb_id"]

    if organism and organism.lower() not in ("unknown", ""):
        species_key = organism
    elif uniprot:
        species_key = f"UNKNOWN::{uniprot}"
    else:
        species_key = f"UNKNOWN::{pdb_id}"

    uniprot_key = uniprot if uniprot else pdb_id

    current = grouped[species_key].get(uniprot_key)
    if current is None or s["resolution"] < current["resolution"]:
        grouped[species_key][uniprot_key] = s

# Flatten — mark paralogs
selected = []
for species_key, uniprot_dict in sorted(grouped.items()):
    is_paralog = len(uniprot_dict) > 1
    for uniprot_key, s in uniprot_dict.items():
        s = dict(s)
        s["is_paralog"]  = is_paralog
        s["species_key"] = species_key
        selected.append(s)

        tag = f"  [paralog {uniprot_key}]" if is_paralog else ""
        print(f"  {s['pdb_id'].upper():<6}  {species_key:<45}  "
              f"{s['resolution']:.2f} Å{tag}", flush=True)

print(f"\n  Total structures selected: {len(selected)}", flush=True)

# =========================================================
# STEP 5 — Save outputs
# =========================================================

# pdb_ids.txt
selected_pdbs = sorted(s["pdb_id"] for s in selected)
with open("data/pdb_ids.txt", "w") as f:
    for pdb in selected_pdbs:
        f.write(pdb + "\n")

# species_map.json
species_map = {}
for s in selected:
    label = s["organism"] if s["organism"].lower() != "unknown" else s["pdb_id"].upper()
    species_map[s["pdb_id"]] = label
with open("data/species_map.json", "w") as f:
    json.dump(species_map, f, indent=2)

# selection_report.tsv
with open("data/selection_report.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["pdb_id", "organism", "resolution_A", "uniprot", "source", "is_paralog"])
    for s in sorted(selected, key=lambda x: x["pdb_id"]):
        writer.writerow([
            s["pdb_id"].upper(),
            s["organism"],
            f"{s['resolution']:.2f}",
            s["uniprot"] or "",
            s["source"],
            "yes" if s["is_paralog"] else "no",
        ])

print(f"\n{'=' * 60}")
print(f"Total structures : {len(selected_pdbs)}")
print(f"Saved -> data/pdb_ids.txt")
print(f"Saved -> data/species_map.json")
print(f"Saved -> data/selection_report.tsv")
print(f"{'=' * 60}\n")
print("Done.\n")
