import os
import json
import re
from Bio.PDB import PDBList

out_dir = "data/raw"
os.makedirs(out_dir, exist_ok=True)

with open("data/pdb_ids.txt") as f:
    pdb_ids = [line.strip().lower() for line in f if line.strip()]

# Load PDB blacklist from config.yaml — nothing hardcoded here
import yaml as _yaml
with open("config.yaml") as _f:
    _config = _yaml.safe_load(_f)
PDB_BLACKLIST = {
    entry["pdb_id"].lower()
    for entry in _config.get("query", {}).get("pdb_blacklist", [])
}
if PDB_BLACKLIST:
    print(f"Skipping blacklisted PDB IDs (from config): {PDB_BLACKLIST}")
pdb_ids = [p for p in pdb_ids if p not in PDB_BLACKLIST]

species_map_path = "data/species_map.json"
with open(species_map_path) as f:
    species_map = json.load(f)

# Load selection report to get resolution + uniprot for dedup
report_path = "data/selection_report.tsv"
resolution_map = {}  # pdb_id -> resolution
uniprot_map    = {}  # pdb_id -> uniprot accession
try:
    with open(report_path) as f:
        for line in f:
            if line.startswith("pdb_id"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                pid = parts[0].lower()
                resolution_map[pid] = float(parts[2])
                uniprot_map[pid]    = parts[3]  # uniprot column
except Exception as e:
    print(f"Warning: could not load report: {e}")

pdbl = PDBList()
successful = []
failed     = []

def extract_organism_from_pdb_header(pdb_path):
    """
    Parse ORGANISM_SCIENTIFIC from PDB file SOURCE record.
    PDB files always have this — last resort for unknown organisms.
    """
    try:
        with open(pdb_path) as f:
            in_source = False
            for line in f:
                record = line[:6].strip()
                if record == "SOURCE":
                    in_source = True
                elif in_source and record not in ("SOURCE", ""):
                    in_source = False
                if in_source:
                    match = re.search(
                        r"ORGANISM_SCIENTIFIC:\s*([^;]+);",
                        line, re.IGNORECASE
                    )
                    if match:
                        name = match.group(1).strip().title()
                        # Normalise to binomial
                        parts = name.split()
                        return " ".join(parts[:2]) if len(parts) >= 2 else name
    except Exception:
        pass
    return None


# =========================================================
# Download all PDB files + patch unknown organisms
# =========================================================

for pdb_id in pdb_ids:
    try:
        print(f"Downloading {pdb_id}...")

        file_path = pdbl.retrieve_pdb_file(
            pdb_id,
            pdir=out_dir,
            file_format="pdb"
        )

        if not file_path or not os.path.exists(file_path):
            print(f"Failed: {pdb_id}")
            failed.append(pdb_id)
            continue

        new_path = os.path.join(out_dir, f"{pdb_id}.pdb")
        os.rename(file_path, new_path)
        successful.append(pdb_id)

        # Patch unknown organisms from PDB header
        current_label = species_map.get(pdb_id, "")
        if not current_label or current_label == pdb_id.upper():
            organism = extract_organism_from_pdb_header(new_path)
            if organism:
                print(f"  Resolved organism: {pdb_id} -> {organism}")
                species_map[pdb_id] = organism

    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        failed.append(pdb_id)

# =========================================================
# Post-download deduplication
#
# After PDB header resolution, structures that were "UNKNOWN"
# may now share the same species name.
# e.g. 3ak8 and 3ak9 both resolve to Salmonella enterica
#      -> keep only the best resolution one
#
# IMPORTANT: only deduplicate structures with the SAME UniProt
# accession (same protein). Structures with different UniProt
# accessions are genuine paralogs and must both be kept,
# even if they share a species name.
#
# This runs AFTER all downloads so we have full organism info.
# =========================================================

print("\n--- Post-download deduplication ---")

# Group by (species, uniprot) -> list of (pdb_id, resolution)
group_key_to_pdbs = {}
for pdb_id in successful:
    species    = species_map.get(pdb_id, pdb_id.upper())
    resolution = resolution_map.get(pdb_id, 999.0)
    uniprot    = uniprot_map.get(pdb_id, "")
    # Key: same species + same UniProt = same protein = dedup candidate
    # Key: same species + different UniProt = paralog = keep both
    key = (species, uniprot)
    if key not in group_key_to_pdbs:
        group_key_to_pdbs[key] = []
    group_key_to_pdbs[key].append((pdb_id, resolution))

# Find groups with multiple structures -> keep best resolution
to_remove = set()
for (species, uniprot), entries in group_key_to_pdbs.items():
    if len(entries) > 1:
        entries.sort(key=lambda x: x[1])
        best = entries[0]
        duplicates = entries[1:]
        print(f"  {species} [{uniprot or 'no UniProt'}]: "
              f"keeping {best[0].upper()} ({best[1]:.2f} Å), "
              f"removing {[d[0].upper() for d in duplicates]}")
        for dup_id, _ in duplicates:
            to_remove.add(dup_id)

if to_remove:
    # Remove duplicate PDB files
    for pdb_id in to_remove:
        pdb_file = os.path.join(out_dir, f"{pdb_id}.pdb")
        if os.path.exists(pdb_file):
            os.remove(pdb_file)
            print(f"  Removed duplicate file: {pdb_file}")
        # Remove from species map
        species_map.pop(pdb_id, None)

    # Rewrite pdb_ids.txt without duplicates
    remaining = [p for p in pdb_ids if p not in to_remove]
    with open("data/pdb_ids.txt", "w") as f:
        for pdb_id in sorted(remaining):
            f.write(pdb_id + "\n")
    print(f"  Removed {len(to_remove)} duplicate(s), "
          f"{len(remaining)} structures remaining")
else:
    print("  No duplicates found after organism resolution")

# =========================================================
# Save updated species map
# =========================================================

with open(species_map_path, "w") as f:
    json.dump(species_map, f, indent=2)

print("\n--- SUMMARY ---")
print(f"Successful : {len(successful)}")
print(f"Failed     : {len(failed)}")
print(f"Removed    : {len(to_remove)} (post-download duplicates)")

with open("data/failed_pdb_ids.txt", "w") as f:
    for pid in failed:
        f.write(pid + "\n")
