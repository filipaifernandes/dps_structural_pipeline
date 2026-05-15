import sys
import requests

input_file = sys.argv[1]
output_file = sys.argv[2]

labels = []

with open(input_file) as f:

    for line in f:

        if line.startswith(">P1;"):

            try:

                header = line.strip()

                if ";" not in header:
                    continue

                pdb_id = header.split(";")[1][:4].lower()

                # =====================================================
                # FETCH ENTRY
                # =====================================================

                url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"

                r = requests.get(url, timeout=30)

                species = "Unknown"

                if r.status_code == 200:

                    data = r.json()

                    entities = (
                        data.get(
                            "rcsb_entry_container_identifiers",
                            {}
                        ).get(
                            "polymer_entity_ids",
                            []
                        )
                    )

                    if entities:

                        entity_id = entities[0]

                        entity_url = (
                            f"https://data.rcsb.org/rest/v1/core/polymer_entity/"
                            f"{pdb_id}/{entity_id}"
                        )

                        er = requests.get(entity_url, timeout=30)

                        if er.status_code == 200:

                            entity_data = er.json()

                            orgs = entity_data.get(
                                "rcsb_entity_source_organism",
                                []
                            )

                            if orgs:

                                species = orgs[0].get(
                                    "scientific_name",
                                    "Unknown"
                                )

                print(f"{pdb_id} -> {species}")

                labels.append(f"{pdb_id}\t{species}")

            except Exception as e:

                print(f"Error for {pdb_id}: {e}")

# =========================================================
# WRITE ITOL FILE
# =========================================================

with open(output_file, "w") as out:

    out.write("LABELS\n")
    out.write("SEPARATOR TAB\n")
    out.write("DATA\n")

    for line in labels:
        out.write(line + "\n")

print(f"\nSaved -> {output_file}\n")
