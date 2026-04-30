from Bio import SeqIO
import subprocess
import os

input_fasta = "data/pdb_sequences.fasta"
ref_fasta = "data/reference/dps_refs.fasta"
output_fasta = "data/pdb_sequences_filtered.fasta"

query_file = "tmp_query.fasta"
db_name = "tmp_db"

# build BLAST database
subprocess.run(
    f"makeblastdb -in {ref_fasta} -dbtype prot -out {db_name}",
    shell=True,
    check=True
)

filtered = []

for record in SeqIO.parse(input_fasta, "fasta"):
    SeqIO.write(record, query_file, "fasta")

    try:
        result = subprocess.check_output(
            f"blastp -query {query_file} -db {db_name} -outfmt '6 pident length qlen'",
            shell=True
        ).decode().strip()

        if result:
            lines = result.split("\n")
            keep = False

            for line in lines:
                pident, length, qlen = line.split()

                pident = float(pident)
                length = float(length)
                qlen = float(qlen)

                coverage = length / qlen

        # 🔥 stricter biological filter
            if pident >= 40 and coverage >= 0.7:
                keep = True
                break

        if keep:
            filtered.append(record)

    except subprocess.CalledProcessError:
        continue

SeqIO.write(filtered, output_fasta, "fasta")

# cleanup
for f in [query_file]:
    if os.path.exists(f):
        os.remove(f)
