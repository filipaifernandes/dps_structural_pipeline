from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import os

pdb_dir = "data/raw/"
output = "data/pdb_sequences.fasta"

parser = PDBParser(QUIET=True)

records = []

valid_residues = {
    "ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL"
}

for file in os.listdir(pdb_dir):
    if not file.endswith(".pdb"):
        continue

    pdb_id = file.replace(".pdb", "")
    structure = parser.get_structure(pdb_id, os.path.join(pdb_dir, file))

    for model in structure:
        for chain in model:
            seq = ""

            for res in chain:
                if res.get_resname() in valid_residues:
                    seq += seq1(res.get_resname())

            # 🔥 length filter (Dps-like)
            if 150 <= len(seq) <= 250:
                records.append(
                    SeqRecord(Seq(seq), id=f"{pdb_id}_{chain.id}", description="")
                )
                break
        break

SeqIO.write(records, output, "fasta")
