#!/usr/bin/env python3
from modeller import *
from modeller.scripts import complete_pdb
from modeller import selection
import os

env = environ()

aln = alignment(env)
aln.append(file="data/tree/structural.ali", align_codes='all')

n = len(aln)
names = [aln[i].name for i in range(n)]

os.makedirs("data/tree", exist_ok=True)

with open("data/tree/distance_matrix.txt", "w") as f:
    f.write("\t" + "\t".join(names) + "\n")
    
    for i in range(n):
        row = [names[i]]
        for j in range(n):
            s1 = selection(aln[i])
            s2 = selection(aln[j])
            rms = s1.superpose(s2)[0]
            row.append(f"{rms:.4f}")
        f.write("\t".join(row) + "\n")

print("Distance matrix written to data/tree/distance_matrix.txt")
