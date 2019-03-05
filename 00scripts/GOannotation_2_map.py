#!/usr/bin/env python
# This script takes a many-many list of gene IDs and
# gene ontology IDs and outputs a map file suitable
# for use with the "


import sys
import os

infile = sys.argv[1]
outfile = os.path.splitext(infile)[0] + "_map.txt"

gene2GO = {}

with open(infile) as ifile:
	for line in ifile:
		geneid, goid = line.split()
		if goid == "":
		    continue
		if geneid not in gene2GO:
		    gene2GO[geneid] = []
		gene2GO[geneid].append(goid)


with open(outfile, "w") as ofile:
    for gid in sorted(gene2GO):
        gos = ", ".join(gene2GO[gid])
        ofile.write("\t".join([gid, gos]) + "\n")

