#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
import personnamenorm as pnn

outfile = open(snakemake.output[0], "w+")
with open(snakemake.input[0], "r") as names_file:
    for names in names_file:
        nameobj = pnn.namenorm(names.strip())
        outfile.write(" ".join(nameobj.name["lastname"]) + "\n")

outfile.close()
