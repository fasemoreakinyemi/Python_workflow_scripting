#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
import personnamenorm as pnn

outfile = open(snakemake.output[0], "w+")
with open(snakemake.input[1], "r") as names_file:
    if snakemake.input[1] == "firstname":
        nametype = "firstname"
    else:
        nametype = "lastname"
    for names in names_file:
        nameobj = pnn.namenorm(names.strip())
        outfile.write(" ".join(nameobj.name[nametype]) + "\n")

outfile.close()
