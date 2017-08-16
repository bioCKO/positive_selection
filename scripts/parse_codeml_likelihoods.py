#!/usr/bin/env python
##parse_codeml_likelihoods.py
##written 12/18/14 by Groves Dixon
ProgramName = 'parse_codeml_likelihoods.py'
LastUpdated = '12/18/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Parse outputs from codeml files.

'''

##Import Modules 
from sys import argv
from sys import exit

null = argv[1]
alt = argv[2]
outfile = argv[3]


def read(fileName):
    geneList = []
    likeList = []
    npList = []
    with open(fileName, 'r') as infile:
        for line in infile:
            line = line.strip("\n").split()
            gene = line[0]
            np = line[4].split(")")[0]
            like = line[5]
            geneList.append(gene)
            likeList.append(like)
            npList.append(np)
    return geneList, likeList, npList
            

nullGenes, nullLikes, nullNps = read(null)
altGenes, altLikes, altNps = read(alt)

with open(outfile, 'w') as out:
    out.write("EST\tnullLike\taltLike\tnullNP\taltNP")
    for i in range(len(nullGenes)):
        nullGene = nullGenes[i].strip(".codeml")
        altGene = altGenes[i].strip(".codeml")
        nullLike = nullLikes[i]
        altLike = altLikes[i]
        nullNp = nullNps[i]
        altNp = altNps[i]
        if nullGene != altGene:
            exit("genes don't match!")
        outstring = "\n{}\t{}\t{}\t{}\t{}".format(nullGene, nullLike, altLike, nullNp, altNp)
        out.write(outstring)
    
        


