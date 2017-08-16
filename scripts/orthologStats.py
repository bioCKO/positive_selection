#!/usr/bin/env python
##Import Modules 
from sys import exit
from sys import argv
import numpy as np

OrthoTable = argv[1]
cut = argv[2]
cut = float(cut)
outfileName = OrthoTable.strip(".txt") + "_rep{}.txt".format(cut)
# print "cut = {}".format(cut)

passing = 0
rejected = 0
totalGenes = 0
with open(OrthoTable, 'r') as infile:
    with open(outfileName, 'w') as out:
        lineNumber = 0
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split("\t")
            if lineNumber == 1:
                header = line
                total = (len(header) - 1)
                minimum = int(np.round(cut*total))
                # print "minium = {}".format(minimum)
                out.write("\t".join(header))
                continue
            else:
                totalGenes += 1
                missing = line.count("NA")
                count = total - missing
                # print line[0] + "\t" + str(count)
                if count >= minimum:
                    passing += 1
                    out.write("\n" + "\t".join(line))
                #skip lines without sufficient representation
                else:
                    rejected += 1
                    continue


print "\nOf {} total genes...".format(totalGenes)
print "{} had at least {} representative species and were kept".format(passing, minimum)
print "{} were rejected because they did not have enough species".format(rejected)
                
        
        
        
        

