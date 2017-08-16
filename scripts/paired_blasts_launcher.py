#!/usr/bin/env python
##Import Modules 
from sys import exit
from sys import argv

FileList = argv[1:]

for query in FileList:
    for database in FileList:
        if query == database:
            continue
        else:
            print "blastp -query {} -db {} -evalue 1e-5 -num_threads 12 -num_alignments 1 -outfmt 5 -out {}-{}.br".format(query, database, query, database)

