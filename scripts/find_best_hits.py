#!/usr/bin/env python
##find_best_hits.py
##written 7/21/15 by Groves Dixon
ProgramName = 'find_best_hits.py'
LastUpdated = '7/21/15'
#Added feature to skip sequences with repeated sequence IDs
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script builds a table of best hits from a blast results file.

'''

AdditionalProgramInfo = '''
Additional Program Information:

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-br', required = True, dest = 'br', help = 'The blast output file')
parser.add_argument('-query', required = True, dest = 'query', help = 'The name of the fasta used as the query file for the blast')
parser.add_argument('-db', required = True, dest = 'db', help = 'The name of the fasta used as the database file for the blast')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
# parser.add_argument('-id', required = False, dest = 'id', default = 0, type = int, help = 'The Position of the Identifying Information in the Sequence File Names for fa1 (split by blank space). Default = 0')
# parser.add_argument('-delim', required = False, dest = 'delim', default = "_", help = 'The delimiting character in the fasta file names that separate the species name from other text')
parser.add_argument('-e', required = False, dest = 'eval', default = 1e-10, help = 'Set an optional theshold to use for e-value. (Default is 1e-10)')
parser.add_argument('-cov', required = False, dest = 'hitPct', default = 75, type = float, help = 'Set an optional threshold for hit percentage (or hit coverage). This is equal to the length of the alignment divided by the length of the hit sequence. Setting this higher will ensure that you are only calling orthologs for entire sequences instead of just portions of them (Default 75)')
parser.add_argument('-pctID', required = False, dest = 'percentID', default = 75, type = float, help = 'Set an optional threshold for percent identity. Here percent identity is the number of identical positions divided by the length of the alignment. Setting this higher will return orthologs with greater sequence similarity (Default 75)')
parser.add_argument('-debug', required = False, dest = 'debug', default = "FALSE", help = 'Set to TRUE to print debugging information')
args = parser.parse_args()

#Assign Arguments
OutfileName = args.out
debug = args.debug
BR = args.br
Ecut = args.eval
Ecut = float(Ecut)
# ID = args.id ##the location of the identifying information in the seq names from fa1 when split by blank space
# Delim = args.delim
HitPctCut = args.hitPct / 100.0
PctIdentityCut = args.percentID / 100.0
NA = 'none'
QueryFasta = args.query
DatabaseFasta = args.db
outfileName = args.out

def get_query_list(QueryFasta):
    """Function to pull the full list of query sequences from a fasta file"""
    print "\nFinding all sequence IDs in the query fasta file..."
    queryList = []
    repeatedList = []
    for seq_record in SeqIO.parse(QueryFasta, "fasta"):
        seqName = seq_record.description.split()[0]
        queryList.append(seqName)
    print "\nFound {} sequences in the query fasta file".format(len(queryList))
    return queryList

def find_top_hits(BR, queryList):
    """ """
    print "\nFinding the top hits in the blast results file..."
    #first set up the dictionary
    resultsList = []
    #now read through the blast results to get the best hits
    result_handle = open(BR) ##open the blast results file
    for blast_record in NCBIXML.parse(result_handle): ##start iterating through the blast records using the parser
        entry = blast_record.query
        entry = entry.split()[0]
        #set up the list of data for this particular sequence, including the query fasta, the database, and the sequence ID
        recordResults = [QueryFasta, DatabaseFasta, entry]
        try:
            topHit = blast_record.alignments[0]
            hitName = topHit.title
            # print "entry {}     hitTitle {}".format(entry, hitName)
            hitName = hitName.split()[1].split()[0]
            alignmentLength = len(topHit.hsps[0].match)
            coverage = float(alignmentLength)/float(topHit.length)
            pctID = float(topHit.hsps[0].identities) / float(alignmentLength)
            # print "entry {}     hitName {}".format(entry, hitName)
            #check that the hit percent is above cutoff, if it isn't then change hitName to 'Pct_Over'
            # if topHit.hsps[0].
            if debug == "TRUE":
                print
                print
                print "Hit Name = {}".format(hitName)
                print "Alignment:"
                print topHit.hsps[0].match                                                 ##prints out the alignment
                print "Alignment Length = {}".format(alignmentLength)                      ##print the length of the alignment
                print "Percent Identity of alignment = {}".format(pctID)                   ##print proportion of identical positions to total alignment length
                print "Evalue =  {}".format(topHit.hsps[0].expect)                         ##prints out the evalue
                print "Hit Sequence Length = {}".format(topHit.length)                     ##prints out the length of the hit sequence
                print "Number Identical Positions = {}".format(topHit.hsps[0].identities)  ##prints out the number of identical positions
                print "Ratio of alignment length to total hit length = {}".format(coverage)
            #check that the hit percentage is higher than the cutoff
            if coverage < HitPctCut:
                # print 'failed hit percent'
                hitName = 'coverage_to_low'
            #look at the evalue for the top hit
            if pctID < PctIdentityCut:
                hitName = 'Percent_ID_to_low'
            if topHit.hsps[0].expect > Ecut:
                hitName = 'Eval_Over'
            if debug == "TRUE":
                print hitName
        except IndexError:
            topHit = 'no hit'
            hitName = 'no_hit'
        # add the top hit name as the final piece of information for this blast record
        recordResults.append(hitName)
        #now add the results for this blast record to the overall results list
        resultsList.append(recordResults)
    return resultsList

def output(resultsList):
    """Outputs the results"""
    with open(outfileName, 'w') as out:
        out.write("QueryFasta\tDatabaseFasta\tquerySeqID\tTopHitSeqID")
        for i in resultsList:
            outString = "\t".join(i)
            out.write("\n{}".format(outString))

    
queryList = get_query_list(QueryFasta)
resultsList = find_top_hits(BR, queryList)
output(resultsList)

 
#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


