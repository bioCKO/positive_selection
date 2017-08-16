#!/usr/bin/env python
##get_multireciprocal_orthos.py
##written 10/21/14 by Groves Dixon
ProgramName = 'get_multireciprocal_orthos.py'
LastUpdated = '6/11/15'
#Added feature to skip sequences with repeated sequence IDs
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This script produces a table of reciprocal best hits
based on sets of top hits from paired blast results files.

'''

AdditionalProgramInfo = '''
Additional Program Information:

To use this script, first get at least two fasta files
that you think contain orthologs.
Blast each fasta against all others using out format 5 
so that biopython can parse the blast output file.
Recommeded that you set up the blasts using paired_blasts_launcher.py.

From the blast results files extract the best hits using find_best_hits.py.
This produces a .hits file from a blast results file, then concatenate all the
results into a single file. The header is like this:
QueryFasta\tDatabaseFasta\tQuerySequenceID\tBestHitID

These, along with the original fasta files are the inputs for this script.
It works by finding all the reciprocal best hits for a single 'anchor' fasta.
These are considered candidate orthologs. This pool of sequences are then compared 
with each other, so ensure that they are all reciprocally best hits for each other.
If they are not, sequences are iteratively removed until the remaining pool is
fully reciprocal in regard to best hits.
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
parser.add_argument('-hits', required = True, dest = 'hits', help = 'The concatenated top hits file produced from the output of find_best_hits.py')
parser.add_argument('-fa', required = True, nargs = "+", dest = 'fa', help = 'The fasta files')
parser.add_argument('-rcut', required = False, default = 0.5, dest = 'rcut', help = 'The cutoff used to remove a potential reciprocal ortholog. This specifies the proportion of other fasta files a given file needs to have reciprocal orthologs with for its sequence to be included in a given orthologous group. Raise this value to increase stringency for including sequences in orthologous groups.')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired name for the output file')
# parser.add_argument('-id', required = False, dest = 'id', default = 0, type = int, help = 'The Position of the Identifying Information in the Sequence File Names for fa1 (split by blank space). Default = 0')
# parser.add_argument('-delim', required = False, dest = 'delim', default = "_", help = 'The delimiting character in the fasta file names that separate the species name from other text')
parser.add_argument('-debug', required = False, dest = 'debug', default = "FALSE", help = 'Set to TRUE to print debugging information')
parser.add_argument('-anchor', required = True, dest = 'anchor', help = 'Enter one of the fasta files to use as the anchor. The sequence IDs for this fasta will be used as rows for the output table. Only reciprocal orthologs including the anchor will be output. This should either be your best transcripome, or the one you are most interested in.')
args = parser.parse_args()

#Assign Arguments
OutfileName = args.out
debug = args.debug
FA = args.fa
# ID = args.id ##the location of the identifying information in the seq names from fa1 when split by blank space
# Delim = args.delim
NA = 'none'
anchorFasta = args.anchor
hitsFile = args.hits
outfileName = args.out
reciprocalCutoff = float(args.rcut)

def show_fastas(FA):
    print "\nLooking for Reciprocal Orhtologs Between the Following {} Fasta Files:".format(len(FA))
    for i in FA: print i

def get_query_list(fa):
    """Function to pull the full list of query sequences from a fasta file"""
    queryList = []
    repeatedList = []
    for seq_record in SeqIO.parse(fa, "fasta"):
        seqName = seq_record.description.split()[0]
        # if seqName in queryList:
        #     repeatedList.append(seqName)
        #     continue
        # else:
        queryList.append(seqName)
    # print "\n{} sequences were repeated at least once in file {}".format(len(repeatedList), fa)
    return queryList

def store_query_lists(FA):
    """Function to iterate through the fasta file list 
    and record the list of the >sequence id tags for each 
    in a dictionary keyed to the fasta file names
    """
    queryNameDict = {}
    for fasta in FA:
        queryNameDict[fasta] = get_query_list(fasta)
    return queryNameDict
        
        

def build_bestHit_dict(hitsFile, queryDict):
    """Reads in the .hits files to produce a dictionary 
    for accessing the best hits results for all pairs of fasta files.
    Results are stored in subdictionaries keyed to the
    pair of fasta files, then keyed with the query sequence name. Like this:
    {spp1.fasta-sp2.fasta : {spp1seq1 : spp2seq1, spp1seq2 : spp2seq2}, spp1.fa-spp3.fa : { etc..}}
    So all top hits for all pairs are stored in a single massive dictionary.
    
    Two of these dictionaries are output. One has all results, including those without top hits.
    The other is the one (reducedHitDict) that is actually used, which has only the entries with filter passing top hits."""
    #first set up the dictionaries
    hitDict = {}
    reducedHitDict = {}
    queryFilesEncountered = {}
    dbFilesEncountered = {}
    #set up some counts for what happened
    noHitCount = 0
    coverageFailCount = 0
    percentIdFailCount = 0
    evalueFailCount = 0
    TotalPossiblePairs = 0
    lineNumber = 0
    #read in the concatentated best hits file
    with open(hitsFile, 'r') as infile:
        for line in infile:
            lineNumber += 1
            #skip the header
            if lineNumber == 1: continue
            TotalPossiblePairs += 1
            line = line.strip("\n").split('\t')
            queryFasta = line[0]
            dbFasta = line[1]
            querySeq = line[2]
            bestHit = line[3]
            fastaPair = "{}-{}".format(queryFasta, dbFasta)
            queryFilesEncountered[queryFasta] = 1
            dbFilesEncountered[dbFasta] = 1
            try:
                hitDict[fastaPair][querySeq] = bestHit
                if bestHit == 'no_hit':
                    noHitCount += 1
                    continue
                if bestHit == 'coverage_to_low':
                    coverageFailCount += 1
                    continue
                if bestHit == 'Percent_ID_to_low':
                    percentIdFailCount += 1
                    continue
                if bestHit == 'Eval_Over':
                    evalueFailCount += 1
                    continue
                reducedHitDict[fastaPair][querySeq] = bestHit
            except KeyError:
                # print "\nFound new fasta pair: {}".format(fastaPair)
                hitDict[fastaPair] = {querySeq : bestHit}
                # print hitDict[fastaPair]
                if bestHit == 'no_hit':
                    noHitCount += 1
                    continue
                if bestHit == 'coverage_to_low':
                    coverageFailCount += 1
                    continue
                if bestHit == 'Percent_ID_to_low':
                    percentIdFailCount += 1
                    continue
                if bestHit == 'Eval_Over':
                    evalueFailCount += 1
                    continue
                reducedHitDict[fastaPair] = {querySeq : bestHit}
    queryFileList = queryFilesEncountered.keys()
    queryFileList.sort()
    dbFileList = dbFilesEncountered.keys()
    dbFileList.sort()
    return hitDict, reducedHitDict, queryFileList, dbFileList




def is_reciprocal(reducedHitDict, queryFile, dbFile, querySeq):
    """Function to that reads a reducedHitDict output from build_bestHit_dict
    for a given query fasta, database fasta, and query sequence and returns 
    whether this pair of fastas had a reciprocal best hit for that particular sequence.
    It returns a boolean result (TRUE or FALSE) as well as the top hit if there is one"""
    result = "FALSE"
    asQuery = "{}-{}".format(queryFile, dbFile)
    asDb = "{}-{}".format(dbFile, queryFile) 
    try:
        bestHit = reducedHitDict[asQuery][querySeq]
    except KeyError:
        bestHit = 'no_best_hit'
    try:
        recipHit = reducedHitDict[asDb][bestHit]
    except KeyError:
        recipHit = 'no_recip_hit'
    if recipHit == querySeq:
        result = "TRUE"
    return result, bestHit
    
 
def pull_reciprocals(reducedHitDict, anchorFasta, FA):
    """Step 1: Make list of fasta files with reciprocal best hits with an anchor fasta
    and a dictionary to keep track of their sequence IDs
    Step 2: Iterate through the members of the group and check that they are reciprocal with each other.
    If the group is not fully reciprocal, then remove the 'least reciprocal' of the group until it is.
    Step 3: Record the reciprocal fasta files and their sequence IDs for each anchor fasta sequence.
    The results are recorded based on the anchor fasta's sequence IDs."""
    #iterate through all entries in the anchor fasta
    resultsDict = {}
    for anchorEntry in queryDict[anchorFasta]:
        if debug != "FALSE":
            print "\n\n--------------"
            print "Looking for orthologs for the following anchor sequence:"
            print anchorEntry
        resultsDict[anchorEntry] = {}
        reciprocalFileList = []
        reciprocalHitDict = {}
        for fastaFile in FA:
            if fastaFile == anchorFasta: continue
            isReciprocal, bestHit = is_reciprocal(reducedHitDict, anchorFasta, fastaFile, anchorEntry)
            #for each fasta file with a reciprocal best hit, append the fasta file and its sequence ID
            if isReciprocal == "TRUE":
                reciprocalFileList.append(fastaFile)
                reciprocalHitDict[fastaFile] = bestHit
            else:
                continue
        #now we have the list of reciprocal best hitting fastas with the anchor
        #next we want to check that they are reciprocal with each other
        if len(reciprocalFileList) > 0:
            if len(reciprocalFileList) == 1: s = ' has'
            else: s = 's have'
            if debug != "FALSE":
                print "{} Fasta file(s) found with reciprocal orthologs with this anchor sequence:".format(len(reciprocalFileList))
                for recipFile in reciprocalFileList:
                    print recipFile
                print "Number of reciprocal files = {}".format(len(reciprocalFileList))
        #if we have no reciprocal files at all, we don't want to shave it down, just skip to the next anchor entry
        else:
            if debug != "FALSE":
                print "No reciprocal files were found for this anchor sequence"
            continue
        fullyReciprocal = 'FALSE'
        #set up scoring method to check if the set of files are fully reciprocal with each other
        #This loop will run until all remaining fasta files in the group are fully reciprocal
        if debug != "FALSE":
            print "Checking Reciprocality Between the Files..."
        #assign a perfect score as the length of the file list. 
        #this score indicates the number of other files each one should be reciprocal with (ie all other fastas including the anchor)
        perfectScore = (len(reciprocalFileList) - 1)
        #assign the acceptable score
        acceptableScore = np.floor(perfectScore*reciprocalCutoff)
        #set up a scoreDict to record how many other members in teh group each fasta is reciprocal with
        scoreDict = {}
        #run through and assign scores for being reciprocal
        #for each file in the current pool of files
        for recipQuery in reciprocalFileList:
            #get the sequence name
            querySeq = reciprocalHitDict[recipQuery]
            scoreDict[recipQuery] = 0
            #now iterate through the other files in the pool to check if they're reciprocal with recipQuery
            for recipDb in reciprocalFileList:
                #skip the current file
                if recipQuery == recipDb: continue
                isReciprocal, bestHit = is_reciprocal(reducedHitDict, recipQuery, recipDb, querySeq)
                if isReciprocal == "TRUE":
                    scoreDict[recipQuery] += 1
                else:
                    continue
        if debug != "FALSE":
            print "Scores:"
            print scoreDict
            print "Perfect Score = {}".format(perfectScore)
            print "Minimum Score = {}".format(acceptableScore)
        keepList = []
        removedList = []
        for z in reciprocalFileList:
            if scoreDict[z] >= acceptableScore:
                keepList.append(z)
            else:
                removedList.append(z)
        #now rewrite the reciprocalFileList as only the keepers before restarting the loop
        reciprocalFileList = []
        for z in keepList: reciprocalFileList.append(z)
        if debug != "FALSE":
            print "{} File(s) removed because score not above acceptable score".format(len(removedList))
            print "{} File(s) had Reciprocal Best Hit with {}:".format(len(reciprocalFileList), anchorEntry)
            print reciprocalFileList
            print "Here are the best hits:"
        #if we don't have any reciprocal files just continue to next sequence
        #should hit this if statement after breaking the while loop because all reciprocal fastas were removed
        if len(reciprocalFileList) == 0:
            if debug != "FALSE":
                print "No set of fully reciprocal orthologs was found for this sequence. Continuing to next anchor sequence."
            continue
        for fasta in reciprocalFileList:
            if debug != "FALSE":
                print reciprocalHitDict[fasta]
            #record the data for this anchor fasta sequence
            resultsDict[anchorEntry][fasta] = reciprocalHitDict[fasta]
    return resultsDict

def output_recips(resultDict, queryDict, FA):
    """Outputs a table of results. For now they are 
    output in terms of the fa1 file entries"""
    print "\nOutputting Results..."
    scoreList = []
    geneList = []
    nonAnchors = []
    for i in FA:
        if i != anchorFasta: nonAnchors.append(i)
    nonAnchors.sort()
    with open(OutfileName, 'w') as out:
        hstring = "\t".join(nonAnchors)
        header = "{}\t{}".format(anchorFasta, hstring)
        out.write(header)
        #iterate through the sequence IDs for hte anchor fasta
        for i in queryDict[anchorFasta]:
            resultList = []
            data = resultDict[i]
            for fasta in nonAnchors:
                try:
                    bestHit = data[fasta]
                except KeyError:
                    bestHit = 'NA'
                resultList.append(bestHit)
            score = 0
            for x in resultList:
                if x != "NA":
                    score += 1
            scoreList.append(score)
            geneList.append(i)
            resultString = "\t".join(resultList)    
            outString = "\n{}\t{}".format(i, resultString)
            out.write(outString)
    with open('representationCounts.txt', 'w') as out:
        out.write('count')
        for i in range(len(scoreList)):
            score = scoreList[i]
            gene = geneList[i]
            out.write("\n{}\t{}".format(gene, score))




#run the functions
show_fastas(FA)
queryDict = store_query_lists(FA)
hitDict, reducedHitDict, queryFileList, dbFileList = build_bestHit_dict(hitsFile, queryDict)
resultDict = pull_reciprocals(reducedHitDict, anchorFasta, FA)
output_recips(resultDict, queryDict, FA)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


