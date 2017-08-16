#!/usr/bin/env python
##output_ortholog_seqs.py
##written 10/15/14 by Groves Dixon
ProgramName = 'output_ortholog_seqs.py'
LastUpdated = '12/14/14'
#revised to output seqs headed by only first 10 characters of species label to match with RAxML style
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Takes a set of orthologs and a list of fasta files and exports the 
sequences for each ortholog as fastas for alignment.

It reads a set of reciprocal orthologs, then reads the paired lists of protein and 
nucleotide fasta files to get ortholog sequences from. It then pulls 
the appropriate sequences for each ortholog from each fasta. 
The sequecnes for each chosen species are then output as fasta files
(one for the nucleotide one for the protein) named for the basal ortholog sequence
name, with sequence names in the fasta named for the species it came from.


'''

AdditionalProgramInfo = '''
Additional Program Information:
Species names must be included in the protein and nucleotide fasta file names.
The species names must match those in the header of the orthologs file so thay can be matched up.
The protein and nucleotide fasta files need to have similar notation, each indicating the species
it came from in the same position in the name. 
For example:
Amillepora_CDS.fasta and 
Amillepora_PRO.fasta would work.

so would 
nuc_amil.fasta and
prot_amil.fasta

This wouldn't work:
Amillepora_Nuc.fasta
amil_prot.fasta

Note that species names included in outputs will be reduced to first 10 characters to make RAxML happy.
The base species must be in the first column.

Note the one file fasta file will be output for each ortholog. The sequences for each species with a sequence in that 
ortholog pool will be included with the species name for the >sequence_id.
The fasta file name will be the name of the sequence ID for the base species in that orthologous group. In the rare
event that pipes "|" are included in sequence IDs, these will be changed to "-" when saving the fasta file names.
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
parser.add_argument('-orthos', required = True, dest = 'orthos', help = 'The table of orthologs. Should be tab delimited and should have species names as the column headings that match the species names included in fasta file names')
parser.add_argument('-prot', required = True, dest = 'prot', nargs="+", help = 'A glob to the protein fastas, for example *PRO.fas')
parser.add_argument('-nucl', required = True, dest = 'nucl', nargs="+", help = 'A glob to the nucleotide fastas, for example *CDS.fas')
parser.add_argument('-spp', required = False, dest = 'spp', default = 0, help = 'Integer designating the pythonic position of the species identifier in the names of the fasta files. The default is 0, which would apply to a fasta file name like this: Amillepora_protein_seqs.fasta. Note that ortholog headings in the ortholog table must match the species identifiers in the fasta file names.')
parser.add_argument('-sep', required = False, dest = 'sep', default = '_', help = 'The delimiter used in the fasta file names. (Default is "_")')
parser.add_argument('-ignore', required = False, dest = 'ignore', default = False, help = 'A species to ignore when outputing. This was added to be able to include')
# parser.add_argument('-as_base', required = False, dest = 'base', default = "none", help = 'The species fasta name that you used as the base.')
args = parser.parse_args()

#Assign Arguments
Orthos = args.orthos
Prots = args.prot
Nucs = args.nucl
Spp = args.spp
Spp = int(Spp)
Sep = args.sep
sppListFromFastas = []
sppListFromOrthologTable = []
IgnoreSpp = args.ignore


def check_files(Nucs, Prots):
    for i in Nucs:
        sppListFromFastas.append(i.split(Sep)[Spp])
    print "\nFound Nucleotide Fastas for the following {} species:".format(len(sppListFromFastas))
    for i in sppListFromFastas:
        print i
    sppList2 = [] 
    for i in Prots:
        sppList2.append(i.split(Sep)[Spp])
    print "\nFound Protein Fastas for the following {} species:".format(len(sppList2))
    for i in sppList2:
        print i
    if IgnoreSpp:
        print("\nIgnoring Species {} for outputs".format(IgnoreSpp))



def read_orthos(Orthos):
    """Function to read in the ortholog table and build 
    a set of nested dictionaries linking each species'
    ortholog to the basal ortholog name (assumed to be the first column of the table)
    Output is like this:
    {BaseOrtholog1 : {spp1 : spp1ortho, spp2 : spp2ortho etc.}, BaseOrtholog2 : {spp1 : spp1ortho, etc.}
    
    It also makes the converse nested dictionary where baseOrthos are values and inidvidual species orthos are keys.
    {Species1 : {sppOrtho1 : base1, sppOrtho2 : base2}, Species2 : {sppOrtho1 : base 1, etc.} etc.}
    """
    with open(Orthos, 'r') as infile:
        sppListFromOrthologTable = []
        base2xDict = {}
        x2baseDict = {}
        orthoList = [] #list of all the base orthologs
        lineNumber = 0
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split('\t')
            if lineNumber == 1:
                print "\nChecking that the header in the ortholog table can be matched to the specified nucleotide and protein fasta files..."
                baseSpecies = line[0].split(Sep)[Spp]
                fastaList = line #grab the list of all fasta files including the base fasta
                #edit them into species names using arguments Sep and Spp to edit fasta file names
                for x in fastaList:
                    speciesName = x.split(Sep)[Spp]
                    sppListFromOrthologTable.append(speciesName)
                    if speciesName not in sppListFromFastas:
                        if speciesName == IgnoreSpp:
                            print(speciesName + "\t" + 'ignoring')
                        else:
                            exit("Failed to Match this entry from header {} from the ortholog header with Fasta species Names".format(x))
                    else:
                        print(speciesName + "\t" + 'check')
                continue
            #for all lines other than the header line
            #assign the base sequence ID, set up the dictionary keyed to this base sequence ID and add it to the ortholist.
            base = line[0]
            base2xDict[base] = {}
            orthoList.append(base)
            for i in range(len(sppListFromOrthologTable)):
                spp = sppListFromOrthologTable[i] #grab the species name
                ortho = line[i]#set the sequence ID for this species
                #attach the species name with its orthologous sequence ID to the base ID in the base2xDict (link the species and its ortholog to the base ortholog ID)
                base2xDict[base][spp] = ortho
                #record the converse dictionary
                #the species may not have been encountered yet, so do it with error exceptions
                try:
                    x2baseDict[spp][ortho] = base
                except KeyError:
                    x2baseDict[spp] = {ortho : base}
    #now we have the ortholist (all base ortholog IDs), base2xDict (links base ortholog IDs to species linked to their specific orhtolog ID), and the converse dictionary that goes {spp: {ortho: base}
    print "\nUsing {} as the base species name".format(baseSpecies)
    return orthoList, base2xDict, x2baseDict, baseSpecies, sppListFromOrthologTable
        
def pull_seqs(orthoList, base2xDict, x2baseDict, fileSet):
    """go through the dictionaries that have the ortholog data and pull out the sequences.
    Generates a single dictionary that holds the orthologous sequences for all species.
    seqDict = {species1 : {baseID1 = species1_orthologous_sequence, baseID2 : species1_orthologous_sequence}}"""
    print "\nPulling sequences from fastas..."
    seqDict = {}
    for fasta in fileSet:
        # print fasta
        #pull the species name from the fasta file name
        species = fasta.split(Sep)[Spp]
        # print species
        #set up the nested dictionary for this sepcies
        seqDict[species] = {}
        #parse the fasta
        fasSeqs = SeqIO.parse(open(fasta), 'fasta')
        #iterate through the seqs
        for seq in fasSeqs:
            try:
                x2baseDict[species]
            except KeyError:
                exit("Species {} is missing for some reason".format(species))
            try:
                #try to pull the base ortholog connected with this sequence id
                SEQID = seq.id.split()[0]
                base = x2baseDict[species][SEQID]
            except KeyError:
                #if you can't find a base ortholog for this one then it is not one of the reciprocal orthologs
                continue
            #if you found a connection with a base ortholog, then record the sequence in the sequence dictionary
            seqDict[species][base] = seq.seq.upper()
    return seqDict

def output(sppListFromOrthologTable, orthoList, seqDict, outSuffix):
    """Output the fasta files for each ortholog. Note that file name is based on the 
    base ortholog sequence name and that pipes switched for dashes if present.
    Iterates through all the base sequence names included in the ortholog table.
    Then iterates through each species included in the species table.
    Fasta files are like this:
    >species1
    ATGCATGC
    >species2
    ATGCAGTAC
    etc.
    """
    for base in orthoList:
        outFileName = base + outSuffix
        #sometimes sequence names have pipe characters in them
        #replace those with "-"s here
        out2 = ""
        for char in outFileName:
            if char != "|":
                out2 += char
            else:
                out2 += "-"
        outFileName = out2
        with open(outFileName, 'w') as out:
            counter = 0
            for species in sppListFromOrthologTable:
                counter += 1
                if counter == 1:
                    carrot = ">"
                else:
                    carrot = "\n>"
                try:
                    out.write("{}{}\n{}".format(carrot, species[0:10], seqDict[species][base])) #output only first 10 characters of species to match RAxML style
                except KeyError:
                    continue
            
                
        
    

check_files(Nucs, Prots)
orthoList, base2xDict, x2baseDict, baseSpecies, sppListFromOrthologTable = read_orthos(Orthos)
# if debug != 'FALSE':
#     print "finished reading orthos"
#     for i in orthoList:
#         print "---------"
#         print i
#         f = base2xDict[i]['Fscutaria']
#         b = x2baseDict['Fscutaria'][f]
#         print "{} == {}?".format(i, b)
#         print b == i

nucSeqDict = pull_seqs(orthoList, base2xDict, x2baseDict, Nucs)
protSeqDict = pull_seqs(orthoList, base2xDict, x2baseDict, Prots)
# for i in nucSeqDict.keys():
#     print nucSeqDict[i]
output(sppListFromOrthologTable, orthoList, nucSeqDict, "_nuc.fasta")
output(sppListFromOrthologTable, orthoList, protSeqDict, "_prot.fasta")

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


