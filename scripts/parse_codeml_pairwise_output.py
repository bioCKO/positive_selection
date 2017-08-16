#!/usr/bin/env python
##parse_codeml_pairwise_output.py
##written 6/26/14 by Groves Dixon
ProgramName = 'parse_codeml_pairwise_output.py'
LastUpdated = '8/5/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Parses a list of codeml output files that were generated using pair-wise
dN/dS estimation (runmode -2). Pairs are set up against one base species
(set as spp1) and all other species (a list file)

'''

AdditionalProgramInfo = '''
Additional Program Information:


'''

##Import Modules 
DEBUG = 'FALSE'
import time
import argparse
from sys import argv
from sys import exit
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-f', required = True, dest = 'files', nargs="+", help = 'A glob to the codeml output files (probably *.codeml)')
parser.add_argument('-spp1', required = True, dest = 'spp1', help = 'The search tag for species 1')
parser.add_argument('-sppList', required = True, dest = 'sppList', help = 'The List of species to pair with species 1')
parser.add_argument('-orthos', required = True, dest = 'orthologTable', help = 'Table of orthologs for reassigning gene names. Note that column headings in this table must contain the species names')
parser.add_argument('-o', required = True, dest = 'out', help = 'The desired output file name')
args = parser.parse_args()

#Assign Arguments
FileList = args.files
Spp1 = args.spp1 #the species we are gathering data for
SppListName = args.sppList
OutfileName = args.out
Orthologs = args.orthologTable
SppList = []
with open(SppListName, 'r') as infile:
    for line in infile:
        SppList.append(line.strip("\n"))


def read_orthologs(Orthologs):
    """Read in the ortholog table and fill out a dictionary"""
    lineNumber = 0
    orthoDict = {}
    with open(Orthologs, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip("\n").split('\t')
            if lineNumber == 1:
                header = line
                sppPos = 'species_not_found'  
                for i in range(len(header)):
                    if Spp1 in header[i]:
                        sppPos = i
                if sppPos == 'species_not_found':
                    print "Could not find column in header that contained species name {}".format(Spp1)
                    print "Please check ortholog table and species name"
                    exit()
                else:
                    print "\nFound column for species of interest called {}".format(header[sppPos])
                continue
            else:
                anchorName = line[0]
                speciesOrtho = line[sppPos]
                orthoDict[anchorName] = speciesOrtho
    return orthoDict

def read_files(FileList, Spp1, SppList):
    '''Function to reads through each file and parses out
    dN and dS estimates for the specified species pair.
    '''
    print "\nLooking for data in {} codeml output files...".format(len(FileList))
    geneList = []
    dNList = []
    dSList = []
    dndsList = []
    speciesList = []
    highDScount = 0
    #iterate through all the species in the species list to find data for each pairing
    for species in SppList:
        if species == Spp1:
            continue #don't look for comparisons of the species to itself
        #iterate through all the codeml output files to gather data
        for file in FileList:
            with open(file, 'r') as infile:
                hit = 0 #this variable specifies whether we have found the species pair yes
                hitCount = 0 #this should never exceed 1
                for line in infile:
                    if hitCount > 1:
                        exit("Found more than one instance of pairing in a file. Something is wrong.")
                    if hit == 0:
                        ##look for your species pair
                        if "("+Spp1+")" in line:
                            if "("+species+")" in line:
                                if "..." in line:
                                    hit = 1
                                    continue
                    elif hit == 1:
                        if "dN/dS=" in line:
                            line = line.strip("\n")
                            noWhite = line.replace(" ", "")
                            dn = noWhite.split("dN=")[1].split("dS")[0]
                            ds = noWhite.split("dN/dS=")[1].split("dN=")[1].split("dS=")[1]
                            dnds = noWhite.split("dN/dS=")[1].split("dN")[0]
                            if DEBUG == 'TRUE':
                                print line
                                print noWhite
                                print "dN = {}".format(dn)
                                print "dS = {}".format(ds)
                                print "dNdS = {}".format(dnds)
                            geneName = file[0:-7] #again here stripping '.codeml' was removing the c from the front of the gene names
                            geneList.append(geneName)
                            dNList.append(dn)
                            dSList.append(ds)
                            dndsList.append(dnds)
                            speciesList.append(species)
                            hit = 0
                            hitCount += 1
                            # print geneName
                            # print species
                            # print dn
    return geneList, dNList, dSList, dndsList, speciesList
                        
def output(OutfileName, geneList, dNList, dSList, dndsList, speciesList, orthoDict):
    """Outputs the data into a table"""
    badValues = []
    lineNums = []
    with open(OutfileName, 'w') as out:
        out.write("contig\tspecies\tdN\tdS\tdNdS")
        for i in range(len(geneList)):
            #########
            ##there is a bug that occurs when the synonymous substitution rate is >99.99
            #these are obviously wrong anyway and they stop the output from uploading into R so skip them
            # fourData = 'TRUE'
            # outList = [geneList[i], speciesList[i], dNList[i], dSList[i]]
            # try:
            #     float(dNList[i])
            #     float(dSList[i])
            # except ValueError:
            #     badValues.append([dNList[i], dSList[i]])
            #     lineNums.append(i)
            #     continue
            # for x in outList:
            #     if x == "":
            #         fourData = 'FALSE'
            # if fourData == 'FALSE':
            #     continue
            ###########
            anchor = geneList[i]
            speciesContig = orthoDict[anchor]
            out.write("\n{}\t{}\t{}\t{}\t{}".format(speciesContig, speciesList[i], dNList[i], dSList[i], dndsList[i]))         
                            
orthoDict = read_orthologs(Orthologs)
geneList, dNList, dSList, dndsList, speciesList = read_files(FileList, Spp1, SppList)
output(OutfileName, geneList, dNList, dSList, dndsList, speciesList, orthoDict)


#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


