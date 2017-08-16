#!/usr/bin/env python
##parse_sites_model_likelihoods.py
##written 7/29/15 by Groves Dixon
ProgramName = 'parse_sites_model_likelihoods.py'
LastUpdated = '7/29/15'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Parses the raw output from Codeml running all 5 sites models M1, M2a, M2b, M7 and M8
and outputs the data for comparisons using R program LRT_for_sites_models.R
'''

AdditionalProgramInfo = '''
Additional Program Information:

Expected input is like this:
TR10006-c0_g1_i1	lnL(ntime:  3  np:  5):	   -334.125176      +0.000000
TR10006-c0_g1_i1	lnL(ntime:  3  np:  6):	   -334.125413      +0.000000
TR10006-c0_g1_i1	lnL(ntime:  3  np:  8):	   -334.125412      +0.000000
TR10006-c0_g1_i1	lnL(ntime:  3  np:  6):	   -334.125344      +0.000000
TR10006-c0_g1_i1	lnL(ntime:  3  np:  8):	   -334.125345      +0.000000
TR10033-c0_g2_i2	lnL(ntime:  3  np:  5):	   -331.822553      +0.000000
TR10033-c0_g2_i2	lnL(ntime:  3  np:  6):	   -331.822700      +0.000000
TR10033-c0_g2_i2	lnL(ntime:  3  np:  8):	   -331.822714      +0.000000
TR10033-c0_g2_i2	lnL(ntime:  3  np:  6):	   -331.822664      +0.000000
TR10033-c0_g2_i2	lnL(ntime:  3  np:  8):	   -331.822682      +0.000000
TR1003-c0_g1_i2	lnL(ntime:  3  np:  5):	   -313.910267      +0.000000
TR1003-c0_g1_i2	lnL(ntime:  3  np:  6):	   -313.885512      +0.000000
TR1003-c0_g1_i2	lnL(ntime:  3  np:  8):	   -313.777036      +0.000000
TR1003-c0_g1_i2	lnL(ntime:  3  np:  6):	   -313.887335      +0.000000
TR1003-c0_g1_i2	lnL(ntime:  3  np:  8):	   -313.777036      +0.000000

We are after the gene name, the number of parameters, and the likelihoods.
Gene name is listed five times with five likelihoods because we ran 5 models.

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the input file')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out

def read_file(InfileName):
    '''Function reads through the raw data file and
    generates a dictionary organizing the data into lists by gene
    Output:
    {gene1 : {NSsites : [0, 1, 2, 7, 8], np : [np0, np1, np2, np7, np8], likelihood : [l0, l1, l2, l7, l8]}}
    
    '''
    lineNumber = 0
    data = {}
    nsSiteList = [0, 1, 2, 7, 8]
    nsIndex = 0
    geneList = []
    with open(InfileName, 'r') as infile:
        for line in infile:
            lineNumber += 1
            line = line.strip('\n').split()
            gene = line[0]
            np = line[4].strip("):")
            l = line[5]
            try:
                data[gene]['NSsites'].append(nsSiteList[nsIndex])
                data[gene]['np'].append(np)
                data[gene]['likelihood'].append(l)
                nsIndex += 1
            #if we have not encountered this gene before, start a new entry in dictionary
            except KeyError:
                geneList.append(gene)
                nsIndex = 0
                data[gene] = {}
                data[gene]['NSsites'] = [nsSiteList[nsIndex]]
                data[gene]['np'] = [np]
                data[gene]['likelihood'] = [l]
    return data, geneList

def output(GeneList, DataDict):
    """output the data for R input"""
    with open(OutfileName, 'w') as out:
        header = "gene\tnp0\tnp1\tnp2\tnp7\tnp8\tl0\tl1\tl2\tl7\tl8"
        out.write(header)
        for i in GeneList:
            data = DataDict[i]
            nps = []
            ls = []
            for x in range(len(data['NSsites'])):
                nps.append(data['np'][x])
                ls.append(data['likelihood'][x])
            outString = i + "\t" + "\t".join(nps) + "\t" + "\t".join(ls)
            out.write("\n{}".format(outString))
        
        
    

DataDict, GeneList = read_file(InfileName)
output(GeneList, DataDict)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


