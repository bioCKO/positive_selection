#!/usr/bin/env python
##clean_seq_definitions.py
##written 12/16/14 by Groves Dixon
ProgramName = 'clean_seq_definitions.py'
LastUpdated = '12/16/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
Cleans the sequence definitions of transcriptomes based on some given arguments.

'''

AdditionalProgramInfo = '''
Additional Program Information:

Use this when working with transcriptomes with lots of junk in their sequence definitions (lines other than the
identifying string). 

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
Start_time = time.time() ##keeps track of how long the script takes to run


##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = True, dest = 'input', help = 'The the input transcriptome to be cleaned')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-delimit', required = False, default = 'use_space', dest = 'delimit', help = 'The delmiter. Default is to use white space')
parser.add_argument('-pos', required = True, dest = 'pos', help = 'The position (starting at 1) of the identifying string within the sequence definition lines split by the delimiter. Example: >seq12345 gene=polymerase GO1;GO2;GO3" should be use pos 1 if you want seq12345')
args = parser.parse_args()

#Assign Arguments
InfileName = args.input
OutfileName = args.out
Delimiter = args.delimit
Pos = args.pos
Index = int(Pos) - 1


def clean_file(InfileName, OutfileName, Delimiter, Index):
    '''Function to read in a file as a list of lists
    '''
    lineNumber = 0
    with open(InfileName, 'r') as infile:
        with open(OutfileName, 'w') as out:
            for line in infile:
                lineNumber += 1
                line = line.strip("\n")
                try:
                    if line[0] == ">":
                        if Delimiter == 'use_space':
                            identifier = line.split()[Index]
                        else:
                            identifier = line.split(Delimiter)[Index]
                        #if it was first then we got the carrot already
                        if Index == 0:
                            line = identifier
                        #otherwise add the carrot
                        else:
                            line = '>' + identifier
                    if lineNumber == 1:
                        out.write(line)
                    else:
                        out.write('\n' + line)
                except IndexError: #this allows skipping of blank lines in case the fasta if formatted that way
                    continue


clean_file(InfileName, OutfileName, Delimiter, Index)

#return time to run
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


