#!/usr/bin/python2.7

import re
import argparse

'''

Script for converting Pindel output to Annovar input

'''

parser = argparse.ArgumentParser(description = "This script takes 3 files and gives an amplification ratio for each region in the amplified file")
parser.add_argument('-a', '--amplifiedfile', help = 'Input file, regions with possible amplification, bed-format (Required)', required = True)
parser.add_argument('-b', '--backgroundfile', help = 'File with background regions likely not background, bed-format (Required)', required = True)
parser.add_argument('-s', '--snpmaniafile', help = 'Variation file from SNPmania (Required)', required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'File with conversion between NC-number and chr', type = str, required = True)
parser.add_argument('-o', '--outputfile', help = 'Output txt file (Required)', required = True)
parser.add_argument('-t', '--tumour', help = 'Tumour type', type = str, required = False)
parser.add_argument('-sample', '--sample', help = 'Sample name', type = str, required = False)

args = parser.parse_args()

# Create dictionary for convertion between nc and chr
nc2chr = {}
with open(args.chr2nc, 'r') as cnfile:
    # Go through the file line by line
    for line in cnfile:
        if line.startswith("chr"):
            line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
            nc2chr[line[0]] = line[1]  # create the dictionary with chr as key and nc as value
    # If the file is not closed, close
    if not cnfile.closed:
        cnfile.close()



amplified = {}
with open(args.amplifiedfile, 'r') as ampfile:
    for line in ampfile:
        if re.match('^chr', line):
            line = line.rstrip('\r\n').split("\t")
            chromosome = nc2chr[line[0]]
            start = line[1]
            end = line[2]
            gene = line[3]

            if not chromosome in amplified:
                amplified[chromosome] = {}
                amplified[chromosome][start] = {}
                amplified[chromosome][start][end] = {}
                amplified[chromosome][start][end]['gene'] = gene
                amplified[chromosome][start][end]['totDepth'] = 0
                amplified[chromosome][start][end]['covBases'] = 0

            else:
                if not start in amplified[chromosome]:
                    amplified[chromosome][start] = {}
                    amplified[chromosome][start][end] = {}
                    amplified[chromosome][start][end]['gene'] = gene
                    amplified[chromosome][start][end]['totDepth'] = 0
                    amplified[chromosome][start][end]['covBases'] = 0
                else:
                    if not end in amplified[chromosome][start]:
                        amplified[chromosome][start][end] = {}
                        amplified[chromosome][start][end]['gene'] = gene
                        amplified[chromosome][start][end]['totDepth'] = 0
                        amplified[chromosome][start][end]['covBases'] = 0

background = {}
with open(args.backgroundfile, 'r') as ampfile:
    for line in ampfile:
        if re.match('^chr', line):
            line = line.rstrip('\r\n').split("\t")
            chromosome = nc2chr[line[0]]
            start = line[1]
            end = line[2]
            gene = line[3]

            if not chromosome in background:
                background[chromosome] = {}
                background[chromosome][start] = {}
                background[chromosome][start][end] = {}
                background[chromosome][start][end]['gene'] = gene
                background[chromosome][start][end]['totDepth'] = 0
                background[chromosome][start][end]['covBases'] = 0

            else:
                if not start in background[chromosome]:
                    background[chromosome][start] = {}
                    background[chromosome][start][end] = {}
                    background[chromosome][start][end]['gene'] = gene
                    background[chromosome][start][end]['totDepth'] = 0
                    background[chromosome][start][end]['covBases'] = 0
                else:
                    if not end in background[chromosome][start]:
                        background[chromosome][start][end] = {}
                        background[chromosome][start][end]['gene'] = gene
                        background[chromosome][start][end]['totDepth'] = 0
                        background[chromosome][start][end]['covBases'] = 0

#             print (chr + "\t" + start + "\t" + end + "\t" + ref + "\t" + var + "\t" + blacklist[chr][start][end][ref][var])

with open(args.snpmaniafile, 'r') as snpfile:

    # Go through the file line by line
    for line in snpfile:

        line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
        rd = int(line[0])
        chromosome = line[3]
        pos = line[4]


        if chromosome in amplified:
            for ampStart in amplified[chromosome]:
                if pos > amplified[chromosome]:
                    for ampEnd in amplified[chromosome][ampStart]:
                        if pos <= ampEnd:
                            amplified[chromosome][ampStart][ampEnd]['totDepth'] += rd
                            amplified[chromosome][ampStart][ampEnd]['covBases'] += 1



        elif chromosome in background:
            for backStart in background[chromosome]:
                if pos > background[chromosome]:
                    for backEnd in background[chromosome][backStart]:
                        if pos <= backEnd:
                            background[chromosome][backStart][backEnd]['totDepth'] += rd
                            background[chromosome][backStart][backEnd]['covBases'] += 1



with open(args.outputfile, 'w') as outfile:
    backTotDepth = 0
    backCovBases = 0

    for backChr in background:
        for backStart in background[backChr]:
            for backEnd in background[backChr][backStart]:
                backTotDepth += background[backChr][backStart][backEnd]['totDepth']
                backCovBases += background[backChr][backStart][backEnd]['covBases']
    backRdBase = 0
    if float(backCovBases > 0):
        backRdBase = float(backTotDepth) / float(backCovBases)

    for ampChr in amplified:
        for ampStart in amplified[ampChr]:
            for ampEnd in amplified[ampChr][ampStart]:
                ampRdBase = 0
                if amplified[ampChr][ampStart][ampEnd]['covBases'] > 0:
                    ampRdBase = float(amplified[ampChr][ampStart][ampEnd]['totDepth']) / float (amplified[ampChr][ampStart][ampEnd]['covBases'])
                amplification = 0
                if backRdBase > 0:
                    amplification = float(ampRdBase) / float(backRdBase)
                else:
                    amplification = "NA"
                # Add info about tumour type
                tumourType = "-"
                if args.tumour:
                    tumourType = (args.tumour)

                # Add info about sample
                sampleName = str(args.snpmaniafile)
                if args.sample:
                    sampleName = str(args.sample)
                outfile.write(sampleName + "\t" + tumourType + "\t" + ampChr + "\t" + ampStart + "\t" + ampEnd + "\t" + amplified[ampChr][ampStart][ampEnd]['gene'] + "\t" + str(amplification) + "\n")
