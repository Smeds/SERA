#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-i', '--infile', help = 'Name of the input file in bed-format with regions to merge', type = str, required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'Name of the input file with chr in first column and NC-annotation in second, if you want to convert chr-annotation to NC-annotation', type = str, required = False)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-r', '--readLength', help = 'Length of the seqeuncing read', type = int, required = True)
parser.add_argument('-if', '--inFormat', choices = ['bed', 'sedd'], help = "Set input format", required = True)
parser.add_argument('-of', '--outFormat', choices = ['bed', 'sedd'], help = "Set output format", required = True)

args = parser.parse_args()
# outfile = open(args.output, 'w')
regions = {}
ncConverter = {}
rl = args.readLength

if args.chr2nc:
    with open(args.chr2nc, 'r') as chr2ncFile:
        for line in chr2ncFile:
            if not re.match('^#', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab
                ncConverter[lineSplit[0]] = lineSplit[1]

# Set in and output format
bedIn = False
seddIn = False
bedOut = False
seddOut = False

if args.inFormat == "bed":
    bedIn = True
else:
    seddIn = True

if args.outFormat == "bed":
    bedOut = True
else:
    seddOut = True

# Open gene file and add gene info to hash
with open(args.infile, 'r') as infile:
    with (open(args.output, mode = 'w'))as outfile:
        if seddOut:
            # Print header
            outfile.write("#Id\tChromosome\tStart\tEnd\n")

        # Go through the infile line by line
        for line in infile:
            # Check that line starts with chr
            if re.match('chr', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                # bed as default input
                chrom = ""
                start = 0
                end = 0

                # If chromosome should be converted to NC-annotation
                if bedIn:
                    chrom = lineSplit[0]
                    if  seddOut:
                        if args.chr2nc:
                            chrom = ncConverter[lineSplit[0]]
                        else:
                            print ("ERROR chr2nc has to be given when input format is bed and output format is sedd!!!")
                    start = int(lineSplit[1]) + 1
                    end = int(lineSplit[2])

                # If the input format is sedd
                if seddIn:
                    chrom = lineSplit[1]
                    # If output format is bed change from NC to chr annotation
                    if bedOut:
                        if re.match("NC", chrom):
                            if re.match("NC_000023", chrom):
                                chrom = "chrX"
                            elif re.match("NC_000024", chrom):
                                chrom = "chrY"
                            else:
                                chrom = re.sub(r'NC_0{4,5}', "chr", chrom)
                                chrom = re.sub(r'\.[0-9]{1,2}', "", chrom)
                    start = int(lineSplit[2])
                    end = int(lineSplit[3])


                if end - start + 1 > 2 * rl:
                    newEnd = start + rl - 1
                    newStart = end - rl + 1
                    if bedOut:
                        start = start - 1
                        newStart = newStart - 1
                        outfile.write(chrom + "\t" + str(start) + "\t" + str(newEnd) + "\t" + lineSplit[3] + "_1" + "\n")
                        outfile.write(chrom + "\t" + str(newStart) + "\t" + str(end) + "\t" + lineSplit[3] + "_2" + "\n")
                    else:
                        outfile.write(lineSplit[3] + "_1" + "\t" + chrom + "\t" + str(start) + "\t" + str(newEnd) + "\n")
                        outfile.write(lineSplit[3] + "_2" + "\t" + chrom + "\t" + str(newStart) + "\t" + str(end) + "\n")

                else:
                    if bedOut:
                        start = start - 1
                        outfile.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + lineSplit[3] + "\n")
                    else:
                        outfile.write(lineSplit[3] + "\t" + chrom + "\t" + str(start) + "\t" + str(end) + "\n")

        # Check if outfile is closed, if not close it
        if not outfile.closed:
            outfile.close()
    if not infile.closed:
            infile.close()
