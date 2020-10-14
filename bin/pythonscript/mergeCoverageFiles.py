#!/usr/bin/python2.7

import argparse
import re
from os import walk

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-d', '--directory', help = 'Name of the directory where the coverage files to merge are', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)

args = parser.parse_args()
# outfile = open(args.output, 'w')
genes = {}
genes['sample'] = []

files =[]

# List all files in the directory given
for (dirpath, dirnames, filenames) in walk (args.directory):
    # Only save the filenames in the given dir not in sub-dir if it exists
    files.extend(filenames)
    break

# Go through all files in the list
for f in files:
    fn = args.directory+f

    # Extract sample name from filename
    idSplit = f.split(".")
    genes['sample'].append(idSplit[0])
    # Open gene file and add gene info to hash
    with open(fn, 'r') as infile:
        for line in infile:
            # Check so the line isn't empty or starts with #
            if not re.match('#', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab
                tup = lineSplit[1]+"#"+lineSplit[2]+"#"+lineSplit[3]+"#"+lineSplit[4]
                if tup in genes:
                    genes[tup].append(lineSplit[5])
                else:
                    genes[tup] = []
                    genes[tup].append(lineSplit[5])

with (open(args.output, mode = 'w'))as outfile:
    outStr = "Gene\tChromosome\tStart\tEnd"
    for sample in genes['sample']:
        outStr += "\t"+sample
    outfile.write(outStr+"\n")

    for region in genes:
        if re.search("#", str(region)):
            infoSplit = region.split("#")
            outStr = '\t'.join(infoSplit)
            for sample in genes[region]:
                outStr+= "\t"+sample
            outfile.write(outStr+"\n")
