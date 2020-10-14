#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-i', '--infile', help = 'Name of the input file in bed-format with regions to merge', type = str, required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'Name of the input file with chr in first column and NC-annotation in second, if you want to convert chr-annotation to NC-annotation', type = str, required = False)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-if', '--inFormat', choices = ['bed', 'sedd'], help = "Set input format", required = True)
parser.add_argument('-of', '--outFormat', choices = ['bed', 'sedd'], help = "Set output format", required = True)


args = parser.parse_args()
# outfile = open(args.output, 'w')
regions = {}
ncConverter = {}
if args.chr2nc:
    with open(args.chr2nc, 'r') as chr2ncFile:
        for line in chr2ncFile:
            if not re.match('^#', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab
                ncConverter[lineSplit[0]] = lineSplit[1]
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
    # Go through the file line by line
    for line in infile:
        # Check that line starts with chr
        if (bedIn and not re.match('^#', line)) or (seddIn and not re.match('^#', line)):
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
                        if not lineSplit[0].startswith("chr"):
                            chrom = ncConverter["chr" + lineSplit[0]]
                        else:
                            chrom = ncConverter[lineSplit[0]]
                    else:
                        print ("ERROR chr2nc has to be given when input format is bed and output format is sedd!!!")
                if not lineSplit[0].startswith("chr"):
                    chrom = "chr" + lineSplit[0]
                start = int(lineSplit[1]) + 1
                end = int(lineSplit[2])

            # If the input format is sedd
            elif seddIn:
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


            # Check if the chromomsome exists as key in dict regions
            # If not add it and then the region
            if not chrom in regions:
                regions[chrom] = {}
                regions[chrom][start] = end
            # If the chrom exists check for overlapping regions
            else:
                overlapping = False
                # Go through all regions on the chromosome
#                for s in regions[chrom]:
                for s in list(regions[chrom]):
                    newStart = s
                    newEnd = regions[chrom][s]
                    # Check if this region overlaps or is adjacent to other region
                    if start <= (regions[chrom][s] + 1) and end >= (s - 1):
                        overlapping = True
                        # If start is less than start for the region set new start position
                        if start < s:
                            newStart = start
                        # If end is greater than end of the region set new end position
                        if end > regions[chrom][s]:
                            newEnd = end

                        # Remove old region and add the new region
                        del regions[chrom][s]
                        regions[chrom][newStart] = newEnd
                # If there were no overlapping regions add the region
                if not overlapping:
                    regions[chrom][start] = end
    # Check if the infile is closed, if not close it
    if not infile.closed:
        infile.close()


with (open(args.output, mode = 'w'))as outfile:
    if not bedOut:
    # Print header
        outfile.write("#Id\tChromosome\tStart\tEnd\n")


    # Go through each gene
    for chrom in regions:
        # Go through all start pos
        for start in regions[chrom]:
            end = regions[chrom][start]
            # If output format is bed set start -1
            if bedOut:
                start = start - 1
                outfile.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + chrom + ":" + str(start) + "-" + str(end) + "\n")
            else:
                outfile.write(chrom + ":" + str(start) + "-" + str(end) + "\t" + chrom + "\t" + str(start) + "\t" + str(end) + "\n")

    # Check if outfile is closed, if not close it
    if not outfile.closed:
        outfile.close()
