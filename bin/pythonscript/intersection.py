#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "")
parser.add_argument('-i1', '--infile1', help = 'Name of the first input file', type = str, required = True)
parser.add_argument('-i2', '--infile2', help = 'Name of the second input file', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-if1', '--inFormat1', choices = ['bed', 'sedd'], help = "Set input format for input file 1", required = True)
parser.add_argument('-if2', '--inFormat2', choices = ['bed', 'sedd'], help = "Set input format for input file 2", required = True)
parser.add_argument('-of', '--outFormat', choices = ['bed', 'sedd'], help = "Set output format", required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'Name of the input file with chr in first column and NC-annotation in second, if you want to convert chr-annotation to NC-annotation', type = str, required = False)


args = parser.parse_args()
# outfile = open(args.output, 'w')
regions = {}
intersection = {}

# Create dict with convertion from chr to NC
ncConverter = {}
if args.chr2nc:
    with open(args.chr2nc, 'r') as chr2ncFile:
        for line in chr2ncFile:
            if not re.match('^#', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab
                ncConverter[lineSplit[0]] = lineSplit[1]

# Set in and output format
bedIn1 = False
seddIn1 = False
bedIn2 = False
seddIn2 = False
bedOut = False
seddOut = False

if args.inFormat1 == "bed":
    bedIn1 = True
else:
    seddIn1 = True

if args.inFormat2 == "bed":
    bedIn2 = True
else:
    seddIn2 = True

if args.outFormat == "bed":
    bedOut = True
else:
    seddOut = True



# Open gene file and add gene info to hash
with open(args.infile1, 'r') as infile1:
    # Go through the file line by line
    for line in infile1:
        # Check that line doesn't start with # or is empty
        if not re.match('^#', line) and not re.match('^$', line) and not re.match("browser", line) and not re.match("track", line):
            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            # bed as default input
            chrom = ""
            start = 0
            end = 0

            # If chromosome should be converted to NC-annotation
            if bedIn1:
                chrom = lineSplit[0]
                if  seddOut:
                    if args.chr2nc:
                        chrom = ncConverter[lineSplit[0]]
                    else:
                        print ("ERROR chr2nc has to be given when input format is bed and output format is sedd!!!")

                start = int(lineSplit[1]) + 1
                end = int(lineSplit[2])

            # If the input format is sedd
            elif seddIn1:
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
            # If the chrom exists check for overlapping regions (make sure overlapping regions are merged)
            else:
                overlapping = False
                # Go through all regions on the chromosome
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
    if not infile1.closed:
        infile1.close()

# Open gene file and add gene info to hash
with open(args.infile2, 'r') as infile2:
    # Go through the file line by line
    for line in infile2:
        # Check that line doesn't start with # or is empty
        if not re.match('^#', line) and not re.match('^$', line) and not re.match("browser", line) and not re.match("track", line):
            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            # bed as default input
            chrom = ""
            start = 0
            end = 0

            # If chromosome should be converted to NC-annotation
            if bedIn1:
                chrom = lineSplit[0]
                if  seddOut:
                    if args.chr2nc:
                        chrom = ncConverter[lineSplit[0]]
                    else:
                        print ("ERROR chr2nc has to be given when input format is bed and output format is sedd!!!")
                start = int(lineSplit[1]) + 1
                end = int(lineSplit[2])

            # If the input format is sedd
            if seddIn2:
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

            # Check if this chromosome was present in the other file
            if chrom in regions:
                # Go through all regions on the chromosome
                for s in list(regions[chrom]):
                    newStart = s
                    newEnd = regions[chrom][s]
                    # Check if this region overlaps with a region in the other file
                    if start <= regions[chrom][s] and end >= s:
                        # If start is greater than start for the region set new start position
                        if start > s:
                            newStart = start
                        # If end is less than end of the region set new end position
                        if end < regions[chrom][s]:
                            newEnd = end

                        # For the intersected region, check if it's chromoome exists in the intersection dict
                        if chrom in intersection:
                            overlapping = False
                            # Go through all regions on the chromosome
                            for sInter in list(intersection[chrom]):
                                newStartInter = sInter
                                newEndInter = intersection[chrom][sInter]
                                # Check if this region overlaps or is adjacent to any other intersection region
                                if int(newStart) <= int(intersection[chrom][sInter]) + 1 and int(newEnd) >= int(sInter) - 1:
                                    overlapping = True
                                    # If newStart is less than startInter for the region set new newStart position
                                    if newStart < sInter:
                                        newStartInter = newStart
                                    # If newEnd is greater than end of intersection region set new newEnd position
                                    if newEnd > intersection[chrom][sInter]:
                                        newEndInter = newEnd

                                    # Remove old region and add the new region
                                    del intersection[chrom][sInter]
                                    intersection[chrom][newStartInter] = newEndInter
                            # If there were no overlapping intersection add the region
                            if not overlapping:
                                intersection[chrom][newStart] = newEnd
                        # If the chromosome wasn't present in the intersection dict add it and it's corresponding region
                        else:  #          Print header
                            intersection[chrom] = {}
                            intersection[chrom][newStart] = newEnd


    # Check if the infile is closed, if not close it
    if not infile2.closed:
        infile2.close()
with (open(args.output, mode = 'w'))as outfile:
    if seddOut:
        # Print header
        outfile.write("#Id\tChromosome\tStart\tEnd\n")


    # Go through each gene
    for chrom in intersection:
        # Go through all start pos
        for start in intersection[chrom]:
            end = intersection[chrom][start]

            if bedOut:
                start = start - 1
                outfile.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + chrom + ":" + str(start) + "-" + str(end) + "\n")
            else:
                outfile.write(chrom + ":" + str(start) + "-" + str(end) + "\t" + chrom + "\t" + str(start) + "\t" + str(end) + "\n")

    # Check if outfile is closed, if not close it
    if not outfile.closed:
        outfile.close()
