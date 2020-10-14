#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a gene file and a variation file from SNPmania and output the mean coverage per gene")
parser.add_argument('-i', '--infile', help = 'Name of the input gene file in bed-format', type = str, required = True)
parser.add_argument('-v', '--variationfile', help = 'File name of variation file from SNPmania', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-g', '--gene', help = 'Set -g to output mean read depth for the whole gene', action = 'store_true')
parser.add_argument('-r', '--region', help = 'Set -r to output mean read depth for every region in the gene file', action = 'store_true')

args = parser.parse_args()
# outfile = open(args.output, 'w')
genes = {}

# Open gene file and add gene info to hash
with open(args.infile, 'r') as infile:
#    print ("Hashing genes!")
    # Go through the file line by line
    for line in infile:
        # Check so the line isn't empty or starts with #
        if re.match('^chr', line):
            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            chrom = lineSplit[0]
            start = int(lineSplit[1]) + 1
            end = int(lineSplit[2])
            gene = lineSplit[3]

            # Check if the tuple exists as key in dict genes
            # If not add it
            # If it does output error message on duplicate genes
            if not gene in genes:
                genes[gene] = {}
                genes[gene][chrom] = {}
                genes[gene][chrom][start] = {}
                genes[gene][chrom][start][end] = {}
                genes[gene][chrom][start][end]['tot_depth'] = 0
                genes[gene][chrom][start][end]['tot_pos'] = 0
            else:
                if not start in genes[gene][chrom]:
                    genes[gene][chrom][start] = {}
                    genes[gene][chrom][start][end] = {}
                    genes[gene][chrom][start][end]['tot_depth'] = 0
                    genes[gene][chrom][start][end]['tot_pos'] = 0
                else:
                    print ("Error gene with coordinates ", gene, " ", chrom, " : ", start, " already exists!")

    # Check if the infile is closed, if not close it
    if not infile.closed:
        infile.close()

with open (args.variationfile) as variationfile:
#    print ("Going through variation file!")
    # Go through the file line by line
    for line in variationfile:

        # Check so the line isn't empty or starts with #
        if not re.match('^#', line) and not re.match('^$', line):
            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab
            chrom = lineSplit[3]

            if re.match("NC", chrom):
                if re.match("NC_000023", chrom):
                    chrom = "chrX"
                elif re.match("NC_000024", chrom):
                    chrom = "chrY"
                else:
                    chrom = re.sub(r'NC_0{4,5}', "chr", chrom)
                    chrom = re.sub(r'\.[0-9]{1,2}', "", chrom)
            # Go through all genes
            for gene in genes:
                # Check if the chromosome exist in dictionary
                if chrom in genes[gene]:
                    # Go through all start positions for that chromosome
                    for start in genes[gene][chrom]:
                        # Check if the startposition is lower than the given position
                        if start <= int(lineSplit[4]):
                            # If so go through all end positions
                            for end in genes[gene][chrom][start]:
                                # Check if the gene end is higher than the position
                                if end >= int(lineSplit[4]):
                                    # If so add pos and read depth to the region
                                    genes[gene][chrom][start][end]['tot_depth'] += int(lineSplit[0])
                                    genes[gene][chrom][start][end]['tot_pos'] += 1

# Open outputfile for writing
with (open(args.output, mode = 'w'))as outfile:
#    print ("Calculating mean read depth!")
    # Print header
    outfile.write("#Type\tName\tChromosome\tStart\tEnd\tMean_read_depth\tTotal_read_depth\tTotal_no_of_positions\n")


    # Go through each gene
    for gene in genes:
        geneTot_depth = 0
        geneTot_pos = 0
        # Set the mean read depth to 0
        meanRD = 0.0
        # Set a super high start pos and an low end pos to be able to store the gene start and end
        s = 1000000000
        e = 0

        # Go through the chromosome
        for chrom in genes[gene]:
            # Go through all start pos
            for start in genes[gene][chrom]:
                # If start pos is lower than s save as the start pos of the gene
                if start < s:
                    s = start
                # Go through each end pos
                for end in genes[gene][chrom][start]:
                    # If the end pos is higher than e save the end pos as the end pos of the gene
                    if end > e:
                        e = end
                    # If mean cov per gene should be printed save info
                    if args.gene:
                        geneTot_depth += genes[gene][chrom][start][end]['tot_depth']
                        geneTot_pos += genes[gene][chrom][start][end]['tot_pos']
                    # If mean read depth per region should be printed, print
                    if args.region:
                        if genes[gene][chrom][start][end]['tot_pos'] > 0:
                            meanRD = float(genes[gene][chrom][start][end]['tot_depth']) / float(genes[gene][chrom][start][end]['tot_pos'])
                        else:
                            meanRD = "-"

                        # Print string
                        outStr = "Region\t" + gene + "\t" + chrom + "\t" + str(start) + "\t" + str(end) + "\t" + str(meanRD) + "\t" + str(genes[gene][chrom][start][end]['tot_depth']) + "\t" + str(genes[gene][chrom][start][end]['tot_pos']) + "\n"
                        outfile.write(str(outStr))
        if args.gene:
            if geneTot_pos > 0:
                meanRD = float(geneTot_depth) / float(geneTot_pos)
            else:
                meanRD = "-"

            # Print string
            outStr = "Gene\t" + gene + "\t" + chrom + "\t" + str(s) + "\t" + str(e) + "\t" + str(meanRD) + "\t" + str(geneTot_depth) + "\t" + str(geneTot_pos) + "\n"
            outfile.write(str(outStr))

    # Check if outfile is closed, if not close it
    if not outfile.closed:
        outfile.close()
