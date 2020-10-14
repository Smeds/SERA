#!/usr/bin/python2.7

"""

Script for filtering and conveting pindel vcf output to annovar input

Elin Falk Sorqvist 20141114

"""

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes vcf from pindel2vcf, filters it and creates an annovar input file")
parser.add_argument('-i', '--infile', help = 'Input file name', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-m', '--minRD', help = 'Minimum read depth for a variant to be included', type = int, default = 1)
parser.add_argument('-v', '--vaf', help = 'Minimum variant allele ratio for a variant to be included', type = float, default = 0)

args = parser.parse_args()
sample = "-"
with open(args.infile, 'r') as infile:
    with (open(args.output, mode = 'w')) as outfile:

        # Go through the file line by line
        for line in infile:
            # Check so the line isn't empty or starts with ##
            if not re.match('^##', line) and not re.match('^$', line):

                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                # If the line start with only 1 # extract the sample name from the last column
                if re.match("^#", line):
                    sample = lineSplit[9]

                # Else extract and print all relevant info
                else:
                    # Split the colummn with reference and variant depth
                    gtype = lineSplit[9].split(":")
                    depth = gtype[1].split(",")
                    alleleFreq = gtype[1]

                    # Calculate tot read depth and the variant allele ratio
                    readDepth = int(depth[0]) + int(depth[1])
                    vRatio = float(depth[1]) / float(readDepth)

                    # Check that the variant allele ratio and the read depth are high enough
                    if vRatio >= args.vaf and readDepth >= args.minRD:
                        # Split the information column
                        infoSplit = lineSplit[7].split(";")
                        end = 0
                        svtype = "-"
                        var = "-"
                        ref = "-"

                        # Go through all information in the information column and extract the relevant parts
                        for info in infoSplit:
                            if re.match("^END", info):
                                tmpSplit = info.split("=")
                                end = tmpSplit[1]
                            elif re.match("^SVTYPE", info):
                                tmpSplit = info.split("=")
                                svtype = tmpSplit[1]

                        # If the variant type is deletion move start one pos forward and set the correct reference
                        start = lineSplit[1]
                        if re.match("DEL", svtype):
                            start = int(lineSplit[1]) + 1
                            ref = lineSplit[3][1:]
                        # If the variant tyoe is insertion set the correct variant
                        elif re.match("INS", svtype):
                            var = lineSplit[4][1:]
                        # If the variant is a combination of deletion and insertion
                        elif re.match("RPL", svtype):
                            start = int(lineSplit[1]) + 1
                            ref = lineSplit[3][1:]
                            var = lineSplit[4][1:]

                        # Check that the chromosome annotation doesn't contain chr or NC notation
                        chrom = lineSplit[0]
                        if re.match("chr", chrom):
                            chrom = chrom.replace("chr", "")
                        elif re.match("NC", chrom):
                            if re.match("NC_000023", chrom):
                                chrom = "chrX"
                            elif re.match("NC_000024", chrom):
                                chrom = "chrY"
                            else:
                                chrom = re.sub(r'NC_0{4,5}', "chr", chrom)
                                chrom = re.sub(r'\.[0-9]{1,2}', "", chrom)

                        # Print in annocar format
                        outfile.write(chrom + "\t" + str(start) + "\t" + end + "\t" + ref + "\t" + var + "\t" + "comments: sample=" + sample + " variantAlleleRatio=" + str(vRatio) + " alleleFreq=" + alleleFreq + " readDepth=" + str(readDepth) + " Tumor_Del=- Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=- Tumor_var_ampliconInfo=- Tumor_ref_ampliconInfo=-\n")

        # Check if outfile is closed, if not close it
    if not infile.closed:
        infile.close()

    # Check if outfile is closed, if not close it
    if not outfile.closed:
        outfile.close()
