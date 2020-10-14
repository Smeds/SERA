#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes an annovar output file and a region file with sequences and output the variants in a vcf-file")
parser.add_argument('-s', '--seqfile', help = 'Name of the input file with region and sequences, tab-delimited', type = str, required = True)
parser.add_argument('-v', '--variantfile', help = 'Name of the variation file from annovar', type = str, required = True)
parser.add_argument('-p', '--pindelfile', help = 'Name of the pindel file from annovar', type = str, required = False)
parser.add_argument('-chr2nc', '--chr2nc', help = 'File with conversion between NC-number and chr', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)

args = parser.parse_args()
# outfile = open(args.output, 'w')

# Create dictionary for convertion between nc and chr
nc2chr = {}
with open(args.chr2nc, 'r') as cnfile:
    # Go through the file line by line
    for line in cnfile:
        if line.startswith("chr"):
            line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
            nc2chr[line[1]] = line[0]  # create the dictionary with nc as key and chr as value

    # If the file is not closed, close
    if not cnfile.closed:
        cnfile.close()

# Create dictionary for the regions included
regions = {}
# Open gene file and add gene info to hash
with open(args.seqfile, 'r') as seqfile:
    # Go through the file line by line
    for line in seqfile:
        # Check so the line isn't empty or starts with #
            line = line.rstrip('\r\n')  # Remove new line character
            line = line.split("\t")  # Split line by tab

            chrom = line[1]
            start = int(line[2])
            end = int(line[3])
            seq = line[4]

            if chrom.startswith("NC"):
                chrom = nc2chr[chrom]
                chrom = chrom[3:]
            elif chrom.startswith('chr'):
                chrom = chrom[3:]


            # Check if the chromosome exist in the region dictionary already
            if chrom in regions:
                if start in regions[chrom]:  # Check if start exists in the dictionary
                    if end > regions[chrom][start]:  # Check if the end is greater than the one already existing
                        regions[chrom][start] = {key: value for key, value in regions[chrom][start].items() if value != regions[chrom][start]}
                        regions[chrom][start][end] = seq
                else:  # If start doesn't exist in dict create it
                    regions[chrom][start] = {}
                    regions[chrom][start][end] = seq  # save sequence
            else:  # If chrom doesn't exist in dict create the whole structure
                regions[chrom] = {}
                regions[chrom][start] = {}
                regions[chrom][start][end] = seq

    # If the file is not closed, close
    if not seqfile.closed:
        seqfile.close()

with open(args.output, 'w') as outputfile:
    outputfile.write("##fileformat=VCFv4.2\n")
    outputfile.write("##reference=hg19\n")
    outputfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    # Open gene file and add gene info to hash
    with open(args.variantfile, 'r') as variantfile:
        # Go through the file line by line
        for line in variantfile:
            # Check so the line isn't empty or starts with #
            if not re.match('^#', line) and not re.match('$', line) and not re.match('^Sample', line):

                line = line.rstrip('\r\n')  # Remove new line character
                line = line.split("\t")  # Split line by tab

                chrom = line[1]
                start = int(line[2])
                end = int(line[3])
                ref = line[4]
                var = line[5]
                refRd = line[10]
                varRd = line[11]
                rd = line[12]
                dbSnp = line[14]
                cosmic = line[17]
                qual = "."
                filter = "."
                info = "."
                id = "NA"

                if dbSnp is "-":
                    if cosmic is "-":
                        id = "."
                    else:
                        cosmic = cosmic.split(";")
                        cosmic = cosmic[0].split("=")
                        id = cosmic[1]
                else:
                    id = dbSnp
                    if not cosmic is "-":
                        cosmic = cosmic.split(";")
                        cosmic = cosmic[0].split("=")
                        id = id + ";" + cosmic[1]

                if chrom.startswith("NC"):
                    chrom = nc2chr[chrom]
                    chrom = chrom[3:]
                elif chrom.startswith('chr'):
                    chrom = chrom[3:]

                if ref is "-":
                    # # INSERTION
                    # print ("Insertion\t"+chrom)
                    if chrom in regions:
                        for s in regions[chrom].keys():  # Go through all start positions on the chromosome in the dictionary
                            if s <= start:  # Check if the region start is lower than the start of the variant
                                for e in regions[chrom][s].keys():  # If so go through the end positions of the regions
                                    if e >= end:  # If the end of the region is greater than the variant end extract seq

                                        rSeq = regions[chrom][s][e]
                                        rSeq = rSeq[(start - s):(start - s + 1)]
                                        outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + rSeq + "\t" + rSeq + var + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")


                elif var is "-":
                    # # DELETION
                    # print("Deletion")
                    if chrom in regions:
                        start = start - 1
                        for s in regions[chrom].keys():  # Go through all start positions on the chromosome in the dictionary
                            if s <= start:  # Check if the region start is lower than the start of the variant
                                for e in regions[chrom][s].keys():  # If so go through the end positions of the regions
                                    if e >= end:  # If the end of the region is greater than the variant end extract seq

                                        rSeq = regions[chrom][s][e]
                                        rSeq = rSeq[(start - s):(start - s + 1)]
                                        outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + rSeq + ref + "\t" + rSeq + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")


                else:
                    # # VARIANT
                    # print ("Variant")
                    outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + ref + "\t" + var + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")
        # If the file is not closed, close
        if not variantfile.closed:
            variantfile.close()

    if args.pindelfile:
        with open(args.pindelfile, 'r') as pindelfile:
            # Go through the file line by line
            for line in pindelfile:
                # Check so the line isn't empty or starts with #
                if not re.match('^#', line) and not re.match('$', line) and not re.match('^Sample', line):

                    line = line.rstrip('\r\n')  # Remove new line character
                    line = line.split("\t")  # Split line by tab

                    chrom = line[1]
                    start = int(line[2])
                    end = int(line[3])
                    ref = line[4]
                    var = line[5]
                    refRd = line[10]
                    varRd = line[11]
                    rd = line[12]
                    dbSnp = line[14]
                    cosmic = line[17]
                    qual = "."
                    filter = "."
                    info = "."
                    id = "NA"

                    if dbSnp is "-":
                        if cosmic is "-":
                            id = "."
                        else:
                            cosmic = cosmic.split(";")
                            cosmic = cosmic[0].split("=")
                            id = cosmic[1]
                    else:
                        id = dbSnp
                        if not cosmic is "-":
                            cosmic = cosmic.split(";")
                            cosmic = cosmic[0].split("=")
                            id = id + ";" + cosmic[1]

                    if chrom.startswith("NC"):
                        chrom = nc2chr[chrom]
                        chrom = chrom[3:]
                    elif chrom.startswith('chr'):
                        chrom = chrom[3:]

                    if ref is "-":
                        # # INSERTION
                        # print ("Insertion\t"+chrom)
                        if chrom in regions:
                            for s in regions[chrom].keys():  # Go through all start positions on the chromosome in the dictionary
                                if s <= start:  # Check if the region start is lower than the start of the variant
                                    for e in regions[chrom][s].keys():  # If so go through the end positions of the regions
                                        if e >= end:  # If the end of the region is greater than the variant end extract seq

                                            rSeq = regions[chrom][s][e]
                                            rSeq = rSeq[(start - s):(start - s + 1)]
                                            outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + rSeq + "\t" + rSeq + var + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")


                    elif var is "-":
                        # # DELETION
                        # print("Deletion")
                        if chrom in regions:
                            start = start - 1
                            for s in regions[chrom].keys():  # Go through all start positions on the chromosome in the dictionary
                                if s <= start:  # Check if the region start is lower than the start of the variant
                                    for e in regions[chrom][s].keys():  # If so go through the end positions of the regions
                                        if e >= end:  # If the end of the region is greater than the variant end extract seq

                                            rSeq = regions[chrom][s][e]
                                            rSeq = rSeq[(start - s):(start - s + 1)]
                                            outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + rSeq + ref + "\t" + rSeq + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")


                    else:
                        if chrom in regions:
                            start = start - 1
                            for s in regions[chrom].keys():  # Go through all start positions on the chromosome in the dictionary
                                if s <= start:  # Check if the region start is lower than the start of the variant
                                    for e in regions[chrom][s].keys():  # If so go through the end positions of the regions
                                        if e >= end:  # If the end of the region is greater than the variant end extract seq

                                            rSeq = regions[chrom][s][e]
                                            rSeq = rSeq[(start - s):(start - s + 1)]
                                            outputfile.write (chrom + "\t" + str(start) + "\t" + id + "\t" + rSeq + ref + "\t" + rSeq + var + "\t" + qual + "\t" + filter + "\t" + info + "\tGT:DP:AD\t0/1:" + rd + ":" + refRd + "," + varRd + "\n")

                # If the file is not closed, close
            if not pindelfile.closed:
                pindelfile.close()
            # If the file is not closed, close
    if not outputfile.closed:
        outputfile.close()
