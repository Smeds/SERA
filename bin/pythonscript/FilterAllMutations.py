#!/usr/bin/python2.7

from variant_functions import *
from print_functions import *
import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a region file, a variation file from SNPmania and an filtered annovaroutput file and outputs the coverage per base in the regions")
parser.add_argument('-i', '--hotspotfile', help = 'Name of the file with hotspot mutations to report in bed-format', type = str, required = True)
parser.add_argument('-v', '--variationfile', help = 'File name of variation file from SNPmania', type = str, required = True)
parser.add_argument('-f', '--annovarfile', help = 'File name of the annovar output file', type = str, required = True)
parser.add_argument('-m', '--multipleBpfile', help = 'File name of the file with annotation of multiple bp variants', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-minRD', '--minRD', help = 'Minimum read depth to be counted as covered', type = str, required = True)
parser.add_argument('-s', '--sample', help = 'Sample ID', type = str, required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'File with conversion between NC-number and chr', type = str, required = True)
parser.add_argument('-g', '--genome1000', help = 'Set the maximum allowed frequency for a variant to be kept, default 0', default = 0, type = float, required = False)
parser.add_argument('-b', '--blacklistfile', help = 'File with blacklisted variants (Required)', required = True)
parser.add_argument('-a', '--ampliconmapped', help = 'Set if data is ampliconmapped. Default: False', action = "store_true")
parser.add_argument('-t', '--transcript', help = 'File with the \"main\" transcript for each gene.', required = True, type = str)
parser.add_argument('-minVaf', '--minVaf', help = 'Minimum variant allele frequency for a variant to be included in the output', type = float, required = True)

args = parser.parse_args()


# Split the given min read depths
minRDs = args.minRD.split(",")

# Create dictionary for convertion between nc and chr
nc2chr = {}
with open(args.chr2nc, 'r') as cnfile:
    # Go through the file line by line
    for line in cnfile:
        if line.startswith("chr"):
            line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
            nc2chr[line[0]] = line[1]  # create the dictionary with nc as key and chr as value
    # If the file is not closed, close
    if not cnfile.closed:
        cnfile.close()

# Create hash with all hotspot positions we want to report
hotspots = {}
intronic = {}
# Open hotspot file and add info about the hotspot positions to the hash
with open(args.hotspotfile, 'r') as hotspotFile:
    # Go through the file line by line
    for line in hotspotFile:

        # Check so the line isn't empty or starts with #
        if not re.match('^#', line) and not re.match('^$', line) and not line.startswith('Chr'):

            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            # Get chromosome
            # chrom = lineSplit[1]
            if lineSplit[0].startswith("chrom"):
                lineSplit[0] = nc2chr[lineSplit[0]]
            elif not lineSplit[0].startswith("NC") and not lineSplit[0].startswith("chr"):
                print(line)
                lineSplit[0] = nc2chr["chr" + lineSplit[1]]

            # Create a hash with the mutations that are ok
            createHotspotHash(lineSplit, hotspots, intronic)

    # Check if the filteredAnnovarfile is closed, if not close it
    if not hotspotFile.closed:
        hotspotFile.close()

# Create a blacklist with the variants (artifacts etc) to remove
blacklist = {}
with open(args.blacklistfile, 'r') as blackfile:
    for line in blackfile:
        if not re.match('^Chr', line):
            createBlacklist(line, blacklist,nc2chr)

    if not blackfile.closed:
        blackfile.close()

# Create a hash with the main transcripts of interest
transcripts = {}
with open(args.transcript, 'r') as transcriptFile:

    for line in transcriptFile:
        # Check so the line isn't empty or starts with # or Gene
        if not re.match('^#', line) and not re.match('^$', line) and not line.startswith('Gene'):
            lineSplit = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
            accNumber = lineSplit[1].split(".")[0]
            if not lineSplit[0] in transcripts:
                transcripts[lineSplit[0]] = accNumber

    if not transcriptFile.closed:
        transcriptFile.close()

# Create a hash with annotation for important multiple bp variants
multipleBp = {}
with open(args.multipleBpfile, 'r') as multipleBpFile:

    for line in multipleBpFile:
        # Check so the line isn't empty or starts with # or Gene
        if not re.match('^#', line) and not re.match('^$', line) and not line.startswith('Chr'):
            lineSplit = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab

            createMultipleBpHash(lineSplit, multipleBp,)
    if not multipleBpFile.closed:
        multipleBpFile.close()

with open(args.output, 'w') as outputFile:
    outputFile.write("Sample\tGene\tVariant_type\tExon\tAA_change\tCDS_change\tAccession_number\tComment\tReport\tFound\tMin_read_depth300\tTotal_read_depth\t",)
    outputFile.write("Reference_read_depth\tVariant_read_depth\tVariant_allele_ratio\tdbSNP_id\tRatio_1000G\tRatio_ESP6500\tClinically_flagged_dbSNP\t",)
    outputFile.write("Cosmic\tClinVar_CLNDB\tClinval_CLINSIG\tReference_plus_amplicons\tReference_minus_amplicons\tVariant_plus_amplicons\tVariant_minus_amplicons\t",)
    outputFile.write("Strands_A_F+F-S+S-\tStrands_G_F+F-S+S-\tStrands_C_F+F-S+S-\tStrands_T_F+F-S+S-\tStrands_Ins\tStrands_Del\tRef_aligned_amplicons\t",)
    outputFile.write("Var_aligned_amplicons\tChr\tStart\tEnd\tReference_base\tVariant_base\tAll_transcripts_annotation\n")

    # Open filteredAnnovar file and add info about the filteredAnnovars to the hash
    with open(args.annovarfile, 'r') as annovarFile:
        # Go through the file line by line
        for line in annovarFile:

            # Check so the line isn't empty or starts with #
            if not re.match('^#', line) and not re.match('^$', line) and not line.startswith('Sample'):

                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                # Get chromosome
                if lineSplit[1].startswith("chr"):
                    lineSplit[1] = nc2chr[lineSplit[1]]
                elif not lineSplit[1].startswith("NC") and not lineSplit[1].startswith("chr"):
                    lineSplit[1] = nc2chr["chr" + lineSplit[1]]

                # Filter annovar output filesCreate a hash with the mutations that are ok
                if not addVariantInfo(lineSplit, minRDs, blacklist, args.genome1000, args.ampliconmapped, hotspots, intronic, args.minVaf):
                    aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(lineSplit, transcripts)  # Get transcript information
                    # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                    found, readLevel = getReadLevel(minRDs, int(lineSplit[12]))  # Get the level of read depth and if it's not analyzable

                    if re.match("-", found):
                        found = "yes"
                    if not re.match("altTranscript", comm):
                        comm = "-"

                    refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(lineSplit, args.ampliconmapped)  # add amplicon information
                    if ampliconCheckNoneHotspot(lineSplit, args.ampliconmapped):
                        outputFile.write (str(lineSplit[0]) + "\t" + str(lineSplit[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comm + "\t4-other\t" + found + "\t" + readLevel + "\t" + str(lineSplit[12]) + "\t" + str(lineSplit[10]) + "\t" + str(lineSplit[11]) + "\t" + str(lineSplit[9]) + "\t" + str(lineSplit[14]) + "\t" + str(lineSplit[13]) + "\t" + str(lineSplit[16]) + "\t" + str(lineSplit[15]) + "\t" + str(lineSplit[17]) + "\t" + str(lineSplit[18]) + "\t" + str(lineSplit[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(lineSplit[20]) + "\t" + str(lineSplit[21]) + "\t" + str(lineSplit[22]) + "\t" + str(lineSplit[23]) + "\t" + str(lineSplit[24]) + "\t" + str(lineSplit[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(lineSplit[1]) + "\t" + str(lineSplit[2]) + "\t" + str(lineSplit[3]) + "\t" + str(lineSplit[4]) + "\t" + str(lineSplit[5]) + "\t" + str(lineSplit[32]) + "\n")


        # Check if the filteredAnnovarfile is closed, if not close it
        if not annovarFile.closed:
            annovarFile.close()
    intronic = {}  # empty intronic hash

    # Open variation file from SNPmania
    with open(args.variationfile, 'r') as variationfile:
    # Go through the file line by line
        for line in variationfile:
            # Check so the line isn't empty or starts with #
            if not re.match('^#', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                chrom = lineSplit[3]
                position = int(lineSplit[4])

                if chrom.startswith("chr"):
                    chrom = nc2chr[chrom]

                # Go through region and add read depth to positions without variant
                if "indel" in hotspots:
                    if chrom in hotspots['indel']:
                        for start in hotspots['indel'][chrom]:
                            if position >= int(start):
                                for end in hotspots['indel'][chrom][start]:
                                    if position <= int(end):
                                        if isinstance(hotspots['indel'][chrom][start][end]['rd'], list):  # if a list already exist only add read depth otherwise create list first
                                            hotspots['indel'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....
                                        else:
                                            rds = ["-"] * (int(end) - int(start) + 1)  # create array to add read depth in
                                            hotspots['indel'][chrom][start][end]['rd'] = rds  # Adding list to rd
                                            hotspots['indel'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....

                # Go through region and add read depth to positions without variant
                if "region" in hotspots:
                    if chrom in hotspots['region']:
                        for start in hotspots['region'][chrom]:
                            if position >= int(start):
                                for end in hotspots['region'][chrom][start]:
                                    if position <= int(end):
                                        if isinstance(hotspots['region'][chrom][start][end]['rd'], list):  # if a list already exist only add read depth otherwise create list first
                                            hotspots['region'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....
                                        else:
                                            rds = ["-"] * (int(end) - int(start) + 1)  # create array to add read depth in
                                            hotspots['region'][chrom][start][end]['rd'] = rds  # Adding list to rd
                                            hotspots['region'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....

                # Go through region_all and add read depth to positions without variant
                if "region_all" in hotspots:
                    if chrom in hotspots['region_all']:
                        for start in hotspots['region_all'][chrom]:
                            if position >= int(start):
                                for end in hotspots['region_all'][chrom][start]:
                                    if position <= int(end):
                                        if isinstance(hotspots['region_all'][chrom][start][end]['rd'], list):  # if a list already exist only add read depth otherwise create list first
                                            hotspots['region_all'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....
                                        else:
                                            rds = ["-"] * (int(end) - int(start) + 1)  # create array to add read depth in
                                            hotspots['region_all'][chrom][start][end]['rd'] = rds  # Adding list to rd
                                            hotspots['region_all'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....

                # Go through hotspot and add read depth to positions without variant
                if "hotspot" in hotspots:
                    if chrom in hotspots['hotspot']:
                        for start in hotspots['hotspot'][chrom]:
                            if position >= int(start):
                                for end in hotspots['hotspot'][chrom][start]:
                                    if position <= int(end):
                                        if isinstance(hotspots['hotspot'][chrom][start][end]['rd'], list):  # if a list already exist only add read depth otherwise create list first
                                            hotspots['hotspot'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....
                                        else:
                                            rds = ["-"] * (int(end) - int(start) + 1)  # create array to add read depth in
                                            hotspots['hotspot'][chrom][start][end]['rd'] = rds  # Adding list to rd
                                            hotspots['hotspot'][chrom][start][end]['rd'][(position - start)] = lineSplit[0]  # Adding total read depth might want to change this to reference read depth....

    if not variationfile.closed:
        variationfile.close()

    # Print hotspot positions
    if "hotspot" in hotspots:
        printHotspots(hotspots, minRDs, args.sample, transcripts, outputFile, args.ampliconmapped, multipleBp)

    if "region" in hotspots:
        printRegion(hotspots, minRDs, args.sample, transcripts, outputFile, args.ampliconmapped, multipleBp)

    if "indel" in hotspots:
        printIndel(hotspots, minRDs, args.sample, transcripts, outputFile, args.ampliconmapped, multipleBp)

    if "region_all" in hotspots:
        printRegionAll(hotspots, minRDs, args.sample, transcripts, outputFile, args.ampliconmapped, multipleBp)

    if not outputFile.closed:
        outputFile.close()
