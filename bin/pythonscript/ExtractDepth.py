#!/usr/bin/python2.7

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a region file, a variation file from SNPmania and an filtered annovaroutput file and outputs the coverage per base in the regions")
parser.add_argument('-r', '--regionfile', help = 'Name of the input region file in bed-format', type = str, required = True)
parser.add_argument('-v', '--variationfile', help = 'File name of variation file from SNPmania', type = str, required = True)
parser.add_argument('-f', '--filteredAnnovarfile', help = 'File name of the filtered annovar output file', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-minRD', '--minRD', help = 'Minimum read depth to be counted as covered', type = str, required = True)
parser.add_argument('-s', '--sample', help = 'Sample ID', type = str, required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'File with conversion between NC-number and chr', type = str, required = True)
parser.add_argument('-header', '--header', help = 'Set if you want header printed. Default: False', action = "store_true")


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

genes = {}
# Open region file and add info about the regions to the hash
with open(args.regionfile, 'r') as regionfile:
#    print ("Hashing genes!")
    # Go through the file line by line
    for line in regionfile:
        # Check so the line isn't empty or starts with #
        if not re.match('^#', line) and not re.match('^$', line):
            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            chrom = lineSplit[0]
            if chrom.startswith("chr"):
                chrom = nc2chr[chrom]

            start = int(lineSplit[1]) + 1  # since it's a bed-file
            end = int(lineSplit[2])
            info = lineSplit[3]

            # Check if the tuple exists as key in dict genes
            # If not add it
            # If it does output error message on duplicate genes
            if not chrom in genes:
                genes[chrom] = {}
                genes[chrom][start] = {}
                genes[chrom][start][end] = info
            else:
                if not start in genes[chrom]:
                    genes[chrom][start] = {}
                    genes[chrom][start][end] = info
                else:
                    print ("Error gene with coordinates ", gene, " ", chrom, " : ", start, " already exists!")

    # Check if the regionfile is closed, if not close it
    if not regionfile.closed:
        regionfile.close()


mutations = {}
# Open filteredAnnovar file and add info about the filteredAnnovars to the hash
with open(args.filteredAnnovarfile, 'r') as filteredAnnovar:
    # Go through the file line by line
    for line in filteredAnnovar:
        # Check so the line isn't empty or starts with #
        if not re.match('^#', line) and not re.match('^$', line) and not line.startswith('Sample'):

            line = line.rstrip('\r\n')  # Remove new line character
            lineSplit = line.split("\t")  # Split line by tab

            # Get chromosome
            chrom = lineSplit[1]
            if chrom.startswith("chr"):
                chrom = nc2chr[chrom]
            elif not chrom.startswith("NC") and not chrom.startswith("chr"):
                chrom = nc2chr["chr" + chrom]

            start = int(lineSplit[2])
            ref = lineSplit[4]
            var = lineSplit[5]

            # Add all known mutations to the mutation dictionary which are NOT from Pindel

            if not chrom in mutations:
                mutations[chrom] = {}
                mutations[chrom][start] = {}
                mutations[chrom][start][ref] = {}
                mutations[chrom][start][ref][var] = {}
                if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):
                    mutations[chrom][start][ref][var]["bwa"] = lineSplit
                else:
                    mutations[chrom][start][ref][var]["pindel"] = lineSplit

            else:
                if not start in mutations[chrom]:
                    mutations[chrom][start] = {}
                    mutations[chrom][start][ref] = {}
                    mutations[chrom][start][ref][var] = {}
                    if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):
                        mutations[chrom][start][ref][var]["bwa"] = lineSplit
                    else:
                        mutations[chrom][start][ref][var]["pindel"] = lineSplit


                else:
                    if not ref in mutations[chrom][start]:
                        mutations[chrom][start][ref] = {}
                        mutations[chrom][start][ref][var] = {}
                        if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):
                            mutations[chrom][start][ref][var]["bwa"] = lineSplit
                        else:
                            mutations[chrom][start][ref][var]["pindel"] = lineSplit


                    else:
                        if not var in mutations[chrom][start][ref]:
                            mutations[chrom][start][ref][var] = {}
                            if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):
                                mutations[chrom][start][ref][var]["bwa"] = lineSplit
                            else:
                                mutations[chrom][start][ref][var]["pindel"] = lineSplit
                        else:
                            if not ("+" in lineSplit[24]) and not ("+" in lineSplit[25]):
                                mutations[chrom][start][ref][var]["bwa"] = lineSplit
                            else:
                                mutations[chrom][start][ref][var]["pindel"] = lineSplit


    # Check if the filteredAnnovarfile is closed, if not close it
    if not filteredAnnovar.closed:
        filteredAnnovar.close()


# Open variation file from SNPmania
with open(args.variationfile, 'r') as variationfile:
    with open(args.output, 'w') as outputfile:

        if args.header:
            commentStr = "#Sample\tGene\tCosmic_id\tReport\tFound"
            for mRD in minRDs:
                commentStr += "\tMin_read_depth" + mRD
                commentStr += "\tTotal_read_depth\tChr\tStart\tEnd\tReference\tVariant\tVariant_read_depth\tReference_read_depth\tVariant_allele_ratio\tRatio_in_1000G\tExonic_type\tCDS_change\tAA_change\tReference_plus_amplicons\tReference_minus_amplicons\tVariant_plus_amplicons\tVariant_minus_amplicons\tTranscript\tProtein\tRef_amplicons\tVar_amplicons\tAll_transcripts"
                outputfile.write (commentStr + "\n")

        # Go through the file line by line
        for line in variationfile:
            # Check so the line isn't empty or starts with #
            if not re.match('^#', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                chrom = lineSplit[3]
                if chrom.startswith("chr"):
                    chrom = nc2chr[chrom]

                if chrom in genes:
                    pos = int(lineSplit[4])
                    for s in genes[chrom]:
                        if pos >= int(s):
                            for e in genes[chrom][s]:
                                if pos <= int(e):
                                    okDepth = ""
                                    depth = lineSplit[0]
                                    found = "no"
                                    flag = 1
                                    for mRD in minRDs:
                                        if int(lineSplit[0]) >= int(mRD):
                                            okDepth += "ok\t"
                                        else:
                                            okDepth += "low\t"
                                            # If the min read depth we are looking at is the lowest set it as not analyzable
                                            if mRD == minRDs[0]:
                                                found = "not analyzable"
                                    nm = "-"
                                    np = "-"
                                    if genes[chrom][s][e].startswith ("KIT"):
                                        nm = "NM_000222.2"
                                        np = "NP_000213.1"
                                    elif genes[chrom][s][e].startswith("PDGFRA"):
                                        nm = "NM_006206.4"
                                        np = "NP_006197.1"
                                    elif genes[chrom][s][e].startswith("EGFR"):
                                        nm = "NM_005228.3"
                                        np = "NP_005219.2"
                                        flag = 3
                                    elif genes[chrom][s][e].startswith("TP53"):
                                        nm = "NM_000546.5"
                                        np = "NP_000537.3"
                                        flag = 3
                                    elif genes[chrom][s][e].startswith("PIK3CA"):
                                        nm = "NM_006218.2"
                                        np = "NP_006209.2"
                                        flag = 3
                                    elif genes[chrom][s][e].startswith("ERBB2"):
                                        nm = "NM_001289937.1"
                                        np = "NP_001276866.1"
                                    elif genes[chrom][s][e].startswith("BRCA1"):
                                        nm = "NM_007294.3"
                                        np = "NP_009225.1"
                                        flag = 3
                                    elif genes[chrom][s][e].startswith("BRCA2"):
                                        nm = "NM_000059.3"
                                        np = "NP_000050.2"
                                        flag = 3
                                    cds = "-"
                                    aa = "-"


                                    # Check if variant in pindels hash
                                    if not "not analyzable" in found and chrom in mutations and pos in mutations[chrom]:
                                        for ref in mutations[chrom][pos]:
                                            for var in mutations[chrom][pos][ref]:
                                                if not "pindel" in mutations[chrom][pos][ref][var]:
                                                    mutationLine = mutations[chrom][pos][ref][var]["bwa"]

                                                    # Save all transcripts without gene name
                                                    allTrans = mutationLine[32]
                                                    allTrans = re.sub('^[A-Z0-9]+:', "", allTrans)
                                                    # Look through all reported transcripts and report the cds and aa change for the right transcript
                                                    for transcript in mutationLine[32:]:
                                                        # If it's not the first transcript add to all transcript column
                                                        if not re.match(transcript, mutationLine[32]):
                                                            transInfo = transcript
                                                            transInfo = re.sub('^[A-Z0-9]+:', "", transInfo)
                                                            allTrans += "|" + transInfo

                                                        nmNumber = nm.split(".")  # split on dot to exclude version number
                                                        if nmNumber[0] in transcript:  # Check that it is the correct transcript
                                                            transcript = transcript.split(":")
                                                            if len(transcript) >= 5:  # Check that both aa and cds info is available
                                                                cds = transcript[3]
                                                                aa = transcript[4]
                                                            elif len(transcript) >= 4:  # If not check if at least cds info is available
                                                                cds = transcript[3]

                                                    # Check read depth for mutation may differ from variation file if it is from pindel or have several variants

                                                    okDepth = ""
                                                    depth = mutationLine[12]
                                                    for mRD in minRDs:
                                                        if int(mutationLine[12]) >= int(mRD):
                                                            okDepth += "ok\t"
                                                        else:
                                                            okDepth += "low\t"
                                                            if mRD == minRDs[0]:
                                                                found = "not analyzable"

                                                    if "not analyzable" in found:
                                                        outStr = args.sample + "\t" + genes[chrom][s][e] + "\t-\t1\t" + found + "\t" + okDepth + lineSplit[0] + "\t" + chrom + "\t" + str(pos) + "\t" + str(pos) + "\t" + lineSplit[5] + "\tN\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + nm + "\t" + np + "\t-\t" + allTrans
                                                        outputfile.write(outStr + "\n")
                                                    else:
                                                        ampliconInfo = mutationLine[28] + "\t" + mutationLine[29] + "\t" + mutationLine[26] + "\t" + mutationLine[27]
                                                        ampliconColumns = mutationLine[31] + "\t" + mutationLine[30]

                                                        outStr = args.sample + "\t" + genes[chrom][s][e] + "\t-\t" + str(flag) + "\tyes\t" + okDepth + depth + "\t" + chrom + "\t" + str(pos) + "\t" + mutationLine[3] + "\t" + ref + "\t" + var + "\t" + mutationLine[11] + "\t" + mutationLine[10] + "\t" + mutationLine[9] + "\t" + mutationLine[13] + "\t" + mutationLine[8] + "\t" + cds + "\t" + aa + "\t" + ampliconInfo + "\t" + nm + "\t" + np + "\t" + ampliconColumns + "\t" + allTrans
                                                        outputfile.write(outStr + "\n")
                                    else:
                                        # Output all positions in KIT, PDGFRA, BRCA1 and BRCA2 independent of mutation or not
                                        if genes[chrom][s][e].startswith("KIT") or genes[chrom][s][e].startswith("PDGFRA") or genes[chrom][s][e].startswith("BRCA1") or genes[chrom][s][e].startswith("BRCA2"):
                                            outStr = args.sample + "\t" + genes[chrom][s][e] + "\t-\t" + str(flag) + "\t" + found + "\t" + okDepth + depth + "\t" + chrom + "\t" + str(pos) + "\t" + str(pos) + "\t" + lineSplit[5] + "\tN\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + nm + "\t" + np + "\t-\t-\t-"
                                            outputfile.write(outStr + "\n")
