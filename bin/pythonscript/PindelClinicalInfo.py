#!/usr/bin/python2.7

import csv
import re
import argparse
from re import split

'''

Script for converting Pindel Annovar output to readable output

'''

def printToFile(args, RDs, chr2nc, ofile, row, cosmicID, tInfo, cdsInfo, aaInfo, stPlus, stMinus, allTrans):
    # Print to file
    flag = 2
    if row[6].startswith("ERBB2"):
        flag = 3
    ofile.write(args['sampleID'] + "\t" + row[6] + "\t" + cosmicID + "\t" + str(flag))  # Go through all read depth and check if the total read depth is above
    count = 1
    for rd in RDs:
        if int(row[12]) >= int(rd):
            if count == 1:
                ofile.write("\tyes\tok")
            else:
                ofile.write("\tok")
        elif count == 1:
            ofile.write("\tnot analyzable\tlow")
        else:
            ofile.write("\tlow")
        count += 1
        varType = "-"
        if re.match("-", row[8]):
            varType = row[7]
        else:
            varType = row[8]

    chr = "chr" + row[1]
    ofile.write("\t" + row[12] + "\t" + chr2nc[chr] + "\t" + row[2] + "\t" + row[3] + "\t" + row[4] + "\t" + row[5] + "\t" + row[11] + "\t" + row[10] + "\t" + row[9] + "\t" + row[13] + "\t" + varType + "\t" + cdsInfo + "\t" + aaInfo + "\t-\t-\t" + str(stPlus) + "\t" + str(stMinus) + "\t" + tInfo + "\t-\t-\t" + str(stPlus) + "|" + str(stMinus) + "\t" + str(allTrans) + "\n")



parser = argparse.ArgumentParser(description = 'Script for converting Pindel Annovar output to readable output')

parser.add_argument('-i', '--inputfile', help = 'Input txt file (Required)', required = True)
parser.add_argument('-o', '--outputfile', help = 'Output txt file (Required)', required = True)
parser.add_argument('-s', '--sampleID', help = 'Sample ID (Required)', required = True)
parser.add_argument('-g', '--genefile', help = 'File with the genes we want to report indels for (Required)', required = True)
parser.add_argument('-minRD', '--minRD', help = 'Minimum read depths to be checked for a position, several RD can be given separated by comma (-minRD 30,50,100) [optional, default is 1]', default = 1)
parser.add_argument('-minVarRatio', '--minVarRatio', help = 'Minimum variant allele ratio for a position to be considered', required = True)
parser.add_argument('-chr2nc', '--chr2nc', help = 'File with chr in first column and NC-number in second column', required = True)

args = vars(parser.parse_args())


# Create list with the read depth to consider
RDs = args['minRD'].rstrip('\r\n').split(",")

# Open gene file and create a list with the genes to consider
genes = {}
with open(args['genefile'], 'r') as gfile:
    # Go through the infile line by line
    for line in gfile:
        line = line.rstrip('\r\n').split("\t")  # Remove new line character

        genes[line[0].lower()] = line[1]
        # genes.append(line.lower())

    # If the gene file is not closed, close
    if not gfile.closed:
        gfile.close()

chr2nc = {}
with open(args['chr2nc'], 'r') as cnfile:
    # Go through the file line by line
    for line in cnfile:
        if line.startswith("chr"):
            line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
            chr2nc[line[0]] = line[1]

    # If the file is not closed, close
    if not cnfile.closed:
        cnfile.close()

with open (args['outputfile'], 'w') as ofile:
    with open(args['inputfile'], 'r') as ifile:
    # Print header
#         ofile.write("#Sample\tGene\tCosmic_id\tReport\tFound")
#         for rd in RDs:
#             ofile.write("\t"+rd)
#         ofile.write("\tTotal_read_depth\tChr\tStart\tEnd\tReference\tVariant\tVariant_read_depth\tReference_read_depth\tVariant_allele_ratio\tExonic_type\tCDS_change\tAA_change\tReference_plus_amplicons\tReference_minus_amplicons\tVariant_plus_amplicons\tVariant_minus_amplicons\tTranscript\tProtein\tRef_amplicons\tVar_amplicons\n")
        # Read every row in inputfile
        for line in ifile:
            line = line.rstrip('\r\n')
            row = line.split("\t")
            done = 0
            if not row[0].startswith("Sample"):  # Check that line does not start with Sample
                if row[6].lower() in genes.keys():
                    if float(row[9]) >= float(args['minVarRatio']):  # Check that the variant allele ratio is above given value
                        # Extract cosmic id
                        cid = row[17].split(";")
                        cID = cid[0].split(",")
                        cosmicID = cID[0].replace("ID=COSM", "")

                        # Extract info about variant strands
                        stPlus = "-"
                        stMinus = "-"
                        # If insertion
                        if re.match("-", row[4]):
                            stInfo = row[24].split("|")
                            stPlus = stInfo[0]
                            stMinus = stInfo[1]
                        else:
                            stInfo = row[25].split("|")
                            stPlus = stInfo[0]
                            stMinus = stInfo[1]

                        allTrans = row[32]
                        allTrans = re.sub('^[A-Z0-9]+:', "", allTrans)
                        for i in range(33, len(row)):
                            transInfo = row[i]
                            transInfo = re.sub('^[A-Z0-9]+:', "", transInfo)
                            allTrans += "|" + transInfo

                        # Extract info about changes
                        transriptInfo = row[32].split(":")
                        # transriptInfo = row[26].split(":")
                        cdsInfo = "-"
                        aaInfo = "-"
                        tInfo = "-"

                        # Check if cds and aa info exists
                        if len(transriptInfo) >= 4:
                            cdsInfo = transriptInfo[3]
                            tInfo = transriptInfo[1]
                            if len(transriptInfo) >= 5:
                                aaInfo = transriptInfo[4]

                        if "exon" in genes[row[6].lower()]:
                            if re.match("exonic", row[7].lower()):  # Check if the gene exist in the list of genes to be considered and that the variant is exonic

                                exonInfo = transriptInfo[2] + ";"

                                if exonInfo in genes[row[6].lower()]:  # Check if the variant is located within the correct exons
                                    if done == 0:
                                        printToFile(args, RDs, chr2nc, ofile, row, cosmicID, tInfo, cdsInfo, aaInfo, stPlus, stMinus, allTrans)
                                        done = 1

                        if re.match("all", genes[row[6].lower()]) or re.match("-", genes[row[6].lower()]):
                            if done == 0:
                                printToFile(args, RDs, chr2nc, ofile, row, cosmicID, tInfo, cdsInfo, aaInfo, stPlus, stMinus, allTrans)
                                done = 1

                        if ":" in genes[row[6].lower()]:
                            reg = genes[row[6].lower()].split(";")
                            regions = {}
                            for r in reg:
                                if ":" in r:
                                    chrInfo = r.split(":")
                                    positions = chrInfo[1].split("-")
                                    if not chrInfo[0] in regions:
                                        regions[chrInfo[0]] = {}
                                        regions[chrInfo[0]][positions[0]] = {}
                                        regions[chrInfo[0]][positions[0]][positions[1]] = 1
                                    else:
                                        if not positions[0] in regions[chrInfo[0]]:
                                            regions[chrInfo[0]][positions[0]] = {}
                                            regions[chrInfo[0]][positions[0]][positions[1]] = 1

                                        else:
                                            if not positions[1] in regions[chrInfo[0]][positions[0]]:
                                                regions[chrInfo[0]][positions[0]][positions[1]] = 1

                                if row[1] in regions:
                                    for start in regions[row[1]]:
                                        if int(row[2]) >= int(start):
                                            for end in regions[row[1]][start]:
                                                if int(row[3]) <= end:
                                                    if done == 0:
                                                        printToFile(args, RDs, chr2nc, ofile, row, cosmicID, tInfo, cdsInfo, aaInfo, stPlus, stMinus, allTrans)
                                                        done = 1






        if not ifile.closed:
            ifile.close()
    if not ofile.closed:
        ofile.close()
