#!/usr/bin/python2.7

import csv
import re
import argparse

'''

Script for converting Pindel Annovar output to readable output


'''

parser = argparse.ArgumentParser(description='Script for converting Pindel Annovar output to readable output')

parser.add_argument('-i', '--inputfile', help='Input txt file (Required)', required=True)
parser.add_argument('-o', '--outputfile', help='Output txt file (Required)', required=True)
parser.add_argument('-s', '--sampleID', help='Sample ID (Required)', required=True)
parser.add_argument('-g', '--genefile', help='File with the genes we want to report indels for (Required)', required=True)
parser.add_argument('-minRD', '--minRD', help='Minimum read depths to be checked for a position, several RD can be given separated by comma (-minRD 30,50,100) [optional, default is 1]', default=1)
parser.add_argument('-minVarRatio', '--minVarRatio', help='Minimum variant allele ratio for a position to be considered', required=True)
parser.add_argument('-chr2nc', '--chr2nc', help='File with chr in first column and NC-number in second column', required=True)

args = vars(parser.parse_args())


# Create list with the read depth to consider
RDs = args['minRD'].rstrip('\r\n').split(",")

# Open gene file and create a list with the genes to consider
genes = []
with open(args['genefile'], 'r') as gfile:
    # Go through the infile line by line
    for line in gfile:
        line = line.rstrip('\r\n')  # Remove new line character
        genes.append(line.lower())

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
            if not row[0].startswith("Sample"):  # Check that line does not start with Sample
                if row[6].lower() in genes and re.match("exonic", row[7].lower()):  # Check if the gene exist in the list of genes to be considered and that the variant is exonic
                    if float(row[9]) >= float(args['minVarRatio']):  # Check that the variant allele ratio is above given value
                        # Extract cosmic id
                        cid = row[17].split(";")
                        cID = cid[0].split(",")
                        cosmicID = cID[0].replace("ID=COSM", "")

                        # Extract info about changes
                        transriptInfo = row[24].split(":")
                        cdsInfo = "-"
                        aaInfo = "-"
                        # Check if cds and aa info exists
                        if len(transriptInfo) >= 4:
                            cdsInfo = transriptInfo[3]
                            if len(transriptInfo) >= 5:
                                aaInfo = transriptInfo[4]

                        if ok == 1:
                            # Print to file
                            ofile.write(args['sampleID'] + "\t" + row[6] + "\t" + cosmicID + "\t3")
                            # Go through all read depth and check if the total read depth is above
                            count = 1
                            for rd in RDs:
                                if int(row[12]) >= int(rd):
                                    if count == 1:
                                        ofile.write("\tyes\tok")
                                    else:
                                        ofile.write("\tok")
                                else:
                                    if count == 1:
                                        ofile.write("\tnot analyzable\tlow")
                                    else:
                                        ofile.write("\tlow")
                                count += 1
                            chr = "chr" + row[1]
                            ofile.write("\t" + row[12] + "\t" + chr2nc[chr] + "\t" + row[2] + "\t" + row[3] + "\t" + row[4] + "\t" + row[5] + "\t" + row[11] + "\t" + row[10] + "\t" + row[9] + "\t" + row[8] + "\t" + cdsInfo + "\t" + aaInfo + "\t-\t-\t-\t-\t" + transriptInfo[1] + "\t-\t-\t-\n")

        if not ifile.closed:
            ifile.close()
    if not ofile.closed:
        ofile.close()
