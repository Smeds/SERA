#!/usr/bin/python2.7

"""

Script for filtering and conveting pindel vcf output to annovar input

Elin Falk Sorqvist 20141114


"""

def extractInfo(information):
    # Extract info
    info = information.split(" ")
    return info[1]

import argparse
import re

# Parse commandline
parser = argparse.ArgumentParser(description = "This script takes vcf from pindel2vcf, filters it and creates an annovar input file")
parser.add_argument('-i', '--infile', help = 'Input file name', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
parser.add_argument('-m', '--minRD', help = 'Minimum read depth for a variant to be included', type = int, default = 1)
parser.add_argument('-r', '--readMethod', help = 'Use the maximum or minimum read depth as reference read depth', type = str, required = True, choices = ("min", "max"))
parser.add_argument('-v', '--vaf', help = 'Minimum variant allele ratio for a variant to be included', type = float, default = 0)
parser.add_argument('-del', '--deletion', help = 'Deletion file from pindel', type = str, required = True)
parser.add_argument('-ins', '--insertion', help = 'Insertion file from pindel', type = str, required = True)


args = parser.parse_args()
sample = "-"

deletions = {}
insertions = {}

with open(args.deletion, 'r') as deletionFile:
    # Go through the file line by line
    for line in deletionFile:

        line = line.rstrip('\r\n')
        line = line.split("\t")  # Remove new line character and split on tab
        # Extract chr, start, end
        chr = extractInfo(line[3])
        start = extractInfo(line[4])
        end = int(line[5]) - 1
        insertInfo = line[2].split(" ")
        insert = insertInfo[2].replace('""', "-")
        insert = insert.replace('"', "")
        readInfo = line[19].split(" ")
        upper = int(readInfo[1])
        lower = int(readInfo[2])
        refRd = upper
        varPlus = int(readInfo[3])
        varMinus = int(readInfo[5])
        # Check if the minimum or maximum reference read depth should be reported and save the correct read depth in the refRd variable
        if re.match('min', args.readMethod.lower()):
            if lower < refRd:
                refRd = lower

        else:
            if lower > refRd:
                refRd = lower
        # Create dictionary
        if not chr in deletions:
            deletions[chr] = {}
            deletions[chr][start] = {}
            deletions[chr][start][end] = {}
            deletions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
        else:
            if not start in deletions[chr]:
                deletions[chr][start] = {}
                deletions[chr][start][end] = {}
                deletions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
            else:
                if not end in deletions[chr][start]:
                    deletions[chr][start][end] = {}
                    deletions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
                else:
                        if not insert in deletions[chr][start][end]:
                            deletions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
                        else:
                            print (str(chr) + "\t" + str(start) + "\t" + str(end) + "\t" + str(insert) + " already exists in deletion dictionary!\n")

    if not deletionFile.closed:
        deletionFile.close()

with open(args.insertion, 'r') as insertionFile:
    # Go through the file line by line
    for line in insertionFile:
        line = line.rstrip('\r\n')
        line = line.split("\t")  # Remove new line character and split on tab
        # Extract chr, start, end
        chr = extractInfo(line[3])
        start = extractInfo(line[4])
        end = int(line[5]) - 1
        insertInfo = line[2].split(" ")
        insert = insertInfo[2].replace('""', "-")
        insert = insert.replace('"', "")
        readInfo = line[19].split(" ")
        upper = int(readInfo[1])
        lower = int(readInfo[2])
        refRd = upper
        varPlus = int(readInfo[3])
        varMinus = int(readInfo[5])

        # Check if the minimum or maximum reference read depth should be reported and save the correct read depth in the refRd variable
        if re.match('min', args.readMethod.lower()):
            if lower < refRd:
                refRd = lower

        else:
            if lower > refRd:
                refRd = lower

        # Create dictionary
        if not chr in insertions:
            insertions[chr] = {}
            insertions[chr][start] = {}
            insertions[chr][start][end] = {}
            insertions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
        else:
            if not start in insertions[chr]:
                insertions[chr][start] = {}
                insertions[chr][start][end] = {}
                insertions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
            else:
                if not end in insertions[chr][start]:
                    insertions[chr][start][end] = {}
                    insertions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
                else:
                        if not insert in insertions[chr][start][end]:
                            insertions[chr][start][end][insert] = {'refRd': refRd, 'varPlus' : varPlus, 'varMinus' : varMinus}
                        else:
                            print (str(chr) + "\t" + str(start) + "\t" + str(end) + "\t" + str(insert) + " already exists in insertion dictionary!\n")

    if not insertionFile.closed:
        insertionFile.close()

with open(args.infile, 'r') as infile:
    with (open(args.output, mode = 'w')) as outfile:

        # Go through the file line by line
        for line in infile:
            # Check so the line isn't empty or starts with ##
            if not re.match('^##', line) and not re.match('^$', line):
                line = line.rstrip('\r\n')  # Remove new line character
                lineSplit = line.split("\t")  # Split line by tab

                # If the line start with only 1 # extract the sample name from the last column
                if re.match("^#CHROM", line):
                    sample = lineSplit[9]

                # Else extract and print all relevant info
                else:
                    chrom = lineSplit[0]
                    startVcf = lineSplit[1]
                    # Split the information column
                    infoSplit = lineSplit[7].split(";")
                    end = 0
                    svtype = "-"
                    varRd = 0
                    refRd = 0
                    varPlusRd = 0
                    varMinusRd = 0
                    var = "-"
                    ref = "-"
                    strandInfo = "-"
                    start = lineSplit[1]

                    # Go through all information in the information column and extract the relevant parts
                    for info in infoSplit:
                        if re.match("^END", info):
                            tmpSplit = info.split("=")
                            end = int(tmpSplit[1])
                        elif re.match("^SVTYPE", info):
                            tmpSplit = info.split("=")
                            svtype = tmpSplit[1]


                    if re.match("DEL", svtype):
                        start = int(lineSplit[1]) + 1
                        ref = lineSplit[3][1:]
                        var = "-"

                        if chrom in deletions and startVcf in deletions[chrom] and end in deletions[chrom][startVcf] and var in deletions[chrom][startVcf][end]:
                            varPlusRd = int(deletions[chrom][startVcf][end][var]['varPlus'])
                            varMinusRd = int(deletions[chrom][startVcf][end][var]['varMinus'])
                            varRd = int(deletions[chrom][startVcf][end][var]['varPlus']) + int(deletions[chrom][startVcf][end][var]['varMinus'])
                            refRd = deletions[chrom][startVcf][end][var]['refRd']
                            strandInfo = "Tumor_Del=+" + str(varPlusRd) + "|-" + str(varMinusRd)

                        else:
                            print (str(chrom) + "\t" + str(startVcf) + "\t" + str(end) + "\t" + str(var) + " DOES NOT exist in deletion dictionary!\n")
                    elif  re.match("RPL", svtype):
                        start = int(lineSplit[1]) + 1
                        ref = lineSplit[3][1:]
                        var = lineSplit[4][1:]
                        if chrom in deletions and startVcf in deletions[chrom] and end in deletions[chrom][startVcf] and var in deletions[chrom][startVcf][end]:
                            varPlusRd = int(deletions[chrom][startVcf][end][var]['varPlus'])
                            varMinusRd = int(deletions[chrom][startVcf][end][var]['varMinus'])
                            varRd = int(deletions[chrom][startVcf][end][var]['varPlus']) + int(deletions[chrom][startVcf][end][var]['varMinus'])
                            refRd = deletions[chrom][startVcf][end][var]['refRd']
                            strandInfo = "Tumor_Del=+" + str(varPlusRd) + "|-" + str(varMinusRd)
                        else:
                            print (str(chrom) + "\t" + str(startVcf) + "\t" + str(end) + "\t" + str(var) + " DOES NOT exist in deletion dictionary!\n")

                    elif re.match("INS", svtype):
                        ref = "-"
                        var = lineSplit[4][1:]
                        if chrom in insertions and startVcf in insertions[chrom] and end in insertions[chrom][startVcf] and var in insertions[chrom][startVcf][end]:
                            varPlusRd = int(insertions[chrom][startVcf][end][var]['varPlus'])
                            varMinusRd = int(insertions[chrom][startVcf][end][var]['varMinus'])
                            varRd = int(insertions[chrom][startVcf][end][var]['varPlus']) + int(insertions[chrom][startVcf][end][var]['varMinus'])
                            refRd = insertions[chrom][startVcf][end][var]['refRd']
                            strandInfo = "Tumor_Ins=+" + str(varPlusRd) + "|-" + str(varMinusRd)
                        else:
                            print (str(chrom) + "\t" + str(startVcf) + "\t" + str(end) + "\t" + str(var) + " DOES NOT exist in insertion dictionary!\n")

                    alleleFreq = str(refRd) + "," + str(varRd)

                    # Calculate tot read depth and the variant allele ratio
                    readDepth = int(varRd) + int(refRd)
                    # vRatio = 0
                    # if readDepth > 0:
                    vRatio = float(varRd) / float(readDepth)

                    # Check that the variant allele ratio and the read depth are high enough
                    if vRatio >= args.vaf and readDepth >= args.minRD:

                        # Check that the chromosome annotation doesn't contain chr or NC notation
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
                        outfile.write(str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(var) + "\t" + "comments: sample=" + sample + " variantAlleleRatio=" + str(vRatio) + " alleleFreq=" + str(alleleFreq) + " readDepth=" + str(readDepth) + " " + strandInfo + " Tumor_var_plusAmplicons=- Tumor_var_minusAmplicons=- Tumor_ref_plusAmplicons=- Tumor_ref_minusAmplicons=-\n")

        # Check if outfile is closed, if not close it
    if not infile.closed:
        infile.close()

    # Check if outfile is closed, if not close it
    if not outfile.closed:
        outfile.close()
