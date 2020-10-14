#!/usr/bin/python2.7

import csv
import re
import argparse

'''

Script for converting Pindel Annovar output to readable output

Viktor Ljungstrom 2013-11-26

'''

parser = argparse.ArgumentParser(description = 'Script for converting Pindel Annovar output to readable output')

parser.add_argument('-i', '--inputfile', help = 'Input txt file (Required)', required = True)
parser.add_argument('-o', '--outputfile', help = 'Output txt file (Required)', required = True)
parser.add_argument('-s', '--sampleID', help = 'Sample ID (Required)', required = True)

args = vars(parser.parse_args())

ifile = open(args['inputfile'], "rb")
reader = csv.reader(ifile, delimiter = '\t')
ofile = open(args['outputfile'], "wb")

# Print header
ofile.write("#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tCosmic\tVariant Type\tTranscripts\n")

# Read every row in inputfile
for row in reader:

	if not row[0].startswith("Chr"):

		info = row[22].replace(":", ",").split(",")
		# print  info[2].split('=')[1]
		sample = args['sampleID']
		variantratio = int(info[2]) / (int(info[1]) + int(info[2]))
		readdepth = (int(info[1]) + int(info[2]))
		transcripts = row[8]
		# transcripts = row[8].split(',')
		variant_type = row[20].split(";")

		# Get variant type
		for i in variant_type:
			if i.startswith("SVTYPE"):
				variant_type_exact = i

		# Print to file
		ofile.write(sample + "\t" + row[0] + "\t" + row[1] + "\t" + row[2] + "\t" + row[3] + "\t" + row[4] + "\t" + row[6] + "\t" + row[5] + "\t" + row[7] + "\t" + str(variantratio) + "\t" + info[1] + "\t" + info[2] + "\t" + str(readdepth) + "\t" + row[9] + "\t" + row[10] + "\t" + row[11] + "\t" + row[12] + "\t" + variant_type_exact.strip("SVTYPE=") + "\t")

		# Iterate through transcripts and print them in separate fields
		# for runner in transcripts:
		# 	ofile.write(runner + "\t")
		ofile.write(transcripts + "\t")
		ofile.write("\n")

ifile.close()
ofile.close()
