#!/usr/bin/python2.7

from __future__ import division
import csv
import re
import argparse

'''

Script for converting Annovar output to readable output

Viktor Ljungstrom 2013-11-19

'''

parser = argparse.ArgumentParser(description='Script for converting Annovar output to readable output')

parser.add_argument('-i','--inputfile', help='Input txt file (Required)', required=True)
parser.add_argument('-o','--outputfile', help='Output txt file (Required)', required=True)

args = vars(parser.parse_args())

ifile  = open(args['inputfile'], "rb")
reader = csv.reader(ifile, delimiter='\t')
ofile  = open(args['outputfile'], "wb")

#Print header
ofile.write("#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tCosmic\tContext\tStrandbias\tCovered\tPower\tTumor_power\tNormal_power\tTranscripts\n")

#Read every row in inputfile
for row in reader:

	if not row[0].startswith("Chr"):

		info = row[13].split(' ')
		#print  info[1].split('=')[1]
		sample = info[1].split('=')[1]
		variantratio = info[17].split('=')[1]
		altreads = info[21].split('=')[1]
		refreads = info[22].split('=')[1]
		readdepth = int(altreads) + int(refreads)
		context = info[47].split('=')[1]
		strandbias = info[39].split('=')[1]
		covered = info[5].split('=')[1]
		power = info[6].split('=')[1]
		tumor_power = info[7].split('=')[1]
		normal_power = info[8].split('=')[1]
		transcripts = row[8].split(',')

		ofile.write(sample+"\t"+row[0]+"\t"+row[1]+"\t"+row[2]+"\t"+row[3]+"\t"+row[4]+"\t"+row[6]+"\t"+row[5]+"\t"+row[7]+"\t"+variantratio+"\t"+altreads+"\t"+refreads+"\t"+str(readdepth)+"\t"+row[9]+"\t"+row[10]+"\t"+row[11]+"\t"+row[12]+"\t"+context+"\t"+strandbias+"\t"+covered+"\t"+power+"\t"+tumor_power+"\t"+normal_power+"\t")

		for runner in transcripts:
			ofile.write(runner+"\t")

		ofile.write("\n")

ifile.close()
ofile.close()
