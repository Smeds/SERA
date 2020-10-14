#!/usr/bin/python2.7

from __future__ import division
import csv
import re
import argparse
import string

'''

Script for converting Pindel output to Annovar input

'''

parser = argparse.ArgumentParser(description = 'Script for converting Pindel output to Annovar input')

parser.add_argument('-i', '--inputfile', help = 'Input txt file (Required)', required = True)
parser.add_argument('-o', '--outputfile', help = 'Output txt file (Required)', required = True)

args = vars(parser.parse_args())

ifile = open(args['inputfile'], "rb")
reader = csv.reader(ifile, delimiter = '\t')
ofile = open(args['outputfile'], "wb")

for row in reader:




	if row[0].startswith("#"):
		ofile.write('\t'.join(row) + "\n")
	else:

		normal_column = row[10].replace(":", ",")
		normal_refreads = normal_column.split(",")[1]
		normal_varreads = normal_column.split(",")[2]
		normal_depth = int(normal_varreads) + int(normal_refreads)

		if normal_depth != 0:

			normal_refratio = int(normal_refreads) / normal_depth

		tumor_column = row[9]
		tumor_reads = tumor_column.split(",")[1]


		if normal_depth > 30 and normal_refratio > 0.98 and int(tumor_reads) > 5:

			del row[9]
			ofile.write('\t'.join(row) + "\n")
