#!/usr/bin/python2.7

from __future__ import division
import csv
import re
import argparse

'''

Script for converting mutect output to Annovar input

Viktor Ljungstrom 2013-11-19

'''

parser = argparse.ArgumentParser(description='Script for converting mutect output to Annovar input')

parser.add_argument('-i','--inputfile', help='Input txt file (Required)', required=True)
parser.add_argument('-o','--outputfile', help='Output txt file (Required)', required=True)

args = vars(parser.parse_args())

ifile  = open(args['inputfile'], "rb")
reader = csv.reader(ifile, delimiter='\t')
ofile  = open(args['outputfile'], "wb")
#writer = csv.writer(ofile, delimiter='\t')

headerlist = []

#Read every row in inputfile
for row in reader:

	#Skip first comment row
	if not row[0].startswith("#"):

		#Find header row and append each column to headerlist
		if row[0].startswith("contig"):
			for headername in row:
				headerlist.append(headername)

		#For all other rows print annovar info first and then add Mutect info
		else:
			ofile.write(row[0]+"\t"+row[1]+"\t"+row[1]+"\t"+row[3]+"\t"+row[4]+"\t"+"comments: ")

			#Write mutect info
			for x in range(5, 51):
				ofile.write(headerlist[x]+"="+row[x]+" ")

			#Add context column
			ofile.write("context="+row[2])

			ofile.write("\n")

ifile.close()
ofile.close()
