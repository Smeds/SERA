#!/usr/bin/python2.7

import re
import argparse

'''

Script for converting Pindel output to Annovar input

'''

parser = argparse.ArgumentParser(description = "This script takes an annovar output file filter it on different parameters.")
parser.add_argument('-i', '--inputfile', help = 'Input txt file (Required)', required = True)
parser.add_argument('-b', '--blacklistfile', help = 'File with blacklisted variants (Required)', required = True)
parser.add_argument('-k', '--intronickeepfile', help = 'File with blacklisted variants')
parser.add_argument('-g', '--genome1000', help = 'Set the maximum allowed frequency for a variant to be kept, default 0', default = 0, type = float, required = False)
parser.add_argument('-o', '--outputfile', help = 'Output txt file (Required)', required = True)
parser.add_argument('-a', '--ampliconmapped', help = 'Set if data is ampliconmapped. Default: False', action = "store_true")

args = parser.parse_args()

blacklist = {}
with open(args.blacklistfile, 'r') as blackfile:
	for line in blackfile:
		if not re.match('^Chr', line):
			line = line.rstrip('\r\n').split("\t")
			chr = line[0]
			start = line[1]
			end = line[2]
			ref = line[3]
			var = line[4]
			gene = line[5]


			if not line[0] in blacklist:
				blacklist[chr] = {}
				blacklist[chr][start] = {}
				blacklist[chr][start][end] = {}
				blacklist[chr][start][end][ref] = {}
				blacklist[chr][start][end][ref][var] = gene
			else:
				if not start in blacklist[chr]:
					blacklist[chr][start] = {}
					blacklist[chr][start][end] = {}
					blacklist[chr][start][end][ref] = {}
					blacklist[chr][start][end][ref][var] = gene
				else:
					if not end in blacklist[chr][start]:
						blacklist[chr][start][end] = {}
						blacklist[chr][start][end][ref] = {}
						blacklist[chr][start][end][ref][var] = gene
					else:
						if not ref in blacklist[chr][start][end]:
							blacklist[chr][start][end][ref] = {}
							blacklist[chr][start][end][ref][var] = gene
						else:
							if not var in blacklist[chr][start][end][ref]:
								blacklist[chr][start][end][ref][var] = gene


# 			print (chr + "\t" + start + "\t" + end + "\t" + ref + "\t" + var + "\t" + blacklist[chr][start][end][ref][var])
# Even if intronic, keep variants in this interval
keeplist = {}
if args.intronickeepfile:
	with open(args.intronickeepfile, 'r') as keepfile:
		for line in keepfile:
			line = line.rstrip('\r\n').split("\t")  # Remove new line character
			regions = line[1].split(";")  # Split regions on semi-colon
			# Go through all regions
			for reg in regions:
				chrInfo = reg.split(":")  # split chr and position info
				positions = chrInfo[1].split("-")  # split position info into start and end
				if not line[0] in keeplist:
					keeplist[line[0]] = {}
					keeplist[line[0]][chrInfo[0]] = {}
					keeplist[line[0]][chrInfo[0]][positions[0]] = positions[1]
				else:
					if not chrInfo[0] in keeplist[line[0]]:
						keeplist[line[0]][chrInfo[0]] = {}
						keeplist[line[0]][chrInfo[0]][positions[0]] = positions[1]
					else:
						if not positions[0] in keeplist[line[0]][chrInfo[0]]:
							keeplist[line[0]][chrInfo[0]][positions[0]] = positions[1]


			keeplist[line[0].lower()] = line[1]
			# If the gene file is not closed, close
		if not keepfile.closed:
			keepfile.close()


with open(args.inputfile, 'r') as infile:
	with open(args.outputfile, 'w') as outfile:
		first = "true"
		# Go through the file line by line
		for line in infile:

			if line.startswith("Sample") and re.match("true", first):
				outfile.write(line)
				first = "false"

			else:

				line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
				chr = line[1]
				start = line[2]
				end = line[3]
				ref = line[4]
				var = line[5]
				gene = line[6]

				done = 0
				# Keep if it's exonic or splicing
				if re.match('^exonic', line[7]) or re.match('splicing', line[7]):
					# Remove all synonymous variants
					if not re.match('^synonymous', line[8]):
						# Keep those without a value in 1000G or a value that is lower than the given input
						if re.match('-', str(line[13])) or float(line[13]) <= args.genome1000:

							varOK = "true"
							# Check if ampliconmapped is set and if so check amplicon status
							if args.ampliconmapped:
								varOK = "false"
								# Set amplicon info
								ref_plus = ref_minus = var_plus = var_minus = 0
								if not re.match('-', line[26]):
									var_plus = int(line[26])
								if not re.match('-', line[27]):
									var_minus = int(line[27])
								if not re.match('-', line[28]):
									ref_plus = int(line[28])
								if not re.match('-', line[29]):
									ref_minus = int(line[29])

								if (ref_plus + ref_minus) >= 3:
									if (var_plus + var_minus) >= 2:
										varOK = "true"
								elif (ref_plus + ref_minus) >= 1:
									if (var_plus + var_minus) >= 1:
										varOK = "true"
								elif (ref_plus + ref_minus) >= 0:
									if (var_plus + var_minus) >= 0:
										varOK = "true"

							# Check if the amplicon info is okey
							if re.match('true', varOK):
								# Remove blacklist variants
								if not chr in blacklist:
									outfile.write("\t".join(line) + "\n")
									done = 1
								else:
									if not start in blacklist[chr]:
										outfile.write ("\t".join(line) + "\n")
										done = 1
									else:
										if not end in blacklist[chr][start]:
											outfile.write ("\t".join(line) + "\n")
											done = 1
										else:
											if not ref in blacklist[chr][start][end]:
												outfile.write ("\t".join(line) + "\n")
												done = 1
											else:
												if not var in blacklist[chr][start][end][ref]:
													outfile.write ("\t".join(line) + "\n")
													done = 1
												else:
													print("Exist in blacklist!\n")
													print ("\t" + "\t".join(line) + "\n")
													done = 1


				# If the variant isn't printed or in blacklist
				if done == 0:
					# Check if the gene is in the
					if gene in keeplist:
 						for regionChr in keeplist[gene]:
							for regionStart in keeplist[gene][regionChr]:
								if regionStart <= start and keeplist[gene][regionChr][regionStart] >= end:
									if done == 0:
										if not re.match('^synonymous', line[8]):
											# Keep those without a value in 1000G or a value that is lower than the given input
											if re.match('-', str(line[13])) or float(line[13]) <= args.genome1000:
												varOK = "true"
												# Check if ampliconmapped is set and if so check amplicon status
												if args.ampliconmapped:
													varOK = "false"
													# Set amplicon info
													ref_plus = ref_minus = var_plus = var_minus = 0
													if not re.match('-', line[26]):
														var_plus = int(line[26])
													if not re.match('-', line[27]):
														var_minus = int(line[27])
													if not re.match('-', line[28]):
														ref_plus = int(line[28])
													if not re.match('-', line[29]):
														ref_minus = int(line[29])

													if (ref_plus + ref_minus) >= 3:
														if (var_plus + var_minus) >= 2:
															varOK = "true"
													elif (ref_plus + ref_minus) >= 1:
														if (var_plus + var_minus) >= 1:
															varOK = "true"
													elif (ref_plus + ref_minus) >= 0:
														if (var_plus + var_minus) >= 0:
															varOK = "true"

												# Check if the amplicon info is okey
												if re.match('true', varOK):
													# Remove blacklist variants
													if not regionChr in blacklist:
														outfile.write("\t".join(line) + "\n")
														done = 1
													else:
														if not start in blacklist[regionChr]:
															outfile.write ("\t".join(line) + "\n")
															done = 1
														else:
															if not end in blacklist[regionChr][start]:
																outfile.write ("\t".join(line) + "\n")
																done = 1
															else:
																if not ref in blacklist[regionChr][start][end]:
																	outfile.write ("\t".join(line) + "\n")
																	done = 1
																else:
																	if not var in blacklist[regionChr][start][end][ref]:
																		outfile.write ("\t".join(line) + "\n")
																		done = 1
																	else:
																		print("Exist in blacklist!\n")
																		print ("\t" + "\t".join(line) + "\n")
																		done = 1
