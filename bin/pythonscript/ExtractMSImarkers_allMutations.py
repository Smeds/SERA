#!/usr/bin/python2.7

import re
import argparse

'''

Script for converting Pindel output to Annovar input

'''

parser = argparse.ArgumentParser(description = "This script takes an annovar output file or a filtered annovar output file and outputs MSI markers.")
parser.add_argument('-i', '--inputfile', help = 'Input annovar ouput file (Required)', required = True)
parser.add_argument('-t1', '--tgfbr2Ratio_1bp', help = 'Minimum variant allele frequency for a 1bp TGFBR2 variant to be extracted', required = False, default = 0.01, type = float)
parser.add_argument('-t2', '--tgfbr2Ratio_2bp', help = 'Minimum variant allele frequency for a two or more bp TGFBR2 variant to be extracted', required = False, default = 0.01, type = float)
parser.add_argument('-a1', '--acvr2aRatio_1bp', help = 'Minimum variant allele frequency for a 1 bp ACVR2A variant to be extracted', required = False, default = 0.01, type = float)
parser.add_argument('-a2', '--acvr2aRatio_2bp', help = 'Minimum variant allele frequency for a two or more bp ACVR2A variant to be extracted', required = False, default = 0.01, type = float)
parser.add_argument('-o', '--outputannovar', help = 'Output in annovar outputformat file (Required)', required = True)

args = parser.parse_args()

with open(args.inputfile, 'r') as infile:
    with open(args.outputannovar, 'w') as outannovar:
        done = 0
        # Go through the file line by line
        for line in infile:
            # If the line is a comment line print it to output file
            if line.startswith("Sample") and done == 0:
                # outannovar.write(line)
                done = 1

            else:
                line = line.rstrip('\r\n').split("\t")  # Remove new line character and split on tab
                if re.match("yes", line[9]):
                    start = int(line[35])
                    end = int(line[36])
                    ref = line[37]
                    var = line[38]
                    gene = line[1]
                    vaf = float(line[14])

                    if ("A" in ref and "A" in var) or "A" in ref or "A" in var:
                        if re.match('^TGFBR2$', gene):
                            if end - start == 0 and len(var) == 1:
                                if vaf >= args.tgfbr2Ratio_1bp:
                                    outannovar.write("\t".join(line) + "\n")
                            else:
                                if vaf >= args.tgfbr2Ratio_2bp:
                                    outannovar.write("\t".join(line) + "\n")
                        elif re.match('^ACVR2A$' , gene):
                            if end - start == 0 and len(var) == 1:
                                if vaf >= args.acvr2aRatio_1bp:
                                    outannovar.write("\t".join(line) + "\n")
                            else:
                                if vaf >= args.acvr2aRatio_2bp:
                                    outannovar.write("\t".join(line) + "\n")
