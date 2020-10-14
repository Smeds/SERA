#!/usr/bin/python2.7

from variant_functions import *
import re

# Print variants, both bwa and pindel on hotspot positions. If now variant found print base and it's total read depth
def printHotspots(hotspots, minRDs, sample, transcripts, OUTPUT, ampliconMapped, multipleBp):
    for chrom in hotspots['hotspot']:
            for start in hotspots['hotspot'][chrom]:
                for end in hotspots['hotspot'][chrom][start]:

                    # If there are mutations from bwa OUTPUT.write them
                    if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['bwa'])):
                        pin = "no"
                        # check if there are also pindel variants present at this position
                        if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):
                            pin = "yes"
                        # Go through all bwa variants
                        for bwaVar in hotspots['hotspot'][chrom][start][end]['bwa']:

                            if re.match("yes", pin):  # If there are pindel variants available extract them
                                pinExist = False;
                                for pindelVar in hotspots['hotspot'][chrom][start][end]['pindel']:  # Go through all pindel variants
                                    # Check if there are any pindel variant which is identical to the bwa variant
                                    if (int(bwaVar[2]) == int(pindelVar[2]) and int(bwaVar[3]) == int(pindelVar[3]) and re.match(bwaVar[4], pindelVar[4]) and re.match(bwaVar[5], pindelVar[5])):
                                        pinExist = True;  # If so set pinExist to True
                                #if not pinExist:  # If no identical pindel variant exists, print the bwa variant
                                if True:
                                    aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                    if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a
                                        exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                        if exist:  # If the variant exists in hash add info
                                            aa = aaTmp
                                            cds = cdsTmp
                                            if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                                accNum = accNumTmp
                                                comm = "altTranscript"

                                    found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                    if re.match("-", found):
                                        found = "yes"
                                        if not ampliconCheckHotspot(bwaVar, ampliconMapped):  # if the amplicon mapping criteria is not ful filled found = no
                                            found = "no"

                                    if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                        exon = hotspots['hotspot'][chrom][start][end]['exon']

                                    refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                                    mutType = checkMutationType(bwaVar, "1-hotspot")  # Get the mutation type

                                    comment = setComment(comm, hotspots['hotspot'][chrom][start][end]['comment'])  # Set the correct comment
                                    OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + " \t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")

                            else:  # If no pindel variants are available write the variant
                                aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                                if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a
                                    exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                    if exist:  # If the variant exists in hash add info
                                        aa = aaTmp
                                        cds = cdsTmp
                                        if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                            accNum = accNumTmp
                                            comm = "altTranscript"

                                found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                    exon = hotspots['hotspot'][chrom][start][end]['exon']
                                if re.match("-", found):
                                    found = "yes"
                                    if not ampliconCheckHotspot(bwaVar, ampliconMapped):  # if the amplicon mapping criteria is not ful filled found = no
                                        found = "no"

                                refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                                mutType = checkMutationType(bwaVar, "1-hotspot")  # Get the mutation type

                                comment = setComment(comm, hotspots['hotspot'][chrom][start][end]['comment'])  # Set the correct comment
                                OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['hotspot'][chrom][start][end]['comment'] + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")


                    # Check if there are pindel variants available
                    if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):

                        for pindelVar in hotspots['hotspot'][chrom][start][end]['pindel']:
                             #only write if not called by normal caller
                             bwaExist = False;
                             if not re.match("no", str(hotspots['hotspot'][chrom][start][end]['bwa'])):
                                for bwaVar in hotspots['hotspot'][chrom][start][end]['bwa']:  # Go through all bwa variants
                                    # Check if there are any pindel variant which is identical to the bwa variant
                                    if (int(bwaVar[2]) == int(pindelVar[2]) and int(bwaVar[3]) == int(pindelVar[3]) and re.match(bwaVar[4], pindelVar[4]) and re.match(bwaVar[5], pindelVar[5])):
                                        bwaExist = True;  # If so set bwaExist to True
                             if not bwaExist:
                                aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(pindelVar, transcripts)

                                if int(pindelVar[3]) - int(pindelVar[2]) > 0 or re.match("-", pindelVar[4]) or (int(pindelVar[3]) - int(pindelVar[2]) > 0 and re.match("-", pindelVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                    exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(pindelVar, multipleBp)
                                    if exist:  # If the variant exists in hash add info
                                        aa = aaTmp
                                        cds = cdsTmp
                                        if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                            accNum = accNumTmp
                                            comm = "altTranscript"

                                # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                                found, readLevel = getReadLevel(minRDs, hotspots['hotspot'][chrom][start][end]['rd'][0])
                                if re.match("-", found):
                                    found = "yes"
                                    if not ampliconCheckHotspot(pindelVar, ampliconMapped):  # if the amplicon mapping criteria is not ful filled found = no
                                       found = "no"
                                if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                    exon = hotspots['hotspot'][chrom][start][end]['exon']

                                refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(pindelVar, ampliconMapped)  # add amplicon information

                                comment = setComment(comm, hotspots['hotspot'][chrom][start][end]['comment'])  # Set the correct comment
                                OUTPUT.write (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + hotspots['hotspot'][chrom][start][end]['comment'] + "\t2-indel\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]) + "\n")


                    if re.match("no", str(hotspots['hotspot'][chrom][start][end]['bwa'])) and re.match("no", str(hotspots['hotspot'][chrom][start][end]['pindel'])):  # If no variant is added OUTPUT.write position
                        for c in range(0, int(end) - int(start) + 1):
                            found = readLevel = rd = "NA"
                            if isinstance(hotspots['hotspot'][chrom][start][end]['rd'], list):  # check that the read depths are within an array, otherwise it's not in design

                                if re.match("-", str(hotspots['hotspot'][chrom][start][end]['rd'][c])):  # if read depth is - the position is not within design
                                    found = "not in design"
                                    readLevel = "-"
                                    rd = "-"
                                else:  # read depth exists
                                    found, readLevel = getReadLevel (minRDs, int(hotspots['hotspot'][chrom][start][end]['rd'][c]))
                                    if re.match("-", found):
                                        found = "no"
                                    rd = getReadLevel (minRDs, int(hotspots['hotspot'][chrom][start][end]['rd'][c]))
                            else:  # position is not within design
                                found = "not in design"
                                readLevel = "-"
                                rd = "-"

                            OUTPUT.write (sample + "\t" + str(hotspots['hotspot'][chrom][start][end]['gene']) + "\t-\t" + str(hotspots['hotspot'][chrom][start][end]['exon']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['aa']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['cds']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['accNum']) + "\t" + str(hotspots['hotspot'][chrom][start][end]['comment']) + "\t1-hotspot\t" + found + "\t" + readLevel + "\t" + str(hotspots['hotspot'][chrom][start][end]['rd'][c]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t-\t-\t-\n")

# OUTPUT.write variants within region, both bwa pindel and the total read depth for positions without variant.
def printRegionAll(hotspots, minRDs, sample, transcripts, OUTPUT, ampliconMapped, multipleBp):
    for chrom in hotspots['region_all']:
        for start in hotspots['region_all'][chrom]:
            for end in hotspots['region_all'][chrom][start]:
                # If there are mutations from bwa OUTPUT.write them
                if not re.match("no", str(hotspots['region_all'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['region_all'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            pinExist = False;
                            for pindelVar in hotspots['region_all'][chrom][start][end]['pindel']:  # Go through all pindel variants
                                    # Check if there are any pindel variant which is identical to the bwa variant
                                if (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                    pinExist = True;  # If so set pinExist to True
                            #if not pinExist:  # If no identical pindel variant exists, print the bwa variant
                            if True:
                                aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                                if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a
                                    exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                    if exist:  # If the variant exists in hash add info
                                        aa = aaTmp
                                        cds = cdsTmp
                                        if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                            accNum = accNumTmp
                                            comm = "altTranscript"

                                found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                if re.match("-", found):
                                    found = "yes"
                                if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                    exon = hotspots['region_all'][chrom][start][end]['exon']

                                refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                                mutType = "3-check"  # Set the mutation type

                                comment = setComment(comm, hotspots['region_all'][chrom][start][end]['comment'])  # Set the correct comment
                                OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")
                        else:  # If no pindel variants are available OUTPUT.write the variant
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                            if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['region_all'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['region_all'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['region_all'][chrom][start][end]['pindel']:
                        #only write if not called by normal caller
                        bwaExist = False;
                        if not re.match("no", str(hotspots['region_all'][chrom][start][end]['bwa'])):
                            for bwaVar in hotspots['region_all'][chrom][start][end]['bwa']:  # Go through all bwa variants
                                # Check if there are any pindel variant which is identical to the bwa variant
                                if (int(bwaVar[2]) == int(pindelVar[2]) and int(bwaVar[3]) == int(pindelVar[3]) and re.match(bwaVar[4], pindelVar[4]) and re.match(bwaVar[5], pindelVar[5])):
                                    bwaExist = True;  # If so set bwaExist to True
                        if not bwaExist:
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(pindelVar, transcripts)

                            if int(pindelVar[3]) - int(pindelVar[2]) > 0 or re.match("-", pindelVar[4]) or (int(pindelVar[3]) - int(pindelVar[2]) > 0 and re.match("-", pindelVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(pindelVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                            found, readLevel = getReadLevel(minRDs, hotspots['region_all'][chrom][start][end]['rd'][0])
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['region_all'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(pindelVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['region_all'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]) + "\n")

                # Check if neither bwa nor pindel variant exists
                if re.match("no", str(hotspots['region_all'][chrom][start][end]['bwa'])) and re.match("no", str(hotspots['region_all'][chrom][start][end]['pindel'])):  # If no variant is added OUTPUT.write position
                    for c in range(0, int(end) - int(start) + 1):
                        found = readLevel = rd = "NA"
                        if isinstance(hotspots['region_all'][chrom][start][end]['rd'], list):  # check that the read depths are within an array, otherwise it's not in design

                            if re.match("-", str(hotspots['region_all'][chrom][start][end]['rd'][c])):  # if read depth is - the position is not within design
                                found = "not in design"
                                readLevel = "-"
                                rd = "-"
                            else:  # read depth exists
                                found, readLevel = getReadLevel (minRDs, int(hotspots['region_all'][chrom][start][end]['rd'][c]))
                                if re.match("-", found):
                                    found = "no"

                        else:  # position is not within design
                            found = "not in design"
                            readLevel = "-"
                            rd = "-"

                        if re.match("low", readLevel) or re.match("-", readLevel) or re.match("ok", readLevel):  # Print all positions, regardless of readLevel
                            mutType = "3-check"  # Set the mutation type
                            OUTPUT.write (sample + "\t" + str(hotspots['region_all'][chrom][start][end]['gene']) + "\t-\t" + str(hotspots['region_all'][chrom][start][end]['exon']) + "\t" + str(hotspots['region_all'][chrom][start][end]['aa']) + "\t" + str(hotspots['region_all'][chrom][start][end]['cds']) + "\t" + str(hotspots['region_all'][chrom][start][end]['accNum']) + "\t" + str(hotspots['region_all'][chrom][start][end]['comment']) + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(hotspots['region_all'][chrom][start][end]['rd'][c]) + "\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t" + str(chrom) + "\t" + str(start + c) + "\t" + str(start + c) + "\t-\t-\t-\n")

# OUTPUT.write variants within region, both bwa and pindel. However it does not OUTPUT.write total read depth for all positions without variant
def printRegion(hotspots, minRDs, sample, transcripts, OUTPUT, ampliconMapped, multipleBp):
    for chrom in hotspots['region']:
        for start in hotspots['region'][chrom]:
            for end in hotspots['region'][chrom][start]:
                # If there are mutations from bwa OUTPUT.write them
                if not re.match("no", str(hotspots['region'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['region'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['region'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            pinExist = False;
                            for pindelVar in hotspots['region'][chrom][start][end]['pindel']:  # Go through all pindel variants
                                    # Check if there are any pindel variant which is identical to the bwa variant
                                if (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                    pinExist = True;  # If so set pinExist to True
                            #if not pinExist:  # If no identical pindel variant exists, print the bwa variant
                            if True:
                                aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information
                                if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                    exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                    if exist:  # If the variant exists in hash add info
                                        aa = aaTmp
                                        cds = cdsTmp
                                        if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                            accNum = accNumTmp
                                            comm = "altTranscript"

                                found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                if re.match("-", found):
                                    found = "yes"
                                if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                    exon = hotspots['region'][chrom][start][end]['exon']

                                refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                                mutType = "3-check"  # Set the mutation type

                                comment = setComment(comm, hotspots['region'][chrom][start][end]['comment'])  # Set the correct comment
                                OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")
                        else:  # If no pindel variants are available OUTPUT.write the variant
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                            if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['region'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['region'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['region'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['region'][chrom][start][end]['pindel']:
                        bwaExist = False;
                        #only write if not called by normal caller
                        if not re.match("no", str(hotspots['region'][chrom][start][end]['bwa'])):
                            for bwaVar in hotspots['region'][chrom][start][end]['bwa']:  # Go through all bwa variants
                                # Check if there are any pindel variant which is identical to the bwa variant
                                if (int(bwaVar[2]) == int(pindelVar[2]) and int(bwaVar[3]) == int(pindelVar[3]) and re.match(bwaVar[4], pindelVar[4]) and re.match(bwaVar[5], pindelVar[5])):
                                    bwaExist = True;  # If so set bwaExist to True
                        if not bwaExist:
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(pindelVar, transcripts)

                            if int(pindelVar[3]) - int(pindelVar[2]) > 0 or re.match("-", pindelVar[4]) or (int(pindelVar[3]) - int(pindelVar[2]) > 0 and re.match("-", pindelVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(pindelVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                            found, readLevel = getReadLevel(minRDs, hotspots['region'][chrom][start][end]['rd'][0])
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['region'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(pindelVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['region'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]) + "\n")


# OUTPUT.write indels within indel region, both bwa and pindel. However it does not OUTPUT.write total read depth for all positions without variant
def printIndel(hotspots, minRDs, sample, transcripts, OUTPUT, ampliconMapped, multipleBp):
    for chrom in hotspots['indel']:
        for start in hotspots['indel'][chrom]:
            for end in hotspots['indel'][chrom][start]:
                # If there are mutations from bwa OUTPUT.write them
                if not re.match("no", str(hotspots['indel'][chrom][start][end]['bwa'])):
                    pin = "no"
                    # check if there are also pindel variants present at this position
                    if not re.match("no", str(hotspots['indel'][chrom][start][end]['pindel'])):
                        pin = "yes"
                    # Go through all bwa variants
                    for bwaVar in hotspots['indel'][chrom][start][end]['bwa']:

                        if re.match("yes", pin):  # If there are pindel variants available extract them
                            for pindelVar in hotspots['indel'][chrom][start][end]['pindel']:
                                # Check that there is no pindel variant which is identical to the bwa variant
                                #if not (bwaVar[2] == pindelVar[2] and bwaVar[3] == pindelVar[3] and bwaVar[4] == pindelVar[4] and bwaVar[5] == pindelVar[5]):
                                if True:
                                    aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                                    if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                        exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                        if exist:  # If the variant exists in hash add info
                                            aa = aaTmp
                                            cds = cdsTmp
                                            if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                                accNum = accNumTmp
                                                comm = "altTranscript"

                                    found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                                    if re.match("-", found):
                                        found = "yes"
                                    if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                        exon = hotspots['indel'][chrom][start][end]['exon']

                                    refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                                    mutType = "3-check"  # Set the mutation type

                                    comment = setComment(comm, hotspots['indel'][chrom][start][end]['comment'])  # Set the correct comment
                                    OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")
                        else:  # If no pindel variants are available OUTPUT.write the variant
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(bwaVar, transcripts)  # Get transcript information

                            if int(bwaVar[3]) - int(bwaVar[2]) > 0 or re.match("-", bwaVar[4]) or (int(bwaVar[3]) - int(bwaVar[2]) > 0 and re.match("-", bwaVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(bwaVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            found, readLevel = getReadLevel(minRDs, bwaVar[12])  # Get the level of read depth and if it's not analyzable
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['indel'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(bwaVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['indel'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(bwaVar[0]) + "\t" + str(bwaVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(bwaVar[12]) + "\t" + str(bwaVar[10]) + "\t" + str(bwaVar[11]) + "\t" + str(bwaVar[9]) + "\t" + str(bwaVar[14]) + "\t" + str(bwaVar[13]) + "\t" + str(bwaVar[16]) + "\t" + str(bwaVar[15]) + "\t" + str(bwaVar[17]) + "\t" + str(bwaVar[18]) + "\t" + str(bwaVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(bwaVar[20]) + "\t" + str(bwaVar[21]) + "\t" + str(bwaVar[22]) + "\t" + str(bwaVar[23]) + "\t" + str(bwaVar[24]) + "\t" + str(bwaVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(bwaVar[1]) + "\t" + str(bwaVar[2]) + "\t" + str(bwaVar[3]) + "\t" + str(bwaVar[4]) + "\t" + str(bwaVar[5]) + "\t" + str(bwaVar[32]) + "\n")


                # Check if there are pindel variants available
                if not re.match("no", str(hotspots['indel'][chrom][start][end]['pindel'])):

                    for pindelVar in hotspots['indel'][chrom][start][end]['pindel']:
                        bwaExist = False;
                        #only write if not called by normal caller
                        if not re.match("no", str(hotspots['indel'][chrom][start][end]['bwa'])):
                            for bwaVar in hotspots['indel'][chrom][start][end]['bwa']:  # Go through all bwa variants
                                # Check if there are any pindel variant which is identical to the bwa variant
                                if (int(bwaVar[2]) == int(pindelVar[2]) and int(bwaVar[3]) == int(pindelVar[3]) and re.match(bwaVar[4], pindelVar[4]) and re.match(bwaVar[5], pindelVar[5])):
                                    bwaExist = True;  # If so set bwaExist to True
                        if not bwaExist:
                            aa, cds, accNum, exon, exonicType, comm = getTranscriptInfo(pindelVar, transcripts)

                            if int(pindelVar[3]) - int(pindelVar[2]) > 0 or re.match("-", pindelVar[4]) or (int(pindelVar[3]) - int(pindelVar[2]) > 0 and re.match("-", pindelVar[4])):  # If it is a multiple bp variant check if there are any match in our multipleBp hash
                                exist, aaTmp, cdsTmp, accNumTmp = getMultipleBpInfo(pindelVar, multipleBp)
                                if exist:  # If the variant exists in hash add info
                                    aa = aaTmp
                                    cds = cdsTmp
                                    if not accNum == accNumTmp:  # if there is a different accession number set comment to altTranscript
                                        accNum = accNumTmp
                                        comm = "altTranscript"

                            # Check the level of reads. Below the lowest => -, between first and second => low, above highest => ok
                            found, readLevel = getReadLevel(minRDs, hotspots['indel'][chrom][start][end]['rd'][0])
                            if re.match("-", found):
                                found = "yes"
                            if re.match("-", exon):  # If no exon info was given for the mutation take it from hotspot input file
                                exon = hotspots['indel'][chrom][start][end]['exon']

                            refPlus, refMinus, varPlus, varMinus, refAll, varAll = getAmpliconInfo(pindelVar, ampliconMapped)  # add amplicon information
                            mutType = "3-check"  # Set the mutation type

                            comment = setComment(comm, hotspots['indel'][chrom][start][end]['comment'])  # Set the correct comment
                            OUTPUT.write (str(pindelVar[0]) + "\t" + str(pindelVar[6]) + "\t" + exonicType + "\t" + exon + "\t" + aa + "\t" + cds + "\t" + accNum + "\t" + comment + "\t" + mutType + "\t" + found + "\t" + readLevel + "\t" + str(pindelVar[12]) + "\t" + str(pindelVar[10]) + "\t" + str(pindelVar[11]) + "\t" + str(pindelVar[9]) + "\t" + str(pindelVar[14]) + "\t" + str(pindelVar[13]) + "\t" + str(pindelVar[16]) + "\t" + str(pindelVar[15]) + "\t" + str(pindelVar[17]) + "\t" + str(pindelVar[18]) + "\t" + str(pindelVar[19]) + "\t" + str(refPlus) + "\t" + str(refMinus) + "\t" + str(varPlus) + "\t" + str(varMinus) + "\t" + str(pindelVar[20]) + "\t" + str(pindelVar[21]) + "\t" + str(pindelVar[22]) + "\t" + str(pindelVar[23]) + "\t" + str(pindelVar[24]) + "\t" + str(pindelVar[25]) + "\t" + str(refAll) + "\t" + str(varAll) + "\t" + str(pindelVar[1]) + "\t" + str(pindelVar[2]) + "\t" + str(pindelVar[3]) + "\t" + str(pindelVar[4]) + "\t" + str(pindelVar[5]) + "\t" + str(pindelVar[32]) + "\n")
