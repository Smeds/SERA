#!/usr/bin/python2.7

"""

Script for creating input file for sera pipeline

Elin Falk Sorqvist 20141114

"""

import argparse
import re
import os
import shutil
import sys
import module_locator


# Parse commandline
parser = argparse.ArgumentParser(description = "This script creates a SERA input file from a file with sample name and barcode")
parser.add_argument('-i', '--infile', help = 'Input file name', type = str, required = True)
parser.add_argument('-p', '--project', help = 'Name of the project', type = str, required = True)
parser.add_argument('-g', '--globals', help = 'Should the HOME or PROJ option be used', choices = ['HOME', 'PROJ', 'MORIARTY'], type = str, required = True)
parser.add_argument('-refDir', '--refDir', help = 'If a directory for reference files is given, they will be copied to the analysis folder', type = str)
parser.add_argument('-n', '--normal', help = 'Name of normal sample to compare with.', type = str)
parser.add_argument('-s', '--software', help = 'Queing system for', type = str, default = "SLURM")
parser.add_argument('-f', '--fileFormat', help = 'File format of reference files. Default: bed', type = str, default = "bed")
parser.add_argument('-platform', '--platform', help = 'Platform for sequencing. Default: Illumina', type = str, default = "Illumina")
parser.add_argument('-readLength', '--readLength', help = 'Read length in sequencing. Default: 150', type = str, default = "150")
parser.add_argument('-matePair', '--matePair', help = 'If sequencing is run with mate pair set true otherwise false. Default: true', type = str, default = "true")
parser.add_argument('-designType', '--designType', help = 'Set the type of design, PCR or MDA. Default: PCR', type = str, default = "PCR")
parser.add_argument('-alignerFlags', '--alignerFlags', help = 'Set the alignment parameters. Default: -hs 15 -act 35 -mhp 100 -ms 7 -p 8 -minp 0.95 -mmp 0.05', type = str, default = "-hs 15 -act 35 -mhp 100 -ms 7 -p 8 -minp 0.95 -mmp 0.05")
parser.add_argument('-pindelFlags', '--pindelFlags', help = 'Set parameters for pindel analysis. Default: -T 5 -x 2 -B 100', type = str, default = "-T 5 -x 2 -B 100")
parser.add_argument('-snpmaniaFlags', '--snpmaniaFlags', help = 'Set the SNPmania parameters. Default: -a 5 -q 20 -e 0 -am', type = str, default = "-a 5 -q 20 -e 0 -am")
parser.add_argument('-tumorNormalFlags', '--tumorNormalFlags', help = 'Set parameters for pre-annovar filtering when both tumor and normal are present. Default: -tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5', type = str, default = "-tminRD 20 -nminRD 20 -tminVarRatio 0.10 -nHomoRefRatio 0.95 -am 5")
parser.add_argument('-annovarFlags', '--annovarFlags', help = 'Set parameters for pre-annovar filtering when only tumor is present. Default: -minRD 20 -minVarRatio 0.01 -am 5', type = str, default = " -minRD 20 -minVarRatio 0.01 -am 5")
parser.add_argument('-annovarPlasmaFlags', '--annovarPlasmaFlags', help = 'Set parameters for pre-annovar filtering when only tumor is present. Default: -minRD 20 -minVarRatio 0.001', type = str, default = " -minRD 20 -minVarRatio 0.001")
parser.add_argument('-pindelAnnovarFlags', '--pindelAnnovarFlags', help = 'Set parameters for pre-annovar filtering of pindel results. Default: -m 20 -v 0.01', type = str, default = "-m 20 -v 0.01")
parser.add_argument('-pindelAnnovarPlasmaFlags', '--pindelAnnovarPlasmaFlags', help = 'Set parameters for pre-annovar filtering of pindel results. Default: -m 20 -v 0.001', type = str, default = "-m 20 -v 0.001")
parser.add_argument('-mutationFlags', '--mutationsFlags', help = 'Set parameters for the filtering of all mutations both hotspots and others. Default: -minRD 30,300 -minVaf 0.01', type = str, default = " -minRD 30,300 -minVaf 0.01")
parser.add_argument('-mutationPlasmaFlags', '--mutationsPlasmaFlags', help = 'Set parameters for the filtering of all mutations both hotspots and others in plasma. Default: -minRD 30,300 -minVaf 0.001', type = str, default = " -minRD 30,300 -minVaf 0.001")
parser.add_argument('-clinicalInfoFile', '--clinicalInfoFile', help = 'File with clinical hotspot and indel filenames per cancer type. If not wanted set to false! Default: clinicalCancerTypeFiles.txt', type = str, default = "clinicalCancerTypeFiles.txt")


args = parser.parse_args()
info = {}

clinicalInfo = {}
# Go through file with info about different clinical files
if not args.clinicalInfoFile.lower() == "false":
    # Using module module_locator to determine the path of the pythonscript directory
    my_path = module_locator.module_path()
    my_pathParts = my_path.split("/")  # Split path on /
    # Get the path of the SERA_script folder
    ciFilePath = ""
    for i in range(1, (len(my_pathParts) - 2)):
        ciFilePath += "/" + my_pathParts[i]
    # Add the path to the file with info about files corresponding to different cancer types
    ciFilePath += "/res/" + args.clinicalInfoFile

    # Go through the cancer type file and add info about the corresponding files to a directory
    with open(ciFilePath, 'r') as cifile:
        for line in cifile:
            if not re.match('$', line):  # If line is not empty start working
                line = line.strip()
                infoParts = line.split("=")
                cancer = infoParts[1].lower()
                clinicalInfo[cancer] = {}

                # Retrieve file info for each cancer type and add to dictionary
                for i in range (0, 4):
                    line = next(cifile)
                    line = line.strip()
                    infoParts = line.split("=")
                    # Check if the file given i - or false, then add false as filename
                    if re.match('-', infoParts[1]) or re.match('false', infoParts[1].lower()):
                        clinicalInfo[cancer][infoParts[0].lower()] = "false"
                    else:  # Add the filename
                        clinicalInfo[cancer][infoParts[0].lower()] = infoParts[1]

        if not cifile.closed:
            cifile.close()

with open(args.infile, 'r') as infile:
    count = 0
    # Go through the file line by line
    for line in infile:
        # Check so the line isn't empty or starts with ##

        if not re.match('^##', line) and not re.match('$', line) and not re.match('^Experiment', line):
            count += 1
            line = line.rstrip('\r\n')  # Remove new line character
            splitPattern = '\t'
            if re.search(";", line):
                splitPattern = ';'
            elif re.search(",", line):
                splitPattern = ','
            lineSplit = line.split(splitPattern)  # Split line by tab
            sample = lineSplit[1]

            info[sample] = {}
            info[sample]['exp'] = lineSplit[0]
            info[sample]['sNummer'] = "S" + str(count)
            info[sample]['barcodeI7'] = lineSplit[2]
            info[sample]['barcodeI5'] = lineSplit[3]
            info[sample]['design'] = lineSplit[4]
            info[sample]['refseq'] = lineSplit[5]
            info[sample]['cutadapt'] = lineSplit[6]
            info[sample]['method'] = lineSplit[7].lower().strip()
            info[sample]['type'] = lineSplit[8].lower().strip()
            info[sample]['tissue'] = lineSplit[9].lower().strip()


            if info[sample]['cutadapt'] == "":
                info[sample]['cutadapt'] = "false"
            if info[sample]['barcodeI7'] == "":
                info[sample]['barcodeI7'] = "false"
            if info[sample]['barcodeI5'] == "":
                info[sample]['barcodeI5'] = "false"

            if len(lineSplit) > 10 and args.normal:
                parser.print_usage()
                print("\nERROR: The normal has to be given either in the file or the commandline - NOT both!\n\n")
                sys.exit()
            else:
                if len(lineSplit) > 10:
                    info[sample]['normal'] = lineSplit[10]
                else:
                    if args.normal:
                        info[sample]['normal'] = args.normal
                    else:
                        parser.print_usage()
                        sys.exit("\nERROR: Either inputfile needs to have 7 columns with the normal given in the last column or the flag -normal has to be used!\n\n")

            # Set the method type used
            if info[sample]['method'] == "halo" or info[sample]['method'] == "haloplex":
                info[sample]['method'] = "haloplex"
            if info[sample]['method'] == "swift" or info[sample]['method'] == "accamp" or info[sample]['method'] == "AccelAmplicon":
                info[sample]['method'] = "swift"

            # If clinical analysis is wanted add file names otherwise add false
            if not args.clinicalInfoFile.lower() == "false":
                # ## COLON
                if info[sample]['tissue'] == "colon" or info[sample]['tissue'] == "kolon" :
                    info[sample]['tissue'] = "colon"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo['colon']['hotspot']):
                        info[sample]['hotspot'] = "false"
                    else:
                        info[sample]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['hotspot']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['amplification']):
                        info[sample]['amplification'] = "false"
                    else:
                        info[sample]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["colon"]['background']):
                        info[sample]['background'] = "false"
                    else:
                        info[sample]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["colon"]['background']

                # ## LUNG
                elif info[sample]['tissue'] == "lung" or info[sample]['tissue'] == "lunga":
                    info[sample]['tissue'] = "lung"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['hotspot']):
                        info[sample]['hotspot'] = "false"
                    else:
                        info[sample]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['hotspot']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['amplification']):
                        info[sample]['amplification'] = "false"
                    else:
                        info[sample]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["lung"]['background']):
                        info[sample]['background'] = "false"
                    else:
                        info[sample]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["lung"]['background']

                # ## GIST
                elif info[sample]['tissue'] == "gist":
                    info[sample]['tissue'] = "gist"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['hotspot']):
                        info[sample]['hotspot'] = "false"
                    else:
                        info[sample]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['hotspot']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['amplification']):
                        info[sample]['amplification'] = "false"
                    else:
                        info[sample]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["gist"]['background']):
                        info[sample]['background'] = "false"
                    else:
                        info[sample]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["gist"]['background']

                # ## OVARIAL
                elif info[sample]['tissue'] == "ovarial" or info[sample]['tissue'] == "ovarian" or info[sample]['tissue'] == "ovary" :
                    info[sample]['tissue'] = "ovarial"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['hotspot']):
                        info[sample]['hotspot'] = "false"
                    else:
                        info[sample]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['hotspot']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['amplification']):
                        info[sample]['amplification'] = "false"
                    else:
                        info[sample]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["ovarial"]['background']):
                        info[sample]['background'] = "false"
                    else:
                        info[sample]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["ovarial"]['background']

                # ## MELANOM
                elif info[sample]['tissue'] == "melanom" or info[sample]['tissue'] == "melanoma":
                    info[sample]['tissue'] = "melanom"
                    # If the given hotspot file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['hotspot']):
                        info[sample]['hotspot'] = "false"
                    else:
                        info[sample]['hotspot'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['hotspot']

                    # If the given amplification file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['amplification']):
                        info[sample]['amplification'] = "false"
                    else:
                        info[sample]['amplification'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['amplification']

                    # If the given background file name is false keep it otherwise add file path
                    if re.match("false", clinicalInfo["melanom"]['background']):
                        info[sample]['background'] = "false"
                    else:
                        info[sample]['background'] = "$FILE_PATH/refFiles/" + clinicalInfo["melanom"]['background']

                else:
                    print ("\nERROR: Unknown cancer type given in input file " + info[sample]['tissue'] + ", so far only lung, colon, GIST, melanoma and ovarial are supported!\n\n")
                    sys.exit()
            else:
                info[sample]['hotspot'] = "false"
                info[sample]['amplification'] = "false"
                info[sample]['background'] = "false"

    if not infile.closed:
        infile.close()

expFolder = ""
count = 1
infoSort = sorted(info)
for sample in infoSort:
    if count == 1:
        expFolder = info[sample]['exp'] + "_" + sample
        count += 1
    else:
        expFolder += "-" + sample

rawPath = "/proj/" + args.project + "/private/" + info[sample]['exp'] + "_rawdata"
filePath = "/proj/" + args.project + "/nobackup/private/" + info[sample]['exp']

# Create folder if it doesn't exist and set permissions
if not os.path.exists(rawPath):
    os.makedirs(rawPath, 0774)

# Create folder if it doesn't exist and set permissions
if not os.path.exists(filePath):
    os.makedirs(filePath, 0774)

# Check if a reference file dir is given, if so copy the files to refFiles in the analysis folder
if args.refDir:
    # List files in the reference file directory
    refDir_files = os.listdir(args.refDir)
    refFilePath = filePath + "/refFiles"  # Setting output folder path
    # Check if the output folder exists, if not create and set permissions
    if not os.path.exists(refFilePath):
        os.makedirs(refFilePath, 0774)
    # Go through all files in the reference file dir
    for file_name in refDir_files:
        full_file_name = os.path.join(args.refDir, file_name)
        # Check that it is a file
        if (os.path.isfile(full_file_name)):
             destFile = refFilePath + "/" + file_name  # Set the destination file name
             os.system("rsync -rlptD " + full_file_name + " " + destFile)  # Copy
             if not re.match((oct(os.stat(destFile).st_mode & 0777)), "0664"):  # Check if the file has permission to read and write for all in the group
                 os.chmod(destFile, 0664)  # If not change the permissions

new_format = re.compile(".+/([0-9]{8})_([A-Za-z0-9-]+)$")
new_format_rerun = re.compile(".+/([0-9]{8})_([A-Za-z0-9-]+)_([0-9]+)")

def extract_date_user_rerun(line):
    result = new_format_rerun.match(line)
    if result:
        return result.group(1),result.group(2),True
    result = new_format.match(line)
    if result:
        return result.group(1),result.group(2),False
    raise Exception("Unable to parse input file for experiment name, user and rerun")

(date, user, rerun) = extract_date_user_rerun(filePath)


output = filePath + "/inputFile"
with (open(output, mode = 'w'))as outfile:
    outfile.write('################################\n### SEQUENCING-SPECIFIC SETTINGS\n################################\n\n# These settings are shared for all samples included in this analysis\n\n')
    outfile.write("SOFTWARE=\"" + args.software + "\";\n")
    outfile.write("UPPNEX_PROJECT_ID=\"" + args.project + "\";\n")
    outfile.write("GLOBALS=\"" + args.globals + "\";\n")
    outfile.write("RAW_PATH=\"" + rawPath + "\";\n")
    outfile.write("FILE_PATH=\"" + filePath + "\";\n")
    outfile.write("\nexport FILE_FORMAT=\"" + args.fileFormat + "\";\n")
    outfile.write("export READS=\"true\";\n")
    outfile.write("export MOLECULES=\"false\";\n")
    outfile.write("\n## Alignment\n")
    outfile.write("export PLATFORM=\"" + args.platform + "\" # SOLiD, Illumina or IonTorrent\n")
    outfile.write("export READ_LENGTH=\"" + args.readLength + "\";\n")
    outfile.write("export MATE_PAIR=\"" + str(args.matePair) + "\"; # true or false\n")
    outfile.write("export DESIGN_TYPE=\"" + args.designType + "\"; # PCR or MDA\n")
    outfile.write("export ALIGNERFLAGS_GENOME=\"" + args.alignerFlags + "\";\n")
    outfile.write("export PINDEL_FLAGS=\"" + args.pindelFlags + "\";\n")
    outfile.write("\n## Filtering\n")
    outfile.write("export CALL_TYPE=\"h.sapiens\";\n")
    outfile.write("export SNPMANIAFLAGS=\"" + args.snpmaniaFlags + "\";\n")
    outfile.write("export TUMOR_NORMAL_FLAGS=\"" + args.tumorNormalFlags + "\";\n")
    outfile.write("export ANNOVAR_FLAGS=\"" + args.annovarFlags + "\";\n")
    outfile.write("export ANNOVAR_PLASMA_FLAGS=\"" + args.annovarPlasmaFlags + "\";\n")
    outfile.write("export PINDEL_ANNOVAR_FLAGS=\"" + args.pindelAnnovarFlags + "\";\n")
    outfile.write("export PINDEL_ANNOVAR_PLASMA_FLAGS=\"" + args.pindelAnnovarPlasmaFlags + "\";\n")
    outfile.write("export MUTATION_FLAGS=\"" + args.mutationsFlags + "\";\n")
    outfile.write("export MUTATION_PLASMA_FLAGS=\"" + args.mutationsPlasmaFlags + "\";\n")




    outfile.write('\n################################\n### SAMPLE SETTINGS\n################################\n\n')

    infoSort = sorted(info)
    for sample in infoSort:
        outfile.write("let COUNT=COUNT+1;\n")
        outfile.write ("SAMPLEID_ARR_[${COUNT}]=\"" + sample + "\";\n")
        outfile.write ("BARCODE_I7_ARR_[${COUNT}]=\"" + info[sample]['barcodeI7'] + "\";\n")
        outfile.write ("BARCODE_I5_ARR_[${COUNT}]=\"" + info[sample]['barcodeI5'] + "\";\n")
        outfile.write ("CUTADAPT_PREFIX_ARR_[${COUNT}]=\"" + info[sample]['cutadapt'] + "\";\n")
        outfile.write ("RAWDATA_PE1_ARR_[${COUNT}]=\"$RAW_PATH/" + sample + "_" + info[sample]['sNummer'] + "_L001_R1_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_PE2_ARR_[${COUNT}]=\"$RAW_PATH/" + sample + "_" + info[sample]['sNummer'] + "_L001_R2_001.fastq.gz" + "\";\n")
        outfile.write ("RAWDATA_INDEX_ARR_[${COUNT}]=\"false\";\n")
        outfile.write ("REFSEQ_ARR_[${COUNT}]=\"" + info[sample]['refseq'] + "\";\n")
        if re.match(info[sample]['method'], "haloplex"):
            outfile.write ("ROIFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Regions.bed" + "\";\n")
        else:
            outfile.write ("ROIFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + ".bed" + "\";\n")
        if re.match(info[sample]['method'], "haloplex"):
            outfile.write ("SELECTIONFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + "_Amplicons.bed" + "\";\n")
        else:
            outfile.write ("SELECTIONFILE_ARR_[${COUNT}]=\"$FILE_PATH/refFiles/" + info[sample]['design'] + ".bed" + "\";\n")
        outfile.write ("NORMAL_SAMPLEID_ARR_[${COUNT}]=\"" + info[sample]['normal'] + "\";\n")
        outfile.write ("HOTSPOTFILE_ARR_[${COUNT}]=\"" + info[sample]['hotspot'] + "\";\n")
        outfile.write ("AMPLIFICATIONFILE_ARR_[${COUNT}]=\"" + info[sample]['amplification'] + "\";\n")
        outfile.write ("BACKGROUNDFILE_ARR_[${COUNT}]=\"" + info[sample]['background'] + "\";\n")
        outfile.write ("METHOD_ARR_[${COUNT}]=\"" + info[sample]['method'] + "\";\n")
        outfile.write ("TYPE_ARR_[${COUNT}]=\"" + info[sample]['type'] + "\";\n")
        outfile.write ("TISSUE_ARR_[${COUNT}]=\"" + info[sample]['tissue'] + "\";\n")
        outfile.write ("\n")
    if not outfile.closed:
        outfile.close()

with open(jsonPath + experiment +".json", mode="w") as json_output:
    for sample in sorted(info):
        json_output.write(str(
            {'experiment.wp': "WP1",
             'experiment.prep': info[sample]['type'].upper(),
             'experiment.method': info[sample]['method'],
             'experiment.user': user,
             'experiment.rerun': rerun,
             'experiment.tissue': info[sample]['tissue'],
             'experiment.sample': sample}) + "\n")
