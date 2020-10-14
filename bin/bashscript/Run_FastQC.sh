#!/bin/bash
#
# Script for running MosaikBuild on either SOLiD or Illumina samples.
#
#SBATCH -p core  -n 1
#SBATCH -t 1:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/FastQC" ]]; then
    mkdir $ROOT_PATH/FastQC;
fi

if [[ "$MATE_PAIR" == "true" ]]; then
    extraLoginfo="mate-pair/pair-ends";
else
    extraLoginfo="single";
fi

fastq_files_r1=($(echo "$RAWDATA_PE1" | tr " " "\n"));

if [[ ${#fastq_files_r1[@]} > 1 ]];
then
    if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R1_001.fastq.gz | head -1)" ];
    then
        READ1=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R1_001.fastq.gz);
    else
        ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read1, please pre-process data!...";
    fi
    if [[ "$MATE_PAIR" == "true" ]]; then
        if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R2_001.fastq.gz | head -1)" ];
        then
            READ2=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R2_001.fastq.gz);
        else
            ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read2, please pre-process data!...";
        fi
    fi
else
    READ1=$RAWDATA_PE1;
    if [[ "$MATE_PAIR" == "true" ]]; then
        READ2=$RAWDATA_PE2;
    fi
fi

SuccessLog "${SAMPLEID}" "Building $extraLoginfo reads Mosaik input file...";

# Run FastQC
# If platform is Illumina
if [[ $PLATFORM = "Illumina" ]]; then

    # Check that the output file doesn't exist or if force is given
    # Do paired end library or not
    if [[ "$MATE_PAIR" == "true" ]]; then
        # Check that input files exist, if not print error message
        if [[ -e ${READ1} && ${READ1} ]]; then
            fastqc -o FastQC ${READ1} ${READ1} > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
        else
            ErrorLog "${SAMPLEID}" "${READ1} and/or ${READ1} do NOT exist!";
        fi

        if [[ ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read1_fastqc.html || ! -e $ROOT_PATH/FastQC/${SAMPLEID}.read2_fastqc.html || ! -z $FORCE ]]; then
            # Check if cutadapted files exist - if so run cutadapt on them
            if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz ]]; then
                fastqc -o FastQC ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_outputCutadapt.txt;
            else
                ErrorLog "${SAMPLEID}" "No cutadapted files exists!"
            fi
        else
            ErrorLog "${SAMPLEID}" "FastQC files already existed for cutadapted output and force was NOT used!";
        fi

    # If single reads are used only care about the first read
    else
        # Check that input file exists, if not print error message
        if [[ -e ${READ2} ]]; then
            fastqc -o FastQC ${READ2} > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_output.txt;
        else
            ErrorLog "${SAMPLEID}" "${READ2} does NOT exist!";
        fi

        # Check if cutadapted files exist - if so run cutadapt on them
        if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]]; then
            fastqc -o FastQC ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz > $ROOT_PATH/FastQC/fastqc_${SAMPLEID}_outputCutadapt.txt;
        else
            ErrorLog "${SAMPLEID}" "No cutadapted file exists!"
        fi
    fi
else
    ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (should be Illumina).";
    exit 1;
fi

# Check if MosaikBuild worked
if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in FastQC...";
else
    SuccessLog "${SAMPLEID}" "Passed FastQC ($extraLoginfo reads)";
fi
