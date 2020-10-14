#!/bin/bash
#
# Script running bwa
#
#SBATCH -p core  -n 3
#SBATCH -t 02:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Bwa" ]]; then
	mkdir $ROOT_PATH/Bwa;
fi

SuccessLog "${SAMPLEID}" "Aligning paired-end reads with bwa mem";

fastq_files_r1=($(echo "$RAWDATA_PE1" | tr " " "\n"));

# Check for cutadapt sequences
if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]]; then
	PE1=${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz;
else
    if [[ ${#fastq_files_r1[@]} > 1 ]];
    then
        if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R1_001.fastq.gz | head -1)" ];
        then
            PE1=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R1_001.fastq.gz);
        else
            ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read1, please pre-process data!...";
        fi
     else
        PE1=$RAWDATA_PE1;
    fi
fi
if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz ]]; then
	PE2=${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz;
else
    if [[ ${#fastq_files_r1[@]} > 1 ]];
    then
        if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R2_001.fastq.gz | head -1)" ];
        then
            PE2=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R2_001.fastq.gz);
        else
            ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read2, please pre-process data!...";
        fi
     else
        PE2=$RAWDATA_PE2;
    fi
fi

SuccessLog "${SAMPLEID}" "Using $PE1 as read1 and $PE2 as read2 as input to bwa...";

# Run bwa
# If platform is Illumina
if [[ $PLATFORM = "Illumina" ]]; then
	# Check that the output file doesn't exist or if force is given
	if [[ ! -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam || ! -z $FORCE ]]; then
	    if [[ ${MATE_PAIR} == "true" ]]; then
    		# Check that input files exist
    		if [[ -e ${PE1} && ${PE2} ]]; then
    			# Get the date when the analysis is run
    			now=$('date' +"%Y%m%d")
    			bwa mem -M -R "@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina" ${GENOME_REF} ${PE1} ${PE2} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;
    			samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
    			samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt;


    			SuccessLog "${SAMPLEID}" "bwa mem -M -R \"@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina\" ${GENOME_REF} ${PE1} ${PE2} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;"
    			SuccessLog "${SAMPLEID}" samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;

    		else
    			ErrorLog "${SAMPLEID}" "${PE1} and/or ${PE2} do NOT exist!";
    		fi
		else
            if [[ -e ${PE1} ]]; then
                # Get the date when the analysis is run
                now=$('date' +"%Y%m%d")
                bwa mem -M -R "@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina" ${GENOME_REF} ${PE1} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;
                samtools index $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
                samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats.txt;


                SuccessLog "${SAMPLEID}" "bwa mem -M -R \"@RG\tID:"$now"_${SAMPLEID}\tSM:${SAMPLEID}\tPL:illumina\" ${GENOME_REF} ${PE1} -t 3 | samtools view -bS /dev/stdin | samtools sort -@ 3 /dev/stdin $ROOT_PATH/Bwa/${SAMPLEID}.sorted;"
                SuccessLog "${SAMPLEID}" samtools flagstat $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam;
            else
                ErrorLog "${SAMPLEID}" "${PE1} do NOT exist!";
            fi
        fi
	else
		ErrorLog "${SAMPLEID}" "$ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam already exists and force was NOT used!";

	fi

else
	ErrorLog "${SAMPLEID}" "Platform $PLATFORM not recognized (only interpreted for Illumina).";
	exit 1;
fi

# Check if bwa worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in bwa alignment...";
else
	SuccessLog "${SAMPLEID}" "Passed bwa alignment";
fi
