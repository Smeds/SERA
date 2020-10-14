#!/bin/bash
#
# Script to run pileup
#
#SBATCH -p core  -n 1
#SBATCH -t 01:00:00
##SBATCH -p core -t 00:15:00 --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Pileup" ]]; then
	mkdir $ROOT_PATH/Pileup;
fi

SuccessLog "${SAMPLEID}" "Starting pileup BAM files...";

# Check that it's set to h.sapiens and that reads are true, otherwise print error message
if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
	if [[ ${READS} == "true" ]]; then
		# Check if the output file already exists or if force is used, otherwise print error message
		if [[ ! -e "$ROOT_PATH/Pileup/${SAMPLEID}.aligned.pileup.gz" || ! -z $FORCE ]]; then
			# Check that the input file exists and run pileup
			if [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam ]]; then
				samtools mpileup -d 10000000 -f $GENOME_FASTA_REF $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | gzip > $ROOT_PATH/Pileup/${SAMPLEID}.pileup.gz;
			else
				ErrorLog "${SAMPLEID}" "$ROOT_PATH/SamBamFiles/${SAMPLEID}.aligned.sorted.bam does not exist!";
			fi
		else
			ErrorLog "${SAMPLEID}" "Output file $ROOT_PATH/Pileup/${SAMPLEID}.aligned.pileup.gz already exists and force was NOT used!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Reads has to be set to true to run the script!";
	fi
else
	ErrorLog "${SAMPLEID}" "Only implemented for CALL_TYPE h.sapiens so far!";
fi
# Check if it worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in reads pileup run.";
else
	SuccessLog "${SAMPLEID}" "Passed reads pileup.";
fi
