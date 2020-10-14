#!/bin/bash
#
# Script to run Pindel
#SBATCH -p core  -n 5
#SBATCH -t 08:00:00
##SBATCH --qos=short -t 00:15:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/PindelOut" ]]; then
	mkdir $ROOT_PATH/PindelOut;
fi

if [[ ! -e $ROOT_PATH/PindelOut/${SAMPLEID}.indels_TD || ! -z $FORCE ]]; then
	if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
		SuccessLog "${SAMPLEID}" "will be analysed in Pindel single sample mode";
		# Check that input file exists
		if [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam ]]; then
			# Check if the pindel config file exists, if so remove
			if [[ -e $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt ]]; then
				rm $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt
			fi
			# Cretae config file
			echo -e "$ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam\t300\t${SAMPLEID}" >> $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt
			echo "$ROOT_PATH_PINDEL/pindel -f ${GENOME_FASTA_REF} -i $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt -J $SERA_PATH/res/excludeChrM.txt -o $ROOT_PATH/PindelOut/${SAMPLEID}.indels ${PINDEL_FLAGS}"
			$ROOT_PATH_PINDEL/pindel -f ${GENOME_FASTA_REF} -i $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt -J $SERA_PATH/res/excludeChrM.txt -o $ROOT_PATH/PindelOut/${SAMPLEID}.indels ${PINDEL_FLAGS}
		else
			ErrorLog "${SAMPLEID}" "Inputfile $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam doesn't exist!";
		fi

	elif [[ ${NORMAL_SAMPLEID} != "false" ]]; then

		SuccessLog "${SAMPLEID}" "will be analysed in Pindel with normal sample:" "${NORMAL_SAMPLEID}";
		# Check that input file exists
		if [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam && $ROOT_PATH/Bwa/${NORMAL_SAMPLEID}.sorted.bam ]]; then
			# Check if the pindel config file exists, if so remove
			if [[ -e $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt ]]; then
				rm $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt
			fi

			echo -e "$ROOT_PATH/Bwa/${NORMAL_SAMPLEID}.sorted.bam\t300\t${NORMAL_SAMPLEID}" >> $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt
			echo -e "$ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam\t300\t${SAMPLEID}" >> $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt

			$ROOT_PATH_PINDEL/pindel -f ${GENOME_FASTA_REF} -i $ROOT_PATH/PindelOut/${SAMPLEID}.config.txt -c ALL  -o $ROOT_PATH/PindelOut/${SAMPLEID}.indels ${PINDEL_FLAGS}
		else
			ErrorLog "${SAMPLEID}" "Inputfile $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam and/or $ROOT_PATH/Bwa/${NORMAL_SAMPLEID}.sorted.bam doesn't exist!";
		fi

	else
			SuccessLog "${SAMPLEID}" "Normal set to false, Pindel will not be started.";
	fi

else

	ErrorLog "${SAMPLEID}" "Pindel output already exists and Force was not used.";

fi
