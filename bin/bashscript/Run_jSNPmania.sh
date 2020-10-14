#!/bin/bash
#
# Script to run jSNPmania
#SBATCH -p core  -n 4
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts jSNPmania ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/SNPmania" ]]; then
	mkdir $ROOT_PATH/SNPmania;
fi


# Start with checking that the reference file exists!

if [[ $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq ]]; then
	if [[ ${READS} == "true" ]]; then
		# Check if the file is ampliconMapped
		if [[ -e $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam ]]; then
			if [[ ! -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations || ! -z $FORCE ]]; then
				# Reheader the ampliconmapped file
				if [[ -e $SERA_PATH/res/Convert2SNPmaniaHeader.sam ]]; then
					samtools reheader $SERA_PATH/res/Convert2SNPmaniaHeader.sam $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam | samtools view /dev/stdin | $ROOT_PATH_JSNPMANIA/JSNPmania.sh -i /dev/stdin -o $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -oi $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -od $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -am ${SNPMANIAFLAGS};
				else
					ErrorLog "$SAMPLEID" "The SNPmania header file doesn't exists!";
				fi
			else
				ErrorLog "$SAMPLEID" "$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations already exists and force was not used!";
			fi
		# If not use the "normal" bam-file
		elif [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam ]]; then
			if [[ ! -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations || ! -z $FORCE ]]; then
				if [[ -e $SERA_PATH/res/Convert2SNPmaniaHeader.sam ]]; then
					samtools reheader $SERA_PATH/res/Convert2SNPmaniaHeader.sam $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | samtools view /dev/stdin | $ROOT_PATH_JSNPMANIA/JSNPmania.sh -i /dev/stdin -o $ROOT_PATH/SNPmania/${SAMPLEID}.variations -oi $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -od $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq ${SNPMANIAFLAGS};
				else
					ErrorLog "$SAMPLEID" "The SNPmania header file doesn't exists!";
				fi
			else
				ErrorLog "$SAMPLEID" "$ROOT_PATH/SNPmania/${SAMPLEID}.variations already exists and force was not used!";
			fi

		else
			ErrorLog "$SAMPLEID" "None of the possible input files could be found ($ROOT_PATH/ampliconMapped/${SAMPLEID}.withAmplicon.bam or $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam)!";
		fi
	else
		ErrorLog "$SAMPLEID" "READS has to be set true to run the analysis!";
	fi
else
	ErrorLog "${SAMPLEID}" "SNPmania reference file $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq does NOT exist!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in JSNPmania";
else
	SuccessLog "$SAMPLEID" "Passed JSNPmania";
fi
