#!/bin/bash
#
# Script creates plot data files seqroi vs seqroi.
#
#SBATCH -p core  -n 1
#SBATCH -t 30:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/plotValues" ]]; then
	mkdir $ROOT_PATH/plotValues;
fi

SuccessLog "${SAMPLEID}" "Creating plot data files (SeqROI vs SeqROI)...";

if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ! -e $ROOT_PATH/plotValues/${SAMPLEID}.seqroi_seqroi.uniq.nusbaum || ! -e $ROOT_PATH/plotValues/${SAMPLEID}.seqroi_seqroi.uniq.stenberg || ! -z $FORCE ]]; then
			if [[ -e $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map ]]; then
				perl $SERA_PATH/bin/perlscript/hitsPerBase2ampregionNormalizedValue.pl -a $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map -i $ROOT_PATH/refFiles/${REFSEQ}.seqroi -c $ROOT_PATH/plotValues/${SAMPLEID}.seqroi_seqroi.uniq.stenberg -n $ROOT_PATH/plotValues/${SAMPLEID}.seqroi_seqroi.uniq.nusbaum;
			else
				ErrorLog "${SAMPLEID}" "$ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map does NOT exist!";
			fi
		else
			ErrorLog "${SAMPLEID}" "Output files already existed and force was not used!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Only supported for call_type h.sapiens so far!";
	fi
else
	ErrorLog "${SAMPLEID}" "READS has to be true to run the analysis!";
fi

# Check if hitsPerBase2ampregionNormalizedValue.pl worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in creating plot values (SeqROI vs SeqROI).";
else
	SuccessLog "${SAMPLEID}" "Passed creating plot values (SeqROI vs SeqROI).";
fi
