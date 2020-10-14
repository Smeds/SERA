#!/bin/bash
#
# Script calculating mean depth on ampregion.
#
#SBATCH -p core  -n 1
#SBATCH -t 15:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Running mean depth calculation...";

# Calculating mean depth amregion
if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ -e $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqregion.uniq.map ]]; then
			depthSeqRegionUniq=`awk 'BEGIN{s=0;}{s+=\$3}END{print FILENAME"\tDEPTH_SEQREGION_UNIQ="(s/NR);}' $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqregion.uniq.map`;
			echo -e "${depthSeqRegionUniq}" >> $ROOT_PATH/TextResults_meanReadDepth.txt;
		else
			ErrorLog "${SAMPLEID}" "Input file $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqregion.uniq.map does NOT exist!";
		fi
		if [[ -e $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map ]]; then
			depthSeqRoiUniq=`awk 'BEGIN{s=0;}{s+=\$3}END{print FILENAME"\tDEPTH_SEQROI_UNIQ="(s/NR);}' $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map`;
			echo -e "${depthSeqRoiUniq}" >> $ROOT_PATH/TextResults_meanReadDepth.txt;
		else
			ErrorLog "${SAMPLEID}" "Input file $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map does NOT exist!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Only CALL_TYPE h.sapiens is supported so far!";
	fi
else
	ErrorLog "${SAMPLEID}" "READS has to be true to run the analysis!";
fi

if
	SuccessLog "${SAMPLEID}" "Passed calculate mean depth in sequenced region and sequenced roi.";
else
	ErrorLog "${SAMPLEID}" "Failed in calculating mean depth.";
fi
