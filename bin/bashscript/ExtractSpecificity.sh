#!/bin/bash
#
# Script to extract specificities from the on- and offtarget files
#
#SBATCH -p core  -n 1
#SBATCH -t 01:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Specificity" ]]; then
	mkdir $ROOT_PATH/Specificity;
fi

SuccessLog "${SAMPLEID}" "Starting counting on- and offtarget bases...";
if [[ ${READS} == "true" ]]; then
	# Check if the ampregion input files exists
	if [[ -e $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget.gz && $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget.gz ]]; then
		gunzip $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget.gz;
		gunzip $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget.gz;

		# Check if the gunzip worked, if so extract the specificities
		if [[ -e $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget && -e $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget ]]; then
			awk -v name=${SAMPLEID} '{if(FILENAME~/offtarget/){offT+=$3} else if(FILENAME~/ontarget/){onT+=$3}} END{print name"\t"onT/(onT+offT)}' $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.o* >> $ROOT_PATH/Specificity/allSamples.ampregion.specificity;
			# gzip the files again to save space
			gzip $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget;
			gzip $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget;
		else
			ErrorLog "${SAMPLEID}" "The gunzip of the ampregion inputfiles failed!";
		fi
	else
		ErrorLog "${SAMPLEID}" "The input files $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget.gz & $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget.gz don't exist!";
	fi

	# Check if the seqregion input files exists
	if [[ -e $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget.gz && $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget.gz ]]; then
		gunzip $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget.gz;
		gunzip $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget.gz;

		# Check if the gunzip worked, if so extract the specificities
		if [[ -e $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget && -e $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget ]]; then
			awk -v name=${SAMPLEID} '{if(FILENAME~/offtarget/){offT+=$3} else if(FILENAME~/ontarget/){onT+=$3}} END{print name"\t"onT/(onT+offT)}' $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.o* >> $ROOT_PATH/Specificity/allSamples.seqregion.specificity;
	        # gzip the files to save space
	        gzip $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget;
	        gzip $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget;
	    else
	    	ErrorLog "${SAMPLEID}" "The gunzip of the seqregion inputfiles failed!";
	    fi
	else
		ErrorLog "${SAMPLEID}" "The input files $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget.gz & $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget.gz don't exist!";
	fi
fi

# Check if specificity calculations worked
if [[ "$?" != "0" ]]; then
        ErrorLog "${SAMPLEID}" "failed in calculating specificity...";
else
        SuccessLog "${SAMPLEID}" "passed calculating specificity...";
fi
