#!/bin/bash
#
# Script to calculate bases on and off target
#
#SBATCH -p core  -n 1
#SBATCH -t 30:00
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

# Check if the call type are set to h.sapiens and that reads true, otherwise print error messages
if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
	if [[ ${READS} == "true" ]]; then
		# Check that the input file exists
		if [[ -e $ROOT_PATH/Pileup/${SAMPLEID}.pileup.gz ]]; then
			# Check if the reference file with ampregion exists
			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.ampregion ]]; then
				zcat $ROOT_PATH/Pileup/${SAMPLEID}.pileup.gz | perl $SERA_PATH/bin/perlscript/pileup2hitsPerBase.pl -i /dev/stdin -o $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget -off $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget -chr2nc $NC2chr -r $ROOT_PATH/refFiles/${REFSEQ}.ampregion;

				# gzip the output files
				gzip -f $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.ontarget;
				gzip -f $ROOT_PATH/Specificity/${SAMPLEID}.ampregion.offtarget;
			else
				ErrorLog "${SAMPLEID}" "The reference file $ROOT_PATH/refFiles/${REFSEQ}.ampregion doesn't exists!";
			fi

			# Check if the reference file seqregion exists
			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqregion ]]; then
				zcat $ROOT_PATH/Pileup/${SAMPLEID}.pileup.gz | perl $SERA_PATH/bin/perlscript/pileup2hitsPerBase.pl -i /dev/stdin -o $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget -off $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget -chr2nc $NC2chr -r $ROOT_PATH/refFiles/${REFSEQ}.seqregion;
				# gzip output files
				gzip -f $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.ontarget;
				gzip -f $ROOT_PATH/Specificity/${SAMPLEID}.seqregion.offtarget;
			else
				ErrorLog "${SAMPLEID}" "The reference file $ROOT_PATH/refFiles/${REFSEQ}.seqregion doesn't exists!";
			fi
		else
			ErrorLog "${SAMPLEID}" "The input file $ROOT_PATH/Pileup/${SAMPLEID}.pileup.gz doesn't exists!";
		fi
	else
		ErrorLog "${SAMPLEID}" "Reads has to be set to true to run the analysis!";
	fi
else
	ErrorLog "${SAMPLEID}" "The analysis is only implemented for h.sapiens so far!";
fi

# Check if the extraction of on- and offtarget bases worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "failed in extracting on- and offtarget bases...";
else
	SuccessLog "${SAMPLEID}" "on- and offtarget bases extracted...";
fi
