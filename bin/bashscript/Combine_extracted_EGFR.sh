#!/bin/bash
#
# Script to run jSNPmania
#SBATCH -p core -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts combining Annovar output ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/Extracted_EGFR ]]; then
	mkdir $ROOT_PATH/Extracted_EGFR;
fi

if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then

		# EGFR T790M - ampliconmapped
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.ampliconmapped.txt ]]; then
			rm $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.ampliconmapped.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*EGFR_T790M.ampliconmapped.txt > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*EGFR_T790M.ampliconmapped.txt > $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.ampliconmapped.txt;

		else
			WarningLog "$SAMPLEID" "No files with ending EGFR_T790M.ampliconmapped.txt were found!";
		fi

		# EGFR T790M
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.txt ]]; then
			rm $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*EGFR_T790M.txt > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*EGFR_T790M.txt > $ROOT_PATH/Extracted_EGFR/EGFR_T790M.allSamples.txt;

		else
			WarningLog "$SAMPLEID" "No files with ending EGFR_T790M.txt were found!";
		fi

		# EGFR G719 - ampliconmapped
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.ampliconmapped.txt ]]; then
			rm $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.ampliconmapped.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*EGFR_G719.ampliconmapped.txt > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*EGFR_G719.ampliconmapped.txt > $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.ampliconmapped.txt;

		else
			WarningLog "$SAMPLEID" "No files with ending EGFR_G719.ampliconmapped.txt were found!";
		fi

		# EGFR G719 - ampliconmapped
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.txt ]]; then
			rm $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*EGFR_G719.txt > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*EGFR_G719.txt > $ROOT_PATH/Extracted_EGFR/EGFR_G719.allSamples.txt;

		else
			WarningLog "$SAMPLEID" "No files with ending EGFR_G719.txt were found!";
                        exit 0;
		fi
	else
		ErrorLog "$SAMPLEID" "Only supported for h.sapiens so far!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true for the analysis to run!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in combining EGFR extractions";
else
	SuccessLog "$SAMPLEID" "Passed combining EGFR extractions";
fi
