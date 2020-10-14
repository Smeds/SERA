#!/bin/bash
#
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts extracting MSI markers ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/Extracted_sampleInfo ]]; then
	mkdir $ROOT_PATH/Extracted_sampleInfo;
fi



if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
			if [[ ${TISSUE} == "colon" ]]; then
				if [[ -e $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv ]]; then
					if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt  || ! -z $FORCE ]]; then
						python $SERA_PATH/bin/pythonscript/ExtractMSImarkers_allMutations.py -i $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05;
						echo "python $SERA_PATH/bin/pythonscript/ExtractMSImarkers_allMutations.py -i $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -o $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt --tgfbr2Ratio_1bp 0.1 --tgfbr2Ratio_2bp 0.05 --acvr2aRatio_1bp 0.1 --acvr2aRatio_2bp 0.05";
					else
						ErrorLog "$SAMPLEID" "$ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}.msiMarkers.txt already exists and force was not used";
					fi
				else
					ErrorLog "$SAMPLEID" "Input file $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv does not exist!";
				fi
			else
				WarningLog "$SAMPLEID" "MSI markers is only extracted for colon cancer!";
                                exit 0;
			fi
		else
			ErrorLog "$SAMPLEID" "Only run for NORMAL_SAMPLEID annovar";
		fi
	else
		ErrorLog "$SAMPLEID" "Only supported for call_type h.sapiens so far!";
	fi
else
	ErrorLog "$SAMPLEID" "Normal_sampleid is false -> annovar is not run!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extracting MSI markers";
else
	SuccessLog "$SAMPLEID" "Passed extracting MSI markers";
fi
