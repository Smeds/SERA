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
if [[ ! -d $ROOT_PATH/MSIanalysis ]]; then
	mkdir $ROOT_PATH/MSIanalysis;
fi

if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then

		# MSI markers - ampliconmapped
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt ]]; then
			rm $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.txt > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Sample\tGene\tVariant_type\tExon\tAA_change\tCDS_change\tAccession_number\tComment\tReport\tFound\tMin_read_depth300\tTotal_read_depth\tReference_read_depth\tVariant_read_depth\tVariant_allele_ratio\tdbSNP_id\tRatio_1000G\tRatio_ESP6500\tClinically_flagged_dbSNP\tCosmic\tClinVar_CLNDB\tClinval_CLINSIG\tReference_plus_amplicons\tReference_minus_amplicons\tVariant_plus_amplicons\tVariant_minus_amplicons\tStrands_A_F+F-S+S-\tStrands_G_F+F-S+S-\tStrands_C_F+F-S+S-\tStrands_T_F+F-S+S-\tStrands_Ins\tStrands_Del\tRef_aligned_amplicons\tVar_aligned_amplicons\tChr\tStart\tEnd\tReference_base\tVariant_base\tAll_transcripts_annotation"} {if($1!~/#Run/){print $0}}' $ROOT_PATH/Extracted_sampleInfo/*.msiMarkers.txt > $ROOT_PATH/MSIanalysis/MSImarkers.allSamples.txt;

		else
			WarningLog "$SAMPLEID" "No files with ending msiMarkers.txt in Extracted_sampleInfo were found!";
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
