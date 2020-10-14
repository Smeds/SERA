#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts extracting info about EGFR positions ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Extracted_sampleInfo" ]]; then
	mkdir $ROOT_PATH/Extracted_sampleInfo;
fi

# Check that READS are true
if [[ ${READS} == "true" ]]; then
	# Check that the call type are set to h.sapiens
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ${TISSUE} == "lung" ]]; then
			# Check if the output contain amplicon information
			if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then

				# Remove output files if they already exists
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.ampliconmapped.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.ampliconmapped.txt;
				fi
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt
				fi

				#EGFR T790M
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.ampliconmapped.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.ampliconmapped.txt
					grep -P "\s55249071\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v cds="c.2369C>T" -v aa="p.T790M" -v r="C" -v v="T" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.ampliconmapped.txt;
				fi


				#EGFR G719
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt;
				fi

				grep -P "\s55241707\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v cds="c.2155G>A" -v aa="p.G719" -v r="G" -v v="A" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt;
				grep -P "\s55241707\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v cds="c.2155G>T" -v aa="p.G719" -v r="G" -v v="T" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt;
				grep -P "\s55241708\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v cds="c.2156G>C" -v aa="p.G719" -v r="G" -v v="C" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt;
				grep -P "\s55241708\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v cds="c.2156G>A" -v aa="p.G719" -v r="G" -v v="A" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.ampliconmapped.txt;


			# If there is not a file with amplicon information, check if there is one without
			elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then

				# Remove output files if they already exists
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.txt;
				fi
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt
				fi

				#EFGR T790M
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.txt;
				fi
				grep -P "\s55249071\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v cds="c.2369C>T" -v aa="p.T790M" -v r="C" -v v="T" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_T790M.txt;

				#EGFR G719
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt;
				fi

				grep -P "\s55241707\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v cds="c.2155G>A" -v aa="p.G719" -v r="G" -v v="A" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt;
				grep -P "\s55241707\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v cds="c.2155G>T" -v aa="p.G719" -v r="G" -v v="T" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt;
				grep -P "\s55241708\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v cds="c.2156G>C" -v aa="p.G719" -v r="G" -v v="C" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt;
				grep -P "\s55241708\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v cds="c.2156G>A" -v aa="p.G719" -v r="G" -v v="A" -f $SERA_PATH/bin/awkscript/extract_T790M.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_EGFR_G719.txt;

			else
				ErrorLog "$SAMPLEID" "Input file $ROOT_PATH/SNPmania/${SAMPLEID}.variations doesn't exist!";
			fi
		else
			WarningLog "$SAMPLEID" "EGFR extraction is only run for lung!";
		fi
	else
		ErrorLog "$SAMPLEID" "The analysis is only supported for CALL_TYPE h.sapiens!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
fi


if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extractin info about EGFR positions";
else
	SuccessLog "$SAMPLEID" "Passed extracting info about EGFR positions";
fi
