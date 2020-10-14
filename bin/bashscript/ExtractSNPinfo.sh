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

SuccessLog $SAMPLEID "Starts extracting info about SNP positions ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Extracted_sampleInfo" ]]; then
	mkdir $ROOT_PATH/Extracted_sampleInfo;
fi

# Check that READS are true
if [[ ${READS} == "true" ]]; then
	# Check that the call type are set to h.sapiens
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ${TISSUE} == "colon" ]]; then
			# Check if the output contain amplicon information
			if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then

				# Remove output files if they already exists
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
				fi

				#EGFR T790M
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tID\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt
					grep -P "\s33025979\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="C" -v v="T" -v name="rs10318" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s111189158\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="A" -v v="G" -v name="rs10749971" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s6404281\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="C" -v v="A" -v name="rs961253" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s128413305\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="G" -v v="T" -v name="rs6983267" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s6388068\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="T" -v v="C" -v name="rs355527" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s128407443\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="A" -v v="G" -v name="rs10505477" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;
					grep -P "\s46459032\s" $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -v tumor=${TISSUE} -v r="T" -v v="C" -v name="rs4464148" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.ampliconmapped.txt;

				fi

			#EFGR T790M
			# Check if the output without amplicon information exists
			elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then

				# Remove output files if they already exists
				if [[ -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt ]]; then
					rm $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
				fi

				#EGFR T790M
				if [[ ! -e $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt ]]; then
					echo -e '#Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tID\tRef_amp\tVar_amp' > $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt
					grep -P "\s33025979\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="C" -v v="T" -v name="rs10318" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s111189158\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="A" -v v="G" -v name="rs10749971" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s6404281\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="C" -v v="A" -v name="rs961253" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s128413305\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="G" -v v="T" -v name="rs6983267" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s6388068\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="T" -v v="C" -v name="rs355527" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s128407443\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="A" -v v="G" -v name="rs10505477" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;
					grep -P "\s46459032\s" $ROOT_PATH/SNPmania/${SAMPLEID}.variations | awk -v fpath=$ROOT_PATH/SNPmania/${SAMPLEID}.variations -v tumor=${TISSUE} -v r="T" -v v="C" -v name="rs4464148" -f $SERA_PATH/bin/awkscript/extract_SNPs.awk /dev/stdin >> $ROOT_PATH/Extracted_sampleInfo/${SAMPLEID}_SNPs.txt;

				fi
			else
				ErrorLog "$SAMPLEID" "Neiher input file $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor $ROOT_PATH/SNPmania/${SAMPLEID}.variations do exist!";
			fi

		else
			WarningLog "$SAMPLEID" "SNP extraction is only run for colon!";
		fi
	else
		ErrorLog "$SAMPLEID" "The analysis is only supported for CALL_TYPE h.sapiens!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
fi


if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extracting info about SNPs";
else
	SuccessLog "$SAMPLEID" "Passed extracting info about SNPs";
fi
