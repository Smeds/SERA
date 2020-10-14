#!/bin/bash
#
# Script runs base mapping on ampregion.
#
#SBATCH -p core  -n 1
#SBATCH -t 30:00
##SBATCH -p core -t 00:15:00 --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/hitsPerBaseFiles" ]]; then
	mkdir $ROOT_PATH/hitsPerBaseFiles;
fi

SuccessLog "${SAMPLEID}" "Running sequenced region base mapping...";

# BaseMapping against ampregion all

# BaseMapping against ampregion unique
if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		if [[ ! -e $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map || ! -z $FORCE ]]; then
			if [[ -e $ROOT_PATH/refFiles/${REFSEQ}.seqroi ]]; then
				if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
					perl $SERA_PATH/bin/perlscript/SNPmania2hitsPerBase.pl -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -i $ROOT_PATH/SNPmania/${SAMPLEID}.variations  -o $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map;
				elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
					perl $SERA_PATH/bin/perlscript/SNPmania2hitsPerBase.pl -r $ROOT_PATH/refFiles/${REFSEQ}.seqroi -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations  -o $ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map;
				else
					ErrorLog "${SAMPLEID}" "None of the possible input files existed ($ROOT_PATH/SNPmania/${SAMPLEID}.variations or $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations)!";
				fi
			else
				ErrorLog "${SAMPLEID}" "ROOT_PATH/refFiles/${REFSEQ}.seqroi.does NOT exist!";
			fi
		else
			ErrorLog "${SAMPLEID}" "$ROOT_PATH/hitsPerBaseFiles/${SAMPLEID}.seqroi.uniq.map already exists and -f was not used!"
		fi
	else
		ErrorLog "${SAMPLEID}" "The analysis is only supported for h.sapiens so far!";
	fi
else
	ErrorLog "${SAMPLEID}" "READS has to be set true to run the analysis!";
fi



if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in baseMapping against sequenced roi (unique aligned reads).";
else
	SuccessLog "${SAMPLEID}" "Passed baseMapping against sequenced roi (unique aligned reads).";
fi
