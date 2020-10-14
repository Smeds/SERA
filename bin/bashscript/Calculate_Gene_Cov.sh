#!/bin/bash
#
# Script to run jSNPmania
#SBATCH -p core  -n 1
#SBATCH -t 15:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts calculating gene coverage ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/GeneCov" ]]; then
	mkdir $ROOT_PATH/GeneCov;
fi


# Start with checking that the reference file exists!
if [[ $FILE_FORMAT == "bed" ]]; then
	if [[ ${READS} == "true" ]]; then
		if [[ ${DESIGN_TYPE} == "PCR" ]]; then
			if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
				python2.7 $SERA_PATH/bin/pythonscript/calculateGeneCov_bed.py -i ${ROIFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o /dev/stdout -g | sort -k3,4 > $ROOT_PATH/GeneCov/${SAMPLEID}.geneCov.txt
			elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
				python2.7 $SERA_PATH/bin/pythonscript/calculateGeneCov_bed.py -i ${ROIFILE }-v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o /dev/stdout | sort -k3,4 > $ROOT_PATH/GeneCov/${SAMPLEID}.geneCov.txt
			else
				 ErrorLog "$SAMPLEID" "No variation file exists!"
			fi
		else
			 ErrorLog "$SAMPLEID" "Only implemented for design type PCR!"
		fi
	else
		 ErrorLog "$SAMPLEID" "Only implemented for reads!"
	fi
else
	 ErrorLog "$SAMPLEID" "Only implemented for bed-format!"
fi

if [[ "$?" != "0" ]]; then
	ErrorLog $SAMPLEID "Failed in Calculating gene coverage";
else
	SuccessLog $SAMPLEID "Passed Calculation gene coverage";
fi
