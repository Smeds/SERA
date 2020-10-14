#!/bin/bash
#
# Script to calculate the amplification ratio
#SBATCH -p core  -n 1
#SBATCH -t 1:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/Amplification" ]]; then
	mkdir $ROOT_PATH/Amplification;
fi


if [[ ${TISSUE} == "lung" || ${TISSUE} == "colon" ]]; then
	if [[ ${AMPLIFICATIONFILE} != "false" && ${BACKGROUNDFILE} != "false" ]]; then
		if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
			python2.7 $SERA_PATH/bin/pythonscript/AmplificationCalculation.py -a ${AMPLIFICATIONFILE} -b ${BACKGROUNDFILE} -s $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -chr2nc $NC2chr -t ${TISSUE} -sample ${SAMPLEID} -o $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt
			sort -k3,3 -k4,4n $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt > $ROOT_PATH/Amplification/${SAMPLEID}.amplification.sorted.txt
			rm $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt

		elif [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
			python2.7 $SERA_PATH/bin/pythonscript/AmplificationCalculation.py -a ${AMPLIFICATIONFILE} -b ${BACKGROUNDFILE} -s $ROOT_PATH/SNPmania/${SAMPLEID}.variations -chr2nc $NC2chr -t ${TISSUE} -sample ${SAMPLEID} -o $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt
			sort -k3,3 -k4,4n $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt > $ROOT_PATH/Amplification/${SAMPLEID}.amplification.sorted.txt
			rm $ROOT_PATH/Amplification/${SAMPLEID}.amplification.txt

		else
			ErrorLog "${SAMPLEID}" "Neither input file $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!";
		fi

	else
		ErrorLog "${SAMPLEID}" "Either Amplifcation and/or background file is set to false!";
	fi
else
	WarningLog "${SAMPLEID}" "Only run for lung and colon!";
fi


if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in calculating amplification ratio";
else
	SuccessLog "$SAMPLEID" "Passed calculating amplification ratio";
fi
