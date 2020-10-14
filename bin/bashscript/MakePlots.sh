#!/bin/bash

#SBATCH -p core  -n 1
#SBATCH -t 00:15:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com
#
# Creating pdf plots.
#

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog ${REFSEQ} "Starts generating plot files...";

# Check if data exist.
if [[ ! -d "$ROOT_PATH/plotValues" ]]; then
	ErrorLog ${REFSEQ} "Failed to locate a plotValues folder.";
	exit 1;
fi



# Making and plotting data.
#echo "perl $SERA_PATH/bin/perlscript/makeGnuplotFile.pl -i $SERA_PATH/plotFile.xml -o $ROOT_PATH/plotPdfs -p $ROOT_PATH/plotValues -plot;"
#perl $SERA_PATH/bin/perlscript/makeGnuplotFile_old.pl -i $SERA_PATH/plotFileOld -o $ROOT_PATH/plotPdfs;
if [[ -d $ROOT_PATH/plotPdfs || ! -z $FORCE ]]; then
	if [[ ${READS} == "true" ]]; then
		if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
			perl $SERA_PATH/bin/perlscript/makeGnuplotFile_old.pl -i $ROOT_PATH/plotPdfs/plotFile.txt -o $ROOT_PATH/plotPdfs;
		fi
	fi
	# Run Gnuplot on all scripts
	cd $ROOT_PATH/plotPdfs;
	find . -name "*.gnuplot" -exec gnuplot '{}' \;

	# Check if makeGnuplotFile.pl worked
	if [[ "$?" != "0" ]]; then
		ErrorLog $REFSEQ "Nusbaum and stenberg plots failed to be generated.";
	else
		SuccessLog $REFSEQ "Nusbaum and stenberg plots successfully created.";
	fi
else
	ErrorLog $REFSEQ "Found previous plots, skipping step.";
fi
