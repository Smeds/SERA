#!/bin/bash
#
# Script to export aligment archives as BAM-files.
#
#SBATCH -p core  -n 1
#SBATCH -t 15:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

if [[ ! -d $ROOT_PATH/plotPdfs ]]; then
	mkdir $ROOT_PATH/plotPdfs;
fi

if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		FILE="$ROOT_PATH/plotPdfs/plotFile.txt";
		# if empty, create file
		if [[ ! -e $FILE || ! -z $FORCE ]]; then
			echo -e "#Plotname type datafile1,title:datafile2,title\n#type=stenberg or nusbaum." > $FILE;

			# options to check for
			types=(
			"seqroi_seqroi.uniq.nusbaum"
			"seqroi_seqroi.uniq.stenberg"
			"seqregion_seqregion.uniq.nusbaum"
			"seqregion_seqregion.uniq.stenberg"
			);

			# if there are no of this reference in file, append
			for type in ${types[@]}; do

				# check that there are files of this type
				if [[ `ls -l $ROOT_PATH/plotValues/ | grep ${type} | wc -l` != 0 ]]; then

					PLOT=`echo $type | awk -F . '{print $NF}'`;

					# temp concat str
					CONCAT="${REFSEQ}.$type\t$PLOT\t$ROOT_PATH/plotValues/${SAMPLEID}.${type},${SAMPLEID}";
					added=${SAMPLEID};

					# loop over all files and add this plot type
					samples=( `ls $ROOT_PATH/plotValues/ | grep $type | awk '{split($0,a,"."); print a[1]}' | uniq | xargs` );
					for SAMPLEID in ${samples[@]}; do
						if [[ $added != $SAMPLEID ]]; then
							CONCAT="$CONCAT:$ROOT_PATH/plotValues/${SAMPLEID}.${type},${SAMPLEID}";
						fi
					done
					echo -e $CONCAT >> $FILE;

				fi

			done
		fi
	fi
fi

if [[ "$?" != "0" ]]; then
	ErrorLog $REFSEQ "Failed in creating gnuplot input file";
else
	SuccessLog $REFSEQ "Inputfile created and can be located at $FILE for further modifications";
fi
