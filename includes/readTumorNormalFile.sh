#!/bin/bash
#
# Include file for reading SERA input file.
#

NUMBEROFSAMPLES="-1";  # Holds the total number of samples 0 based.

TUMOR="N/A";      # Holds the id of the tumor sample.
NORMAL="N/A";     # Holds the id of the normal sample.
MINCOV="N/A";	  # The min cov in normal sample to keep position
MINALLELEDIFF="N/A"; # The min difference in reference allele ratio to keep the position
PLATFORM="N/A";
#RAWDATA="N/A";       # Raw data files from the seq. run.
#REFSEQ="N/A";        # Reference sequence file.
#ALIGNERFLAGS="N/A":  # Flags set for aligner.
#ROIFILE="N/A";       # Holds path to the ROI file.
#SELECTIONFILE="N/A"; # Holds path to the selection file.

function readTumorNormalFile {

	input=$1; # read input file
	count=0;

	while read line
	do
		
		firstSign=`echo $line | awk '{print substr(\$0,0,1);}'`;

		if [[ "$firstSign" != "#" ]]; then
			
			TUMOR[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$1}'`;
	          NORMAL[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$2}'`;
			MINCOV[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$3}'`;
			MINALLELEDIFF[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$4}'`;
			PLATFORM[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$5}'`;
			#RAWDATA[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$3}'`;
	                #REFSEQ[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$4}'`;
			#ALIGNERFLAGS[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$5}'`;
			#ROIFILE[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$6}'`;
			#SELECTIONFILE[${count}]=`echo $line | awk 'BEGIN{FS="#"}{print \$7}'`;
			
			let count=count+1;
		fi

	done < $input

	# Set correct number of samples.
	let count=count-1;
	NUMBEROFSAMPLES=$count; 
}