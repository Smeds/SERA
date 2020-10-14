#!/bin/bash
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core -n 2
#SBATCH -t 1:00:00
#module load bioinfo-tools blast/2.2.25
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

# calculate numbers or ppn:s.
#NUMPROC=`expr $(cat /proc/cpuinfo | grep processor | wc -l) - 1`;
NUMPROC="4";

# Check if BLAST analysis is possible.
if [[ $PLATFORM == "SOLiD" ]]; then
	ErrorLog "${SAMPLEID}" "${0##*/}" "Analysis is not supported for SOLiD data.";
	exit 0;
fi
if [[ $MATE_PAIR != "true" ]]; then
	ErrorLog "${SAMPLEID}" "${0##*/}" "Mate-pair data needed for this analysis.";
	exit 0;
fi

# Check if BlastDB exist.
if [[ ! -f "$ROOT_PATH/BLAST/Databases/${SAMPLEID}.ReadsDB.nhr" ]] && [[ ! -f "$ROOT_PATH/BLAST/Databases/${SAMPLEID}.ReadsDB.nal" ]]; then
	ErrorLog "${SAMPLEID}" "${0##*/}" "No BlastDB found for sample ID.";
	exit 0;
fi

# Check if folder exist.
if [[ ! -d "$ROOT_PATH/BLAST/Runs" ]]; then
	mkdir $ROOT_PATH/BLAST/Runs;
fi

SuccessLog "${SAMPLEID}" "${0##*/}" "Starting Blast analysis.";

# Running BlastAll
cd $ROOT_PATH/BLAST/Runs;

if [[ ${READS} == "true" ]]; then
	zcat $ROOT_PATH/BLAST/Databases/${SAMPLEID}.fragments.fasta.gz | blastall -p blastn -a ${NUMPROC} -b 1000000 -d $ROOT_PATH/BLAST/Databases/${SAMPLEID}.ReadsDB -i /dev/stdin -m 8 | gzip - > $ROOT_PATH/BLAST/Runs/${SAMPLEID}.blast.gz;

	SuccessLog "${SAMPLEID}" "${0##*/}" "Filtering BLAST Results with max mm = ${BLAST_MM} bp, and minimal match length of ${BLAST_ML} bp.";
	zcat $ROOT_PATH/BLAST/Runs/${SAMPLEID}.blast.gz | awk -v ml=${BLAST_ML} -v mm=${BLAST_MM} '{if(($4>=ml)&&($5<=mm)&&($6==0)&&($9==1)) {print $1;}}' | sort | uniq -c | awk '{print $2 "\t" $1;}' | sort -k2nr | gzip - > $ROOT_PATH/BLAST/Runs/${SAMPLEID}.blast.filtered.gz;

	# Creating hits per fragment file
	echo -e "FragmentID\tReadHits" > $ROOT_PATH/BLAST/${SAMPLEID}.blast.hits.txt;
	for i in $(zcat $ROOT_PATH/BLAST/Databases/${SAMPLEID}.fragments.fasta.gz | awk '{if(substr($0,1,1)==">") { split($0,s,">"); print s[2];}}'); do

		fragId=$i;

		hits=$(zcat $ROOT_PATH/BLAST/Runs/${SAMPLEID}.blast.filtered.gz | grep "${fragId}\b" | cut -f 2);
		if [[ "$hits" != "" ]]; then
			echo -e "$fragId\t$hits" >> $ROOT_PATH/BLAST/${SAMPLEID}.blast.hits.txt;
		else
			echo -e "$fragId\t0" >> $ROOT_PATH/BLAST/${SAMPLEID}.blast.hits.txt;
		fi
	done;
fi
# Creates a list of all mapped reads.
#zcat $ROOT_PATH/BLAST/Runs/${SAMPLEID}.blast.gz | awk -v ml=${BLAST_ML} -v mm=${BLAST_MM} '{if(($4>=ml)&&($5<=mm)&&($6==0)&&($9==1)) {print $2;}}' | sort | uniq -c | awk 'BEGIN{print "ReadID\tFragmentHits";}{print $2 "\t" $1;}' | sort -k2nr | gzip - > $ROOT_PATH/BLAST/Runs/${SAMPLEID}.mapped.reads.gz;

SuccessLog "${SAMPLEID}" "${0##*/}" "Blast analysis done.";
