#!/bin/bash -l

#SBATCH -p core  -n 4
#SBATCH -t 02:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

fastq_files_r1=($(echo "$RAWDATA_PE1" | tr " " "\n"));

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/seqdata" ]]; then
    mkdir $ROOT_PATH/seqdata;
fi

if [[ ${#fastq_files_r1[@]} > 1 ]];
then
    pre_filename=$(basename ${fastq_files_r1[0]} | sed -e 's/_L[0-9]\+_R1_001\.fastq\.gz//')
    echo "zcat $RAWDATA_PE2 | pigz -p 4 > ${ROOT_PATH}/seqdata/${pre_filename}_L000_R1_001.fastq.gz";
    zcat $RAWDATA_PE1 | pigz -p 4 > "$ROOT_PATH/seqdata/${pre_filename}_L000_R1_001.fastq.gz";
    echo "zcat $RAWDATA_PE2 | pigz -p 4 > $ROOT_PATH/seqdata/${pre_filename}_L000_R2_001.fastq.gz";
    zcat $RAWDATA_PE2 | pigz -p 4 > "$ROOT_PATH/seqdata/${pre_filename}_L000_R2_001.fastq.gz";
else
    echo " $RAWDATA_PE1 $ROOT_PATH/seqdata/";
    cp $RAWDATA_PE1 "$ROOT_PATH/seqdata/";
    echo " $RAWDATA_PE1 $ROOT_PATH/seqdata/";
    cp $RAWDATA_PE2 "$ROOT_PATH/seqdata/";
fi


if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed pre-process fastq.";
else
    SuccessLog "${SAMPLEID}" "Passed pre-process fastq.";
fi
