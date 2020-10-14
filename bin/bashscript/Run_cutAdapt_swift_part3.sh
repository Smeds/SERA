#!/bin/bash -l

#SBATCH -p core  -n 6
#SBATCH -t 01:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Start running cutadapt";

export TRIM_LOG="${ROOT_PATH}/seqdata/${SAMPLEID}.trim.log"
PREFIX="${SNIC_TMP}/${SAMPLEID}";

cputhreads=12;

ptrim()
{
    fqt1=$1
    fqt2=$2
    tprefix=${fqt1%%_R1_001.fastq}

    cutadaptFile3prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_3ptrim.fa";

    #3’trim
    cutadapt \
        -a file:$cutadaptFile3prim \
        -o ${tprefix}_R2_primertrimd.fq -p ${tprefix}_R1_primertrimd.fq \
        $fqt2 $fqt1 --minimum-length 40 -e 0.12 >> $TRIM_LOG;


}

r1reformat()
{
    tr '~' '\n' < $1 > ${1}_R1_001.fastq
}

r2reformat()
{
    tr '~' '\n' < $1 > ${1}_R2_001.fastq
}

export -f ptrim
export -f r1reformat
export -f r2reformat

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/seqdata" ]]; then
    mkdir $ROOT_PATH/seqdata;
fi

if [[ $PLATFORM = "Illumina" ]]; then
    # If MATE_PAIR is set to true in the input file
    if [[ "$MATE_PAIR" == "true" ]]; then
        if [[ ${METHOD} == "swift" ]]; then
            # Check that input files exist
            if [[ -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp2.fastq.gz && -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp2.fastq.gz ]]; then
                # Check that output file doesn't exist then run cutAdapt, if it does print error message
                if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then

                    # Set parameters for data to continue trimming
                    PE1_G_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp2.fastq.gz";
                    PE2_G_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp2.fastq.gz";
                    gunzip -f ${PE1_G_T};
                    gunzip -f ${PE2_G_T};
                    PE1_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp2.fastq";
                    PE2_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp2.fastq";

                    # convert fastq format to one line per record for splitting
                    paste - - - - < $PE1_T | tr '\t' '~' > ${PE1_T}.tmp1;
                    paste - - - - < $PE2_T | tr '\t' '~' > ${PE2_T}.tmp1;

                    # get number of fastq records in sample before converting back to fastq format
                    l=$(wc -l ${PE1_T}.tmp1 | awk '{print $1}')
                    chunklinecnt=$(( $l / $cputhreads ))

                    # split re-formatted fastq files into chunks
                    split -d -l $chunklinecnt ${PE1_T}.tmp1 ${PREFIX}_r1split
                    split -d -l $chunklinecnt ${PE2_T}.tmp1 ${PREFIX}_r2split

                    # 20160722 NOTE: adding the sample-specific prefix to the temporary fastq
                    #                files should fix the sample concatenation bug
                    # convert each chunk back to fastq format
                    parallel r1reformat ::: ${PREFIX}_r1split*
                    parallel r2reformat ::: ${PREFIX}_r2split*

                    ls ${PREFIX}_r1*.fastq > ${SNIC_TMP}/r1infiles
                    ls ${PREFIX}_r2*.fastq > ${SNIC_TMP}/r2infiles

                    # run parallel on paired chunks of fastq files with ptrim() function
                    parallel --xapply ptrim {1} {2} ::: $(cat ${SNIC_TMP}/r1infiles) ::: $(cat ${SNIC_TMP}/r2infiles)

                    # concatenate primer-trimmed fastq chunks
                    cat ${PREFIX}*_R1_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz;
                    cat ${PREFIX}*_R2_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz

                    # remove all intermediate files
                    rm ${PREFIX}_r[1,2]split*
                    rm ${PE1_T}
                    rm ${PE2_T}
                    rm ${PE1_T}.tmp1
                    rm ${PE2_T}.tmp1

                else
                    ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
                fi
            else
                ErrorLog "${SAMPLEID}" "One or both of the input files (${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz & ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz) do not exist!";
            fi
        else
            ErrorLog "This step only has to be run for METHOD swift!";
        fi
     else
        ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
    fi
fi


if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in cutadapt step 2...";
else
    SuccessLog "${SAMPLEID}" "Passed cutadapt step 2";
fi
