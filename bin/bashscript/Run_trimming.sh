#!/bin/bash -l

#SBATCH -p core -n 6
#SBATCH --nodes=1
#SBATCH -t 12:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

export TRIM_LOG="${ROOT_PATH}/seqdata/${SAMPLEID}.trim.log"

if [[ $GLOBALS == "MORIARTY" ]] || [[ $GLOBALS == "MARVIN" ]] ; then
    SNIC_TMP="${ROOT_PATH}/tmp";
    if [[ ! -d ${SNIC_TMP} ]]; then
        mkdir ${SNIC_TMP};
    fi
fi

PREFIX="${SNIC_TMP}/${SAMPLEID}";

cputhreads=12;

fastq_files_r1=($(echo "$RAWDATA_PE1" | tr " " "\n"));

if [[ ${#fastq_files_r1[@]} > 1 ]];
then
    if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R1_001.fastq.gz | head -1)" ];
    then
        READ1=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R1_001.fastq.gz);
    else
        ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read1, please pre-process data!...";
    fi
    if [ -n "$(find ${ROOT_PATH}/seqdata -name ${SAMPLEID}_S*_L000_R2_001.fastq.gz | head -1)" ];
    then
        READ2=$(ls ${ROOT_PATH}/seqdata/${SAMPLEID}_S*_L000_R2_001.fastq.gz);
    else
        ErrorLog "${SAMPLEID}" "Multiple lanes for sample ${SAMPLEID}, read2, please pre-process data!...";
    fi
else
    READ1=$RAWDATA_PE1;
    READ2=$RAWDATA_PE2;
fi

ptrim()
{
    fqt1=$1
    fqt2=$2
    tprefix=${fqt1%%_R1_001.fastq}

    cutadaptFile5prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_5ptrim.fa";
    cutadaptFile3prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_3ptrim.fa";

    #5’trim
    cutadapt \
        -g file:$cutadaptFile5prim \
        -o ${tprefix}_tmpR1.fq -p ${tprefix}_tmpR2.fq \
    $fqt1 $fqt2 --minimum-length 40 -e 0.12 >> $TRIM_LOG;

    #5’trim
    cutadapt \
        -g file:$cutadaptFile5prim \
        -o ${tprefix}_5ptmpR2.fq -p ${tprefix}_5ptmpR1.fq \
        ${tprefix}_tmpR2.fq ${tprefix}_tmpR1.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

        #3' trim
    cutadapt \
        -a file:$cutadaptFile3prim \
        -o ${tprefix}_tmp3R1.fq -p ${tprefix}_tmp3R2.fq \
        ${tprefix}_5ptmpR1.fq ${tprefix}_5ptmpR2.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

    #3’trim
    cutadapt \
        -a file:$cutadaptFile3prim \
        -o ${tprefix}_R2_primertrimd.fq -p ${tprefix}_R1_primertrimd.fq \
        ${tprefix}_tmp3R2.fq ${tprefix}_tmp3R1.fq --minimum-length 40 -e 0.12 >> $TRIM_LOG;

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


        if [[ ${METHOD} == "haloplex" ]]; then
                # Get sequencing tags
                . $SERA_PATH/config/sequencingTags.sh;
            # Check that output file doesn't exist then run cutAdapt, if it does print error message
            if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then

                # Run cutadapt
                cutadapt -a $tTag -A `perl $SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 ${READ1} ${READ2} > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;

            else
            ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
            fi

        elif [[ ${METHOD} == "swift" && ( ${CUTADAPT_PREFIX} == "cp288_masterfile_191114" || ${CUTADAPT_PREFIX} == "Accel-Amplicon-Plus_Lung_Cancer_masterfile" || ${CUTADAPT_PREFIX} == "18-2132_EGFR_MID_Masterfile_mod20191002") ]]; then
            # Check that output file doesn't exist then run trimmomatic, if it does print error message
            if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz || ! -z $FORCE ]]; then

                ILLUMINA_ADAPTER_TRIMMOMATIC="${ROOT_PATH}/refFiles/TruSeq3-PE-2.fa";
                PE1_G_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz";
                PE2_G_T="${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz";

                java -Xmx32g -Xms16g -jar ${TRIMMOMATIC_JAR} PE \
                    -threads 8 -trimlog $TRIM_LOG \
                    ${READ1} ${READ2} ${PE1_G_T} ${PREFIX}_unpaired_R1.fq.gz \
                    ${PE2_G_T} ${PREFIX}_unpaired_R2.fq.gz \
                    ILLUMINACLIP:${ILLUMINA_ADAPTER_TRIMMOMATIC}:2:30:10 \
                    MINLEN:30

            else
                ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz already exists and force was NOT used!";
            fi
        elif [[ ${METHOD} == "swift" ]]; then 
            # Check that output file doesn't exist then run cutAdapt, if it does print error message
            if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz || ! -z $FORCE ]]; then

                ILLUMINA_ADAPTER_TRIMMOMATIC="${ROOT_PATH}/refFiles/TruSeq3-PE-2.fa";
                PE1_G_T="${PREFIX}_R1_trimd.fq.gz";
                PE2_G_T="${PREFIX}_R2_trimd.fq.gz";

                java -Xmx32g -Xms16g -jar ${TRIMMOMATIC_JAR} PE \
                    -threads 8 -trimlog $TRIM_LOG \
                    ${READ1} ${READ2} ${PE1_G_T} ${PREFIX}_unpaired_R1.fq.gz \
                    ${PE2_G_T} ${PREFIX}_unpaired_R2.fq.gz \
                    ILLUMINACLIP:${ILLUMINA_ADAPTER_TRIMMOMATIC}:2:30:10 \
                    MINLEN:30

                gunzip -f ${PE1_G_T};
                gunzip -f ${PE2_G_T};
                PE1_T="${PREFIX}_R1_trimd.fq";
                PE2_T="${PREFIX}_R2_trimd.fq";

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

                ls ${PREFIX}_r1*.fastq > ${PREFIX}_r1infiles
                ls ${PREFIX}_r2*.fastq > ${PREFIX}_r2infiles

                # run parallel on paired chunks of fastq files with ptrim() function
                parallel --xapply ptrim {1} {2} ::: $(cat ${PREFIX}_r1infiles) ::: $(cat ${PREFIX}_r2infiles)

                # concatenate primer-trimmed fastq chunks
                cat ${PREFIX}_*_R1_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz;
                cat ${PREFIX}_*_R2_primertrimd.fq | gzip -f > ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz

                # remove all intermediate files
                rm ${PREFIX}_*
            else
                ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.tmp.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.tmp.fastq.gz already exists and force was NOT used!";
            fi

        else
            ErrorLog "${SAMPLEID}" "Only implemented for method haloplex and swift so far!!!";
        fi
    else
        ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
    fi
fi


if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in cutadapt (haloplex and swift) or trimmomatic (swift+breast)...";
else
    SuccessLog "${SAMPLEID}" "Passed cutadapt (haloplex and swift) or trimmomatic (swift+breast)";
fi
