#!/bin/bash -l
#SBATCH -p core  -n 1
#SBATCH -t 02:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog "${SAMPLEID}" "Start running cutadapt";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/seqdata" ]]; then
	mkdir $ROOT_PATH/seqdata;
fi

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

if [[ $PLATFORM = "Illumina" ]]; then

	# get sequencing tags
	. $SERA_PATH/config/sequencingTags.sh;

	# If MATE_PAIR is set to true in the input file
	if [[ "$MATE_PAIR" == "true" ]]; then
		# Check that output file doesn't exist then run cutAdapt, if it does print error message
		if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz && ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz || ! -z $FORCE ]]; then

		    if [[ ${METHOD} == "haloplex" ]]; then
    			cutadapt -a $tTag -A `perl $SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $READ1 $READ2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;

	       		SuccessLog "${SAMPLEID}" "cutadapt -a $tTag -A `perl $SERA_PATH/bin/perlscript/reverseComplement.pl $fTag` -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz --minimum-length 1 $READ1 $READ2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;"
            elif [[ ${METHOD} == "swift" ]]; then
                if [[ ${CUTADAPT_PREFIX} != "false" ]]; then

                    TMP1_PE1="$SNIC_TMP/pe1.tmp1.fastq.gz";
                    TMP1_PE2="$SNIC_TMP/pe2.tmp1.fastq.gz";

                    TMP2_PE1="$SNIC_TMP/pe1.tmp2.fastq.gz";
                    TMP2_PE2="$SNIC_TMP/pe2.tmp2.fastq.gz";

                    TEMP1_PE1="$SNIC_TMP/pe1.temp1.fastq.gz";
                    TEMP1_PE2="$SNIC_TMP/pe2.temp1.fastq.gz";

                    TEMP2_PE1="$SNIC_TMP/pe1.temp2.fastq.gz";
                    TEMP2_PE2="$SNIC_TMP/pe2.temp2.fastq.gz";

                    TEMP3_PE1="$SNIC_TMP/pe1.temp3.fastq.gz";
                    TEMP3_PE2="$SNIC_TMP/pe2.temp3.fastq.gz";

                    cutadaptAdapterSeq="${ROOT_PATH}/refFiles/TruSeq3-PE-2.fa";
                    cutadaptFile5prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_5ptrim.fa";
                    cutadaptFile3prim="${ROOT_PATH}/refFiles/${CUTADAPT_PREFIX}_3ptrim.fa";

                    cutadapt -g file:${cutadaptAdapterSeq} -e 0.12 -o $TMP1_PE1 -p $TMP1_PE2 --minimum-length 1  $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:${cutadaptAdapterSeq} -e 0.12 -o $TMP1_PE1 -p $TMP1_PE2 --minimum-length 1  $RAWDATA_PE1 $RAWDATA_PE2 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";

                    cutadapt -g file:${cutadaptAdapterSeq} -e 0.12 -o $TMP2_PE2 -p $TMP2_PE1 --minimum-length 1  $TMP1_PE2 $TMP1_PE1 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:${cutadaptAdapterSeq} -e 0.12 -o $TMP2_PE2 -p $TMP2_PE1 --minimum-length 1  $TMP1_PE2 $TMP1_PE1 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";

                    cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 $TMP1_PE1 $TMP1_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "cutadapt -g file:$cutadaptFile5prim -o $TEMP1_PE1 -p $TEMP1_PE2 $TMP_PE1 $TMP_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";

                    cutadapt -g file:$cutadaptFile5prim -o $TEMP2_PE2 -p $TEMP2_PE1 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile5prim -o $TEMP2_PE2 -p $TEMP2_PE2 $TEMP1_PE2 $TEMP1_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";

                    cutadapt -a file:$cutadaptFile3prim -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE1 $TEMP2_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile3prim -o $TEMP3_PE1 -p $TEMP3_PE2 $TEMP2_PE1 $TEMP2_PE2 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";

                    cutadapt -a file:$cutadaptFile3prim -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
                    SuccessLog "${SAMPLEID}" "cutadapt -g file:$cutadaptFile3prim -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz -p ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz $TEMP3_PE2 $TEMP3_PE1 --minimum-length 40 -e 0.12 >> ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log";
                else
                    ErrorLog "${SAMPLEID}" "CUTADAPT_PREFIX has to be set for swift in order to be able to run cutadapt!";
                fi
            else
                ErrorLog "${SAMPLEID}" "Only implemented for METHOD haloplex and swift so far!";
            fi

		else
			ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz and ${ROOT_PATH}/seqdata/${SAMPLEID}.read2.fastq.gz already exists and force was NOT used!";
		fi
	else
		if [[ ! -e ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz ]]; then

            cutadapt -a $tTag -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz --minimum-length 1 $RAWDATA_PE1 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;
            SuccessLog "${SAMPLEID}" "cutadapt -a $tTag -o ${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz --minimum-length 1 $RAWDATA_PE1 > ${ROOT_PATH}/seqdata/${SAMPLEID}.cutadapt.log;"

        else
            ErrorLog "${SAMPLEID}" "${ROOT_PATH}/seqdata/${SAMPLEID}.read1.fastq.gz already exists and force was NOT used!";
        fi
	fi
fi


if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in cutadapt...";
else
	SuccessLog "${SAMPLEID}" "Passed cutadapt";
fi
