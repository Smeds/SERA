#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts Annovar ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/AnnovarOutput ]]; then
    mkdir $ROOT_PATH/AnnovarOutput;
fi
if [[ ! -d $ROOT_PATH/FilteredAnnovarOutput ]]; then
    mkdir $ROOT_PATH/FilteredAnnovarOutput;
fi

ANNOVARFILE=""
PINDELANNOVARFILE=""

if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
    if [[ ${READS} == "true" ]]; then
        if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
            if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
                if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput";

                elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.ampliconmapped.singleSample.annovarOutput";

                fi

                if [[ -e $ANNOVARFILE && $PINDELANNOVARFILE ]]; then
                    if [[ ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ]]; then
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -a
                            echo "cat  $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -a"

                            if [[ ${METHOD} == "swift" ]]; then
                                cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -a
                                echo "cat  $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -a"

                        else
                            cat  $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -k $KEEPFILE -a
                            echo "cat  $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -k $KEEPFILE -a"
                        fi
                    else
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -a
                            echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -a"
                        else
                            cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -k $KEEPFILE -a
                            echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -k $KEEPFILE -a"
                        fi
                    fi
                else
                    ErrorLog "$SAMPLEID" "The annovar output file $ANNOVARFILE and/or pindel annovar output file $PINDELANNOVARFILE do NOT exist!!!";
                fi

            else
                if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then

                    if [[ ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ]]; then
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02
                        else
                            cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -k $KEEPFILE
                        fi
                    else
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01
                        else
                            cat $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -k $KEEPFILE
                        fi
                    fi

                elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput && $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then

                    if [[ ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ]]; then
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02
                        else
                            cat $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.02 -k $KEEPFILE
                        fi
                    else
                        # Check if there are particular regions we want to keep
                        if [[ $KEEPFILE == "false" ]]; then
                            cat $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01
                        else
                            cat $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput | python $SERA_PATH/bin/pythonscript/FilterAnnovarOutput.py -i /dev/stdin -o $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.filtered.annovarOutput -b $BLACKLIST_FILE -g 0.01 -k $KEEPFILE
                        fi
                    fi
                else
                    ErrorLog "$SAMPLEID" "Using NORMAL_SAMPLEID=normal_sampleId none of the possible annovar input files exist ($ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput or $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput)!";
                fi
            fi
        else
            ErrorLog "$SAMPLEID" "Only supported for call_type h.sapiens so far!";
        fi
    else
        ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
    fi

else
    ErrorLog "$SAMPLEID" "Normal_sampleid is false -> annovar is not run!";
fi

if [[ "$?" != "0" ]]; then
    ErrorLog "$SAMPLEID" "Failed in running Annovar";
else
    SuccessLog "$SAMPLEID" "Passed running Annovar";
fi
