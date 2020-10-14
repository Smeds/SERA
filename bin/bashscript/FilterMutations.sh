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

if [[ ! -d $ROOT_PATH/FilteredMutations ]]; then
    mkdir $ROOT_PATH/FilteredMutations;
fi

ANNOVARFILE=""
PINDELANNOVARFILE=""
SNPMANIAFILE=""
ampMapped=""

if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
	if [[ ${READS} == "true" ]]; then
        if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
            if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
                if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput && -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput";
                    SNPMANIAFILE="$ROOT_PATH/SNPmania/${SAMPLEID}.variations";

                elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput && -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
                    ANNOVARFILE="$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput";
                    PINDELANNOVARFILE="$ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput";
                    SNPMANIAFILE="$ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations";
                    ampMapped=" -a ";

                else
                    ErrorLog "$SAMPLEID" "All inputfiles needed can not be found!";
                fi

                if [[ -e $ANNOVARFILE && -e $PINDELANNOVARFILE && -e $SNPMANIAFILE ]]; then
                    if [[ (! -e $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv) || ($FORCE == "true") ]]; then
                        if [[ ${TYPE} == "ffpe" ]]; then
                            if [[ ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ]]; then
                                echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_FLAGS}"
                                cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_FLAGS}
                            else
                                echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python2.7 $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_FLAGS}"
                                cat $ANNOVARFILE $PINDELANNOVARFILE | python2.7 $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_FLAGS}
                            fi
                        elif [[ ${TYPE} == "plasma" ]]; then
                            if [[ ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ]]; then
                                echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_PLASMA_FLAGS}"
                                cat $ANNOVARFILE $PINDELANNOVARFILE | python $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.02 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_PLASMA_FLAGS}
                            else
                                echo "cat $ANNOVARFILE $PINDELANNOVARFILE | python2.7 $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_PLASMA_FLAGS}"
                                cat $ANNOVARFILE $PINDELANNOVARFILE | python2.7 $SERA_PATH/bin/pythonscript/FilterAllMutations.py -i ${HOTSPOTFILE} -o $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -f /dev/stdin -b $BLACKLIST_FILE -s ${SAMPLEID} -chr2nc $NC2chr -g 0.01 ${ampMapped} -t ${MAIN_TRANSCRIPTS} -v ${SNPMANIAFILE} -m $ROOT_PATH/refFiles/MultipleBp_variations.csv ${MUTATION_PLASMA_FLAGS}
                            fi
                        else
                            ErrorLog "$SAMPLEID" "Only implemented for TYPE ffpe and plasma so far!";
                        fi
                    else
                        ErrorLog "$SAMPLEID" "Output file $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv already exists and Force was not used!";
                    fi
                else
                    ErrorLog "$SAMPLEID" "The annovar output file $ANNOVARFILE and/or pindel annovar output file $PINDELANNOVARFILE and/or SNPmania varition file $SNPMANIAFILE do NOT exist!!!";
                fi
            else
                    ErrorLog "$SAMPLEID" "Only implemented for NORMAL_SAMPLEID annovar so far!";
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
