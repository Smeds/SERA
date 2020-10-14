#!/bin/bash
#
# Script to convert Pindel output to VCF, filter to normal and annotate with annovar
#SBATCH -p core  -n 1
#SBATCH -t 1:00:00
##SBATCH --mail-user=elin.falk_sorqvist@igp.uu.se
##SBATCH --mail-type=END
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

##SBATCH --qos=short -t 00:15:00

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/ClinicalPositions" ]]; then
	mkdir $ROOT_PATH/ClinicalPositions;
fi

if [[ ${TYPE} == "ffpe" ]]; then
    if [[ ${INDELFILE} != "false" ]]; then
    	if [[ -e $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt ]]; then
    		if [[ -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
    			python2.7 $SERA_PATH/bin/pythonscript/PindelClinicalInfo.py -i $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -o /dev/stdout -s $SAMPLEID -g ${INDELFILE} -chr2nc $NC2chr ${PINDEL_CLINICAL_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt
    		else
    			ErrorLog "${SAMPLEID}" "The input file $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput does not exist!";
    		fi

    		if [[ ${REGIONFILE} != "false" ]]; then
    			# Check if the clinical position output file exists, if so add info to it. Ohterwise create it
    			if [[ -e $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt ]]; then
    				if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
    					python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat /dev/stdin >> $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;

    				elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
    					python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat /dev/stdin >> $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;
    				else
    					ErrorLog "${SAMPLEID}" "Neither the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!";
    				fi
    			else
    				ErrorLog "${SAMPLEID}" "The pindel analysis has not worked and $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt has not been created -> Not running region analysis!";
    			fi

    		else
    			SuccessLog "${SAMPLEID}" "REGIONFILE is set to false -> no region analysis is run!";
    		fi

    	else
    		ErrorLog "${SAMPLEID}" "$ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt does NOT exist!!!";
    	fi
    else
    	if [[ -e $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt ]]; then
    		cp $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt
    		SuccessLog "${SAMPLEID}" "INDELFILE set to false, no indels to add!";

    		if [[ ${REGIONFILE} != "false" ]]; then
    			if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
    				python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt ;

    			elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
    				python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;
    			else
    				ErrorLog "${SAMPLEID}" "Neither the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!";
    			fi
    		else
    			SuccessLog "${SAMPLEID}" "REGIONFILE is set to false -> no region analysis is run!";
    		fi
    	else
    		ErrorLog "$SAMPLEID" "$ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt does NOT exist!!!";
    	fi
    fi
elif [[ ${TYPE} == "plasma" ]]; then
    if [[ ${INDELFILE} != "false" ]]; then
        if [[ -e $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt ]]; then
            if [[ -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput ]]; then
                python2.7 $SERA_PATH/bin/pythonscript/PindelClinicalInfo.py -i $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -o /dev/stdout -s $SAMPLEID -g ${INDELFILE} -chr2nc $NC2chr ${PINDEL_CLINICAL_PLASMA_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt
            else
                ErrorLog "${SAMPLEID}" "The input file $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput does not exist!";
            fi

            if [[ ${REGIONFILE} != "false" ]]; then
                # Check if the clinical position output file exists, if so add info to it. Ohterwise create it
                if [[ -e $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt ]]; then
                    if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
                        python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat /dev/stdin >> $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;

                    elif [[  -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
                        python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat /dev/stdin >> $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;
                    else
                        ErrorLog "${SAMPLEID}" "Neither the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!";
                    fi
                else
                    ErrorLog "${SAMPLEID}" "The pindel analysis has not worked and $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt has not been created -> Not running region analysis!";
                fi

            else
                SuccessLog "${SAMPLEID}" "REGIONFILE is set to false -> no region analysis is run!";
            fi

        else
            ErrorLog "${SAMPLEID}" "$ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt does NOT exist!!!";
        fi
    else
        if [[ -e $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt ]]; then
            cp $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt
            SuccessLog "${SAMPLEID}" "INDELFILE set to false, no indels to add!";

            if [[ ${REGIONFILE} != "false" ]]; then
                if [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations ]]; then
                    python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt ;

                elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput && $ROOT_PATH/SNPmania/${SAMPLEID}.variations ]]; then
                    python2.7 $SERA_PATH/bin/pythonscript/ExtractDepth.py -r ${REGIONFILE} -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -o /dev/stdout -s $SAMPLEID -chr2nc $NC2chr -f $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ${REGION_CLINICAL_FLAGS} | cat $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt /dev/stdin > $ROOT_PATH/ClinicalPositions/${SAMPLEID}.clinicalPositions.txt;
                else
                    ErrorLog "${SAMPLEID}" "Neither the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations nor the combination of $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput & $ROOT_PATH/SNPmania/${SAMPLEID}.variations exist!";
                fi
            else
                SuccessLog "${SAMPLEID}" "REGIONFILE is set to false -> no region analysis is run!";
            fi
        else
            ErrorLog "$SAMPLEID" "$ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt does NOT exist!!!";
        fi
    fi
else
    ErrorLog "${SAMPLEID}" "Only implemented for TYPE ffpe and plasma so far!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in merging clinical info";
else
	SuccessLog "$SAMPLEID" "Passed merging clinical info";
fi
