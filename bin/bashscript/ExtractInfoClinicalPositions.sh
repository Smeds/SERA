#!/bin/bash
#
# Script to run jSNPmania
##SBATCH --qos=short
#SBATCH -p core  -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts extrating info about clinical positions ...";

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/ClinicalVariants" ]]; then
	mkdir $ROOT_PATH/ClinicalVariants;
fi

if [[ ${HOTSPOTFILE}!="false" ]]; then
	# Check that READS are true
	if [[ ${READS} == "true" ]]; then
		# Check that the call type are set to h.sapiens
		if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
            if [[ ${TYPE} == "ffpe" ]]; then
                if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions ]]; then
    				echo "perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt"
    				perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt  -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt

    			elif  [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.deletions ]]; then
    				perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt  -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt

    			else
    				ErrorLog "$SAMPLEID" "All files needed to run the analysis is not available!";
    			fi
            elif [[ ${TYPE} == "plasma" ]]; then
                if [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions ]]; then
                    echo "perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_PLASMA_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt"
                    perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.ampliconmapped.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt  -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_PLASMA_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt

                elif  [[ -e $ROOT_PATH/SNPmania/${SAMPLEID}.variations && -e $ROOT_PATH/SNPmania/${SAMPLEID}.insertions && -e $ROOT_PATH/SNPmania/${SAMPLEID}.deletions ]]; then
                    echo "perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt  -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_PLASMA_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt"
                    perl $SERA_PATH/bin/perlscript/ClinicalPositions.pl -v $ROOT_PATH/SNPmania/${SAMPLEID}.variations -i $ROOT_PATH/SNPmania/${SAMPLEID}.insertions -d $ROOT_PATH/SNPmania/${SAMPLEID}.deletions -o $ROOT_PATH/ClinicalVariants/${SAMPLEID}.clinicalVariants.txt  -s ${SAMPLEID} -c ${HOTSPOTFILE} ${CLINICAL_PLASMA_FLAGS} >> $ROOT_PATH/ClinicalVariants/extractClinicalVariants_output.txt

                else
                    ErrorLog "$SAMPLEID" "All files needed to run the analysis is not available!";
                fi
            else
                WarningLog "$SAMPLEID" "Only implemented for TYPE ffpe and plasma so far!!!"
            fi
		else
			ErrorLog "$SAMPLEID" "The analysis is only supported for CALL_TYPE h.sapiens!";
		fi

	else
		ErrorLog "$SAMPLEID" "READS has to be true to run the analysis!";
	fi
else
	SuccessLog "$SAMPLEID" "HOTSPOTFILE set to false, no clinical positions to extract!";

fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in extracting info for clinical positions";
else
	SuccessLog "$SAMPLEID" "Passed extracting info for clinical positions";
fi
