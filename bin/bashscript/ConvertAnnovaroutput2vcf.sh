#!/bin/bash
#
# Script to convert Annovar output to VCF
#SBATCH -p core  -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts converting Annovar output to vcf...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/vcfOutput ]]; then
	mkdir $ROOT_PATH/vcfOutput;
fi
ANNOVARFILE=""
if [[ ${READS} == "true" ]]; then
    if [[ (! -e $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf) || ($FORCE == "true") ]]; then
        # Check which output file from Annovar that exists
    	if [[ -e $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv ]]; then
            python ${SERA_PATH}/bin/pythonscript/Annovar2vcf_filteredMutations.py -m $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.all.vcf

            awk '{if(/^#/){print($0)}else{if(!/-,-$/){print($0)}}}' $ROOT_PATH/vcfOutput/${SAMPLEID}.all.vcf > $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf

            if [[ ${METHOD} == "swift" && ( ${TISSUE} == "ovarial" || ${TISSUE} == "breast" ) ]]; then
                awk 'BEGIN{FS="\t"}{if($1!~/^Sample/ && $15>=0.05){print $0}}' $ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv | python ${SERA_PATH}/bin/pythonscript/Annovar2vcf_filteredMutations.py -m /dev/stdin -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}_vaf0.05.vcf
            fi
        else
            if [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput ]]; then
                ANNOVARFILE=$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput

    	    elif [[ -e $ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput ]]; then
                ANNOVARFILE=$ROOT_PATH/FilteredAnnovarOutput/${SAMPLEID}.singleSample.filtered.annovarOutput

            elif [[ -e $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput ]]; then
                ANNOVARFILE=$ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput

    	    else
    		    ErrorLog "$SAMPLEID" "No output file from Annovar was found!";
    	    fi

            if [[ -e $ANNOVARFILE ]]; then
                    python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v ${ANNOVARFILE} -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}.all.vcf

                    awk '{if(/^#/){print($0)}else{if(!/-,-$/){print($0)}}}' $ROOT_PATH/vcfOutput/${SAMPLEID}.all.vcf > $ROOT_PATH/vcfOutput/${SAMPLEID}.vcf

                if [[ ${METHOD} == "swift" && ${TISSUE} == "ovarial" ]]; then
                    awk 'BEGIN{FS="\t"}{if($1!~/^Sample/ && $10>=0.05){print $0}}' ${ANNOVARFILE} | python ${SERA_PATH}/bin/pythonscript/Annovar2vcf.py -v /dev/stdin -s $ROOT_PATH/refFiles/${REFSEQ}.ampregion.SNPseq -chr2nc $NC2chr -o $ROOT_PATH/vcfOutput/${SAMPLEID}._vaf0.05.vcf
                fi
            else
                ErrorLog "$SAMPLEID" "Neither filtered mutations output ($ROOT_PATH/FilteredMutations/${SAMPLEID}.filteredMutations.tsv) nor the output file from Annovar ($ANNOVARFILE) was not found!";
            fi
        fi
    else
        ErrorLog "$SAMPLEID" "$ROOT_PATH/vcfOutput/${SAMPLEID}.vcf already exists and Force was not used!";
    fi
else
	ErrorLog "$SAMPLEID" "READS has to be true for the analysis to run!";

fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in converting Annovar output to vcf";
else
	SuccessLog "$SAMPLEID" "Passed converting Annovar output to vcf";
fi
