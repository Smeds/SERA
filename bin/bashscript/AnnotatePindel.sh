#!/bin/bash
#
# Script to convert Pindel output to VCF, filter to normal and annotate with annovar
#SBATCH -p core  -n 1
#SBATCH -t 1:00:00
##SBATCH --qos=short -t 00:15:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh

# Check if the directory exists, if not create it
if [[ ! -d "$ROOT_PATH/PindelAnnovar" ]]; then
	mkdir $ROOT_PATH/PindelAnnovar;
fi

if [[ ! -d "$ROOT_PATH/PindelOut/Unused" ]]; then
	mkdir $ROOT_PATH/PindelOut/Unused;
fi

if [[ ! -d "$ROOT_PATH/PindelAnnovarOutput" ]]; then
	mkdir $ROOT_PATH/PindelAnnovarOutput;
fi

if [[ ! -e $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.annovaroutput.txt || ! -z $FORCE ]]; then

	if [[ ${NORMAL_SAMPLEID} != "false" ]]; then

		if [[ -e $ROOT_PATH/PindelOut/${SAMPLEID}.indels_BP ]]; then

			#Move unwanted analysisfiles
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_BP $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_CloseEndMapped $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_INT $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_INT_final $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_INV $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_LI $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_RP $ROOT_PATH/PindelOut/Unused
			mv $ROOT_PATH/PindelOut/${SAMPLEID}.indels_TD $ROOT_PATH/PindelOut/Unused

		fi


			#Convert Pindel deletions to VCF
			grep chr $ROOT_PATH/PindelOut/${SAMPLEID}.indels_D > $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_D
			grep chr $ROOT_PATH/PindelOut/${SAMPLEID}.indels_SI > $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_SI

			$ROOT_PATH_PINDEL/pindel2vcf -P $ROOT_PATH/PindelOut/${SAMPLEID}.indels -r $GENOME_FASTA_REF -R Hg19 -d 20090401 -v $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.vcf -is 3

			# Check if ffpe or plasma crieteria should be used
            if [[ ${METHOD} == "haloplex" ]]; then
                # Filter indels
    			python2.7 $SERA_PATH/bin/pythonscript/FilterPindelVCF_onlyTumor.py -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.vcf -o $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput -del $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_D -ins $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_SI -r min ${PINDEL_ANNOVAR_FLAGS}

    			#Annotate with annovar
    			perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,

    			perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -s -am

            elif [[ ${METHOD} == "swift" ]]; then
                if [[ ${TYPE} == "ffpe" ]]; then
                    #Filter indels
                    python2.7 $SERA_PATH/bin/pythonscript/FilterPindelVCF_onlyTumor.py -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.vcf -o $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput -del $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_D -ins $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_SI -r min ${PINDEL_ANNOVAR_FLAGS}

                    #Annotate with annovar
                    perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,

                    perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -s
                elif [[ ${TYPE} == "plasma" ]]; then
                    #Filter indels
                    python2.7 $SERA_PATH/bin/pythonscript/FilterPindelVCF_onlyTumor.py -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.vcf -o $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput -del $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_D -ins $ROOT_PATH/PindelOut/${SAMPLEID}.indels.onlyChr_SI -r min ${PINDEL_ANNOVAR_PLASMA_FLAGS}

                else
                    ErrorLog "${SAMPLEID}" "Only implemented for ffpe and plasma so far!!!"
                fi

                if [[ -e $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput ]]; then

                    #Annotate with annovar
                    perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,

                    perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/PindelAnnovarOutput/${SAMPLEID}.pindel.singleSample.annovarOutput -s
                else
                    ErrorLog "${SAMPLEID}" "Pindel annovar input file $ROOT_PATH/PindelAnnovar/${SAMPLEID}.pindel.filtered.annovarInput does NOT exist!!!"
                fi
            else
                ErrorLog "${SAMPLEID}" "Only implemented for METHOD haloplex and swift so far!!!"
            fi
	else

			SuccessLog "${SAMPLEID}" "Normal set to false, Pindel will not be Annotated.";

	fi

else

	ErrorLog "${SAMPLEID}" "Pindel output already exists and Force was not used.";

fi
