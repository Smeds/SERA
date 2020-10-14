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


if [[ ${NORMAL_SAMPLEID} != "false" ]]; then
	if [[ ${READS} == "true" ]]; then
			if [[ ${NORMAL_SAMPLEID} == "annovar" ]]; then
				if [[ -e $ROOT_PATH/Annovar/${SAMPLEID}.annovarInput ]]; then
					perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Annovar/${SAMPLEID}.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,
					perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Annovar/${SAMPLEID}.annovarInput.hg19_multianno.txt -o $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.annovarOutput -s

				elif [[ -e $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput ]]; then
					perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,

					perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput.hg19_multianno.txt -o $ROOT_PATH/AnnovarOutput/${SAMPLEID}.singleSample.ampliconmapped.annovarOutput -s -am

				else
					ErrorLog "$SAMPLEID" "Using NORMAL_SAMPLEID=annovar none of the possible annovar input files exist ($ROOT_PATH/Annovar/${SAMPLEID}.annovarInput or $ROOT_PATH/Annovar/${SAMPLEID}.ampliconmapped.annovarInput)!";
				fi
			else
				if [[ -e $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput ]]; then
					perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,

					perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput.hg19_multianno.txt -o $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.annovarOutput -tn

				elif [[ -e $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput ]]; then
					perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,
					perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput.hg19_multianno.txt -o $ROOT_PATH/AnnovarOutput/${SAMPLEID}_${NORMAL_SAMPLEID}.tumorNormalSample.ampliconmapped.annovarOutput -tn -am

				else
					ErrorLog "$SAMPLEID" "Using NORMAL_SAMPLEID=normal_sampleId none of the possible annovar input files exist ($ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.annovarInput or $ROOT_PATH/Annovar/${SAMPLEID}_${NORMAL_SAMPLEID}.ampliconmapped.annovarInput)!";
				fi
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
