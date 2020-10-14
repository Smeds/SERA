#!/bin/bash
#
# Script to run jSNPmania
#SBATCH -p core -n 1
#SBATCH -t 30:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

SuccessLog $SAMPLEID "Starts combining Annovar output ...";

# Check if the directory exists, if not create it
if [[ ! -d $ROOT_PATH/AnnovarOutput ]]; then
	mkdir $ROOT_PATH/AnnovarOutput;
fi

if [[ ${READS} == "true" ]]; then
	if [[ ${CALL_TYPE} == "h.sapiens" ]]; then
		# Check for singleSamples without ampliconmapping
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/AnnovarOutput/${REFSEQ}_sinlgeSamples.annovarOutput.txt ]]; then
			rm $ROOT_PATH/AnnovarOutput/${REFSEQ}_singleSamples.annovarOutput.txt;
		fi
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/AnnovarOutput/*singleSample.annovarOutput > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tStrands_A(F+|F-|S+|S-)\tStrands_G(F+|F-|S+|S-)\tStrands_C(F+|F-|S+|S-)\tStrands_T(F+|F-|S+|S-)\tStrands_Ins(F+|F-|S+|S-)\tStrands_Del(F+|F-|S+|S-)\tTranscripts"} {if($1!~/Sample/){print $0}}' $ROOT_PATH/AnnovarOutput/*singleSample.annovarOutput > $ROOT_PATH/AnnovarOutput/${REFSEQ}_singleSamples.annovarOutput.txt;

		else
			ErrorLog "$SAMPLEID" "No files with ending singleSample.annovarOutput were found!";
		fi

		# Check for single samples with ampliconmapping
		# Check if a combined file exist for ampliconmapped, then remove it
		if [[ -e $ROOT_PATH/AnnovarOutput/${REFSEQ}_singleSamples.ampliconmapped.annovarOutput.txt ]]; then
			rm $ROOT_PATH/AnnovarOutput/${REFSEQ}_singleSamples.ampliconmapped.annovarOutput.txt;
		fi
		# If there still exist input files - create combined file
		# If there still exist input files - create combined file
		ls -1 $ROOT_PATH/AnnovarOutput/*singleSample.ampliconmapped.annovarOutput > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tStrands_A(F+|F-|S+|S-)\tStrands_G(F+|F-|S+|S-)\tStrands_C(F+|F-|S+|S-)\tStrands_T(F+|F-|S+|S-)\tStrands_Ins(F+|F-|S+|S-)\tStrands_Del(F+|F-|S+|S-)\t#variant_+_amplicons\t#variant_-_amplicons\t#reference_+_amplicons\t#reference_-_amplicons\tVariant_ampliconinfo\tReference_ampliconinfo\tTranscripts"} {if($1!~/Sample/){print $0}}' $ROOT_PATH/AnnovarOutput/*singleSample.ampliconmapped.annovarOutput > $ROOT_PATH/AnnovarOutput/${REFSEQ}_singleSamples.ampliconmapped.annovarOutput.txt;
		else
			ErrorLog "$SAMPLEID" "No files with ending singleSample.ampliconmapped.annovarOutput were found!";
		fi

		# Check for tumor-normal samples without ampliconmapping
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.annovarOutput.txt ]]; then
			rm $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.annovarOutput.txt;
		fi
		# If there still exist input files
		ls -1 $ROOT_PATH/AnnovarOutput/*tumorNormalSample.annovarOutput > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tTumor_Strands_A(F+|F-|S+|S-)\tTumor_Strands_G(F+|F-|S+|S-)\tTumor_Strands_C(F+|F-|S+|S-)\tTumor_Strands_T(F+|F-|S+|S-)\tTumor_Strands_Ins(F+|F-|S+|S-)\tTumor_Strands_Del(F+|F-|S+|S-)\tNormal_Strands_A(F+|F-|S+|S-)\tNormal_Strands_G(F+|F-|S+|S-)\tNormal_Strands_C(F+|F-|S+|S-)\tNormal_Strands_T(F+|F-|S+|S-)\tTranscripts"} {if($1!~/Sample/){print $0}}' $ROOT_PATH/AnnovarOutput/*tumorNormalSample.annovarOutput > $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.annovarOutput.txt;
		else
			ErrorLog "$SAMPLEID" "No files with ending tumorNormalSample.annovarOutput were found!";
		fi

		# Check for tumor-normal samples with ampliconmapping
		# Check if a combined file exist, then remove it
		if [[ -e $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.ampliconmapped.annovarOutput.txt ]]; then
			rm $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.ampliconmapped.annovarOutput.txt;
		fi
		# If there still exist input files
		ls -1 $ROOT_PATH/AnnovarOutput/*tumorNormalSample.ampliconmapped.annovarOutput > /dev/null 2>&1
		if [[  "$?" = "0" ]]; then
			awk 'BEGIN{print "#Sample\tChr\tStart\tEnd\tReference_base\tVariant_base\tGene\tType\tExonic_type\tVariant_allele_ratio\t#reference_alleles\t#_variant_alleles\tRead_depth\tRatio_in_1000Genome\tdbSNP_id\tClinically_flagged_dbSNP\tESP_6500\tCosmic\tTumor_Strands_A(F+|F-|S+|S-)\tTumor_Strands_G(F+|F-|S+|S-)\tTumor_Strands_C(F+|F-|S+|S-)\tTumor_Strands_T(F+|F-|S+|S-)\tTumor_Strands_Ins(F+|F-|S+|S-)\tTumor_Strands_Del(F+|F-|S+|S-)\tNormal_Strands_A(F+|F-|S+|S-)\tNormal_Strands_G(F+|F-|S+|S-)\tNormal_Strands_C(F+|F-|S+|S-)\tNormal_Strands_T(F+|F-|S+|S-)\t#Tumor_variant_+_amplicons\t#Tumor_variant_-_amplicons\t#Tumor_reference_+_amplicons\t#Tumor_reference_-_amplicons\t#Normal_reference_+_amplicons\t#Normal_reference_-_amplicons\tTumor_Variant_ampliconinfo\tTumor_Reference_ampliconinfo\tReference_ampliconinfo\tTranscripts"} {if($1!~/Sample/){print $0}}' $ROOT_PATH/AnnovarOutput/*tumorNormalSample.ampliconmapped.annovarOutput > $ROOT_PATH/AnnovarOutput/${REFSEQ}_tumorNormalSamples.ampliconmapped.annovarOutput.txt;
		else
			ErrorLog "$SAMPLEID" "No files with ending tumorNormalSample.ampliconmapped.annovarOutput were found!";
		fi
	else
		ErrorLog "$SAMPLEID" "Only supported for h.sapiens so far!";
	fi
else
	ErrorLog "$SAMPLEID" "READS has to be true for the analysis to run!";
fi

if [[ "$?" != "0" ]]; then
	ErrorLog "$SAMPLEID" "Failed in combining Annovar output";
else
	SuccessLog "$SAMPLEID" "Passed combining Annovar output";
fi
