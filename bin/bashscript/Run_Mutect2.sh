#!/bin/bash
#
# Script running bwa
#
#SBATCH -p core -n 8
#SBATCH --nodes=1
#SBATCH -t 02:00:00
##SBATCH --qos=short
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

if [[ $GLOBALS == "MORIARTY" ]]; then
    SNIC_TMP="${ROOT_PATH}/tmp";
    if [[ ! -d ${SNIC_TMP} ]]; then
        mkdir ${SNIC_TMP};
    fi
fi

declare -a StringArray=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
 
if [[ $PLATFORM = "Illumina" ]]; then
    # If MATE_PAIR is set to true in the input file
    if [[ "$MATE_PAIR" == "true" ]]; then
        if [[ ! -d $ROOT_PATH/Mutect2 ]]; then
            mkdir $ROOT_PATH/Mutect2;
        fi

        BAMOUT=""
        #Run Mutect2 and filtering for each chromosome
        for chrom in ${StringArray[@]}; do
            if [[ ${METHOD} == "haloplex" ]]; then
                singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4-4.1.7.0--py38_0.simg gatk --java-options '-Xmx4g' Mutect2 -L chr$chrom -I $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam -bamout $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.bam -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf -DF NotDuplicateReadFilter --disable-adaptive-pruning && singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4.1.4.1.simg java -Xmx4g -jar /gatk/gatk-package-4.1.4.1-local.jar FilterMutectCalls -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -V $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf &&  bgzip $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf && tabix $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf.gz &
           elif [[ ${METHOD} == "swift" ]]; then
                singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4-4.1.7.0--py38_0.simg gatk --java-options '-Xmx4g' Mutect2 -L chr$chrom -I $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam -bamout $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.bam -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf && singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4.1.4.1.simg java -Xmx4g -jar /gatk/gatk-package-4.1.4.1-local.jar FilterMutectCalls -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -V $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf &&  bgzip $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf && tabix $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf.gz &
           else
                ErrorLog "${SAMPLEID}" "Only implemented for METHOD haloplex and swift so far!!!"
           fi

            BAMOUT="${BAMOUT} $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.bam"
        done
        wait        

        samtools merge $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.bam $BAMOUT && samtools index $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.bam 
        module load bcftools/1.3.1
        bcftools concat --allow-overlaps --remove-duplicates $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.*.vcf.gz > $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.vcf.gz
        bgzip -d $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.vcf.gz
        rm $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.*.vcf.*
        rm $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.*.ba*
        rm $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.*.raw.vcf
#        if [[ ${METHOD} == "haloplex" ]]; then
#            python $SERA_PATH/bin/pythonscript/FilterMutect2.py $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput 20 0.01
#            #Annotate with annovar
#            perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,
#            perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.singleSample.annovarOutput -s -am
#    
#        elif [[ ${METHOD} == "swift" ]]; then
#            if [[ ${TYPE} == "ffpe" ]]; then
#                python $SERA_PATH/bin/pythonscript/FilterMutect2.py $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput 20 0.01
#                #Annotate with annovar
#                perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,
#                perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.singleSample.annovarOutput -s
#            elif [[ ${TYPE} == "plasma" ]]; then
#                python $SERA_PATH/bin/pythonscript/FilterMutect2.py $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.vcf $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput 20 0.001
#                #Annotate with annovar
#                perl $ROOT_PATH_ANNOVAR/table_annovar.pl $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput $ANNOVAR_MODULES/annovar_humandb/ -protocol refGene,1000g2015aug_eur,snp138,snp138NonFlagged,esp6500siv2_ea,cosmic70,clinvar_20150629 -operation g,f,f,f,f,f,f -nastring "-" -otherinfo -buildver hg19 -remove -arg '-splicing_threshold 5',,,,,,
#                perl $SERA_PATH/bin/perlscript/createAnnovarOutput.pl -i $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.filtered.annovarInput.hg19_multianno.txt -o $ROOT_PATH/Mutect2/${SAMPLEID}.mutect2.sorted.singleSample.annovarOutput -s
#            else
#                ErrorLog "${SAMPLEID}" "Only implemented for ffpe and plasma so far!!!"
#            fi
#        else
#            ErrorLog "${SAMPLEID}" "Only implemented for METHOD haloplex and swift so far!!!"
#        fi

        SuccessLog "${SAMPLEID}" "singularity exec -B /data -B /projects /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4-4.1.7.0--py38_0.simg gatk --java-options '-Xmx4g' Mutect2 -L chr$chrom -I $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam -bamout $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.bam -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf -DF NotDuplicateReadFilter && singularity exec -B /data /projects/wp2/nobackup/Twist_Myeloid/Containers/gatk4.1.4.1.simg java -Xmx4g -jar /gatk/gatk-package-4.1.4.1-local.jar FilterMutectCalls -R /data/ref_genomes/hg19/genome_fasta/hg19.with.mt.fasta -V $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.raw.vcf -O $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf &&  bgzip $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf && tabix $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.${chrom}.vcf.gz &"
        SuccessLog "${SAMPLEID}" "samtools merge $ROOT_PATH/Mutect2/${SAMPLEID}.sorted.bam $ROOT_PATH/Mutect2/${SAMPLEID}.*.sorted.bam"
    else
        ErrorLog "${SAMPLEID}" "Only implemented for paired-end sequencing!";
    fi
fi

if [[ "$?" != "0" ]]; then
    ErrorLog "${SAMPLEID}" "Failed in Mutect2...";
else
    SuccessLog "${SAMPLEID}" "Passed Mutect2";
fi


