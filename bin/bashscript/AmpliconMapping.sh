#!/bin/bash
#
# Script creats BlastDB on mate-pair reads.
#
#SBATCH -p core  -n 8
#SBATCH -t 02:00:00
#SBATCH --mail-type=FAIL --mail-user=bioinfo-clinical-genomics-uu@googlegroups.com

##SBATCH --qos=short

. $SERA_PATH/includes/load_modules.sh

# Include functions
. $SERA_PATH/includes/logging.sh;

if [[ ${METHOD} == "haloplex" ]]; then
    if [[ ! -d $ROOT_PATH/AmpliconCoverage ]]; then
    	mkdir $ROOT_PATH/AmpliconCoverage
    fi

    if [[ ! -d $ROOT_PATH/AmpliconMapped ]]; then
    	mkdir $ROOT_PATH/AmpliconMapped
    fi

    # Check that platform is set to Illumina.
    if [[ ${PLATFORM} == "Illumina" ]]; then
    	# Check if the selection file exists in bed-format, if not create
    	if [[ ! -e $ROOT_PATH/refFiles/$REFSEQ.selection.bed || ! -z $FORCE ]]; then
    		awk -f ${SERA_PATH}/bin/awkscript/sedd2bed.awk -v name="$REFSEQ" -v desc="Selection file" $NC2chr $ROOT_PATH/refFiles/$REFSEQ.selection > $ROOT_PATH/refFiles/$REFSEQ.selection.bed;
    	else
    		SuccessLog "${SAMPLEID}" "$ROOT_PATH/refFiles/$REFSEQ.selection.bed already exists and force was not used!";
    	fi
    	# Check that reads are true and molecules are false
    	if [[ ${READS} == "true" && ${MOLECULES} == "false" ]]; then
            # Check if the output file exists
            if [[ ! -e $ROOT_PATH/AmpliconMapped/${SAMPLEID}.aligned.withAmplicon.bam || ! -z $FORCE ]]; then
                # Check if the aligned file from bwa exists
                if [[ -e $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam ]]; then
                    # Querysort bam-file
                    samtools view -h -b -F 0x100 $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | samtools sort -n -@ 8 /dev/stdin $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted
                    # Run ampliconmapping
                    # java -Xmx8g -jar ${SERA_PATH}/bin/java/GenomeAnalysisTKLite_molecules.jar -T MapReadToAmpliconsIlluminaReadPair -R $GENOME_FASTA_REF -I $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam -o $ROOT_PATH/AmpliconCoverage/${SAMPLEID}.amplicon.bed -fragments $ROOT_PATH/refFiles/$REFSEQ.selection.bed -ampAnReads $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam -U ALL -nonunique -allowPotentiallyMisencodedQuals --downsample_to_coverage 90000 -molBarCode 0
                    java -Xmx8g -jar ${SERA_PATH}/bin/java/GenomeAnalysisTKLite_molecules.jar -T MapReadToAmpliconsIlluminaReadPair -R $GENOME_FASTA_REF -I $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam -o $ROOT_PATH/AmpliconCoverage/${SAMPLEID}.amplicon.bed -fragments $ROOT_PATH/refFiles/$REFSEQ.selection.bed -ampAnReads $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -U ALL -nonunique -allowPotentiallyMisencodedQuals --downsample_to_coverage 90000 -molBarCode 0
                    samtools sort -@ 8 $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon;
                    samtools index $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;

                    samtools flagstat $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam > $ROOT_PATH/Bwa/${SAMPLEID}.alignmentStats_noDuplicateReads.txt;
                    # Remove the querysorted file
                    rm $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam;
                    rm $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam;

                    SuccessLog "${SAMPLEID}" "samtools view -h -b -F 0x100 $ROOT_PATH/Bwa/${SAMPLEID}.sorted.bam | samtools sort -n -@ 8 /dev/stdin $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted";
                    SuccessLog "${SAMPLEID}" "java -Xmx8g -jar ${SERA_PATH}/bin/java/GenomeAnalysisTKLite_molecules.jar -T MapReadToAmpliconsIlluminaReadPair -R $GENOME_FASTA_REF -I $ROOT_PATH/AmpliconMapped/${SAMPLEID}.querysorted.bam -o $ROOT_PATH/AmpliconCoverage/${SAMPLEID}.amplicon.bed -fragments $ROOT_PATH/refFiles/$REFSEQ.selection.bed -ampAnReads $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam -U ALL -nonunique -allowPotentiallyMisencodedQuals --downsample_to_coverage 90000 -molBarCode 0";
                    SuccessLog "${SAMPLEID}" "samtools sort -@ 8 $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.unsorted.bam $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon;";
                    SuccessLog "${SAMPLEID}" "samtools index $ROOT_PATH/AmpliconMapped/${SAMPLEID}.withAmplicon.bam;";

                    # To create reference file dictionary
                    #java -Xmx2g -jar /sw/apps/bioinfo/picard/1.141/nestor/picard.jar CreateSequenceDictionary R=/proj/a2013225/private/reference_genomes/BWA.ref.seqs/BWA_0.7.10_hg19_refseqs/hg19.with.mt.fasta O=/proj/a2013225/private/reference_genomes/BWA.ref.seqs/BWA_0.7.10_hg19_refseqs/hg19.with.mt.fasta.dict

                else
                    ErrorLog "${SAMPLEID}" "$ROOT_PATH/AmpliconMapped/${SAMPLEID}.sorted.bam does not exist!";
                fi
            else
                ErrorLog "${SAMPLEID}" "$ROOT_PATH/AmpliconMapped/${SAMPLEID}.aligned.withAmplicon.bam already exists and force was not used!";
            fi
        else
            ErrorLog "${SAMPLEID}" "READS has to be true and MOLECULES false in the input file for the analysis to run!";
        fi
    else
    	ErrorLog ${SAMPLEID} "So far only supported for Illumina!";
    fi
elif [[ ${METHOD} == "swift" ]]; then
    SuccessLog "${SAMPLEID}" "Ampliconmapping is not needed for swift samples - skipping!";

else
    ErrorLog "${SAMPLEID}" "halo and swift are the only methods supported so far!";
fi
# Check if readsVSmolecules worked
if [[ "$?" != "0" ]]; then
	ErrorLog "${SAMPLEID}" "Failed in ampliconMapping...";
else
	SuccessLog "${SAMPLEID}" "Passed ampliconMapping...";
fi
