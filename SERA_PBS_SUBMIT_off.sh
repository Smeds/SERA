#!/bin/bash -l

# SERA Analysis Software


##############
### PREPARE ENVIRONMENT
##############

# Get script path
export SERA_PATH=$(dirname `readlink -f $0`)
. $SERA_PATH/includes/prepareEnvironment.sh

##############
### PROGRAM GLOBALS
##############

# program title
title="SlurmSERA"

# array structure:
# id brief_explanation default_status associated_script dependency_on_id run_after_samples(0: do not wait, 1: wait and run all samples, 2: wait and run once)
steps=(
1 "Create ampregion, ampROI, seqregion and seqROI Files" off "Create_region_files_python.sh" false 2 \
2 "Create SNPseq file" off "Create_SNPseq_file.sh" "1" 2 \
5 "Pre-process fastq-files" off "Prepare_fastq_files.sh" false 0 \
9 "Remove adapter and primer (haloplex and swift) or adapter (swift+breast) from reads" off "Run_trimming.sh" "5" 0 \
10 "Remove adapter sequences from reads" off "Run_cutAdapt_swift_all.sh" "5" 0 \
11 "Run FastQC" off "Run_FastQC.sh" "5:9" 0 \
12 "Align with BWA against Genome" off "Run_Bwa_sam.sh" "5:9" 0 \
13 "Remove adapter sequences in Swift runs" off "Run_primerclip.sh" "12" 0 \
14 "Amplicon mapping" off "AmpliconMapping.sh" "13" 0 \
15 "Run jSNPmania" off "Run_jSNPmania.sh" "2:13:14" 0 \
16 "SNPmania output to Annovar input" off "jSNPmania2AnnovarInput.sh" "15" 0 \
17 "Run Annovar" off "Run_Annovar.sh" "16" 0 \
18 "Combine Annovar output to one file" off "Combine_annovarOutput_to_one_file.sh" "17" 2 \
20 "Run Pindel" off "Run_Pindel.sh" "13" 0 \
21 "Annotate Pindel with Annovar" off "AnnotatePindel.sh" "20" 0 \
22 "Run Mutect2" off "Run_Mutect2.sh" "13" 0 \
25 "Output all mutations (hotspots and others)" off "FilterMutations.sh" "17:21" 0 \
26 "Convert Annovar output to vcf-format" off "ConvertAnnovaroutput2vcf.sh" "25" 0 \
30 "Extract MSI markers" off "ExtractMSI.sh" "17:21:25" 0 \
31 "Combine extracted MSI markers to one file" off "Combine_extracted_MSI.sh" "30" 2 \
32 "Extract EGFR information" off "Extract_T790M.sh" "15" 0 \
33 "Combine extracted EGFR information to one file" off "Combine_extracted_EGFR.sh" "32" 2 \
34 "Extract SNP information" off "ExtractSNPinfo.sh" "15" 0 \
35 "Combine extracted SNP information to one file" off "Combine_extracted_SNPs.sh" "34" 2 \
40 "Calculate amplification ratio" off "CalculateAmplificationRatio.sh" "15" 0 \
50 "Calculate fraction covered" off "FractionCovered.sh" "15" 0 \
51 "Calculate mean cov per gene" off "Calculate_Gene_Cov.sh" "15" 0 \
52 "Calculate mean cov per region (exon)" off "Calculate_Region_Cov.sh" "15" 0 \
60 "Create hits per base file from SNPmania file SeqRegion" off "SNPmania2hitsPerBase_seqregion.sh" "15" 0 \
61 "Create hits per base file from SNPmania file SeqRoi" off "SNPmania2hitsPerBase_seqroi.sh" "15" 0 \
62 "Produce plot-values SeqRegion Against SeqRegion" off "Create_plotValues_seqregion_against_seqregion.sh" "60" 0 \
63 "Produce plot-values SeqRoi Against SeqRoi" off "Create_plotValues_seqroi_against_seqroi.sh" "61" 0 \
64 "Create stenberg nusbaum input file" off "makeGnuplotInputFilePost.sh" "62:63" 2 \
65 "Make stenberg nusbaum plots" off "MakePlots.sh" "64" 2 \
66 "Calculate Mean Sequencing Depth" off "Calculate_mean_depth.sh" "60:61" 0 \
70 "Run Pileup" off "Run_PileUp.sh" "4" 0 \
71 "Separate on- and offtarget bases from pileup" off "PileUp2onOffTarget.sh" "70" 0 \
72 "Extract specificities" off "ExtractSpecificity.sh" "71" 0 \
73 "Calculate variant allele ratio per strand" off "StrandVarAlleleRatio.sh" "11" 0 \
100 "Copy results to OUTBOX" off "Move_results2outbox.sh" "17:21:22:25:26:31:33:35:40" 2 \
101 "Copy results to storage" off "Move_results2storage.sh" "17:21:22:25:26:31:33:35:40" 2);

let NUMBOFSTEPS=${#steps[@]}-1;

# adjust steps based on input file
. $SERA_PATH/includes/adjustSteps.sh



##############
### START PROGRAM
##############

# Start interface
if [[ -z "$INTERFACE" ]]; then
    . $SERA_PATH/includes/commandLineInterface.sh
else
    . $SERA_PATH/includes/dialogInterface.sh
fi

# Run steps for samples
. $SERA_PATH/includes/submitScripts.sh

#submit successfull exit status
#exit 0;
