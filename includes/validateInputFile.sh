#!/bin/bash
#
# Include file for reading SERA input file.
#

COUNT=-1;

# Read input file.
. $INPUT_FILE

# Get number of samples
count=${#SAMPLEID_ARR_[@]};
let NUMBEROFSAMPLES=count-1;

# error handling
ERROR=0;
for i in $(seq 0 $NUMBEROFSAMPLES); do

	# check that required files are available for all input samples
	if [[ ! -e ${ROIFILE_ARR_[${i}]} ]]; then fileioerror ${ROIFILE_ARR_[${i}]};  ERROR=1; fi;
	if [[ ! -e ${SELECTIONFILE_ARR_[${i}]} ]]; then fileioerror ${SELECTIONFILE_ARR_[${i}]}; ERROR=1; fi;
#	if [[ $CNV_FILE && ! -e ${CNV_FILE} ]]; then fileioerror $CNV_FILE; ERROR=1; fi;
	
	# check that optional files are available
#	if [[ ${CONTROL_REGIONS_ARR_[${i}]} != "false" && ! -e ${CONTROL_REGIONS_ARR_[${i}]}.selection ]]; then fileioerror ${CONTROL_REGIONS_ARR_[${i}]}; ERROR=1; fi;
#	if [[ ${HAPMAP_SNP_REF_ARR_[${i}]} != "false" && ! -e ${HAPMAP_SNP_REF_ARR_[${i}]} ]]; then fileioerror ${HAPMAP_SNP_REF_ARR_[${i}]}; ERROR=1; fi;

	# check that no required options are null
	if [[ ! ${RAWDATA_PE1_ARR_[${i}]} ]]; then settingserror "Read file have to be specified for ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
	if [[ ! ${ROIFILE_ARR_[${i}]} ]]; then settingserror "ROI file have to be specified for ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
	if [[ ! ${SELECTIONFILE_ARR_[${i}]} ]]; then settingserror "Selection file have to be specified for ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
	if [[ ! ${SAMPLEID_ARR_[${i}]} ]]; then settingserror "A id have to be set for sample  ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
	if [[ ! ${REFSEQ_ARR_[${i}]} ]]; then settingserror "A design reference have to be set for sample ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
#	if [[ ! $NUMBER_OF_MM ]]; then settingserror "Number of mismatches have to be specified"; ERROR=1; fi;
	if [[ ! $READ_LENGTH ]]; then settingserror "Read length have to be specified"; ERROR=1; fi;

	# check that options are allowed
	if [[ $PLATFORM != "SOLiD" && $PLATFORM != "Illumina" ]]; then settingserror "Platform have to be specified as SOLiD or Illumina"; ERROR=1; fi;
	if [[ $DESIGN_TYPE != "PCR" && $DESIGN_TYPE != "MDA" ]]; then settingserror "Design type have to be specified as PCR or MDA"; ERROR=1; fi;
	if [[ $MATE_PAIR != "true" && $MATE_PAIR != "false" ]]; then settingserror "Paired end reads have to be specified as true or false"; ERROR=1; fi;
	if [[ $READS != "true" && $READS != "false" ]]; then settingserror "Analysis of reads have to be specified as true or false"; ERROR=1; fi;
	if [[ $MOLECULES != "true" && $MOLECULES != "false" ]]; then settingserror "Analysis of molecules have to be specified as true or false"; ERROR=1; fi;

	# check that variables are set and not empty
	if [[ -z "$ALIGNERFLAGS_GENOME" ]]; then settingserror "Mosaik aligner flags must be specified for genome"; ERROR=1; fi;
	if [[ -z "$SNPMANIAFLAGS" ]]; then settingserror "jSNPmaniaFlags must be set"; ERROR=1; fi;
	if [[ -z "$CALL_TYPE" ]]; then settingserror "Call type must be specified"; ERROR=1; fi;
	if [[ -z "$PINDEL_FLAGS" ]]; then settingserror "Pindel alignment flags have to be set (or false)"; ERROR=1; fi;
	if [[ -z "$TUMOR_NORMAL_FLAGS" ]]; then settingserror "Tumor-normal flags have to be set (or false)"; ERROR=1; fi;
	if [[ -z "$ANNOVAR_FLAGS" ]]; then settingserror "Annovar flags has to be set (or false)"; ERROR=1; fi;
	if [[ -z "$PINDEL_ANNOVAR_FLAGS" ]]; then settingserror "Pindel annovar flags has to be set (or false)"; ERROR=1; fi;
	if [[ -z "$MUTATION_FLAGS" ]]; then settingserror "Mutation flags has to be set (or false)"; ERROR=1; fi;
	if [[ -z "$MUTATION_PLASMA_FLAGS" ]]; then settingserror "Mutation plasma flags has to be set (or false)"; ERROR=1; fi;
	if [[ -z "${BARCODE_I7_ARR_[${i}]}" ]]; then settingserror "Sequencing tag must be specified or set to false for sample ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
    if [[ -z "${BARCODE_I5_ARR_[${i}]}" ]]; then settingserror "Sequencing tag must be specified or set to false for sample ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;


	# more complex
	if [[ -z "$SOFTWARE" ]]; then settingserror "Software must be specified."; ERROR=1; fi;
	if [[ $SOFTWARE = "SLURM" && ! $UPPNEX_PROJECT_ID ]]; then settingserror "Uppnex project id must be specified when running SLURM"; ERROR=1; fi;
	if [[ $MATE_PAIR = "true" && ${RAWDATA_PE2_ARR_[${i}]} = "false" ]]; then settingserror "A second read file have to be specified for paired ends in sample ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;
	if [[ $NORMAL_SAMPLEID != "false" && $NORMAL_SAMPLEID != "annovar" && $TUMOR_NORMAL_FLAGS = "false" ]]; then settingserror "If Normal sample id isn't set to false or annovar the TUMOR_NORMAL_FLAGS have to be set (not false)"; ERROR=1; fi
    if [[ ${METHOD_ARR_[${i}]} == "swift" && ${CUTADAPT_PREFIX_ARR_[${i}]} == "false" ]]; then settingserror "If method is swift then cutadapt_prefix also has to be set (not false) for sample ${SAMPLEID_ARR_[${i}]}"; ERROR=1; fi;

	# check for genefile
#	if [[ ! -z $NORMAL_SAMPLEID && ! -e $GENEFILE && $GENEFILE != "false" && ! -z $GENEFILE ]]; then fileioerror "Gene file ${GENEFIGILE_ARR_[{$i}]} do not exist, set to false or specify correct location"; ERROR=1; fi;
	
	# check for sample information file if specified
#	if [[ -n $SAMPLES_INFORMATION && $SAMPLES_INFORMATION != "false" && ! -e $SAMPLES_INFORMATION ]]; then fileioerror "File $SAMPLES_INFORMATION was not found, please provide new or set variable to false"; ERROR=1; fi;

	# check whether reads are in valid format and readjust if necessarry
	function checkData {
	
		if [[ $PLATFORM = "SOLiD" ]]; then
			if [[ -e ${1} ]]; then settingserror "SOLiD data should point to file stub w/o qual.gz or csfasta.gz for ${SAMPLEID_ARR_[${i}]}"; ERROR=1;
			elif [[ -e "${1}.csfasta" || -e "${1}.qual" ]]; then settingserror "SOLiD data expected in gzip format for ${SAMPLEID_ARR_[${i}]}."; ERROR=1;
			elif [[ ! -e "${1}.csfasta.gz" || ! -e "${1}.qual.gz" ]]; then fileioerror ${1}.csfasta.gz; fileioerror ${1}.qual.gz; ERROR=1;
			fi
		elif [[ $PLATFORM = "Illumina" ]]; then
			if [[ ! -e ${1} ]]; then fileioerror ${1}; ERROR=1;
			elif [[ -e ${1} && "`file -ib ${1} | awk '{print $1}'`" != "application/x-gzip;" ]]; then settingserror "Illumina data expected in gzip format for ${SAMPLEID_ARR_[${i}]}"; ERROR=1;
			fi
		fi
	}
	checkData ${RAWDATA_PE1_ARR_[${i}]};
	checkData ${RAWDATA_PE2_ARR_[${i}]};

done;

# more cheks for checking array length and amount of samples
if [[ $count = "0" ]]; then settingserror "No samples have been specified"; ERROR=1; fi;
if [[ "$COUNT" != "$NUMBEROFSAMPLES" ]]; then settingserror "The array in inputfile is not of correct length."; ERROR=1; fi

if [[ $ERROR = 1 ]]; then
	kill -SIGINT $$;
fi
