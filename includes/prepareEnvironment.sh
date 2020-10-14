#!/bin/bash -i

# Read input options
ROOT_PATH_INPUT=""; INPUT_FILE=""; INTERFACE="" CNV_FILE=""; FORCE="";
while getopts "h:p:i:u:f" OPTION; do
	case $OPTION in
		h) usage; exit 1 ;;
		p) ROOT_PATH_INPUT=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		u) INTERFACE="gui" ;;
#		c) CNV_FILE=$OPTARG ;;
		f) export FORCE="true" ;;
		?) usage; exit ;;
     esac
done

# validate input
. $SERA_PATH/includes/validateInput.sh

# Go to output folder
cd $ROOT_PATH_INPUT;
export ROOT_PATH=$(pwd);

# Include and validate  inputfile data
. $SERA_PATH/includes/validateInputFile.sh;
COUNT=-1;

# Read input file.
. $INPUT_FILE

# Get number of samples
count=${#SAMPLEID_ARR_[@]};
let NUMBEROFSAMPLES=count-1;


## Include server-specific globals
#if [[ "$SOFTWARE" = "SLURM" ]]; then
#	. $SERA_PATH/config/globalsUppnex.sh;
#else
#	. $SERA_PATH/config/globalsSmaug.sh;
#fi

# create files for logging
if [[ ! -d "$ROOT_PATH/slurmOutput" && $SOFTWARE = "SLURM" ]]; then
	mkdir slurmOutput;
fi
#if [[ ! -e "$LOG_FILE" ]]; then
#	echo -e "<?xml version='1.0'?>\n<document info=\"$title\">\n</document>" > $LOG_FILE;
#fi

# Create log file for pipeline
rpath=($(echo $ROOT_PATH | tr "/" " "));
rtpath=${rpath[ ${#rpath[@]}-1]};
PIPELINE_LOG="$ROOT_PATH/PipelineLog_"$rtpath".txt";
if [[ -e $ROOT_PATH/PipelineLog_"$rtpath".txt ]]; then
    PIPELINE_LOG="$ROOT_PATH/PipelineLog_"$rtpath"_2.txt"
fi
export PIPELINE_LOG;

echo -e "#####General#####" > $PIPELINE_LOG;
echo -e "Analysis_date="`date +%Y-%m-%d` >> $PIPELINE_LOG;
echo -e "Analysis_start_time="`date +%H:%M:%S` >> $PIPELINE_LOG;
echo -e "Number_of_samples="$count >> $PIPELINE_LOG;
echo -e "\n#####Programs#####" >> $PIPELINE_LOG;
echo -e "SERA_version="$SERA_PATH >> $PIPELINE_LOG;
echo -e "SNPmania_version="$ROOT_PATH_JSNPMANIA >> $PIPELINE_LOG;
echo -e "Annovar_version="$ROOT_PATH_ANNOVAR >> $PIPELINE_LOG;
echo -e "Pindel_version="$ROOT_PATH_PINDEL >> $PIPELINE_LOG;
echo -e "\n#####Files#####" >> $PIPELINE_LOG;
echo -e "Genome_used="$GENOME_REF >> $PIPELINE_LOG;
echo -e "Fasta_reference="$GENOME_FASTA_REF >> $PIPELINE_LOG;
echo -e "Blacklist_file="$BLACKLIST_FILE >> $PIPELINE_LOG;
echo -e "Convert_NC_to_chr="$NC2chr >> $PIPELINE_LOG;
echo -e "BlastDB="$BLAST_DB >> $PIPELINE_LOG;
echo -e "Main_transcripts="$MAIN_TRANSCRIPTS >> $PIPELINE_LOG;
echo -e "Illumina_adapter_file_for_trimmomatic="$ILLUMINA_ADAPTER_TRIMMOMATIC >> $PIPELINE_LOG;
echo -e "Trimmomatic_jarfile="$TRIMMOMATIC_JAR >> $PIPELINE_LOG;

echo -e "="$ >> $PIPELINE_LOG;

echo -e "\n#####Modules#####" >> $PIPELINE_LOG;
moduleList=$( (module list) 2>&1 );
for mod in $moduleList
do
    if [[ ${mod} != *\)* && ${mod} != *Currently* && ${mod} != *Loaded* && ${mod} != *Module* ]]; then
        echo -e "Loaded_module="$mod >> $PIPELINE_LOG;
    fi
done

echo -e "\n#####Samples#####" >> $PIPELINE_LOG;

# inputfile information
let ADJ_SAMPLENUMBER=NUMBEROFSAMPLES+1;
MP="single end"; if [[ $MATE_PAIR = "true" ]]; then MP="paired end"; fi;
if [[ $SOFTWARE = "SLURM" ]]; then INFO="Running on $SOFTWARE with project $UPPNEX_PROJECT_ID. ";
else INFO="Running in pure bash environment. "; fi
INFO=${INFO}"Analyzing $ADJ_SAMPLENUMBER samples from $DESIGN_TYPE design with $READ_LENGTH bp $PLATFORM $MP sequencing data. ";
if [[ ! -z $FORCE ]]; then INFO=${INFO}"All previous data will be overwritten for associated steps. "; fi
if [[ $ALIGN_NONUNIQUE == "true" ]]; then INFO=${INFO}"Analyzing both unique and nonunique. "; else INFO=${INFO}"Analyzing only unique reads. "; fi
if [[ ! -z $CALL_TYPE ]]; then INFO=${INFO}"SNP calling on $CALL_TYPE."; fi

#prepare variable array for concatenation for running post submit scripts
POST="";
