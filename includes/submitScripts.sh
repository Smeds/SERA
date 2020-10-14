#!/bin/bash

function runStep {
	
	# if submitting to slurm cluster
	if [[ "$SOFTWARE" = "SLURM" ]]; then
		
		#transform dependency list to sbatch parameter
		DEP="";
		if [[ ${3} != "false" ]]; then
			DEP_ARR=( `echo ${3} | awk 'BEGIN{FS=":"}{for (i=1; i<=NF; i++) printf $i" "}'` );
			
			# if elements found in array, add the correct dependency
			for depstep in ${DEP_ARR[@]}; do
				if [[ ! -z "${RUN_STEP[${depstep}]}" ]]; then
					if [[ ! -z ${DEP} ]]; then
						DEP_TMP=$DEP
						DEP="$DEP_TMP:`echo ${RUN_STEP[${depstep}]} | awk '{print $4}'`";
						DEP_TMP=""
					else
						DEP=":`echo ${RUN_STEP[${depstep}]} | awk '{print $4}'`";
					fi
				fi
				
				# if it is a POST step, add the correct dependency from the array with all stepts and dependencies to the POST array
				if [[ ${4} = 1 || ${4} = 2 ]]; then
					if [[ ! -z "${ALL_STEP[${depstep}]}" ]]; then
						for dependency in ${ALL_STEP[${depstep}]}; do
							POST+=":"$dependency;
						done;
					fi;
				fi;
			done;

			# if it is a POST step, add the correct dependency based on POST array
			if [[ ${4} = 1 || ${4} = 2 ]]&&[[ ! -z $POST ]]; then

				DEP="-d afterok${POST}";
#				echo "Running post steps with dependency $DEP";	# TMP REMOVE ME!

			# if there are jobs in queue that this step depends upon, add dependency
			elif [[ ! -z $DEP ]]; then
				DEP="-d afterok${DEP}";
			fi;
		fi
		
		# submit scripts to SLURM
		SUBMIT_TEXT=$(sbatch -p core -J ${SAMPLEID}-S${1} -e $ROOT_PATH/slurmOutput/${SAMPLEID}-S${1}.error -o $ROOT_PATH/slurmOutput/${SAMPLEID}-S${1}.output $SLURM_MISC_SETTINGS $DEP $SERA_PATH/bin/bashscript/${2});
		SUBMIT_TEXT_ARR=($SUBMIT_TEXT);
		
		RUN_STEP[${1}]=$SUBMIT_TEXT
				
		if [[ ! -z ALL_STEP[${1}] ]]; then
			ALL_STEP[${1}]+=" "${SUBMIT_TEXT_ARR[3]}
		else
			ALL_STEP[${1}]=${SUBMIT_TEXT_ARR[3]}
		fi
	
	# if not on SLURM, submit with regular bash script
	else
		source ${SERA_PATH}/bin/bashscript/${2};
	fi;
}


# SUBMIT SCRIPTS
if [[ "$SOFTWARE" = "SLURM" ]]; then echo "Submitting to cluster..."; else echo "Running scripts..."; fi
for i in $(seq 0 $NUMBEROFSAMPLES); do
	
	#export sample-specific variables
	export SAMPLEID=${SAMPLEID_ARR_[${i}]} BARCODE_I7=${BARCODE_I7_ARR_[${i}]} BARCODE_I5=${BARCODE_I5_ARR_[${i}]} RAWDATA_PE1=${RAWDATA_PE1_ARR_[${i}]} RAWDATA_PE2=${RAWDATA_PE2_ARR_[${i}]} RAWDATA_INDEX=${RAWDATA_INDEX_ARR_[${i}]} REFSEQ=${REFSEQ_ARR_[${i}]} ROIFILE=${ROIFILE_ARR_[${i}]} SELECTIONFILE=${SELECTIONFILE_ARR_[${i}]} NORMAL_SAMPLEID=${NORMAL_SAMPLEID_ARR_[${i}]}  HOTSPOTFILE=${HOTSPOTFILE_ARR_[${i}]} INDELFILE=${INDELFILE_ARR_[${i}]} TISSUE=${TISSUE_ARR_[${i}]} KEEPFILE=${KEEPFILE_ARR_[${i}]} REGIONFILE=${REGIONFILE_ARR_[${i}]} AMPLIFICATIONFILE=${AMPLIFICATIONFILE_ARR_[${i}]} BACKGROUNDFILE=${BACKGROUNDFILE_ARR_[${i}]} TYPE=${TYPE_ARR_[${i}]} CUTADAPT_PREFIX=${CUTADAPT_PREFIX_ARR_[${i}]} METHOD=${METHOD_ARR_[${i}]};

    # Add info about samples to pipeline log
    echo -e "Sample="$SAMPLEID >> $PIPELINE_LOG;
    echo -e "\tTissue="$TISSUE >> $PIPELINE_LOG;
    echo -e "\tType="$TYPE >> $PIPELINE_LOG;
    echo -e "\tMethod="$METHOD >> $PIPELINE_LOG;
    echo -e "\tRefseq="$REFSEQ >> $PIPELINE_LOG;
    hotspots=($(echo $HOTSPOTFILE | tr "/" " "));
    echo -e "\tHotspots"=${hotspots[${#hotspots[@]}-1]} >> $PIPELINE_LOG;
    amp=($(echo $AMPLIFICATIONFILE | tr "/" " "));
    echo -e "\tAmplification_regions"=${amp[ ${#amp[@]}-1]} >> $PIPELINE_LOG;
    backg=($(echo $BACKGROUNDFILE | tr "/" " "));
    echo -e "\tBackground_regions"=${backg[ ${#backg[@]}-1]} >> $PIPELINE_LOG;
    
	# loop over all available steps
	for j in $(seq 0 6 $NUMBOFSTEPS); do

		# get current step in loop
		STEP=${steps[${j}]}; DEP=${steps[${j}+4]}; SCRIPT=${steps[${j}+3]}; PRE_POST=${steps[${j}+5]};

		# run step if on and not a post-step
		if [[ ${steps[${j}+2]} = "on" && $PRE_POST = 0 ]]; then
			runStep $STEP $SCRIPT $DEP $PRE_POST;

		fi

	done;
done;

# SUBMIT POST SCRIPTS FOR EACH SAMPLE
for i in $(seq 0 $NUMBEROFSAMPLES); do

	#export sample-specific variables
	export SAMPLEID=${SAMPLEID_ARR_[${i}]} BARCODE_I7=${BARCODE_I7_ARR_[${i}]} BARCODE_I5=${BARCODE_I5_ARR_[${i}]} RAWDATA_PE1=${RAWDATA_PE1_ARR_[${i}]} RAWDATA_PE2=${RAWDATA_PE2_ARR_[${i}]} RAWDATA_INDEX=${RAWDATA_INDEX_ARR_[${i}]} REFSEQ=${REFSEQ_ARR_[${i}]} ROIFILE=${ROIFILE_ARR_[${i}]} SELECTIONFILE=${SELECTIONFILE_ARR_[${i}]} NORMAL_SAMPLEID=${NORMAL_SAMPLEID_ARR_[${i}]}  HOTSPOTFILE=${HOTSPOTFILE_ARR_[${i}]} INDELFILE=${INDELFILE_ARR_[${i}]} TISSUE=${TISSUE_ARR_[${i}]} KEEPFILE=${KEEPFILE_ARR_[${i}]} REGIONFILE=${REGIONFILE_ARR_[${i}]} AMPLIFICATIONFILE=${AMPLIFICATIONFILE_ARR_[${i}]} BACKGROUNDFILE=${BACKGROUNDFILE_ARR_[${i}]} TYPE=${TYPE_ARR_[${i}]} CUTADAPT_PREFIX=${CUTADAPT_PREFIX_ARR_[${i}]} METHOD=${METHOD_ARR_[${i}]};

	for j in $(seq 0 6 $NUMBOFSTEPS); do
	
		if [[ ${steps[${j}+2]} = "on" && ${steps[${j}+5]} = 1 ]]; then
			# get current step in loop
			STEP=${steps[${j}]}; DEP=${steps[${j}+4]}; SCRIPT=${steps[${j}+3]}; PRE_POST=${steps[${j}+5]};
			runStep $STEP $SCRIPT $DEP  $PRE_POST;
		fi
	done;
done;

# SUBMIT POST SCRIPTS
for i in $(seq 0 6 $NUMBOFSTEPS); do

	if [[ ${steps[${i}+2]} = "on" && ${steps[${i}+5]} = 2 ]]; then
		# get current step in loop
		STEP=${steps[${i}]}; SCRIPT=${steps[${i}+3]}; DEP=${steps[${i}+4]};  PRE_POST=${steps[${i}+5]};
		runStep $STEP $SCRIPT $DEP $PRE_POST
		

	fi
done;
