#/bin/bash

# remove steps in step_array based on user input
function removeSteps {
	for i in $(seq 0 6 $NUMBOFSTEPS); do
		if [[ ${steps[${i}]} = ${1} ]]; then

			# set to unaviable if command line ui
			if [[ -z "$INTERFACE" ]]; then
				steps[$i+2]="na";

			# remove form list if dialog gui
			else
				steps=( "${steps[@]:0:$i}" "${steps[@]:$i+6}" );
			fi
		fi
	done
}

#if [[ $INPUT_DESIGN != "MDA" ]]; then
#	removeSteps 50; removeSteps 51; removeSteps 52;
#fi
if [[ $SNPMANIAFLAGS = "false" || $CALL_TYPE = "false" ]]; then
	removeSteps 21; removeSteps 22; removeSteps 23; removeSteps 26;
fi

# if mosaikaligner region flag is missing, adjust many steps
#if [[ -z $ALIGNERFLAGS_REGION || $ALIGNERFLAGS_REGION = "false" ]]; then
#	removeSteps 1; removeSteps 3; removeSteps 6;
#fi

#if [[ ! $CNV_FILE ]]; then
#	removeSteps 61;
#fi
