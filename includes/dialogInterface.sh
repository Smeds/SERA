#!/bin/bash

# generate list for dialog
SELECTARR=();
for i in $(seq 0 6 $NUMBOFSTEPS); do
	if [[ ! -z ${steps[${i}]} ]]; then
		SELECTARR+=( "${steps[${i}]}" "${steps[${i}+1]}" "${steps[${i}+2]}" );
	fi
done

# run dialog
dialog_steps=`dialog --stdout --separate-output --title "$title" --checklist "\n$INFO\n\nViewing current settings." 0 0 0 "${SELECTARR[@]}"`;

# change original array
for i in $(seq 0 6 $NUMBOFSTEPS); do
	# set all to off
	steps[$i+2]="off";
	
	# if step id in dialog output, set to on
	for j in ${dialog_steps[@]}; do
		if [[ "${steps[${i}]}" = "$j" ]]; then
			steps[$i+2]="on";
		fi
	done;
done
