#!/bin/bash

# renders list for choosing steps
function renderList {
	clear;
	echo -e $title;
	echo -e "Viewing current settings.\n";

	echo -e "$INFO\n";

	# print the list of available steps
	for i in $(seq 0 6 $NUMBOFSTEPS); do
		echo -n "${steps[$i]} ";
		echo -n "${steps[${i}+1]} ";
		if [[ ${steps[${i}+2]} = "on" ]]; then echo "[ON]";
		elif [[ ${steps[${i}+2]} = "off" ]]; then echo "[OFF]";
		else echo "[UNAVAILABLE]";
		fi
	done;

	# prompt for input
	echo -e "\nChange status by typing step #;"
	echo -n "> ";
	read inp;
	
	# react according to input
	case $inp in
		
		# if start, dont render new list
		"start") ;;

		# change status for corresponding entry
		(*)
			for i in $(seq 0 6 $NUMBOFSTEPS); do
				case ${steps[$i]} in
					"$inp") 
						if [[ ${steps[$i+2]} == "on" ]]; then
							steps[$i+2]="off";
						elif [[ ${steps[$i+2]} == "off" ]]; then
							steps[$i+2]="on";
						fi
					;;
				esac;
			done

			# iterate list if input chosen
			renderList;
		;;	
	esac;
}

renderList;
