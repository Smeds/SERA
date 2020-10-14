#!/bin/bash

. $SERA_PATH/includes/errorHandling.sh;

# check if valid input
if [[ -z $ROOT_PATH_INPUT ]] || [[ -z $INPUT_FILE ]]; then
	usage; exit 1;
elif [[ ! -d $ROOT_PATH_INPUT ]]; then
	fileioerror $ROOT_PATH_INPUT; kill -SIGINT $$;
elif [[ ! -e $INPUT_FILE ]]; then
	fileioerror $INPUT_FILE; kill -SIGINT $$;
elif [[ ( ! -z $INTERFACE ) && ( "$INTERFACE" != "gui" ) ]]; then
	usage; exit 1;
# if gui set, check that dialog is found in bin folders set in system path
elif [[ $INTERFACE && `echo $PATH | awk 'BEGIN{FS=":"}{for (i=1; i<=NF; i++) print $i}' | xargs ls | grep -w -P '^dialog$' | wc -l` = 0 ]]; then
	settingserror "Dialog not found in system, required for graphical user interface."; kill -SIGINT $$;
fi


