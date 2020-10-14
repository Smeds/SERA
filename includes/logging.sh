#!/bin/bash
#
# Include file with logging functions.
#
#. config/globals.sh # Include globals.

function SuccessLog {
	
	# Print a success post to log file
#	entry="\t<entry>\n\t\t<date>`date +%y-%m-%d`</date>\n\t\t<time>`date +%H:%M:%S`</time>\n\t\t<status>OK</status>\n\t\t<host>${HOSTNAME}</host>\n\t\t<sample>${1}</sample>\n\t\t<text>${2}</text>\n\t</entry>\n</document>";
#	eval `sed -i "s:</document>:$entry:" $LOG_FILE`;
#	sed -i "s|</document>|$entry|g" $LOG_FILE;
#	sed -i "s|</document>|$entry|" $LOG_FILE;
	
	if [[ "$SOFTWARE" != "SLURM" ]]; then
		echo "[${1}] SUCCESS: $2";
	fi

}  

function ErrorLog {

    # Print a error post to log file
#	entry="\t<entry>\n\t\t<date>`date +%y-%m-%d`</date>\n\t\t<time>`date +%H:%M:%S`</time>\n\t\t<status>ERROR</status>\n\t\t<host>${HOSTNAME}</host>\n\t\t<sample>${1}</sample>\n\t\t<text>${2}</text>\n\t</entry>\n</document>";
#	eval `sed -i "s:</document>:$entry:" $LOG_FILE`;
#	sed -i "s|</document>|$entry|" $LOG_FILE;
#	sed -i "s|</document>|$entry|g" $LOG_FILE;

	if [[ "$SOFTWARE" != "SLURM" ]]; then
		echo "[${1}] ERROR: $2";
	fi
	exit 1;
}

function WarningLog {

    # Print a warning post to log file
#	entry="\t<entry>\n\t\t<date>`date +%y-%m-%d`</date>\n\t\t<time>`date +%H:%M:%S`</time>\n\t\t<status>WARNING</status>\n\t\t<host>${HOSTNAME}</host>\n\t\t<sample>${1}</sample>\n\t\t<text>${2}</text>\n\t</entry>\n</document>";
#	eval `sed -i "s:</document>:$entry:" $LOG_FILE`;
#	sed -i "s|</document>|$entry|g" $LOG_FILE;
#	sed -i "s|</document>|$entry|" $LOG_FILE;

	if [[ "$SOFTWARE" != "SLURM" ]]; then
		echo "[${1}] WARNING: $2";
	fi
}
