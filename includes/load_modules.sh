#!/bin/bash -i

# Include home, proj or moriarty specific globals
if [[ $GLOBALS = "HOME" ]]; then
    . $SERA_PATH/config/globalsHome.sh;

elif [[ $GLOBALS = "MORIARTY" ]]; then
    . $SERA_PATH/config/globalsMoriarty.sh;
    
elif [[ $GLOBALS = "MARVIN" ]]; then
    . $SERA_PATH/config/globalsMarvin.sh;
else
    . $SERA_PATH/config/globalsProj.sh;
fi
