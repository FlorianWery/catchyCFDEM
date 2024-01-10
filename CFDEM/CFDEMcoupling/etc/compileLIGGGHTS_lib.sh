#!/bin/bash

#===================================================================#
# compile routine for LIGGGHTS libraries, part of CFDEMproject 
# Christoph Goniva - March. 2014, DCS Computing GmbH
#===================================================================#
 
#- include functions
source $CFDEM_PROJECT_DIR/etc/functions.sh

NOW="$(date +"%Y-%m-%d-%H:%M")"
logDir="log"


cd $CFDEM_PROJECT_DIR/etc
mkdir -p $logDir

#================================================================================#
# compile src
#================================================================================#
    whitelist="$CFDEM_PROJECT_DIR/etc/library-liggghts-list.txt"
    echo ""
    echo "Compiling sub-libraries of LIGGGHTS now..."
    echo "Please provide the libraries to be compiled in the $CWD/$whitelist file."
    echo "Libraries must be in the $CFDEM_LIGGGHTS_SRC_DIR/../lib directory."

    if [ ! -f "$CWD/$whitelist" ];then
        echo "$whitelist does not exist in $CWD. Nothing will be done."
        NLINES=0
        COUNT=0
    else
        NLINES=`wc -l < $CWD/$whitelist`
        COUNT=0
    fi

    logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"

    while [ $COUNT -lt $NLINES ]
    do
            let COUNT++  
            LINE=`head -n $COUNT $CWD/$whitelist | tail -1`
  
            # white lines
            if [[ "$LINE" == "" ]]; then
                echo "compile $LINE"
                continue
            # comments
            elif [[ "$LINE" == \#* ]]; then
                continue
             # paths
            elif [[ "$LINE" == */dir ]]; then
                echo "will change path..."
                LINE=$(echo "${LINE%????}")
                path="$CFDEM_LIGGGHTS_SRC_DIR/../lib"
                cd $path
                echo $PWD
                #continue
            fi

            #--------------------------------------------------------------------------------#
            #- define variables
            #logpath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")/$logDir"
            logfileName="log_compile$LINE""lib"
            headerText="$logfileName""-$NOW"
            libVarMakefileName="CFDEM_$LINE""LIB_MAKEFILENAME"
            makeFileName="Makefile.${!libVarMakefileName}"
            libVarName="CFDEM_$LINE""LIB_PATH"
            libraryPath="${!libVarName}"
            #--------------------------------------------------------------------------------#

            compileLMPlib $logpath $logfileName $headerText $makeFileName $libraryPath
    done
