#!/bin/bash

source $CFDEM_PROJECT_DIR/etc/functions.sh

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cleanCFDEMcase $casePath
rm -f */*/*.png */*.png
