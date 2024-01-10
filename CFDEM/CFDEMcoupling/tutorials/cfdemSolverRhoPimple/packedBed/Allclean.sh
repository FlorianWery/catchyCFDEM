#!/bin/bash

source $CFDEM_PROJECT_DIR/etc/functions.sh

casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"

cleanCFDEMcase $casePath
rm -rf CFD/0
rm -f CFD/constant/mechanism/*.gas CFD/constant/mechanism/*.surf
rm -f */*/*.png */*.png *.png
rm -rf __pycache__ *.pyc