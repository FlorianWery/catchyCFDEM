#!/bin/bash -l

source $WM_PROJECT_DIR/bin/tools/CleanFunctions
source $WM_PROJECT_DIR/bin/tools/RunFunctions

PBS_NP=6

# clean
rm -r DEM/post/*.vtk CFD/evolutionTemperature.csv

# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init
cd ..

# run CFDDEM
cd CFD
cleanCase
blockMesh > log.blockMesh
topoSet > log.topoSet
createPatch -overwrite > log.createPatch
decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemPimple -parallel > log.cfdemPimple 2>&1
