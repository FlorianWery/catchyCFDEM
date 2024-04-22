#!/bin/bash

export PBS_NP=4
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !


# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > ../log.liggghts_init
cd ..

# run CFDDEM
cd DEM
cp in.liggghts_runfirst in.liggghts_run

cd ../CFD
cp -r 0.orig 0
blockMesh > log.blockMesh

decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemCatalytic -parallel > log.cfdemCatalytic 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate.py
