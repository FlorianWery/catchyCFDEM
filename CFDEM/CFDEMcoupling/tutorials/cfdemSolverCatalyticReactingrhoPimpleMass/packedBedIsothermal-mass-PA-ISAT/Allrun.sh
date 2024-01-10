#!/bin/bash

export PBS_NP=4
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !

# run CFDDEM
cd DEM
cp in.liggghts_runfirst in.liggghts_run

cd ../CFD
cp -r 0.orig 0
blockMesh > log.blockMesh

reconstructPar -noLagrangian -latestTime > log.reconstructPar
decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimpleMass -parallel > log.cfdemSolverCatalyticReactingRhoPimpleMass 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar
rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
#postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate.py
