#!/bin/bash

export PBS_NP=4
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !



./Allclean.sh

# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > ../log.liggghts_init
cd ..

# run CFDDEM
cd DEM
cp in.liggghts_runfirst in.liggghts_run

cd ../CFD
sed -i 's/endTime         0.15;/endTime         0.10;/g' system/controlDict
cp -r 0.orig 0
blockMesh > log.blockMesh

decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimpleMass -parallel > log.cfdemSolverCatalyticReactingRhoPimpleMass 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar
#rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
#postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate.py

# run CFDDEM
cd DEM
cp in.liggghts_runrestart in.liggghts_run

cd ../CFD

sed -i 's/endTime         0.10;/endTime         0.15;/g' system/controlDict
decomposePar -latestTime -force > log.decomposePar

mpirun -np $PBS_NP cfdemSolverCatalyticReactingRhoPimpleMass -parallel > log.cfdemSolverCatalyticReactingRhoPimpleMassrestart 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar
#rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
#postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate_end.py
