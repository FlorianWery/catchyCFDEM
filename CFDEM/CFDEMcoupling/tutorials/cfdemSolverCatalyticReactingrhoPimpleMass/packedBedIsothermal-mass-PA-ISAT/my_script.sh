#!/bin/bash -l
#
#PBS -N test-mass-PAISAT
#PBS -l nodes=1:ppn=4
#PBS -l walltime=71:59:59
#PBS -m abe

source $HOME/setup_CFDEM-8_env-kirlia-foss-2020b.sh

cd $PBS_O_WORKDIR

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
mpirun -np $PBS_NP mycfdemSolverCatalyticReactingRhoPimpleMass -parallel > log.mycfdemSolverCatalyticReactingRhoPimpleMass 2>&1
reconstructPar -noLagrangian -latestTime > log.reconstructPar
#rm -r processor*

# post-processing
postProcess -func "components(U)" -latestTime > log.postProcess
#postProcess -func sample -latestTime >> log.postProcess
calcMixingCup -fields "(voidfraction T p Ux Uy Uz CH4 O2 C2H4 C2H6 CO CO2)" -nPoints 20 -latestTime >> log.postProcess
cd ..

# compare cantera
python3 validate.py
