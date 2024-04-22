#!/bin/bash

export PBS_NP=4 
#! manually check decomposeParDict and in.liggghts_run to be consistent with number of processors !

# run liggghts init
cd DEM
mpirun -np $PBS_NP $CFDEM_LIGGGHTS_BIN_DIR/liggghts -in in.liggghts_init > ../log.liggghts_init
cd ..

# run CFDDEM
cd CFD
blockMesh > log.blockMesh
decomposePar -force > log.decomposePar
mpirun -np $PBS_NP cfdemPimple -parallel > log.cfdemPimple 2>&1
reconstructPar -latestTime -noLagrangian > log.reconstructPar
rm -r processor*

# post-process
cd octave 
octave --no-gui totalPressureDrop.m
eog cfdemSolverPiso_ErgunTestMPI.png
cd ../..
