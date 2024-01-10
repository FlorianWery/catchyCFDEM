# Installation instructions
All commands to be run from this location (catchyCFDEM/CFDEM)

## This can omitted if LIGGGHTS files are already present. Continue with normal installation instructions (README in catchyCFDEM/.).

## Requirements

sudo apt-get install build-essential cmake openmpi-bin libopenmpi-dev python-dev libvtk6-dev python-numpy

## Download LIGGGHTS from Github

git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git LIGGGHTS

## Installation CFDEMcoupling

source CFDEMcoupling/etc/bashrc
cfdemCompCFDEMall
