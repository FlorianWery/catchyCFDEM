# README for catchyCFDEM
Original: Florian Wéry, January 2024

## About catchyCFDEM
catchyCFDEM (CATCHY = CATalytic CHemistrY) is an extension of the catchyFOAM framework.<br/><br/>

## Installation prerequisites
### OpenFOAM
catchyCFDEM is compatible with OpenFOAM-8 (OpenFOAM Foundation). For download instructions, see [https://openfoam.org/version/8].

Do not forget to install all prerequisites before starting the compilation of catchyCFDEM ([https://www.cfdem.com/media/CFDEM/docu/CFDEMcoupling_Manual.html#install-4]).

### (optional) Cantera
catchyCFDEM comes with a utility to convert Cantera mechanism files into OpenFOAM format, <span style="font-family:Courier;">canteraToFoam</span> (as an alternative to the <span style="font-family:Courier;">chemkinToFoam</span> utility available by default). To install <span style="font-family:Courier;">canteraToFoam</span>, Cantera needs to be installed.

Get the Cantera source code from Github as reported in catchyFOAM and switch to the correct branch via
```
    git clone https://github.com/lavdwall/cantera.git --branch gas_transport_access
```

Detailed installation instructions can be found [here](https://cantera.org/install/compiling-install.html). 
For Ubuntu, the following commands should do the trick
```
    sudo apt-get install g++ gfortran scons libboost-dev python3 python3-dev python3-setuptools python3-numpy python3-ruamel.yaml python3-pip
    pip3 install Cython
    
    cd cantera
    scons build prefix='$HOME/.local' python_cmd=/usr/bin/python3
    scons install
    cd ..
```

## Install catchyCFDEM
Get the catchyCFDEM source code from Github 
```
    git clone https://github.com/FlorianWery/catchyCFDEM.git
```
By default, it is assumed that the catchyCFDEM repository is located in *$HOME/OpenFOAM*. If this is not the case, adjust *catchyCFDEM/etc/bashrc* accordingly. Also, if Cantera is installed in a non-default location, it can be set in *catchyCFDEM/etc/bashrc*.

Navigate to the catchyCFDEM repository and run
```
    ./Allwmake
```
This will install all libraries, solvers and utilities of catchyFOAM.

To use catchyCFDEM, the catchyCFDEM environment should be set by sourcing the *catchyCFDEM/etc/bashrc* and *catchyCFDEM/CFDEM/CFDEMcoupling/etc/bashrc* files. This can be done by adding the following line at the end of the user's *.bashrc* file (to be adjusted accordingly if catchyCFDEM is not located in *$HOME/OpenFOAM*)
```
    source $HOME/OpenFOAM/catchyCFDEM/etc/bashrc
    source $HOME/OpenFOAM/catchyCFDEM/CFDEM/CFDEMcoupling/etc/bashrc
```
To install the CFD-DEM part of the code, run the following command:
```
    cfdemCompCFDEMall
```
## Uninstall catchyCFDEM
To clean the installation of catchyFOAM, run
```
    ./Allwclean
```
