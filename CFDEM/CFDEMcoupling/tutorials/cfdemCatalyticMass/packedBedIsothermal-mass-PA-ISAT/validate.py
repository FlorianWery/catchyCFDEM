import cantera as ct
import numpy as np

import os
import sys
catchy_etc = os.environ.get('CATCHY_ETC') 
sys.path.append(os.path.join(catchy_etc, "validationTools"))

from reactormodels_thermo import surfacePFR_epsB

#######################################################################
# Input Parameters
#######################################################################

T0 = 1098.15 #1048.15  # Temperature 
P0 = 1.1e5 #1.01325e5   # Pressure
Y0 = {'CH4':0.66667, 'O2':0.33333} # Composition

#epsB = 0.45             # Voidage
epsC = 0.27              # Catalyst porosity
velocity = 5           # Superficial inlet velocity
cat_area_per_vol = 18.29e6 # Catalyst surface area per unit volume

ID = 0.005
L = 0.04

# input file containing the surface reaction mechanism
cti_file = 'CFD/constant/mechanism/ocm_polimi31_srlasic.cti'
gas = ct.Solution(cti_file, 'gas')
surf = ct.Interface(cti_file,'surface1', [gas])

# void fraction profile
latestTime=os.popen("foamListTimes -case CFD | tail -1").read().split('\n')[0]
ofdata = np.genfromtxt('CFD/postProcessing/mixingCup/{}/mixingCupResults.csv'.format(latestTime), dtype=float, delimiter = ',', names=True)
epsB_func = ct.Func1(lambda z: np.interp(z, ofdata['distance'], ofdata['voidfraction']))

#####################################################################

pfr = surfacePFR_epsB(gas, surf, 1.0, ID, epsB_func, cat_area_per_vol, epsC = epsC, isothermal=True)
pfr.gas.TPY = T0, P0, Y0
pfr.mdot = pfr.A*pfr.gas.density*velocity
y0 = np.hstack((pfr.gas.T, pfr.gas.P, pfr.gas.Y))
states, surfstates = pfr.solve(y0, L, 1000)
states.write_csv('pfr_gas.csv', cols=('extra','T','density','P','Y'))

#####################################################################

import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = ['Times']
matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'

matplotlib.rcParams.update({'font.size': 13})

latestTime=os.popen("foamListTimes -case CFD | tail -1").read().split('\n')[0]
pfrdata = np.genfromtxt('pfr_gas.csv', dtype=float, delimiter = ',', names=True)
ofdataG = np.genfromtxt('CFD/postProcessing/mixingCup/{}/mixingCupResults.csv'.format(latestTime), dtype=float, delimiter = ',', names=True)

fig, [[ax1,ax2,ax3],[ax4,ax5,ax6]] = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
ax1.set_xlabel('Axial position [m]')
ax1.set_ylabel('Voidage [-]')
ax1.plot(pfrdata['z'], pfrdata['epsB'], color='crimson', label='Cantera')
ax1.plot(ofdataG['distance'], ofdataG['voidfraction'], '--', linewidth=2, color='crimson', label='CFDEM')
handles, labels = ax1.get_legend_handles_labels()
ax2.set_xlabel('Axial position [m]')
ax2.set_ylabel('Velocity [m/s]')
ax2.plot(pfrdata['z'], pfrdata['vel'], color='crimson')
ax2.plot(ofdataG['distance'], ofdataG['Uz'], '--', linewidth=2, color='crimson', label='Uz')
ax2.plot(ofdataG['distance'], ofdataG['Ux'], '-.', linewidth=2, color='darkorange', label='Ux')
ax2.plot(ofdataG['distance'], ofdataG['Uy'], ':', linewidth=2, color='gold', label='Uy')
ax2.legend(loc=0, frameon=False)
ax3.set_xlabel('Axial position [m]')
ax3.set_ylabel('CH4 mass fraction [-]')
ax3.plot(pfrdata['z'], pfrdata['Y_CH4'], color='crimson', label='Cantera')
ax3.plot(ofdataG['distance'], ofdataG['CH4'], '--', linewidth=2, color='crimson', label='CFDEM')
ax4.set_xlabel('Axial position [m]')
ax4.set_ylabel('O2 mass fraction [-]')
ax4.plot(pfrdata['z'], pfrdata['Y_O2'], color='crimson', label='Cantera')
ax4.plot(ofdataG['distance'], ofdataG['O2'], '--', linewidth=2, color='crimson', label='CFDEM')
ax5.set_xlabel('Axial position [m]')
ax5.set_ylabel('Mass fraction [-]')
ax5.plot(pfrdata['z'], pfrdata['Y_C2H4'], color='limegreen', label='C2H4')
ax5.plot(ofdataG['distance'], ofdataG['C2H4'], '--', linewidth=2, color='limegreen')
ax5.plot(pfrdata['z'], pfrdata['Y_C2H6'], color='dodgerblue', label='C2H6')
ax5.plot(ofdataG['distance'], ofdataG['C2H6'], '--', linewidth=2, color='dodgerblue')
ax5.plot(pfrdata['z'], pfrdata['Y_CO'], color='darkorange', label='CO')
ax5.plot(ofdataG['distance'], ofdataG['CO'], '--', linewidth=2, color='darkorange')
ax5.plot(pfrdata['z'], pfrdata['Y_CO2'], color='crimson', label='CO2')
ax5.plot(ofdataG['distance'], ofdataG['CO2'], '--', linewidth=2, color='crimson')
ax5.legend(loc=0, frameon=False)
ax6.set_xlabel('Axial position [m]')
ax6.set_ylabel('Pressure [Pa]')
ax6.plot(pfrdata['z'], pfrdata['P'], color='crimson', label='Cantera')
ax6.plot(ofdataG['distance'], ofdataG['p'], '--', linewidth=2, color='crimson', label='CFDEM')
ax6.set_ylim([1.09e5,1.11e5])

plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False, fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')
