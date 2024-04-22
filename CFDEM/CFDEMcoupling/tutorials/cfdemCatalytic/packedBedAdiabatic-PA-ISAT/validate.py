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

separate_figs = True

T0 = 998.15 # Temperature 
P0 = 1.1e5 #1.01325e5   # Pressure
X0 = {'CH4':4., 'O2':1., 'N2':20.}  # Composition

#epsB = 0.45             # Voidage
epsC = 0.27              # Catalyst porosity
velocity = 7.5           # Superficial inlet velocity
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

pfr = surfacePFR_epsB(gas, surf, 1.0, ID, epsB_func, cat_area_per_vol, epsC = epsC, isothermal=False)
pfr.gas.TPX = T0, P0, X0
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
ax6.set_ylabel('Temperature [K]')
ax6.plot(pfrdata['z'], pfrdata['T'], color='crimson', label='Cantera')
ax6.plot(ofdataG['distance'], ofdataG['T'], '--', linewidth=2, color='crimson', label='CFDEM')

plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False, fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')

if separate_figs:
	fig, ax = plt.subplots(figsize=(5,5))
	ax.set_xlabel('Axial position [m]')
	ax.set_xlim([0,0.03])
	ax.set_ylim([970,1130])
	ax.set_ylabel('Temperature [K]')
	ax.plot(pfrdata['z'], pfrdata['T'], color='blue')
	ax.plot(ofdataG['distance'], ofdataG['T'], '--', linewidth=2, color='blue')
	plt.xticks([0,0.01,0.02,0.03])
	plt.yticks([975,1000,1025,1050,1075,1100,1125])
	plt.tight_layout()
	plt.savefig('results_T_adi.png')
	plt.close()

	fig, ax = plt.subplots(figsize=(5,5))
	ax.set_xlabel('Axial position [m]')
	ax.set_xlim([0,0.03])
	ax.set_ylim([0.082,0.102])
	ax.set_ylabel('$CH_{4}$ mass fraction [-]')
	ax.plot(pfrdata['z'], pfrdata['Y_CH4'], color='blue')
	ax.plot(ofdataG['distance'], ofdataG['CH4'], '--', linewidth=2, color='blue')
	plt.xticks([0,0.01,0.02,0.03])
	plt.yticks([0.085,0.09,0.095,0.1])
	plt.tight_layout()
	plt.savefig('results_CH4_adi.png')
	plt.close()

	fig, ax = plt.subplots(figsize=(5,5))
	ax.set_xlabel('Axial position [m]')
	ax.set_xlim([0,0.03])
	ax.set_ylim([0.028,0.052])
	ax.set_ylabel('$O_{2}$ mass fraction [-]')
	ax.plot(pfrdata['z'], pfrdata['Y_O2'], color='blue')
	ax.plot(ofdataG['distance'], ofdataG['O2'], '--', linewidth=2, color='blue')
	plt.xticks([0,0.01,0.02,0.03])
	plt.yticks([0.03,0.035,0.04,0.045,0.05])
	plt.tight_layout()
	plt.savefig('results_O2_adi.png')
	plt.close()

	fig, ax = plt.subplots(figsize=(5,5))
	ax.set_xlabel('Axial position [m]')
	ax.set_xlim([0,0.03])
	ax.set_ylabel('Mass fraction [-]')
	ax.set_ylim(bottom = -0.005, top = 0.008)
	ax.plot(pfrdata['z'], pfrdata['Y_CO'], color='green', label='$C$$O$')
	ax.plot(ofdataG['distance'], ofdataG['CO'], '--', linewidth=2, color='green')
	ax.plot(pfrdata['z'], pfrdata['Y_CO2'], color='blue', label='$C$$O_2$')
	ax.plot(ofdataG['distance'], ofdataG['CO2'], '--', linewidth=2, color='blue')
	ax.legend(loc=0, frameon=False)
	plt.xticks([0,0.01,0.02,0.03])
	plt.tight_layout()
	plt.savefig('results_COx_adi.png')
	plt.close()

	fig, ax = plt.subplots(figsize=(5,5))
	ax.set_xlabel('Axial position [m]')
	ax.set_xlim([0,0.03])
	ax.set_ylabel('Mass fraction [-]')
	ax.set_ylim(bottom = -0.0002, top = 0.0052)
	ax.plot(pfrdata['z'], pfrdata['Y_C2H4'], color='green', label='$C_{2}$$H_{4}$')
	ax.plot(ofdataG['distance'], ofdataG['C2H4'], '--', linewidth=2, color='green')
	ax.plot(pfrdata['z'], pfrdata['Y_C2H6'], color='blue', label='$C_{2}$$H_{6}$')
	ax.plot(ofdataG['distance'], ofdataG['C2H6'], '--', linewidth=2, color='blue')
	ax.legend(loc=0, frameon=False)
	plt.xticks([0,0.01,0.02,0.03])
	plt.tight_layout()
	plt.savefig('results_C2_adi.png')
	plt.close()

	fig, ax = plt.subplots(figsize=(10,10))
	ax.plot(pfrdata['z'], pfrdata['Y_CO2'], color='blue', label='Cantera')
	ax.plot(ofdataG['distance'], ofdataG['CO2'], '--', linewidth=2, color='blue', label='catchyCFDEM')
	ax.legend(loc=0, frameon=False, ncol = 2, bbox_to_anchor=(1,-0.05))
	plt.savefig('Legend_adi.png')
