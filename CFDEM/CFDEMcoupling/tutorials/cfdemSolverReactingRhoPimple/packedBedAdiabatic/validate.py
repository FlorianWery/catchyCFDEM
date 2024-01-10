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

T0 = 1123.0  # Temperature 
P0 = 20.e5   # Pressure
Y0 = {'CH4':0.66666667, 'O2':0.333333333} # Composition

#epsB = 1.0             # Voidage
epsC = 0.0              # Catalyst porosity
velocity = 4.0          # Superficial inlet velocity
cat_area_per_vol = 0.0  # Catalyst surface area per unit volume

ID = 0.005
area = ID**2*np.pi/4.
L = 0.02

# input file containing the surface reaction mechanism
cti_file = 'CFD/constant/mechanism/ocm_polimi31_srla2o3.cti'
gas = ct.Solution(cti_file, 'gas')
surf = ct.Interface(cti_file,'surface1', [gas])

# bed voidage
latestTime=os.popen("foamListTimes -case CFD | tail -1").read().split('\n')[0]
ofdata = np.genfromtxt('CFD/postProcessing/mixingCup/{}/mixingCupResults.csv'.format(latestTime), dtype=float, delimiter = ',', names=True)
epsB_func = ct.Func1(lambda z: np.interp(z, ofdata['distance'], ofdata['voidfraction']))

#######################################################################
#
#   PFR
#
#######################################################################

#pfr = surfacePFR(gas, surf, 1.0, ID, epsB, cat_area_per_vol, epsC = epsC, isothermal=False)
pfr = surfacePFR_epsB(gas, surf, 1.0, ID, epsB_func, cat_area_per_vol, epsC = epsC, isothermal=False)
pfr.gas.TPY = T0, P0, Y0
pfr.mdot = pfr.A*pfr.gas.density*velocity
y0 = np.hstack((pfr.gas.T, pfr.gas.P, pfr.gas.Y))
states, surfstates = pfr.solve(y0, L, 100)
states.write_csv('pfr_gas.csv', cols=('extra','T','density','P','Y'))

#######################################################################
#
#   CASCADE OF CSTR
#
#######################################################################

# resolution
NReactors = 25
rlen = L/(NReactors-1)

# import the gas model and set the initial conditions
gas.TPY = T0, P0, Y0
surf.TP = T0, P0

epsB = epsB_func(0)

cat_area = cat_area_per_vol * area * rlen * (1.-epsB)
mdot = velocity * gas.density * area

TDY = gas.TDY
cov = surf.coverages

gas.TDY = TDY
r = ct.IdealGasReactor(gas, energy='on')
r.volume = area * rlen * (epsB + (1.-epsB)*epsC)

upstream = ct.Reservoir(gas, name='upstream')
downstream = ct.Reservoir(gas, name='downstream')

rsurf = ct.ReactorSurface(surf, r, A=cat_area)

m = ct.MassFlowController(upstream, r, mdot=mdot)
v = ct.PressureController(r, downstream, master=m, K=1e-5)

sim = ct.ReactorNet([r])
sim.max_err_test_fails = 12

sim.rtol = 1.0e-9
sim.atol = 1.0e-21

states = ct.SolutionArray(r.thermo, 1, extra={'z': [0.0], 'epsB':epsB, 'vel':mdot/gas.density/area})
for n in range(NReactors):
    dist = n * rlen
    epsB = epsB_func(dist)
    cat_area = cat_area_per_vol * area * rlen * (1.-epsB)
    r.volume = area * rlen * (epsB + (1.-epsB)*epsC)
    rsurf.area = cat_area

    gas.TDY = r.thermo.TDY
    upstream.syncState()
    sim.reinitialize()
    sim.advance_to_steady_state()

    states.append(r.thermo.state, z=dist, epsB=epsB, vel=mdot/r.thermo.density/area/epsB)

states.write_csv('pfr_gas_cascade.csv', cols=('extra','T','density','P','Y'))

#######################################################################
#
#   PLOT COMPARISON
#
#######################################################################

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
pfrdataC = np.genfromtxt('pfr_gas_cascade.csv', dtype=float, delimiter = ',', names=True)
ofdataG = np.genfromtxt('CFD/postProcessing/mixingCup/{}/mixingCupResults.csv'.format(latestTime), dtype=float, delimiter = ',', names=True)

fig, [[ax1,ax2,ax3],[ax4,ax5,ax6]] = plt.subplots(nrows=2, ncols=3, figsize=(12,8))
ax1.set_xlabel('Axial position [m]')
ax1.set_ylabel('Voidage [-]')
ax1.plot(pfrdata['z'], pfrdata['epsB'], color='crimson', label='Cantera PFR')
ax1.plot(pfrdataC['z'], pfrdataC['epsB'], '--', linewidth=2, color='crimson', label='Cantera 25 CSTR')
ax1.plot(ofdataG['distance'], ofdataG['voidfraction'], ':', linewidth=2, color='crimson', label='CFDEM')
handles, labels = ax1.get_legend_handles_labels()
ax2.set_xlabel('Axial position [m]')
ax2.set_ylabel('Velocity [m/s]')
ax2.plot(pfrdata['z'], pfrdata['vel'], color='crimson')
ax2.plot(pfrdataC['z'], pfrdataC['vel'], '--', linewidth=2, color='crimson')
ax2.plot(ofdataG['distance'], ofdataG['Uz'], ':', linewidth=2, color='crimson', label='Uz')
ax2.plot(ofdataG['distance'], ofdataG['Ux'], '-.', linewidth=2, color='darkorange', label='Ux')
ax2.plot(ofdataG['distance'], ofdataG['Uy'], '--', linewidth=2, color='gold', label='Uy')
ax2.legend(loc=0, frameon=False)
ax3.set_xlabel('Axial position [m]')
ax3.set_ylabel('CH4 mass fraction [-]')
ax3.plot(pfrdata['z'], pfrdata['Y_CH4'], color='crimson', label='Cantera')
ax3.plot(pfrdataC['z'], pfrdataC['Y_CH4'], '--', linewidth=2, color='crimson')
ax3.plot(ofdataG['distance'], ofdataG['CH4'], ':', linewidth=2, color='crimson', label='CFDEM')
ax4.set_xlabel('Axial position [m]')
ax4.set_ylabel('O2 mass fraction [-]')
ax4.plot(pfrdata['z'], pfrdata['Y_O2'], color='crimson', label='Cantera')
ax4.plot(pfrdataC['z'], pfrdataC['Y_O2'], '--', linewidth=2, color='crimson')
ax4.plot(ofdataG['distance'], ofdataG['O2'], ':', linewidth=2, color='crimson', label='CFDEM')
ax5.set_xlabel('Axial position [m]')
ax5.set_ylabel('Mass fraction [-]')
ax5.plot(pfrdata['z'], pfrdata['Y_C2H4'], color='limegreen', label='C2H4')
ax5.plot(pfrdataC['z'], pfrdataC['Y_C2H4'], '--', linewidth=2, color='limegreen')
ax5.plot(ofdataG['distance'], ofdataG['C2H4'], ':', linewidth=2, color='limegreen')
ax5.plot(pfrdata['z'], pfrdata['Y_C2H6'], color='dodgerblue', label='C2H6')
ax5.plot(pfrdataC['z'], pfrdataC['Y_C2H6'], '--', linewidth=2, color='dodgerblue')
ax5.plot(ofdataG['distance'], ofdataG['C2H6'], ':', linewidth=2, color='dodgerblue')
ax5.plot(pfrdata['z'], pfrdata['Y_CO'], color='darkorange', label='CO')
ax5.plot(pfrdataC['z'], pfrdataC['Y_CO'], '--', linewidth=2, color='darkorange')
ax5.plot(ofdataG['distance'], ofdataG['CO'], ':', linewidth=2, color='darkorange')
ax5.plot(pfrdata['z'], pfrdata['Y_CO2'], color='crimson', label='CO2')
ax5.plot(pfrdataC['z'], pfrdataC['Y_CO2'], '--', linewidth=2, color='crimson')
ax5.plot(ofdataG['distance'], ofdataG['CO2'], ':', linewidth=2, color='crimson')
ax5.legend(loc=0, frameon=False)
ax6.set_xlabel('Axial position [m]')
ax6.set_ylabel('Temperature [K]')
ax6.plot(pfrdata['z'], pfrdata['T'], color='crimson', label='Cantera')
ax6.plot(pfrdataC['z'], pfrdataC['T'], '--', linewidth=2, color='crimson')
ax6.plot(ofdataG['distance'], ofdataG['T'], ':', linewidth=2, color='crimson', label='CFDEM')

plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False, fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')
