import cantera as ct
import numpy as np
import scipy.integrate

import os
import sys


def get_net_lh_production_rates(gas):
    conc = gas.concentrations
    ich4, io2 = gas.species_index('CH4'), gas.species_index('O2')
    ico, ico2 = gas.species_index('CO'), gas.species_index('CO2')
    ic2h6, ic2h4 = gas.species_index('C2H6'), gas.species_index('C2H4')
    ih2o, ih2 = gas.species_index('H2O'), gas.species_index('H2')
    
    rtemp = gas.T*ct.gas_constant
    k0s = np.array([3.83e7,1.58e4,2.19e2,1.14e6,1.23e8])
    Eas = np.array([1.58e8,1.23e8,9e7,1.23e8,1.56e8])
    ks = k0s*np.exp(-Eas/rtemp)
    k_co2s = np.array([8.2,8.2,8.2,8.2,8.2])
    Ea_co2s = np.array([-5.6e7,-5.6e7,-5.6e7,-5.6e7,-5.6e7])
    KCO2s = k_co2s*np.exp(-Ea_co2s/rtemp)
    KO2 = 1.18e-1*np.exp(7.9e7/rtemp)
    
    r1 = ks[0]*conc[ich4]*conc[io2]**0.5/(1.+KO2*conc[io2]+KCO2s[0]*conc[ico2])**2
    r2 = ks[1]*conc[ich4]**0.5*conc[io2]**0.5/(1.+KCO2s[1]*conc[ico2])
    r3 = ks[2]*conc[ich4]**0.5*conc[io2]**0.5/(1.+KCO2s[2]*conc[ico2])
    r4 = ks[3]*conc[ico]/(1.+KCO2s[3]*conc[ico2])
    r5 = ks[4]*conc[ic2h6]*conc[io2]**0.5/(1.+KCO2s[4]*conc[ico2])
    
    sdot = np.zeros(gas.n_species)
    sdot[ich4] = - 2*r1 - r2 - r3
    sdot[io2] = - 0.5*r1 - 1.5*r2 - r3 - 0.5*r4 - 0.5*r5
    sdot[ico] = r3 - r4
    sdot[ico2] = r2 + r4
    sdot[ih2o] = r1 + r2 + r3 + r5
    sdot[ih2] = r2 + r3
    sdot[ic2h6] = r1 - r5
    sdot[ic2h4] = r5
    
    return sdot

class ocmSurfacePFR:
    def __init__(self, gas, mdot, diam, epsB, rhocat, epsC = 0, isothermal = False, pressuredrop = False):
        self.gas = gas
        self.mdot = mdot
        self.diam = diam
        self.A = np.pi*diam**2/4
        self.epsB = epsB
        self.epsC = epsC
        self.rhocat = rhocat
        self.isothermal = isothermal
        self.deltap = pressuredrop

    def __call__(self, z, y):
        """the ODE function, y' = f(z,y) """
        # State vector is [T, p, Y_1, Y_2, ..., Y_K]
        #self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TPY = y[0], y[1], y[2:]

        wdot = self.gas.net_production_rates
        sdot = get_net_lh_production_rates(self.gas)
        rdot = (self.epsB + (1.0-self.epsB)*self.epsC) * wdot + (1.0 - self.epsB) * self.rhocat * sdot
        dYdz = rdot * self.gas.molecular_weights / self.mdot * self.A
        
        if (self.isothermal):
            dTdz = 0.0
        else:
            dTdz = (-np.dot(self.gas.partial_molar_enthalpies, rdot)*self.A) / (self.mdot*self.gas.cp_mass)  
        
        if (self.deltap):
            rho = self.gas.density
            u = self.mdot/rho/self.A
            mu = self.gas.viscosity
            dpdz = - 150*mu*u*(1-self.epsB)**2/self.epsB**3 - 1.75*rho*u**2*(1-self.epsB)/self.epsB**3
        else:
            dpdz = 0.0

        return np.hstack((dTdz, dpdz, dYdz))
    
    def state(self, y):
        self.gas.set_unnormalized_mass_fractions(y[2:])
        self.gas.TP = y[0], y[1]
        return self.gas.state
        return self.gas.state

    def solve(self, y0, length, nsteps):
        dz = length/nsteps

        # Set up solver
        solver = scipy.integrate.ode(self)
        solver.set_integrator('vode', method='bdf', with_jacobian=True, nsteps=10000, atol=1e-15, rtol=1e-9, order=5, min_step = 1e-15)
        solver.set_initial_value(y0, 0.0)

        # Integrate the equations
        restime = 0.0
        u = self.mdot/self.gas.density/self.A
        states = ct.SolutionArray(self.gas, 1, extra={'z': [0.0], 'tau':[restime], 'vel':[u], 'MW':[self.gas.mean_molecular_weight], 'epsB':self.epsB})
        while solver.successful() and solver.t < length:
            solver.integrate(solver.t + dz)
            self.state(solver.y)
            u = self.mdot/self.gas.density/self.A
            restime=restime+dz/u
            states.append(self.state(solver.y), z=solver.t, tau=restime, vel=u, MW=self.gas.mean_molecular_weight, epsB=self.epsB)
        
        if states.z[-1] < length:
            print("PFR simulation did not converge up to specified length")

        return states

#######################################################################
# Input Parameters
#######################################################################

T0 = 973.0   # Temperature 
P0 = 1.1e5   # Pressure
Y0 = {'CH4':0.0977562, 'O2':0.048746, 'N2':0.853498} # Composition

epsB = 0.4             # Voidage
epsC = 0.0             # Catalyst porosity
velocity = 1.5         # Interstitial inlet velocity
cat_mass_per_vol = 1700. # Catalyst surface area per unit reactor(!) volume

ID = 0.005
L = 0.01

# input file containing the surface reaction mechanism
cti_file = 'constant/mechanismLH/ocm_lh.cti'
gas = ct.Solution(cti_file, 'gas')

# output files
output_filename_gas = 'pfr_gas.csv'

#####################################################################

pfr = ocmSurfacePFR(gas, 1.0, ID, epsB, cat_mass_per_vol/(1.0-epsB), epsC = epsC, isothermal=False)
pfr.gas.TPY = T0, P0, Y0
pfr.mdot = pfr.A*pfr.gas.density*velocity
y0 = np.hstack((pfr.gas.T, pfr.gas.P, pfr.gas.Y))
states = pfr.solve(y0, L, 1000)

states.write_csv('pfr_gas.csv', cols=('extra','T','density','P','Y'))
print("Gas phase results saved to '{0}'".format(output_filename_gas))

#####################################################################

import matplotlib
import matplotlib.pyplot as plt
import glob

matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['font.family'] = "serif"

latestTime = os.popen("foamListTimes | tail -1").read().split('\n')[0]
pfrdata = np.genfromtxt('pfr_gas.csv', dtype=float, delimiter = ',', names=True)
offiles = glob.glob('postProcessing/sample/{}/*.csv'.format(latestTime))
ofdata = np.genfromtxt(offiles[0], dtype=float, delimiter = ',', names=True)

fig, [[ax1,ax2],[ax3,ax4]] = plt.subplots(nrows=2, ncols=2, figsize=(8,8))
ax1.set_xlabel('Axial position [m]')
ax1.set_ylabel('CH4 mass fraction [-]')
ax1.plot(pfrdata['z'], pfrdata['Y_CH4'], color='crimson', label='Cantera')
ax1.plot(ofdata['y'], ofdata['CH4'], '--', color='crimson', label='catchyFOAM')
#ax1.plot(ofdata['y'], ofdata['CH4catalyst'], ':', color='crimson', label='catchyFOAM(cat)')
handles, labels = ax1.get_legend_handles_labels()
ax2.set_xlabel('Axial position [m]')
ax2.set_ylabel('O2 mass fraction [-]')
ax2.plot(pfrdata['z'], pfrdata['Y_O2'], color='crimson', label='Cantera')
ax2.plot(ofdata['y'], ofdata['O2'], '--', color='crimson', label='catchyFOAM')
ax3.set_xlabel('Axial position [m]')
ax3.set_ylabel('Temperature [K]')
ax3.plot(pfrdata['z'], pfrdata['T'], color='crimson', label='Cantera')
ax3.plot(ofdata['y'], ofdata['T'], '--', color='crimson', label='catchyFOAM')
#ax3.set_ylabel('Mass fraction N2 [-]')
#ax3.plot(pfrdata['z'], pfrdata['Y_N2'], color='crimson', label='Cantera')
#ax3.plot(ofdata['y'], ofdata['N2'], '--', color='crimson', label='catchyFOAM')
ax4.set_xlabel('Axial position [m]')
ax4.set_ylabel('Mass fraction [-]')
ax4.plot(pfrdata['z'], pfrdata['Y_C2H4'], color='navy', label='C2H4')
ax4.plot(ofdata['y'], ofdata['C2H4'], '--', color='navy')
ax4.plot(pfrdata['z'], pfrdata['Y_C2H6'], color='dodgerblue', label='C2H6')
ax4.plot(ofdata['y'], ofdata['C2H6'], '--', color='dodgerblue')
ax4.plot(pfrdata['z'], pfrdata['Y_CO'], color='darkorange', label='CO')
ax4.plot(ofdata['y'], ofdata['CO'], '--', color='darkorange')
ax4.plot(pfrdata['z'], pfrdata['Y_CO2'], color='crimson', label='CO2')
ax4.plot(ofdata['y'], ofdata['CO2'], '--', color='crimson')
ax4.legend(loc=0, frameon=False)
plt.figlegend(handles, labels, loc='upper center', ncol=3, frameon=False)
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('results_comparison.png')
#plt.show()
print("Plot saved to 'results_comparison.png'")

#####################################################################
