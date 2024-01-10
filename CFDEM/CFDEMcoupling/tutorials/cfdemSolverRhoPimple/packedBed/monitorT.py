import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

matplotlib.rcParams['font.family'] = "serif"
matplotlib.rcParams['font.sans-serif'] = ['Times']
matplotlib.rcParams['mathtext.fontset'] = 'dejavuserif'

matplotlib.rcParams.update({'font.size': 13})

offileV = 'CFD/postProcessing/volAverage/0/volFieldValue.dat'
ofdataV = np.genfromtxt(offileV, dtype=float, delimiter = '\t', names=['time','T','Ts'])
offileO = 'CFD/postProcessing/outletAverage/0/surfaceFieldValue.dat'
ofdataO = np.genfromtxt(offileO, dtype=float, delimiter = '\t', names=['time','T'])

fig, ax = plt.subplots(figsize=(6,4))
ax.set_xlabel('Time [s]')
ax.set_ylabel('Average temperature [K]')
ax.plot(ofdataO['time'][1:], ofdataO['T'][1:], linewidth=2, color='crimson', label='gas outlet')
ax.plot(ofdataV['time'][1:], ofdataV['T'][1:], linewidth=2, color='darkorange', label='gas bed')
ax.plot(ofdataV['time'][1:], ofdataV['Ts'][1:], linewidth=2, color='dodgerblue', label='particles bed')
#ax.set_ylim([300,500])
ax.legend(loc=0, frameon=False)
plt.tight_layout()
plt.savefig('T_evolution.png')
plt.show()
