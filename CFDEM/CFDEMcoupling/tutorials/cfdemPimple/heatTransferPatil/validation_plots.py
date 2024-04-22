import sys
import glob
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['font.family'] = "serif"

validation1=np.genfromtxt('patil_u1.2.csv', dtype=float, delimiter = ',')
validation2=np.genfromtxt('patil_u1.54.csv', dtype=float, delimiter = ',')
validation3=np.genfromtxt('patil_u1.71.csv', dtype=float, delimiter = ',')
cfdemcase1=np.genfromtxt('hw350_u1.20/CFD/evolutionTemperature.csv', dtype=float, delimiter = ',')
cfdemcase2=np.genfromtxt('hw350_u1.54/CFD/evolutionTemperature.csv', dtype=float, delimiter = ',')
cfdemcase3=np.genfromtxt('hw350_u1.71/CFD/evolutionTemperature.csv', dtype=float, delimiter = ',')

fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(8,4))
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Mean particle temperature [°C]')
ax1.set_title('75 g particles (1 mm, 2500 kg/m$^3$)', fontsize=10)
ax1.plot(validation1[:,0], validation1[:,1], 'sk', fillstyle='none', label=r'$u_{bg}$ = 1.2 m/s')
ax1.plot(cfdemcase1[:,0], cfdemcase1[:,1]-273, '--k', linewidth=2)
ax1.plot(validation2[:,0], validation2[:,1], 'sr', fillstyle='none', label=r'$u_{bg}$ = 1.54 m/s')
ax1.plot(cfdemcase2[:,0], cfdemcase2[:,1]-273, '--r', linewidth=2)
ax1.plot(validation3[:,0], validation3[:,1], 'sb', fillstyle='none', label=r'$u_{bg}$ = 1.71 m/s')
ax1.plot(cfdemcase3[:,0], cfdemcase3[:,1]-273, '--b', linewidth=2)
ax1.set_xlim([0,10])
ax1.set_ylim([55,90])
ax1.legend(loc='upper right', frameon=False)
plt.tight_layout()
plt.show()
