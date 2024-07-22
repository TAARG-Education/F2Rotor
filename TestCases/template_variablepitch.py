
from ... import bemt
from ... import aero
from ... import cad_v2 as cad
import os
import numpy as np
import math
import matplotlib.pyplot as plt

## USER INTERFACE ##
# operating condition
z = 1
rpm = 2000 
J_vec = np.arange(0., 2.6, 0.02)

# Geometry
D = 1.66         # blade diameter
nbl = 2          # blade number
Rhub = 0.1       # hub final section 
dx = 0.01        # radius increment along blade
c = 0.14         # chord distribution
beta = 30.       # beta distribution
pitch = 0.       # pitch
r_R_known = [0.12048,0.25, 0.5, 0.75, 0.9999999999999999]       # r/R station
c_known = [0.083,0.1079,0.1411,0.0913,0.0166]                   # chord distribution
beta_known = [70,55,35,25,0.1]                                  # beta distribuction
Airfoil = ['4412', '4412', '4412', '4412','4412']
kind = 'cubic'

# initialize pitch vector
pitchd_vec = [0., 5., 10., 15., 20]

### CAD ###
geom = cad.Geometry(D/2, Rhub, nbl, rpm, r_R_known, c_known, beta_known,0,Airfoil,pitch, kind)

### AERODYNAMICS ###
### Instantiate the aerodynamics class ###
aero_method = 2
aero_params2 = {                            # For aero method 2
    'alpha_0_lift': np.deg2rad(-4.3865),      # In deg inside the bracket
    'Cl_alpha':  6.4582,                   # 1/rad
    'Cl_alpha_stall': -0.034,
    'Cl_max':   1.5209,
    'Cl_min':    -0.932,
    'Cl_incr_to_stall':0.284,
    'Cd_min': 0.0059,
    'Cl_at_cd_min': 0.5296,
    'dCd_dCl2':   0.00804,
    'Mach_crit': 0.8,
    'Re_scaling_exp': -0.15
}
aero = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params2, M_corr=True, Rey_corr=True, z = z)


# initialize matrices to store coefficients
ct_matrix = np.full((len(J_vec),len(pitchd_vec)),np.nan)
cp_matrix = np.full((len(J_vec),len(pitchd_vec)),np.nan)
eta_matrix = np.full((len(J_vec),len(pitchd_vec)),np.nan)
# BEMT ##
for m,pitchd in enumerate(pitchd_vec):
    geom.pitch = np.array(pitchd)
    for i,J in enumerate(J_vec):
        ct_tvor, cp_tvor = bemt.BEMT_tvor(z,J,dx,geom,aero)
     
        ct_matrix[i,m] = ct_tvor
        cp_matrix[i,m] = cp_tvor
        eta_matrix[i,m] = J*ct_tvor/cp_tvor
        if (ct_matrix[:,m] < 0).any(): break
ind = np.where(ct_matrix < 0)[0]
indx = np.where(ct_matrix < 0)[0]


# OUTPUT ##
var = [ct_matrix, cp_matrix, eta_matrix]
label = [r'$\bf C_T$', r'$\bf C_P$', r'$\bf \eta$']
fig_lab = ['ct_var-pitch', 'cp_var-pitch', 'eta_var-pitch']
ls = ['--','-','-.']
values = np.linspace(0, 1, len(pitchd_vec))
cmap = plt.get_cmap('jet')
colors = [cmap(value) for value in values]

# for loop on variables
for i in range(3):
    fig = plt.figure()
    max_y = 1.1*np.max(var[i][~np.isnan(var[i])]) if i!=2 else 1.
    for m, pitch in enumerate(pitchd_vec):
        max = np.max(var[i][~np.isnan(var[i][:,m]),m])
        plt.plot(J_vec, var[i][:,m], 'k', label=r'$\phi={:.0f}°$'.format(pitch))
        plt.xlabel(r'$\bf J$', fontsize=12)  
        plt.ylabel(label[i], fontsize=12)
        plt.grid(True, which='both')
        plt.minorticks_on()
        if i==0: plt.text(J_vec[np.where(np.round(var[i][:,m],2) == 0.06)[0][0]], 0.06, r'$\phi={:.0f}°$'.format(pitch), fontsize=8, rotation=-60)
        if i==1: plt.text(J_vec[np.where(var[i][:,m] == max)[0]], 1.05*max, r'$\phi={:.0f}°$'.format(pitch), fontsize=8)
        if i==2: plt.text(1.02*J_vec[np.where(np.round(var[i][:,m],1) == 0.7)[0][-1]], 0.4, r'$\phi={:.0f}°$'.format(pitch), fontsize=8, rotation=-90, backgroundcolor='white')
        plt.grid(True)
        plt.ylim([0, max_y])
        plt.xlim([0, np.max(J_vec)])
    plt.savefig(fig_lab[i])
plt.show()
