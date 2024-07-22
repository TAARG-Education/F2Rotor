from modules import bemt
from modules import aero
from modules import cad_v2 as cad
from ambiance import Atmosphere as atm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt



## USER INTERFACE ##
# operating condition
z = 1 #km
rpm = 2000
J_vec = np.arange(0, 2, 0.01)      # advance ratio
mu = atm(z*1000).dynamic_viscosity       # dynamic viscosity
rho = atm(z*1000).density[0]       # density
# Geometry
D = 1.66         # blade diameter
nbl = 2          # blade number
Rhub = 0.1       # hub final section 
dx = 0.01        # radius increment along blade
pitch = 0.       # pitch
r_R_known = [0.12048,0.25, 0.5, 0.75, 0.9999999999999999]    # 0.12048 #r/R station
c_known = [0.083,0.1079,0.1411,0.0913,0.0166]                # chord distribution
beta_known = [70,55,35,25,0.1]                               # beta distribuction
Airfoil = ['4412', '4412', '4412', '4412','4412']            # airfoil

### CAD ### (ARBI Xrotor)
geom = cad.Geometry(D/2, Rhub, nbl, rpm, r_R_known, c_known, beta_known,0,Airfoil,pitch)

### Instantiate the aerodynamics class ###
aero_method = 2
M_corr = True; Rey_corr = True
if aero_method not in (1, 2, 3):
    raise ValueError("The variable 'aero_method' must be 1, 2 or 3.")
else: print(f"You are currently using aerodynamic model {aero_method}.")

aero_params1 = {                            # For aero method 1
    'alpha_0_lift': np.deg2rad(-2.34),      # In deg inside the bracket
    'Cl_alpha':  6.28,                      # 1/rad
    'Cl_max':   1.5,
    'Cl_min':    -1.5,
    'Cd': 0.02
}
aero_params2 = {                                # For aero method 2
    'alpha_0_lift': np.deg2rad(-4.3865),        # In deg inside the bracket
    'Cl_alpha':  6.4582,                        # 1/rad
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

if aero_method == 1:
    aero = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params1, M_corr=M_corr, Rey_corr=Rey_corr, z = z)

elif aero_method == 3:
    alpha_vec_xFoil = np.arange(-25,25,1)
    aero = aero.Aerodynamics(aero_method = aero_method, 
                            alpha_vec_xFoil=alpha_vec_xFoil,                # AoAs to perform xFoil analysis
                            M_corr = M_corr, Rey_corr = Rey_corr, z = z) 
    
    aero.f_cl, aero.f_cd = aero.aero_database_generator(Re_ref = 1e+6, airfoils_folder = 'airfoil_input', r_R = r_R_known)

else: aero = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params2, M_corr=M_corr, Rey_corr=Rey_corr, z = z)

# initialize variables
ct_matrix = np.full((len(J_vec),3),np.nan)
cp_matrix = np.full((len(J_vec),3),np.nan)
eta_matrix = np.full((len(J_vec),3),np.nan)

## BEMT ##
for i,J in enumerate(J_vec):
        ct_timp, cp_timp= bemt.BEMT_timp(z,J,dx,geom,aero, False, False) # Momentum theory

        ct_tvorpd, cp_tvorpd = bemt.BEMT_tvorpd(z,J,dx,geom,aero, True, False) # Vortical theory and small disturbance	

        ct_tvor, cp_tvor = bemt.BEMT_tvor(z,J,dx,geom,aero, False, True) # vortical theory 

        # Calculate efficiency
        eta_timp = ct_timp*J/cp_timp
        eta_tvorpd = ct_tvorpd*J/cp_tvorpd
        eta_tvor = ct_tvor*J/cp_tvor
       # Fill matrices 
        ct_matrix[i,:] = [ct_timp, ct_tvorpd, ct_tvor]
        cp_matrix[i,:] = [cp_timp, cp_tvorpd, cp_tvor]
        eta_matrix[i,:] = [eta_timp, eta_tvorpd, eta_tvor]
       # Stop iterations if the power coefficient is negative
        if (cp_matrix[:,2] < 0).any(): break
        
        print(f'simulation for J = {J} done.')
        
ind = np.where(ct_matrix[:,2] < 0)[0][0]
indx = np.where(ct_matrix[:,0] < 0)[0][0]
eta_matrix[indx+1:,0] = np.nan
eta_matrix[ind+1:,:] = np.nan

##### OUTPUT #################

# Plot CT vs J
plt.figure()
plt.plot(J_vec, ct_matrix[:, 0], label=r'$T_{i}$', color='black', linestyle='--')
plt.plot(J_vec, ct_matrix[:, 1], label=r'$T_{v_{PD}}$', color='black', linestyle='-')
plt.plot(J_vec, ct_matrix[:, 2], 'k-', label=r'$VT$')
plt.xlabel(r'$\bf J$'); plt.ylabel(r'$\bf C_T$')  
plt.grid(True, which='both')
plt.minorticks_on()
plt.ylim(0)
plt.xlim(0)
plt.subplots_adjust(hspace=0.5)
plt.legend()
plt.tight_layout()

# Plot CP vs J
plt.figure()
plt.plot(J_vec, cp_matrix[:, 0], label=r'$T_{i}$', color='black', linestyle='--')
plt.plot(J_vec, cp_matrix[:, 1], label=r'$T_{v_{PD}}$', color='black', linestyle='-')
plt.plot(J_vec, cp_matrix[:, 2], 'k-')
plt.xlabel(r'$\bf J$'); plt.ylabel(r'$\bf C_P$')  
plt.grid(True, which='both')
plt.minorticks_on()
plt.ylim(0)
plt.xlim(0)
plt.subplots_adjust(hspace=0.5)
plt.tight_layout()
plt.savefig('cp_xrot-comp', bbox_inches='tight')

# Plot Eta vs J
plt.figure()
plt.plot(J_vec, eta_matrix[:, 0], label=r'$T_{i}$', color='black', linestyle='--')
plt.plot(J_vec, eta_matrix[:, 1], label=r'$T_{v_{PD}}$', color='black', linestyle='-')
plt.plot(J_vec, eta_matrix[:, 2], 'k-')
plt.xlabel(r'$\bf J$'); plt.ylabel(r'$\bf \eta$')
plt.grid(True, which='both')
plt.minorticks_on()
plt.ylim(0, 1)
plt.xlim(0)
plt.subplots_adjust(hspace=0.5)
plt.tight_layout()
plt.savefig('eta_xrot-comp', bbox_inches='tight')
plt.show()
        
