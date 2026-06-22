from modules import bemt
from modules import aero
from modules import cad_v2 as cad
from ambiance import Atmosphere
import os
import numpy as np
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt

results_aero_testcases_folder = 'results_aero'      # Results folder name

## Manage folders (Delete the previous results folder if it already exists, then create it)
if not os.path.exists(results_aero_testcases_folder): os.mkdir(results_aero_testcases_folder)
else: shutil.rmtree(results_aero_testcases_folder); os.mkdir(results_aero_testcases_folder)


## Graphics options
font_size = 18
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman, Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman'
mpl.rcParams['font.size'] = font_size
plt.rc('legend',fontsize=font_size) # using a size in points
# Set ggplot styles and update Matplotlib with them.
ggplot_styles = {
    'axes.grid': True,
    'axes.grid.which': 'both',
    'xtick.major.bottom': True,
    'xtick.minor.bottom': False,
    'ytick.major.left': True,
    'ytick.minor.left': False,
}
plt.rcParams.update(ggplot_styles)
linestyles_list = [':', '--', '-']



## Operating condition
z = 1                                                           # Altitude [km]
rpm = 2000                                                      # Revolutions per minute [1/min]
J_vec = np.arange(0, 2, 0.01)                                   # Advance ratio [-]

# Geometry
D = 1.66                                                        # Blade diameter [m]
nbl = 2                                                         # Number of blades [-]
Rhub = 0.1                                                      # non-dim hub final section [-]
dx = 0.01                                                       # radius increment along blade [-]
pitch = 0.                                                      # pitch [deg]

r_R_known = [0.12048,0.25, 0.5, 0.75, 0.99999999999999999]      # non-dim known radial stations [-]
c_R_known = [0.083,0.1079,0.1411,0.0913,0.0166]                 # non-dim known chords at radial stations [-]
beta_known = [70,55,35,25,0.1]                                  # known pitch at radial stations [-]
Airfoil = ['4412', '4412', '4412', '4412','4412']               # known airfoils at radial stations

'''
Aerodynamic parameters for model 1 and 2:

DEFAULT VALUES:
    aero_params1 = {                            # For aero method 1
        'alpha_0_lift': np.deg2rad(0),      # In deg inside the bracket
        'Cl_alpha':  6.28,                   # 1/rad
        'Cl_max':   1.5,
        'Cl_min':    -1.5,
        'Cd': 0.02
    }
    aero_params2 = {                            # For aero method 2
        'alpha_0_lift': np.deg2rad(0),      # In deg inside the bracket
        'Cl_alpha':  6.28,                   # 1/rad
        'Cl_alpha_stall': 0.1,
        'Cl_max':   1.50,
        'Cl_min':    -0.5,
        'Cl_incr_to_stall':    0.1,
        'Cd_min': 0.013,
        'Cl_at_cd_min': 0.5,
        'dCd_dCl2':   0.004,
        'Mach_crit': 0.8,
        'Re_scaling_exp': -0.4
    }

ESTIMATED VALUE:
    aero_params1 = {                            # For aero method 1
        'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
        'Cl_alpha':  6.3216,                   # 1/rad
        'Cl_max':   1.5789,
        'Cl_min':    -0.9522,
        'Cd': 0.00582
    }
    aero_params2 = {                            # For aero method 2
        'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
        'Cl_alpha':  6.3216,                   # 1/rad
        'Cl_alpha_stall': -2.7215,
        'Cl_max':   1.5789,
        'Cl_min':    -0.9522,
        'Cl_incr_to_stall':    0.2602,
        'Cd_min': 0.00582,
        'Cl_at_cd_min': 0.5522,
        'dCd_dCl2':   0.01145,
        'Mach_crit': 0.8,
        'Re_scaling_exp': -0.4
    }
'''
### CAD ### (ARBI Xrotor)
geom = cad.Geometry(D/2, Rhub, nbl, rpm, r_R_known, c_R_known, beta_known, 0., Airfoil,pitch)  # Generate the blade geometry

# Initialization
ct_matrix = np.full((3, 1000), np.nan); 
cp_matrix = np.full((3, 1000), np.nan); 
eta_matrix = np.full((3, 1000), np.nan)



########################## SECTION 1 - EXTRACT AERODYNAMICS #############################################
'''
In this sections, the extraction of aerodynamic curves is presented (Lift curve, Drag polar)
'''

aero_params1 = {                            # For aero method 1
    'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
    'Cl_alpha':  6.3216,                   # 1/rad
    'Cl_max':   1.5789,
    'Cl_min':    -0.9522,
    'Cd': 0.00582
}
aero_params2 = {                            # For aero method 2
    'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
    'Cl_alpha':  6.3216,                   # 1/rad
    'Cl_alpha_stall': -2.7215,
    'Cl_max':   1.5789,
    'Cl_min':    -0.9522,
    'Cl_incr_to_stall':    0.2602,
    'Cd_min': 0.00582,
    'Cl_at_cd_min': 0.5522,
    'dCd_dCl2':   0.01145,
    'Mach_crit': 0.8,
    'Re_scaling_exp': 0
}

# Instantiate the aero class for the 3 methods
aero_1 = aero.Aerodynamics(aero_method=1, aero_params=aero_params1, M_corr=False, Rey_corr=False, z = z)
aero_2 = aero.Aerodynamics(aero_method=2, aero_params=aero_params2, M_corr=False, Rey_corr=False, z = z)

alpha_vec_xFoil = np.arange(-25,25,1)       # AoAs for xFoil simulations [deg]
aero_3 = aero.Aerodynamics(aero_method = 3, 
                        alpha_vec_xFoil=alpha_vec_xFoil,        
                        M_corr = False, Rey_corr = False, z = z)
aero_3.f_cl, aero_3.f_cd = aero_3.aero_database_generator(Re_ref = 1e+6, airfoils_folder = 'airfoil_input', r_R = r_R_known)   # Generate aerodynamic databases of Cl and Cd with Xfoil
alpha_vec = np.arange(-25, 25, 1)           # AoAs for aerodynamic curves evaluation [deg]

# Initialize Cl and Cd vectors as empty lists
cl_1 = []; cl_2 = []; cl_3 = []
cd_1 = []; cd_2 = []; cd_3 = []
for i, alpha in enumerate(alpha_vec):
    alpha = np.deg2rad(alpha)
    
    cl, cd = aero_1.clcd1(Re_ref =1e+6 , Vr = 100 ,chord = 0.1, AoA=alpha)                        # Compute Cl and Cd with method 1
    cl_1.append(cl); cd_1.append(cd)
    
    cl,cd = aero_2.clcd2(Re_ref = 1e+6, Vr = 100 ,chord = 0.1, AoA = alpha)                       # Compute Cl and Cd with method 2
    cl_2.append(cl); cd_2.append(cd)

    cl, cd = aero_3.clcd3(AoA = np.rad2deg(alpha), Re_ref = 1e+6, Vr = 100 ,chord = 0.1, x = 0.5) # Compute Cl and Cd with method 3
    cl_3.append(cl); cd_3.append(cd)


## Graphics

# Lift curve
plt.figure(figsize=(9,6))
plt.plot(alpha_vec, cl_1, color = 'black', linestyle = ':', label = 'Method 1')
plt.plot(alpha_vec, cl_2, color = 'black', linestyle = '--', label = 'Method 2')
plt.plot(alpha_vec, cl_3, color = 'black', linestyle = '-', label = 'Method 3')
plt.xlabel(r'$\bf \alpha [°]$'); plt.ylabel(r'$\bf Cl$') 
plt.minorticks_on()
plt.grid(True, which = 'both')
plt.legend()
plt.savefig('Lift curve comparison')

# Drag polar
plt.figure(figsize=(9,6))
plt.plot(cl_1, cd_1, color = 'black', linestyle = ':', label = 'Method 1')
plt.plot(cl_2, cd_2, color = 'black', linestyle = '--', label = 'Method 2')
plt.plot(cl_3, cd_3, color = 'black', linestyle = '-', label = 'Method 3')
plt.xlabel(r'$\bf Cl$'); plt.ylabel(r'$\bf Cd$')
plt.minorticks_on()
plt.grid(True, which = 'both')
plt.legend()
plt.savefig('Drag polar comparison')

##################### SECTION 2 - EXTRACT PROPELLER PERFORMANCE ####################################

## Instantiate the aerodynamics class ###
for k in [3,2,1]:

    #Initialize Ct, Cp, eta as empty lists
    ct_vec = []; cp_vec = []; eta_vec = []
    aero_method = k
    M_corr = True; Rey_corr = True
    if aero_method not in (1, 2, 3):              # Control on the aero_method variable
        raise ValueError("The variable 'aero_method' must be 1, 2 or 3.")
    else: print(f"You are currently using aerodynamic model {aero_method}.")
    
    aero_params1 = {                              # For aero method 1
        'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
        'Cl_alpha':  6.3216,                      # 1/rad
        'Cl_max':   1.5789,
        'Cl_min':    -0.9522,
        'Cd': 0.00582
    }
    aero_params2 = {                              # For aero method 2
        'alpha_0_lift': np.deg2rad(-4.4107),      # In deg inside the bracket
        'Cl_alpha':  6.3216,                      # 1/rad
        'Cl_alpha_stall': -2.7215,
        'Cl_max':   1.5789,
        'Cl_min':    -0.9522,
        'Cl_incr_to_stall':    0.2602,
        'Cd_min': 0.00582,
        'Cl_at_cd_min': 0.5522,
        'dCd_dCl2':   0.01145,
        'Mach_crit': 0.8,
        'Re_scaling_exp': 0
    }

    # Instantiate aero class
    if aero_method == 1:
        aero_k = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params1, M_corr=M_corr, Rey_corr=Rey_corr, z = z)

    elif aero_method == 3:
        alpha_vec_xFoil = np.arange(-25,25,1)
        aero_k = aero.Aerodynamics(aero_method = aero_method, 
                                alpha_vec_xFoil=alpha_vec_xFoil,        # AoAs to perform xFoil analysis
                                M_corr = M_corr, Rey_corr = Rey_corr, z = z) 
        aero_k.f_cl, aero_k.f_cd = aero_k.aero_database_generator(Re_ref = 1e+6, airfoils_folder = 'airfoil_input', r_R = r_R_known) # Generate aerodynamic databases of Cl and Cd with Xfoil

    else: aero_k = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params2, M_corr=M_corr, Rey_corr=Rey_corr, z = z)


    ## BEMT ##
    for i,J in enumerate(J_vec):

        ct_tvor, cp_tvor= bemt.BEMT_tvor(z,J,dx,geom,aero_k,curvature=False, thickness=False)       # Compute Ct and Cp (ONLY with theory tvor!!)

        eta_tvor = ct_tvor*J/cp_tvor
        ct_vec.append(ct_tvor); cp_vec.append(cp_tvor); eta_vec.append(eta_tvor)

        if (cp_tvor < 0).any(): break
    
    ct_vec = np.array(ct_vec); cp_vec = np.array(cp_vec); eta_vec = np.array(eta_vec)       # From lists to arrays
    ind = np.where(ct_vec < 0)[0][0]                                                        # Find index where Ct < 0 for efficiency graphics
    eta_vec[ind+1:] = np.nan
    ct_matrix[k-1, :len(ct_vec)] = ct_vec                                                   # Insert Ct in the Ct matrix
    cp_matrix[k-1, :len(cp_vec)] = cp_vec
    eta_matrix[k-1, :len(eta_vec)] = eta_vec
    if k == 3:
        fig_ct, axs_ct = plt.subplots(figsize = (9,6))
        fig_cp, axs_cp = plt.subplots(figsize = (9,6))
        fig_eta, axs_eta = plt.subplots(figsize = (9,6))

    ## Plot Ct
    axs_ct.plot(J_vec[:len(cp_vec)], ct_vec[:len(cp_vec)], color = 'black', linestyle = linestyles_list[k-1], label = f'Method {k}')
    axs_ct.set_xlim(left = 0); axs_ct.set_ylim(bottom = 0, top = 0.12)
    axs_ct.set_xlabel(r'$\bf J$'); axs_ct.set_ylabel(r'$\bf C_T$') 

    ## Plot Cp
    axs_cp.plot(J_vec[:len(cp_vec)], cp_vec[:len(cp_vec)], color = 'black',linestyle = linestyles_list[k-1], label = f'Method {k}')
    axs_cp.set_xlim(left = 0); axs_cp.set_ylim(bottom = 0, top = 0.08)
    axs_cp.set_xlabel(r'$\bf J$'); axs_cp.set_ylabel(r'$\bf C_P$')

    ## Plot eta
    axs_eta.plot(J_vec[:len(cp_vec)], eta_vec[:len(cp_vec)], color = 'black',linestyle = linestyles_list[k-1], label = f'Method {k}')
    axs_eta.set_xlim(left = 0); axs_eta.set_ylim([0, 1])
    axs_eta.set_xlabel(r'$\bf J$'); axs_eta.set_ylabel(r'$\eta$')

# Graphics options
axs_ct.minorticks_on()
axs_ct.grid(True, which = 'both')
axs_cp.minorticks_on()
axs_cp.grid(True, which = 'both')
axs_eta.minorticks_on()
axs_eta.grid(True, which = 'both')
axs_ct.legend()
axs_cp.legend()
axs_eta.legend()

# Save figures
fig_ct.savefig(os.path.join(results_aero_testcases_folder, 'ct_aero_testcase'))
fig_cp.savefig(os.path.join(results_aero_testcases_folder, 'cp_aero_testcase'))
fig_eta.savefig(os.path.join(results_aero_testcases_folder, 'eta_aero_testcase'))

# Show graphics
plt.show()