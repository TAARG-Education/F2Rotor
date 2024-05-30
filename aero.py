#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#x File: aero.py                                                                                                                                                x
#x Authors: Ciro Cuozzo, Daniele Trincone                                                                                                                       x
#x Date: 28/05/2024                                                                                                                                             x
#x Version: 1.03                                                                                                                                                x
#x                                                                                                                                                              x
#x This code calculates the aerodynamic parameters of a chosen station on the blade through a Python class.                                                     x
#x Within the class there are functions capable of calculating the Cl_alpha, the alpha_zl, the Cl and the Cd for any                                            x
#x chosen station.                                                                                                                                              x
#x Three methods are implemented for calculating the coefficients:                                                                                              x
#x                                                                                                                                                              x
#x - Method '1': The aerodynamics is given by Cl = Cl_alpha * (alpha - alpha_0_lift), limited by Cl_max and Cl_min and constant Cd;                             x
#x   With this method the class must be instantiated by providing a dictionary with the following parameters:                                                   x
#x   'alpha_0_lift' (rad), 'Cl_alpha' (1/rad), 'Cl_max', 'Cl_min', 'Cd'                                                                                         x
#x                                                                                                                                                              x
#x - Method '2': The aerodynamics is reconstructed through synthetic parameters.                                                                                x
#x   With this method the class must be instantiated by providing a dictionary with the following parameters:                                                   x
#x   'alpha_0_lift' (rad), 'Cl_alpha' (1/rad), 'Cl_alpha_stall' (1/rad), 'Cl_max', 'Cl_min' 'Cl_incr_to_stall', 'Cd_min', 'Cl_at_cd_min',                       x
#x   'dCd_dCl2', 'Mach_crit', 'Re_scaling_exp'                                                                                                                  x
#x                                                                                                                                                              x
#x - Method '3': The aerodynamics is computed by interpolating on a response surface generated using xFoil.                                                     x
#x   The independant variables are the section and the AoA, the dependant ones are Cl and Cd.                                                                   x
#x                                                                                                                                                              x
#x   Input to instantiate the class:                                                                                                                            x
#x   - aero_method: Int (1, 2 or 3), one of the 3 methods above specified;                                                                                      x
#x   - M_corr: Bool, if 'True' applies compressibility correction;                                                                                              x
#x   - Rey_corr: Bool, if 'True' applies Reynolds number correction;                                                                                            x
#x   - z: float, altitude in kilometers [km];                                                                                                                        x
#x   - aero_params (Optional): Dictionary, must be provided for methods '1' and '2'. Dictionary with the necessary parameters specified above;                  x
#x   - alpha_vec_xFoil (Optional): List, must be provided for method '3'. List of desired AoAs to perform xFoil analysis.                                       x
#x                                                                                                                                                              x
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import math
from ambiance import Atmosphere as atm
from scipy.interpolate import interp2d

os.system('cls')

class Aerodynamics():
    def __init__(self, aero_method, z, M_corr = True, Rey_corr = True, aero_params = None, alpha_vec_xFoil = None):
        self.aero_method = aero_method
        self.M_corr = M_corr
        self.Rey_corr = Rey_corr
        self.a_sound = atm(z*1000).speed_of_sound[0]
        self.mu = atm(z*1000).dynamic_viscosity[0]
        self.rho = atm(z*1000).density[0]
        if self.aero_method == 1 or self.aero_method == 2:
            self.aero_params = aero_params
        else:
            self.alpha_vec = alpha_vec_xFoil

    def aero_database_generator(self, Re_ref, airfoils_folder, r_R):
        '''
        This function generates a response surface Cl = Cl (r/R, alpha) using xFoil at the known sections using the AoAs assigned
        from the user in the class initialization.

        Input: 
        - Re_ref: Reference Reynolds number (used for xFoil calculations);
        - airfoils_folder: The folder in which the airfoils are stored;
        - r_R: the sections assigned from the user.
        
        Output:
        - f_cl: The interpolation function for 2D lift coefficient
        - f_cd: The interpolation function for 2D drag coefficient
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''
        if not os.path.exists(airfoils_folder): # Control on the folder containing airfoils
            print("airfoils_folder does not exist.")
            return None
        
        
        file_coord = [file for file in os.listdir(airfoils_folder) if file.endswith('.txt')] # Count .txt files in the airfoils folder
    
        Cl_database = np.full((len(file_coord), len(self.alpha_vec)), np.nan) # Generate the Cl database for interpolation
        Cd_database = np.full((len(file_coord), len(self.alpha_vec)), np.nan) # Generate the Cd database for interpolation
        for i in range(len(file_coord)):
            print(f'Currently generating aerodynamic database for airfoil {i+1}...')
            if os.path.exists(f'airfoil_polar{i+1}.dat'): os.remove(f'airfoil_polar{i+1}.dat')
            for alpha in self.alpha_vec:
                self.xfoil_run(file_coord[i], Re_ref, alpha, i) # Run xFoil for each AoA in a for loop

            polar_i = np.loadtxt(f'airfoil_polar{i+1}.dat', skiprows=12)
            Cl_database[i,:] = np.interp(self.alpha_vec, polar_i[:,0], polar_i[:,1]) # Constant extrapolation
            Cd_database[i,:] = np.interp(self.alpha_vec, polar_i[:,0], polar_i[:,2]) # Constant extrapolation
            
        f_cl = interp2d(self.alpha_vec, r_R, Cl_database, kind='linear') # Create the interpolation function for Cl
        f_cd = interp2d(self.alpha_vec, r_R, Cd_database, kind='linear') # Create the interpolation function for Cd
        self.f_cl = f_cl # Assign the Cl function to the class
        self.f_cd = f_cd # Assign the Cd function to the class
        return f_cl, f_cd 
    
    def eval_lift_properties(self, x, Vr):
        '''
        This function evaluates Cl_alpha and alpha_zero lift of the airfoil at section 'x' and Mach number 'M'
        The evaluation is performed interpolating on the grid defined in 'aero_database generator'.

        Input:
        - x: the query r/R;
        - Vr: section effective velocity.
        
        Output:
        - Cl_alpha (1/rad)
        - alpha_0_lift (rad)
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''
        
        # Points to interpolate
        Cl_alpha = (self.f_cl(2, x)[0] - self.f_cl(-2, x)[0])/(2-(-2))   # Evaluate Cl alpha based on alpha = +-2°
        # Interpolation on 0
        Cl = []
        alpha_0_lift_vec_interp = np.linspace(-6,6,10) 
        for alpha in alpha_0_lift_vec_interp:
            Cl.append(self.f_cl(alpha, x)[0])  # Evaluate Cl at the predefined AoAs = +-2°
    
        alpha_0_lift = np.interp(0, Cl, alpha_0_lift_vec_interp)   # Evaluate alpha_0_lift by interpolating [deg]
        M = Vr/self.a_sound
        if self.M_corr:
            Cl_alpha /= (1 - (M)**2)**0.5                              # Mach number correction
        return np.rad2deg(Cl_alpha), np.deg2rad(alpha_0_lift)      # Cl_alpha in 1/rad, alpha_0_lift in rad    
    
    def clcd1(self, AoA, Re_ref, Vr, chord):
        '''
        This function evaluates lift coefficient and drag coefficient using the aerodynamics method '1'. 
        Further details in the aerodynamics class.

        Input:
        - AoA: Section angle of attack (rad);
        - Re_ref: Reference Reynolds number;
        - Vr: section effective velocity;
        - chord: section chord.
        
        Output:
        - lift_coeff: Section lift coefficient
        - drag_coeff: Section drag coefficient
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''
        ni = self.mu/self.rho
        M = Vr/self.a_sound
        Re =Vr*chord/ni
        lift_coeff = min(self.aero_params['Cl_alpha']*(AoA - self.aero_params['alpha_0_lift'])/math.sqrt(1-M**2), self.aero_params['Cl_max'])    # Positive Stall at Cl = Cl_max    
        lift_coeff = max(lift_coeff, self.aero_params['Cl_min'])       # Negative Stall at Cl = Cl_min
        drag_coeff = self.aero_params['Cd']                                        # Default Cd
        if Re < 1e+5:
            f = -0.4                                    # Empirical factor for Reynolds number correction
        elif Re >= 1e+5 and Re < 1e+6:
            f = -1
        else:
            f = -0.15
        if self.Rey_corr:
            drag_coeff *= (Re/Re_ref)**f                # Empirical Reynolds number correction
        if self.M_corr:
            lift_coeff /=  math.sqrt(1-M**2)            # Prandtl Glauert correction
        return lift_coeff, drag_coeff
    
    def clcd2(self, AoA, Re_ref, Vr, chord):
        '''
        This function evaluates lift coefficient and drag coefficient using the aerodynamics method '2'. 
        Further details in the aerodynamics class.
        
        Reference: xRotor documentation
        
        Input:
        - AoA: Section angle of attack (rad);
        - Re_ref: Reference Reynolds number;
        - Vr: section effective velocity;
        - chord: section chord.
        
        Output:
        - lift_coeff: Section lift coefficient
        - drag_coeff: Section drag coefficient
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''
        # atmosphere = Atmosphere(0)
        ni = self.mu/self.rho
        M = Vr/self.a_sound
        Re =Vr*chord/ni
        
        CDMSTALL  =  0.1000
        CDMFACTOR = 10.0
        CLMFACTOR =  0.25
        MEXP      =  3.0
        CDMDD     =  0.0020
        PG = 1/math.sqrt(1-M**2)
        lift_coeff = self.aero_params['Cl_alpha'] * (AoA - self.aero_params['alpha_0_lift'])
        if self.M_corr:
            lift_coeff *= PG

        # --- Effective CLmax is limited by Mach effects ---

        DMSTALL = (CDMSTALL / CDMFACTOR) ** (1.0 / MEXP)
        CLMAXM = max(0.0, (self.aero_params['Mach_crit'] + DMSTALL - M) / CLMFACTOR) + self.aero_params['Cl_at_cd_min']
        CLMAX = min(self.aero_params['Cl_max'], CLMAXM)
        CLMINM = min(0.0, -(self.aero_params['Mach_crit'] + DMSTALL - M) / CLMFACTOR) + self.aero_params['Cl_at_cd_min']
        CLMIN = max(self.aero_params['Cl_min'], CLMINM)

        # --- CL limiter function (turns on after +-stall) ---

        ECMAX = math.exp(min(200.0, float((lift_coeff - CLMAX) / self.aero_params['Cl_incr_to_stall'])))
        ECMIN = math.exp(min(200.0, float((CLMIN - lift_coeff) / self.aero_params['Cl_incr_to_stall'])))
        CLLIM = self.aero_params['Cl_incr_to_stall'] * math.log((1.0 + ECMAX) / (1.0 + ECMIN))

        # --- Subtract off a (nearly unity) fraction of the limited CL function ---

        FSTALL = self.aero_params['Cl_alpha_stall'] / self.aero_params['Cl_alpha']
        lift_coeff = lift_coeff - (1.0 - FSTALL) * CLLIM

        # ------------------------- Drag coefficient --------------------------#
        if Re <= 0:
            RCORR = 1.0
        else:
            RCORR = (Re / Re_ref) ** self.aero_params['Re_scaling_exp']
        # --- Drag parabolic in the linear lift range ---
        drag_coeff = (self.aero_params['Cd_min'] + self.aero_params['dCd_dCl2'] * (lift_coeff - self.aero_params['Cl_at_cd_min']) ** 2)
        if self.Rey_corr:
            drag_coeff*= RCORR

        # --- Post Stall Drag ---

        FSTALL = self.aero_params['Cl_alpha_stall'] / self.aero_params['Cl_alpha']
        DCDX = (1.0 - FSTALL) * CLLIM / (PG * self.aero_params['Cl_alpha'])
        DCD = 2.0 * DCDX ** 2

        # --- Compressibility Drag ---

        DMDD = (CDMDD / CDMFACTOR) ** (1.0 / MEXP)
        CRITMACH = self.aero_params['Mach_crit'] - CLMFACTOR * abs(lift_coeff - self.aero_params['Cl_at_cd_min']) - DMDD       
        if M < CRITMACH:
            CDC = 0.0
        else:
            CDC = CDMFACTOR * (M - CRITMACH) ** MEXP

        FAC = 1.0
        # FAC = PG     

        # Total drag terms
        drag_coeff = FAC * drag_coeff + DCD + CDC
        
        return lift_coeff, drag_coeff

    def clcd3(self, AoA, Re_ref, Vr, chord, x):
        '''
        This function evaluates lift coefficient and drag coefficient using the aerodynamics method '3'. 
        Further details in the aerodynamics class.

        Input:
        - AoA: Section angle of attack (rad);
        - Re_ref: Reference Reynolds number;
        - Vr: section effective velocity;
        - chord: section chord.
        
        Output:
        - lift_coeff: Section lift coefficient
        - drag_coeff: Section drag coefficient
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''
        ni = self.mu/self.rho
        M = Vr/self.a_sound
        Re =Vr*chord/ni
        
        # Points to interpolate
        lift_coeff = self.f_cl(AoA, x)[0]
        drag_coeff = self.f_cd(AoA, x)[0]

        if Re < 1e+5:
            f = -0.4                                    # Empirical factor for Reynolds number correction
        elif Re >= 1e+5 and Re < 1e+6:
            f = -1
        else:
            f = -0.15
        if self.Rey_corr:
            drag_coeff *= (Re/Re_ref)**f    
        if self.M_corr:
            lift_coeff /=  math.sqrt(1-M**2) 
        return lift_coeff, drag_coeff
    
    def xfoil_run(self, file_coord, Re, alpha,  i, n_iter = 50, airfoil_folder = 'airfoil_input'):
        '''Function to run xFoil
        'xfoil.exe' must be in the same folder as the main
        
        Reference: xFoil documentation
        
        Input:
        - file_coord: Coordinate file;
        - Re: Reynolds number used xFoil;
        - alpha: Angle of Attack used in xFoil;
        - i: index of the airfoil to save correct polars;
        - n_iter: Number of iterations in xFoil;
        - airfoil_folder: folder in which the airfoils are stored.
        
        Output:
        - Generation of the file containing Cl, Cd, Cm... at each angle of attack given to xFoil.
        
        Authors: Ciro Cuozzo, Daniele Trincone
        Date: 30/05/2024
        Version: 1.04
        '''

        xfoil_path = 'xfoil.exe'
        # Write xfoil command file
        XfoilCmd_fileName = 'Xfoil_cmd.inp'
        with open(XfoilCmd_fileName, 'w') as f:
            f.write('PLOP\n') 
            f.write('G\n\n')  
            # Load coordinates
            f.write(f'load {airfoil_folder}\{file_coord}\n')
            # Set the number of panels
            f.write('ppar\n')  
            f.write('N\n')
            f.write('160\n\n\n')
            # Close TE
            f.write('gdes\n')  
            f.write('tgap\n')  
            f.write('0\n')     
            f.write('0.2\n')   
            f.write('exec gset\n\n')  
            # Set Re and M
            f.write('oper\n')
            f.write('Re\n')
            f.write(f'{Re}\n')
            f.write('M\n') 
            f.write('0\n') # Incompressible analysis
            f.write('v\n') # Viscous analysis
            # Change number of iterations
            f.write('iter\n')
            f.write('%d\n' % n_iter)
            # Polar accumulation
            f.write('pacc\n')
            f.write(f'airfoil_polar{i+1}.dat\n\n')
            f.write('a\n')
            f.write(f'{alpha}\n')
            f.write('pacc\n\n')
            f.write('quit\n')
            
        # Execute xfoil
        cmd = f'"{xfoil_path}" < Xfoil_cmd.inp > xfoil.out'
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('Xfoil execution failed:', e)
            return [float('NaN')]

####################### THE FOLLOWING PART IS JUST AN EXAMPLE OF HOW TO USE IT. IT MUST BE COMMENTED ########################
'''
How to Instantiate the class: 
aero_method = 3
M_corr = True; Rey_corr = True
if aero_method not in (1, 2, 3):
    raise ValueError("The variable 'aero_method' must be 1, 2 or 3.")
else: print(f"You are currently using aerodynamic model {aero_method}.")

aero_params1 = {                            # For aero method 1
    'alpha_0_lift': np.deg2rad(-2.34),      # In deg inside the bracket
    'Cl_alpha':  6.28,                   # 1/rad
    'Cl_max':   1.5,
    'Cl_min':    -1.5,
    'Cd': 0.02
}
aero_params2 = {                            # For aero method 2
    'alpha_0_lift': np.deg2rad(-2.34),      # In deg inside the bracket
    'Cl_alpha':  6.160,                   # 1/rad
    'Cl_alpha_stall': -0.624,
    'Cl_max':   1.40,
    'Cl_min':    -1.12,
    'Cl_incr_to_stall':    0.136,
    'Cd_min': 0.0054,
    'Cl_at_cd_min': 0.318,
    'dCd_dCl2':   0.0125,
    'Mach_crit': 0.8,
    'Re_scaling_exp': -1.4
}

if aero_method == 1:
    aero = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params1, M_corr=M_corr, Rey_corr=Rey_corr)

elif aero_method == 3:
    alpha_vec_xFoil = np.arange(-25,25,1)
    aero = aero.Aerodynamics(aero_method = aero_method, 
                            alpha_vec_xFoil=alpha_vec_xFoil,        # AoAs to perform xFoil analysis
                            M_corr = M_corr, Rey_corr = Rey_corr) 
    
    aero.f_cl, aero.f_cd = aero.aero_database_generator(Re_ref = 1e+6, airfoils_folder = 'airfoil_input', r_R = r_R_known)

else: aero = aero.Aerodynamics(aero_method=aero_method, aero_params=aero_params2, M_corr=M_corr, Rey_corr=Rey_corr)


## --------------------- HOW TO COMPUTE Cl alpha and alpha zero lift --------------------- ##

    if aero.aero_method == 1 or aero.aero_method == 2:
        cla = aero.aero_params['Cl_alpha']
        if aero.M_corr:
            cla/=(1-(Vr/a_sound)**2)**0.5
        bo = aero.aero_params['alpha_0_lift']
    else: 
        cla, bo = aero.eval_lift_properties(x = x, M = Vr/a_sound)

## --------------------- HOW TO COMPUTE Cl and Cd --------------------- ##
    if aero.aero_method == 1:
        cl, cd = aero.clcd1(Re_ref =1e+6 ,Re =Vr*chord/ni ,M = Vr/a_sound,AoA=alpha)
    elif aero.aero_method == 2:
        cl,cd = aero.clcd2(Re_ref = 1e+6, Re = Vr*chord/ni, M = Vr/a_sound, AoA = alpha)
    else: 
        cl, cd = aero.clcd3(AoA = np.rad2deg(alpha), Re_ref = 1e+6, Re = Vr*chord/ni, M = Vr/a_sound, x = x)
'''
