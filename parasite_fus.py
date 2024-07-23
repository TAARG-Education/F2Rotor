#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#x Simple estimation of drag coefficient f_Af and equivalent wetted area f_f of the fuselage in forward flight in case the following geometrical data and basic parameters   x
#x are known (all parameters listed are defined as objects of the DataImport class):                                                                                         x
#x  - Length L and cross-sectional dimensions h and d of the fuselage                                                                                                        x
#x  - Helicopter type (conventional or compound)                                                                                                                             x
#x  - Reference area S (rotor disk area)                                                                                                                                     x
#x  - Asymptotic Mach number M_inf                                                                                                                                           x
#x  - Asymptotic Reynolds number R_inf based on fuselage length                                                                                                              x
#x                                                                                                                                                                           x
#x Procedure:                                                                                                                                                                x
#x  1. For compound helicopters: calculation of wing-fuselage interference factor R_WB from experimental data stored in the values.json file; for conventional helicopters:  x
#x     R_WB = 1 is assumed (inter_fact function)                                                                                                                             x
#x  2. Application of the NASA method for calculating the fuselage drag coefficient (fuselage function):                                                                     x
#x       - Calculation of the form factor FF                                                                                                                                 x
#x       - Calculation of the average coefficient of friction using the formula for the average Cf of a flat plate in the case of fully turbulent flow                       x
#x       - Calculation of the wetted area of the fuselage S_wet                                                                                                              x
#x       - Calculation of the drag coefficient f_Af of the fuselage                                                                                                          x
#x       - Calculation of equivalent wetted area (or parasite drag area) f_f of the fuselage                                                                                 x
#x                                                                                                                                                                           x
#x References:                                                                                                                                                               x
#x  - Carlo de Nicola, (2018-2019), Appunti per un corso di Aerodinamica degli Aeromobili, Appendix E: http://wpage.unina.it/denicola/AdA/DOWNLOAD/Appunti_AdA_2018_2019.pdf x                                                                    
#x  - Renato Tognaccini, (2023-2024), Lezioni di Aerodinamica dell'Ala Rotante, section 7.5: http://wpage.unina.it/rtogna/ar2019.pdf                                         x
#x                                                                                                                                                                           x
#x Notes: for usage examples refer to test_cases_parasite_fus.py                                                                                                             x                                                                                       
#x                                                                                                                                                                           x
#x Author: Iole Paolillo                                                                                                                                                     x
#x Latest update: 19/07/2024                                                                                                                                                 x
#x Version: 1.0                                                                                                                                                              x                                                                                                          
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Import necessary modules
import json                                          # needed to read the values.json file 
import numpy as np                                   
from scipy.interpolate import RectBivariateSpline    # needed to interpolate experimental data
import math as ma

# Definition of the parameter class to be provided as input
class DataImport:                                
    def __init__(self, Re_inf, M_inf, config, S, helicopter_type):
        self.Re_inf = Re_inf
        self.M_inf = M_inf
        self.config = config                         # geometric dimensions of the fuselage
        self.S = S                                   # reference area 
        self.R_WB = None                             # initialize R_WB to None
        self.helicopter_type = helicopter_type       # helicopter type (conventional or compound)

def inter_fact(data_import):
    '''Function that evaluates the fuselage and wing correlation factor R_WB. The correlation factor is interpolated from experimental curves (data are present within the values.json file)
	     - Input (DataImport instance): 
           1) geometric dimensions of the fuselage (parameters h, d, L)
           2) S: reference area (rotor disk area)
           3) M_oo: asymptotic Mach number 
           4) R_oo: asymptotic Reynolds number
           5) Helicopter type (conventional or compound) 
	     - Output: R_WB (fuselage-wing correlation factor)
    '''
    if data_import.helicopter_type == "conventional":
       data_import.R_WB = 1
    else:
        # The values.json file is opened in read mode (“r”) and its contents are loaded into the config_data dictionary
        data = "values.json"
        with open(data, "r") as file:
             config_data = json.load(file)

        # Data extraction from the dictionary (Reynolds numbers and fuselage-wing correlation factors):
        #  - extraction of Reynolds number and interference factor corresponding to M_inf1 = 0.25
        Re1 = config_data["Re_inf1"]           
        Rwb1 = config_data["R_wb1"]
        #  - extraction of Reynolds number and interference factor corresponding to M_inf2 = 0.4
        Re2 = config_data["Re_inf2"]
        Rwb2 = config_data["R_wb2"]
            
        # Determination of the range of Reynolds numbers over which to do the interpolation
        min_max = min(max(Re1), max(Re2))
        max_min = max(min(Re1), min(Re2))
        Re = np.linspace(max_min, min_max, 1000)  # array of Reynolds numbers equally spaced between max_min and min_max
            
        # Interpolation of Rwb1 and Rwb2 values on the new Re values. These interpolated values are combined into a 2D array (Rwb_2D) transposed with .T 
        Rwb1 = np.interp(Re, Re1, Rwb1)
        Rwb2 = np.interp(Re, Re2, Rwb2)
        Rwb_2D = np.array([Rwb1, Rwb2]).T

        # Exceptions handling (checks whether the Mach number (M_inf) and Reynolds number (Re_inf) fall within valid ranges. If this does not happen, an error message is printed)
        if data_import.M_inf > 0.4:
            print('Error! Mach number exceeds the maximum value. ')
            return None  # Return None in case of error

        if data_import.Re_inf < max_min:
            print('Error! Reynolds number too low.')
            return None  # Return None in case of error

        if data_import.Re_inf > min_max:
           print('Error! Reynolds number too high.')
           return None  # Return None in case of error

        M_vec = np.array([0.25, 0.4])

        # Bivariate interpolation:
        # Re and M_vec are used as input coordinates
        # Rwb_2D is the corresponding output data matrix
        # bbox specifies the bounds of the grid
        # kx=1 and ky=1 indicate that a linear polynomial is used for interpolation
        spl = RectBivariateSpline(Re, M_vec, Rwb_2D, bbox=[min(Re), max(Re), min(M_vec), max(M_vec)], kx=1, ky=1)
            
        # Calculation of the interpolated value of R_WB for the specified values of Re_inf and M_inf
        R_WB_arr = spl(data_import.Re_inf, data_import.M_inf)   

        # The interpolated value is extracted as a scalar number (item) and provided in output
        data_import.R_WB = R_WB_arr.item()
    return data_import.R_WB
    
def fuselage(data_import):
    '''This is the method to calculate the parasite drag area of the fuselage (the numbering of the formulas refers to C. de Nicola)
        - Input (DataImport instance): 
            1) geometric dimensions of the fuselage (parameters h, d, L)
            2) S: reference area (rotor disk area)
            3) M_oo: asymptotic Mach number 
            4) R_oo: asymptotic Reynolds number
            5) R_WB: wing-fuselage interference factor (call to inter_fact function)

	    - Output: 
            1) f_Af (drag coefficient of the fuselage)
	        2) f_f (equivalent wetted area f=Cf*Swet)
        '''
    # Extraction of geometric and aerodynamic parameters
    h = data_import.config["ParasiteArea"][0]["h"]
    d = data_import.config["ParasiteArea"][0]["d"]
    L = data_import.config["ParasiteArea"][0]["L"]
    S = data_import.S
    M_oo = data_import.M_inf
    Re_oo = data_import.Re_inf
        
    # Shape factor calculation: formula (E.11)
    FR = L / (h * d) ** 0.5                   # fuselage slenderness (fuselage length L / square root of the product between the two fuselage cross dimensions)
    FF = 1 + 60 * 1 / FR ** 3 + 0.0025 * FR   # fuselage shape factor

    # Wet area of the fuselage: formula (E.10)
    Swet = 0.75 * 3.14 * (h * d) ** 0.5 * L

    # Wing-fuselage interference factor extraction (call inter_fact here)
    IF = inter_fact(data_import)   

    # Calculation the Reynolds number of the fuselage (minimum between Re_L1 and Re_L2 is chosen): formula (E.7)
    Re_L1 = Re_oo
    K1 = 37.587 + 4.615 * M_oo + 2.949 * M_oo ** 2 + 4.132 * M_oo ** 3
    K = 0.000045                              # allowable roughness of the component (dimensionally it is a length)
    Re_L2 = K1 * (L / K) ** 1.0489
    Re_L = min([Re_L1, Re_L2])

    # Calculation of mean fuselage friction coefficient (hp: totally turbulent flow): formula (E.6)
    r = 0.89
    gamma = 1.4                                                           # specific heat ratio (air)
    t = (1 + r * (gamma - 1) / 2 * M_oo ** 2) ** (-1)
    F = 1 + 0.03916 * M_oo ** 2 * t                                       # t and F are factors that account for the effects of compressibility   
    Cf = t * F ** 2 * 0.430 / (ma.log10(Re_L * t ** 1.67 * F)) ** 2.56    # mean fuselage friction coefficient

    # Calculation of output parameters 
    f_Af = (Cf * Swet / S) * FF * IF                                      # fuselage profile drag coefficient: formula (E.5)
    f_f = f_Af * S                                                        # equivalent wetted area
    f_Af = round(f_Af, 5)                                                 # rounding to five decimal digits 
    f_f = round(f_f, 4)                                                   # rounding to four decimal digits

    return f_Af, f_f 


