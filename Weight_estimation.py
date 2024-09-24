# Author: Liberti Marco
# Rotary Wing Aerodynamics course, Prof. Renato Tognaccini
# University of Naples Federico II
# Academic Year 2023-2024
# Date: 17/09/2024
#
# Version: 1.0.0
# import numpy as np: Import NumPy for advanced mathematical operations and array handling.

# Mass Conversion Function (convmass):
# Converts masses between pounds (lbm) and kilograms (kg). The function uses a dictionary of conversion factors and returns the converted masses.

# Main Function (estimation_weight_prouty):
# Parameters: Takes a helicopter dictionary containing the necessary parameters for weight calculations.

# Weight Calculations:
# Each variable represents a specific weight of a part of the helicopter, calculated with an empirical formula based on the provided parameters.
# The formulas are used to estimate the weights of the main rotor, fuselage, tail, landing gear, propulsion system, etc.

# Mass Conversion:
# Matrix_masses is an array containing the weights of all components. It is converted from pounds to kilograms.
# M_empty_kg represents the estimated total weight of the helicopter in kilograms.

# Output:
# Returns M_empty_kg (total weight) and Matrix_masses_kg (array of individual component masses).

import numpy as np

def convmass(masses, from_unit, to_unit):
    """
    Converts masses between units of measurement (lbm to kg and vice versa).
    Parameters:
        masses (array or float): Masses to convert.
        from_unit (str): Original unit of measure (e.g. 'lbm').
        to_unit (str): Destination unit of measure (e.g. 'kg').
    Returns:
        Array or float of converted masses.
    """
    # Conversion factors between lbm and kg
    conversion_factors = {
        ('lbm', 'kg'): 0.453592,
        ('kg', 'lbm'): 2.20462
    }
    
    if from_unit == to_unit:
        return masses
    return masses * conversion_factors[(from_unit, to_unit)]

def estimation_weight_prouty(helicopter):
    """
    Estimates the weights of a helicopter based on provided parameters.
    Parameters:
        helicopter (dict): Dictionary containing the necessary parameters for calculation.
    Returns:
        M_empty_kg (float): Estimated total weight of the helicopter in kilograms.
        Matrix_masses_kg (array): Array of component weights in kilograms.
    """

    # Helicopter parameters (extracted from the 'helicopter' dictionary)
    b = helicopter['b']  # Main rotor blade width (m)
    c = helicopter['c']  # Main rotor blade chord (width) (m)
    R = helicopter['R']  # Main rotor radius (m)
    Omega = helicopter['Omega']  # Main rotor angular velocity (rad/s)
    g = helicopter['g']  # Gravitational acceleration (m/s²)
    A_H = helicopter['A_H']  # Horizontal rotor area (m²)
    AR_H = helicopter['AR_H']  # Horizontal rotor aspect ratio
    A_V = helicopter['A_V']  # Vertical rotor area (m²)
    AR_V = helicopter['AR_V']  # Vertical rotor aspect ratio
    n_tailgearboxes = helicopter['n_tailgearboxes']  # Number of tail gearboxes
    R_T = helicopter['R_T']  # Tail rotor radius (m)
    transm_hp_rating = helicopter['transm_hp_rating']  # Transmission system power (hp)
    Omega_T = helicopter['Omega_T']  # Tail rotor angular velocity (rad/s)
    GW = helicopter['GW']  # Gross weight of the helicopter (lbm)
    L_F = helicopter['L_F']  # Fuselage length (m)
    S_wetF = helicopter['S_wetF']  # Fuselage wetted area (m²)
    n_wheellegs = helicopter['n_wheellegs']  # Number of landing gear legs (wheels)
    N_eng = helicopter['N_eng']  # Number of engines
    Installed_wt_pereng = helicopter['Installed_wt_pereng']  # Installed weight per engine (lbm)
    Cap_In_Gal = helicopter['Cap_In_Gal']  # Fuel tank capacity (gallons)
    N_tanks = helicopter['N_tanks']  # Number of tanks
    RPM_eng = helicopter['RPM_eng']  # Engine speed (rpm)
    tail_hp_rating = helicopter['tail_hp_rating']  # Tail engine power (hp)
    N_gearboxes = helicopter['N_gearboxes']  # Number of gearboxes
    # Note: 'b', 'c', 'Omega', and 'R' are reused and do not require further definition

    # Calculation of main rotor blade weight
    W_bM = 0.026 * b**0.66 * c * R**1.3 * (Omega * R)**0.67
    
    # Calculation of main rotor blade moment of inertia
    J = 0.0311 * W_bM * R**2
    
    # Calculation of main rotor hub and joint weight
    W_hM = 0.0037 * b**0.28 * R**1.5 * (Omega * R)**0.43 * (0.67 * W_bM + g * J / R**2)**0.55
    
    # Total main rotor weight
    W_rotor = W_bM + W_hM
    
    # Calculation of horizontal rotor weight
    W_H = 0.72 * A_H**1.2 * AR_H**0.32
    
    # Calculation of vertical rotor weight
    W_V = 1.05 * A_V**0.94 * AR_V**0.53 * n_tailgearboxes**0.71
    
    # Calculation of tail rotor weight
    W_T = 1.4 * R_T**0.09 * (transm_hp_rating / Omega)**0.9
    
    # Total tail weight
    W_tail = W_H + W_V + W_T
    
    # Calculation of fuselage (body) weight
    W_body = 6.9 * (GW / 1000)**0.49 * L_F**0.61 * S_wetF**0.25
    
    # Calculation of landing gear weight
    W_landinggear = 40 * (GW / 1000)**0.67 * n_wheellegs**0.54
    
    # Calculation of engine installation weight
    W_eng = N_eng * Installed_wt_pereng
    
    # Calculation of propulsion subsystem weight
    W_Pss = 2 * W_eng**0.59 * N_eng**0.79
    
    # Total propulsion weight
    W_prop = W_eng + W_Pss
    
    # Calculation of fuel system weight
    W_FS = 0.43 * Cap_In_Gal**0.77 * N_tanks**0.59
    
    # Calculation of transmission system weight
    W_DS = 13.6 * transm_hp_rating**0.82 * (RPM_eng / 1000)**0.037 * \
           ((tail_hp_rating / transm_hp_rating) * (Omega / Omega_T))**0.068 * \
           N_gearboxes**0.066 / Omega**0.64
    
    # Total transmission weight
    W_trasm = W_FS + W_DS
    
    # Calculation of cockpit controls weight
    W_CC = 11.5 * (GW / 1000)**0.40
    
    # Calculation of system control weight
    W_SC = 36 * b * c**2.2 * (Omega * R / 1000)**3.2
    
    # Total control system weight
    W_controls = W_CC + W_SC
    
    # Auxiliary power unit weight
    W_APP = 150
    W_aux = W_APP
    
    # Calculation of instrument weight
    W_I = 3.5 * (GW / 1000)**1.3
    W_instr = W_I
    
    # Calculation of hydraulic system weight
    W_hyd = 37 * b**0.63 * c**1.3 * (Omega * R / 1000)**2.1
    
    # Calculation of electrical system weight (excluding hydraulic weight)
    W_EL = 9.6 * transm_hp_rating**0.65 / (GW / 1000)**0.4 - W_hyd
    W_hpe = W_hyd + W_EL
    
    # Avionics weight
    W_av = 50
    
    # Calculation of furnishing and equipment weight
    W_FE = 6 * (GW / 1000)**1.3
    W_equip = W_FE
    
    # Weight of air conditioning and anti-icing system
    W_AC_AI = 8 * (GW / 1000)
    W_antiicing = W_AC_AI
    
    # Creating an array with the weights of each component
    Matrix_masses = np.array([W_rotor, W_tail, W_body, W_landinggear, W_prop,
                              W_trasm, W_controls, W_instr, W_hpe, W_av,
                              W_equip, W_antiicing])
    
    # Conversion of masses from pounds (lbm) to kilograms (kg)
    Matrix_masses_kg = convmass(Matrix_masses, 'lbm', 'kg')
    
    # Calculation of total helicopter weight and conversion to kilograms
    M_empty_kg = convmass(np.sum(Matrix_masses), 'lbm', 'kg')
    
    return M_empty_kg, Matrix_masses_kg

