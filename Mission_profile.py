# MISSION PROFILE
#
# Description: # This code provides the mission profile of a helicopter, 
#                calculating various flight parameters for an assigned mission such as fuel consumption, 
#                flight time, and performance metrics for different flight phases, including climb,
#                cruise, descent, and hovering. Functions for different flight phases are imported from F2Rotor library
#
# References:  Renato Tognaccini - Lezioni per il corso di Aerodinamica dell'Ala Rotante - Eliche rotori ed aeromotori
#              con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024
#              Di Giorgio, G. (2024), Lezioni integrative dell’insegnamento di Aerodinamica dell'Ala Rotante      
#
# Authors:     Gabriele Giangreco
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic Year 2023-2024.
#
# Date: 01/20/2025.

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd 
from matplotlib import cm                          
from scipy.constants import pi, g
from ambiance import Atmosphere 

from Drone_functions import RigidFW2
from Drone_functions import RotorSolving
from Climb import Helicopter_properties
from BEMT_Hover import BEMT_Hover

def Mission_profile(Helicopter, SFC, Pd, height_c, MTOW, Voo, Voo_d, height_h, time_h,\
                     V_infty, dist_f, Vinf_vec, N, r_segn, c, theta, W, height_f, R, A, f_A, Clalpha, Cd_mean,\
                          Omega, N_rotors, F_Prandtl, height_d, height_f2, dist_f2):
    
    """
    Function to compute the mission profile of a helicopter including climb, hover, forward flight, and descent phases.

    Args:
        Helicopter: Helicopter object containing properties and methods related to the helicopter.
            - An object that contains the helicopter's attributes and relevant functions for the mission profile.
        SFC: Specific fuel consumption [kg/kWh].
            - A measure of fuel efficiency in terms of fuel consumed per unit of power produced over time.
        Pd: Power required at hover [W].
            - The power required to maintain hover for the helicopter.
        height_c: Climb altitude [m].
            - The altitude at which the helicopter will begin climbing.
        MTOW: Maximum take-off weight [N].
            - The maximum allowable weight for take-off, including fuel, passengers, and cargo.
        Voo: Climb airspeed [m/s].
            - The airspeed at which the helicopter will climb.
        Voo_d: Descent airspeed [m/s].
            - The airspeed during the descent phase of the mission.
        height_h: Hover altitude [m].
            - The altitude at which the helicopter performs the hover phase of the mission.
        time_h: Time spent in hover [s].
            - The total time spent in the hover phase.
        V_infty: Free-stream velocity [m/s].
            - The velocity of the airflow relative to the helicopter in forward flight.
        dist_f: Forward flight distance [m].
            - The total distance to be covered during the forward flight phase.
        Vinf_vec: Velocity vector for forward flight [m/s].
            - A vector containing the various speeds at which the helicopter will travel during forward flight.
        N: Number of blades.
            - The number of blades on the rotor.
        r_segn: Normalized vector of radial positions (r/R).
            - The normalized positions of the rotor blades in terms of the rotor radius.
        c: Chord length of the rotor blade [m].
            - The chord length of the rotor blade used for aerodynamic calculations.
        theta: Pitch angle of the blade [rad].
            - The pitch angle of the rotor blade relative to the plane of rotation.
        W: Helicopter weight [N].
            - The weight of the helicopter in newtons.
        height_f: Forward flight altitude [m].
            - The altitude at which the helicopter will perform the forward flight phase.
        R: Rotor radius [m].
            - The radius of the rotor in meters.
        A: Rotor area [m^2].
            - The area swept by the rotor blades.
        f_A: Equivalent flat plate area [m^2].
            - The equivalent flat plate area for drag calculations.
        Clalpha: Lift curve slope [rad^-1].
            - The slope of the lift curve with respect to the angle of attack.
        Cd_mean: Average drag coefficient.
            - The mean drag coefficient for the rotor blades during the flight.
        Omega: Angular velocity of the rotor [rad/s].
            - The angular velocity of the rotor during flight.
        N_rotors: Number of rotors.
            - The total number of rotors on the helicopter.
        F_Prandtl: Prandtl's tip loss correction factor.
            - A factor used to account for the loss of lift at the rotor tips.
        height_d: Descent altitude [m].
            - The altitude at which the helicopter will begin its descent phase.
        height_f2: Second forward flight altitude [m].
            - The altitude at which the helicopter will perform a second forward flight phase.
        dist_f2: Second forward flight distance [m].
            - The distance to be covered during the second forward flight phase.

    Returns:
        climb_data: Dictionary containing data for the climb phase.
            - Contains information on climb altitude, speed, distance, time, power rating, fuel flow, fuel used, and gross weight.
        hover_data: Dictionary containing data for the hover phase.
            - Contains information on hover altitude, time, power rating, fuel flow, fuel used, and gross weight.
        forward_data: Dictionary containing data for the first forward flight phase.
            - Contains information on forward flight altitude, speed, distance, time, power rating, fuel flow, fuel used, and gross weight.
        Descend_data: Dictionary containing data for the descent phase.
            - Contains information on descent altitude, speed, distance, time, fuel flow, fuel used, and gross weight.
        forward_data2: Dictionary containing data for the second forward flight phase.
            - Contains information on second forward flight altitude, speed, distance, time, power rating, fuel flow, fuel used, and gross weight.
        fuel_data: Dictionary containing total fuel consumption and maximum takeoff weight.
            - Contains information on total fuel consumption and MTOW after all flight phases.
        gamma_deg: Climb angle in degrees.
            - The angle of climb during the climb phase.
        gamma_deg_d: Descent angle in degrees.
            - The angle of descent during the descent phase.
    """

    atmosphere = Atmosphere(height_h)               # Initialize atmospheric properties at the hover altitude
    rho_inf = atmosphere.density[0]                 # Air density at hover altitude (kg/m^3)
    f = f_A * (4 * A)                               # Equivalent flat plate area for drag calculations (m^2)
    meters_to_feet = 3.28                           # Conversion factor from meters to feet
    meters_to_nautical_miles = 0.00054              # Conversion factor from meters to nautical miles
    meters_per_second_to_knots = 1.944              # Conversion factor from meters per second to knots


                                                                    
    # Climb

    P_parasite = 0.5 * rho_inf * f * Voo**3                                                   # Parasite power in climb
    vi = np.sqrt(-0.5 * Voo**2 + 0.5 * np.sqrt(Voo**4 + 4 * (W / (2 * rho_inf * A))**2))      # Induced velocity in climb
    P_induced = 1.2 * W * vi                                                                  # Induced power in climb
    P_profile = (1/8) * rho_inf * Helicopter.Cd0 * N * c[0] * R * (Omega*R)**3 * (1 + 4.7 * (Voo / (Omega*R))**2)  # Profile power
    Pn_c = P_parasite + P_induced + P_profile                                                 # Required power during climb

    result_c = Helicopter_properties(Pd, Pn_c, MTOW, Voo)      # Get helicopter properties for climb phase
    properties, climb = result_c                               # Retrieve data from result climb vector
                                                                             
    ROC_ft_min, gamma_deg = climb()                            # Climb function to get Rate of Climb (ft/min) and climb angle (deg)
    ROC = (ROC_ft_min) / (3.28 * 60)                           # Convert Rate of Climb to m/s
    gamma = np.deg2rad(gamma_deg)                              # Convert climb angle to radians
    TAS_c = (ROC**2 + Voo**2)**0.5                             # True Airspeed during climb (m/s)
    power_rate_c = Pn_c / Pd * 100                             # Power rate in climb (%)
    time_c = height_c / TAS_c                                  # Time required for climb (s)
    dist_c = Voo * time_c                                      # Horizontal distance covered during climb (m)
    F_flow_c = SFC * (Pn_c / 1000)                             # Fuel flow during climb (kg/h)
    F_used_c = (time_c / 3600) * F_flow_c                      # Fuel used during climb (kg)
    GW_c = W / 9.81 - F_used_c                                 # Remaining gross weight after climb (kg)
                      
   
    # Hover

    V_infty = 0 
    result_h = BEMT_Hover(Helicopter, V_infty)                 # Retrieve the function for hover data with BEMT       

    (Tc, Pc_i, Pc0, Qc, dTcdr_bar_PR, dPcdr_bar, FM, _) = result_h

    Omega_h = np.sqrt(W / (rho_inf * R**2 * A * Tc))           # Hover angular velocity [rad/s].
    lambda_i = np.sqrt(Tc / 2)                                 # Axial interference factor
    w_h = lambda_i * Omega_h * R                               # Induced velocity in hover [m/s]

    P_ih = Pc_i * rho_inf * Omega_h**3 * R**3 * A              # Power in hovering (W)
    Pn_h = P_ih + Pc0                                          # Required power in hovering (W)
    
    power_rate_h = Pn_h/Pd * 100                               # Power rate in hovering (%)
    F_flow_h = SFC * (Pn_h/1000)                               # Fuel flow in hovering (kg/h)
    F_used_h = (time_h/3600) * F_flow_h                        # Fuel used in hovering (Kg)

    # Forward flight

    W = GW_c * g
    result_f = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height_f, R, A, P_ih, f, Pc0, Clalpha, Cd_mean, Omega_h, N_rotors, F_Prandtl)

    (P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
    Omega_new, alpha_e_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
    T_vec, Y_vec, Q_vec, dT_dr, dQ_dr, psi) = result_f  
    
    Vinf = Vinf_vec[-1]                                      # Assigned speed from speed vector (m/s)
    time_f = dist_f / Vinf_vec[-1]                           # Cruising time (s)
    Pn_f = P[-1]                                             # Required power (W)
    power_rate_f = Pn_f/Pd * 100                             # Power rate for cruise (%)
    F_flow_f = SFC * (Pn_f/1000)                             # Fuel flow for cruise (kg/h)
    F_used_f = (time_f/3600) * F_flow_f                      # Fuel used for cruise (kg)
    GW_f = GW_c - F_used_f                                   # Remaining gross weight after cruise (Kg)

    # Descend

    W = GW_f * g
    Pn_d = 0.5 * Pn_h                                         # Power required during descent, assumed half of hovering power
    result_d = Helicopter_properties(Pd, Pn_d, MTOW, Voo_d)   # Retrieve the function from the helicopter properties
    properties, descend = result_d                            # Retrieve data from result descend vector
    ROD_ft_min, gamma_deg_d = descend()                       # Call the descent function to get Rate of Descent and angle
    ROD = w_h / 2                                             # Rate of Descent assumed as half the induced velocity in hover

    power_rate_d = Pn_d / Pd * 100                            # Power rate during descent (%)
    TAS_d = (ROD**2 + Voo**2)**0.5                            # True Airspeed during descent (m/s)
    time_d = (height_f - height_d) / TAS_d                    # Time required for descent (s)
    dist_d = Voo_d * time_d                                   # Horizontal distance covered during descent (m)
    F_flow_d = SFC * (Pn_d / 1000)                            # Fuel flow during descent (kg/h)
    F_used_d = (time_d / 3600) * F_flow_d                     # Fuel used during descent (kg)
    GW_d = GW_f - F_used_d                                    # Remaining gross weight after descent (kg)
    GW_h = GW_d - F_used_h                                    # Remaining gross weight after the next hover phase (kg)

    # Forward flight 2

    W = GW_h * g
    result_f2 = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height_f2, R, A, P_ih, f, Pc0, Clalpha, Cd_mean, Omega_h, N_rotors, F_Prandtl)

    (P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
    Omega_new, alpha_e_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
    T_vec, Y_vec, Q_vec, dT_dr, dQ_dr, psi) = result_f2 
    
    Vinf = Vinf_vec[-2]                                       # Assigned speed from speed vector (m/s)
    time_f2 = dist_f2 / Vinf_vec[-1]                          # Second Cruising time (t)
    Pn_f2 = P[-1]                                             # Required power (W)
    power_rate_f2 = Pn_f2/Pd * 100                            # Power rate for cruise (%)
    F_flow_f2 = SFC * (Pn_f/1000)                             # Fuel flow for second cruise (kg/h)
    F_used_f2 = (time_f/3600) * F_flow_f2                     # Fuel used for second cruise (kg)
    GW_f2 = GW_h - F_used_f2                                  # Remaining gross weight after second cruise (Kg)

    # Output vectors
 
    climb_data = {
            "Altitude": height_c * meters_to_feet,             
            "Speed ": TAS_c,                                               
            "Distance": dist_c * meters_to_nautical_miles,
            "Time": time_c,
            "Power Rating": "MCP",
            "Fuel Flow": "-",
            "Fuel Used": F_used_c,
            "Gross Weight": GW_c,
        }
    
    forward_data = {
        "Altitude": height_f * meters_to_feet,                             
        "Speed ": Vinf * meters_per_second_to_knots,
        "Distance": dist_f * meters_to_nautical_miles,
        "Time": time_f,
        "Power Rating": power_rate_f,
        "Fuel Flow": F_flow_f,
        "Fuel Used": F_used_f,
        "Gross Weight": GW_f,
    }

    Descend_data = {
        "Altitude": height_d * meters_to_feet,
        "Speed ": TAS_d,
        "Distance": dist_d * meters_to_nautical_miles,
        "Time": time_d,
        "Power Rating": "-",
        "Fuel Flow": F_flow_d,
        "Fuel Used": F_used_d,
        "Gross Weight": GW_d,
    }

    hover_data = {
        "Altitude": height_h * meters_to_feet,
        "Speed ": "-",
        "Distance": "-",
        "Time": time_h,
        "Power Rating": power_rate_h,
        "Fuel Flow": F_flow_h,
        "Fuel Used": F_used_h,
        "Gross Weight": GW_h,
    }

    forward_data2 = {
        "Altitude": height_f2 * meters_to_feet,
        "Speed ": Vinf * meters_per_second_to_knots,
        "Distance": dist_f2 * meters_to_nautical_miles,
        "Time": time_f2,
        "Power Rating": power_rate_f2,
        "Fuel Flow": F_flow_f2,
        "Fuel Used": F_used_f2,
        "Gross Weight": GW_f2,
    }

    fuel_data = {
            "MTOW": MTOW/g,
            "Total fuel Consuption": F_used_c + F_used_f + F_used_d + F_used_h,
        }

    return climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d

