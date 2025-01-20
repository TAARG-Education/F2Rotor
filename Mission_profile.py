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
from matplotlib.table import Table
from scipy.constants import pi, g
from ambiance import Atmosphere 

from Drone_functions import RigidFW2
from Drone_functions import RotorSolving
from Climb import Helicopter_properties
from BEMT_Hover import BEMT_Hover

def Mission_profile(Helicopter, SFC, Pd, height_c, MTOW, Voo, Voo_d, height_h, time_h,\
                     V_infty, dist_f, Vinf_vec, N, r_segn, c, theta, W, height_f, R, A, f_A, Clalpha, Cd_mean,\
                          Omega, N_rotors, F_Prandtl, height_d, height_f2, dist_f2):

    atmosphere = Atmosphere(height_h)                                                                     # Initialize atmospheric properties at the hover altitude
    rho_inf = atmosphere.density[0]                                                                       # Air density at hover altitude (kg/m^3)
    f = f_A * (4 * A)                                                                                     # Equivalent flat plate area for drag calculations (m^2)
                                                                    
    # Climb

    P_parasite = 0.5 * rho_inf * f * Voo**3                                                                         # Parasite power in climb
    vi = np.sqrt(-0.5 * Voo**2 + 0.5 * np.sqrt(Voo**4 + 4 * (W / (2 * rho_inf * A))**2))                            # Induced velocity in climb
    P_induced = 1.2 * W * vi                                                                                        # Induced power in climb
    P_profile = (1/8) * rho_inf * Helicopter.Cd0 * N * c[0] * R * (Omega*R)**3 * (1 + 4.7 * (Voo / (Omega*R))**2)   # Profile power required for rotor drag
    Pn_c = P_parasite + P_induced + P_profile                                                                       # Required power during climb

    result_c = Helicopter_properties(Pd, Pn_c, MTOW, Voo)                                                           # Get helicopter properties for climb phase
    properties, climb = result_c                                                                                    # Retrieve data from result climb vector
                                                                             
    ROC_ft_min, gamma_deg = climb()                                                       # Call the climb function to get Rate of Climb (ft/min) and climb angle (deg)
    ROC = (ROC_ft_min) / (3.28 * 60)                                                      # Convert Rate of Climb to m/s
    gamma = np.deg2rad(gamma_deg)                                                         # Convert climb angle to radians
    TAS_c = (ROC**2 + Voo**2)**0.5                                                        # True Airspeed during climb
    power_rate_c = Pn_c / Pd * 100                                                        # Power rate in climb
    time_c = height_c / TAS_c                                                             # Time required for climb (s)
    dist_c = Voo * time_c                                                                 # Horizontal distance covered during climb (m)
    F_flow_c = SFC * (Pn_c / 1000)                                                        # Fuel flow during climb (kg/s)
    F_used_c = (time_c / 3600) * F_flow_c                                                 # Fuel used during climb (kg)
    GW_c = W / 9.81 - F_used_c                                                            # Remaining gross weight after climb (kg)
                      
   
    # Hover

    V_infty = 0 
    result_h = BEMT_Hover(Helicopter, V_infty)                                             # # Retrieve the function for hover data with BEMT       

    (Tc, Pc_i, Pc0, Qc, dTcdr_bar_PR, dPcdr_bar, FM, _) = result_h

    Omega_h = np.sqrt(W / (rho_inf * R**2 * A * Tc))                                       # Hover angular velocity [rad/s].
    lambda_i = np.sqrt(Tc / 2)                                                             # Axial interference factor
    w_h = lambda_i * Omega_h * R                                                           # Induced velocity in hover [m/s]

    P_ih = Pc_i * rho_inf * Omega_h**3 * R**3 * A                                          # Power in hovering
    Pn_h = P_ih + Pc0                                                                      # Required power in hovering
    
    power_rate_h = Pn_h/Pd * 100                                                           # Power rate in hovering
    F_flow_h = SFC * (Pn_h/1000)                                                           # Fuel flow in hovering
    F_used_h = (time_h/3600) * F_flow_h                                                    # Fuel used in hovering

    # Forward flight

    W = GW_c * g
    result_f = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height_f, R, A, P_ih, f, Pc0, Clalpha, Cd_mean, Omega_h, N_rotors, F_Prandtl)

    (P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
    Omega_new, alpha_e_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
    T_vec, Y_vec, Q_vec, dT_dr, dQ_dr, psi) = result_f  
    
    Vinf = Vinf_vec[-1]                                                                   # Assigned speed from speed vector
    time_f = dist_f / Vinf_vec[-1]                                                        # Cruising time
    Pn_f = P[-1]                                                                          # Required power
    power_rate_f = Pn_f/Pd * 100                                                          # Power rate for cruise
    F_flow_f = SFC * (Pn_f/1000)                                                          # Fuel flow for cruise
    F_used_f = (time_f/3600) * F_flow_f                                                   # Fuel used for cruise
    GW_f = GW_c - F_used_f                                                                # Remaining gross weight after cruise (Kg)

    # Descend

    W = GW_f * g
    Pn_d = 0.5 * Pn_h                                                                     # Power required during descent, assumed to be half of hovering power
    result_d = Helicopter_properties(Pd, Pn_d, MTOW, Voo_d)                               # Retrieve the function from the helicopter properties
    properties, descend = result_d                                                        # Retrieve data from result descend vector
    ROD_ft_min, gamma_deg_d = descend()                                                   # Call the descent function to get Rate of Descent and angle
    ROD = w_h / 2                                                                         # Rate of Descent assumed as half the induced velocity in hover

    power_rate_d = Pn_d / Pd * 100                                                        # Power rate during descent
    TAS_d = (ROD**2 + Voo**2)**0.5                                                        # True Airspeed during descent
    time_d = (height_f - height_d) / TAS_d                                                # Time required for descent (s)
    dist_d = Voo_d * time_d                                                               # Horizontal distance covered during descent (m)
    F_flow_d = SFC * (Pn_d / 1000)                                                        # Fuel flow during descent (kg/s)
    F_used_d = (time_d / 3600) * F_flow_d                                                 # Fuel used during descent (kg)
    GW_d = GW_f - F_used_d                                                                # Remaining gross weight after descent (kg)
    GW_h = GW_d - F_used_h                                                                # Remaining gross weight after the next hover phase (kg)

    # Forward flight 2

    W = GW_h * g
    result_f2 = RigidFW2(Vinf_vec, N, r_segn, c, theta, W, w_h, height_f2, R, A, P_ih, f, Pc0, Clalpha, Cd_mean, Omega_h, N_rotors, F_Prandtl)

    (P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
    Omega_new, alpha_e_vec, dTc, U_P, dH_dpsi_0, dQ_dpsi_0,
    T_vec, Y_vec, Q_vec, dT_dr, dQ_dr, psi) = result_f2 
    
    Vinf = Vinf_vec[-2]                                                                   # Assigned speed from speed vector
    time_f2 = dist_f2 / Vinf_vec[-1]                                                      # Second Cruising time
    Pn_f2 = P[-1]                                                                          # Required power
    power_rate_f2 = Pn_f2/Pd * 100                                                          # Power rate for cruise
    F_flow_f2 = SFC * (Pn_f/1000)                                                          # Fuel flow for second cruise
    F_used_f2 = (time_f/3600) * F_flow_f2                                                   # Fuel used for second cruise
    GW_f2 = GW_h - F_used_f2                                                                # Remaining gross weight after second cruise (Kg)

    # Output vectors
 
    climb_data = {
            "Altitude": "To " + f"{round(height_c * 3.28)}",
            "Speed ": round(TAS_c,1),
            "Distance": round(dist_c * 0.00054,1),
            "Time": round(time_c),
            "Power Rating": "MCP",
            "Fuel Flow": "-",
            "Fuel Used": round(F_used_c,1),
            "Gross Weight": round(GW_c),
        }
    
    forward_data = {
        "Altitude": round(height_f * 3.28),
        "Speed ": round(Vinf * 1.944,1),
        "Distance": round (dist_f * 0.00054,1),
        "Time": round(time_f),
        "Power Rating": f"{round(power_rate_f)}%",
        "Fuel Flow": round(F_flow_f,1),
        "Fuel Used": round(F_used_f,1),
        "Gross Weight": round(GW_f,1),
    }

    Descend_data = {
        "Altitude": "To " + f"{round(height_d * 3.28)}",
        "Speed ": round(TAS_d,1),
        "Distance": round (dist_d * 0.00054,1),
        "Time": round(time_d),
        "Power Rating": "-",
        "Fuel Flow": round(F_flow_d,1),
        "Fuel Used": round(F_used_d,1),
        "Gross Weight": round(GW_d,1),
    }

    hover_data = {
        "Altitude": round(height_h * 3.28),
        "Speed ": "-",
        "Distance": "-",
        "Time": time_h,
        "Power Rating": f"{round(power_rate_h)}%",
        "Fuel Flow": round(F_flow_h,1),
        "Fuel Used": round(F_used_h,1),
        "Gross Weight": round(GW_h,1),
    }

    forward_data2 = {
        "Altitude": round(height_f2 * 3.28),
        "Speed ": round(Vinf * 1.944,1),
        "Distance": round (dist_f2 * 0.00054,1),
        "Time": round(time_f2),
        "Power Rating": f"{round(power_rate_f2)}%",
        "Fuel Flow": round(F_flow_f2,1),
        "Fuel Used": round(F_used_f2,1),
        "Gross Weight": round(GW_f2,1),
    }

    fuel_data = {
            "MTOW": round(MTOW/9.81),
            "Total fuel Consuption": round(F_used_c + F_used_f + F_used_d + F_used_h),
        }

    return climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d

# Graphics

def create_table_figure(climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d):
    
    headers = ["Segment", "Altitude \n(ft)", "Speed \n(knots - TAS)", "Distance \n(nm)", "Time \n(s)", 
               "Power Rating \n(%)", "Fuel Flow \n(kg/h)", "Fuel Used \n(kg)", "Gross Weight \n(kg)"]
    
    rows = [
        ["Warm-up", "-", "-", "-", "-", "-", "-", "-", "-"],
        ["Climb", climb_data["Altitude"], climb_data["Speed "], climb_data["Distance"], climb_data["Time"], 
         climb_data["Power Rating"], climb_data["Fuel Flow"], climb_data["Fuel Used"], climb_data["Gross Weight"]],
        ["Cruise", forward_data["Altitude"], forward_data["Speed "], forward_data["Distance"], forward_data["Time"], 
         forward_data["Power Rating"], forward_data["Fuel Flow"], forward_data["Fuel Used"], forward_data["Gross Weight"]],
        ["Descend", Descend_data["Altitude"], Descend_data["Speed "], Descend_data["Distance"], Descend_data["Time"], 
         Descend_data["Power Rating"], Descend_data["Fuel Flow"], Descend_data["Fuel Used"], Descend_data["Gross Weight"]],
        ["Hover", hover_data["Altitude"], hover_data["Speed "], hover_data["Distance"], hover_data["Time"], 
         hover_data["Power Rating"], hover_data["Fuel Flow"], hover_data["Fuel Used"], hover_data["Gross Weight"]],
        ["Cruise", forward_data2["Altitude"], forward_data2["Speed "], forward_data2["Distance"], forward_data2["Time"], 
         forward_data2["Power Rating"], forward_data2["Fuel Flow"], forward_data2["Fuel Used"], forward_data2["Gross Weight"]],
        ["Final descent", "-", "-", "-", "-", "-", "-", "-", "-"]
    ]
    
    fig, ax = plt.subplots(figsize=(12, 6)) 
    ax.axis('tight')
    ax.axis('off')
    
    # Table creation
    table = Table(ax, bbox=[0, 0.2, 1, 0.7])  
    n_cols = len(headers)
    n_rows = len(rows) + 1  

    header_fontsize = 15
    cell_fontsize = 15
    footer_fontsize = 10

    for col in range(n_cols):
        cell = table.add_cell(0, col, width=1/n_cols, height=0.1, text=headers[col], loc='center', facecolor='lightgrey')
        cell.set_fontsize(header_fontsize)

    for row_idx, row in enumerate(rows):
        for col_idx, cell_value in enumerate(row):
            if row_idx == 1 and col_idx == 1:  
                cell_value_with_note = f"{str(cell_value)}\nγ = {round(gamma_deg)}°"
            elif row_idx == 3 and col_idx == 1:
                cell_value_with_note = f"{str(cell_value)}\nγ = - {round(gamma_deg_d)}°"

            else:
                cell_value_with_note = str(cell_value)
        
            cell = table.add_cell(row_idx + 1, col_idx, width=1/n_cols, height=0.1, text=str(cell_value), loc='center')
            cell = table.add_cell(row_idx + 1, col_idx, width=1/n_cols, height=0.1, text=cell_value_with_note, loc='center')
            cell.set_fontsize(cell_fontsize)

    ax.add_table(table)

    footer_text = (
        "Take-off GW: " + str(fuel_data["MTOW"]) + " kg\n"
        "Total fuel consumption: " + str(fuel_data["Total fuel Consuption"]) + " kg\n"
        "Total flight time: " + str(round((climb_data["Time"] + forward_data["Time"] + hover_data["Time"] + Descend_data["Time"])/60)) + " minutes"
    )
    plt.figtext(0.5, 0.02, footer_text, wrap=True, horizontalalignment='center', fontsize=footer_fontsize)

    plt.show()
