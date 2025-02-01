# MISSION PROFILE TESTCASE
#
# Description: This is the test case function in order obtain an example of a mission profile for a light twin engine
#              helicopter following climb, cruise, descent and hover for a fixed height and an assigned initial weight
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

import math
import numpy as np
import matplotlib.pyplot as plt                     
from matplotlib import cm
from matplotlib.table import Table                  
from ambiance import Atmosphere                    
from scipy.constants import pi, g
from Helicopter_Class import Helicopter
from Mission_profile import Mission_profile
from Drone_functions import RigidFW2
from Climb import Helicopter_properties
from BEMT_Hover import BEMT_Hover
from Drone_functions import Calculate_alpha_hover

# Baisc Input

MTOW = 40000                                                    # Maximum Take-Off Weight (N).
W = 40000                                                       # Initial Weight (N).
c_r = 0.45                                                      # Root chord [m].
c_t = 0.35                                                      # Tip chord [m].
r_segn = np.linspace(0.1, 1, 19)                                # Dimensionless disk rotor radius.
c_tw = (c_t-c_r)/(r_segn[-1]-r_segn[0])                         # Twist chord [m].
Vinf_vec = np.linspace(0, 65, 10)                               # Asymptotic air speed [m/s].

Omega = 27                                                      # Main rotor angular velocity.
R = 5.5                                                         # Disk rotor radius [m].
R_hub = 0.1 * R                                                 # Hub rotor radius [m].
c = c_r + c_tw*(r_segn-r_segn[0])                               # Chord law [m].
N = 4                                                           # Number of blades.
theta0 = 25                                                     # Initial pitch angle [deg].
theta_end = 15                                                  # Tip pitch angle [deg].
theta_tw = (theta_end-theta0)/(r_segn[-1]-r_segn[0])            # Pitch law [deg].
Cla = 2 * pi                                                    # Lift slope curve [rad^-1].
Cd0 = 0.012                                                     # Parasitic drag coefficient.
Omega_r_tr = 80                                                 # Tail rotor tip speed (m/s)
R_tr = 1.5                                                      # Tail rotor radius (m)
R_hub_tr = 0.3                                                  # Tail rotor hub radius (m)
c_tr = 0.25                                                     # Tail rotor chord length (m)
N_tr = 2                                                        # Number of blades of the tail rotor
l_tr = 3                                                        # Distance between the main rotor and tail rotor (m)

Helicopter = Helicopter(Omega, R, R_hub, c, N, theta0, theta_tw, Cla, Cd0, 
                        Omega_r_tr, R_tr, R_hub_tr, c_tr, N_tr, l_tr)

# Input variables
SFC = 0.7                                       # Specific Fuel Consumption (kg/kWh)
Pd = 500000 * 2                                 # Available power (W)
height_c = 3298.78                              # Climb altitude (m)
Voo = 20                                        # Speed during climb (m/s)
Voo_d = 50                                      # Speed during descent (m/s)
height_h = 451.10                               # Hover altitude (m)
time_h = 120                                    # Hover time (s)
V_infty = 0                                     # Speed during hovering (m/s)
dist_f = 75200                                  # Forward flight distance (m)
height_f = 3298.78                              # Forward flight altitude (m)
A = math.pi * R**2                              # Rotor area
f_A = 0.00001                                   # Friction coefficient
F_Prandtl = np.full(19, 0.98)                   # Prandtl factor
height_d = 451.1                                # Final altitude (m)
height_f2 = 451.1                               # Second Forward flight altitude (m)
dist_f2 = 70000                                 # Second Forward flight distance (m)

theta = Helicopter.theta0 + Helicopter.theta_tw*(r_segn-r_segn[0])                          # Pitch law [deg].
theta = np.deg2rad(theta)     

# Mission_profile function and table printing
climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d = Mission_profile(
    Helicopter=Helicopter,
    SFC=SFC,
    Pd=Pd,
    height_c=height_c,
    MTOW=MTOW,
    Voo=Voo,
    Voo_d=Voo_d,
    height_h=height_h,
    time_h=time_h,
    V_infty=V_infty,
    dist_f=dist_f,
    Vinf_vec=Vinf_vec,
    N=N,
    r_segn=r_segn,
    c=c,
    theta=theta,
    W=W,
    height_f=height_f,
    R=R,
    A=A,
    f_A=f_A,
    Clalpha=Cla,
    Cd_mean=Cd0,
    Omega=Omega,
    N_rotors=N,
    F_Prandtl=F_Prandtl,
    height_d=height_d,
    height_f2=height_f2,
    dist_f2=dist_f2
)

# Graphics

def create_table_figure(climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d):
    
    headers = ["Segment", "Altitude \n(ft)", "Speed \n(knots - TAS)", "Distance \n(nm)", "Time \n(s)", 
               "Power Rating \n(%)", "Fuel Flow \n(kg/h)", "Fuel Used \n(kg)", "Gross Weight \n(kg)"]

    rows = [
        ["Warm-up", "-", "-", "-", "-", "-", "-", "-", "-"],
        ["Climb", f"To {round(float(climb_data['Altitude']))}", climb_data["Speed "], climb_data["Distance"], climb_data["Time"], 
         climb_data["Power Rating"], climb_data["Fuel Flow"], climb_data["Fuel Used"], climb_data["Gross Weight"]],
        ["Cruise", forward_data['Altitude'], forward_data["Speed "], forward_data["Distance"], forward_data["Time"], 
         f"{round(float(forward_data["Power Rating"]))}%", forward_data["Fuel Flow"], forward_data["Fuel Used"], forward_data["Gross Weight"]],
        ["Descend",f"To {round(float(Descend_data['Altitude']))}", Descend_data["Speed "], Descend_data["Distance"], Descend_data["Time"], 
         Descend_data["Power Rating"], Descend_data["Fuel Flow"], Descend_data["Fuel Used"], Descend_data["Gross Weight"]],
        ["Hover", hover_data["Altitude"], hover_data["Speed "], hover_data["Distance"], hover_data["Time"], 
         f"{round(float(hover_data["Power Rating"]))}%", hover_data["Fuel Flow"], hover_data["Fuel Used"], hover_data["Gross Weight"]],
        ["Cruise", forward_data2["Altitude"], forward_data2["Speed "], forward_data2["Distance"], forward_data2["Time"], 
         f"{round(float(forward_data2["Power Rating"]))}%", forward_data2["Fuel Flow"], forward_data2["Fuel Used"], forward_data2["Gross Weight"]],
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
            
            try:
                cell_value = round(float(cell_value), 1)
            except (ValueError, TypeError):
                pass

            if col_idx == 1 and row_idx in [2, 4, 5]:
                try:
                    cell_value = int(float(cell_value))  
                except (ValueError, TypeError):
                    pass 

            if row_idx == 1 and col_idx == 1:
                cell_value_with_note = f"{cell_value}\nγ = {round(gamma_deg)}°"
            elif row_idx == 3 and col_idx == 1:
                cell_value_with_note = f"{cell_value}\nγ = - {round(gamma_deg_d)}°"
            else:
                cell_value_with_note = str(cell_value)

            cell = table.add_cell(row_idx + 1, col_idx, width=1/n_cols, height=0.1, text=cell_value_with_note, loc='center')
            cell.set_fontsize(cell_fontsize)


        ax.add_table(table)

    footer_text = (
        "Take-off GW: " + str(round(fuel_data["MTOW"])) + " kg\n"
        "Total fuel consumption: " + str(round(fuel_data["Total fuel Consuption"])) + " kg\n"
        "Total flight time: " + str(round((climb_data["Time"] + forward_data["Time"] + hover_data["Time"] + Descend_data["Time"])/60)) + " minutes"
    )
    plt.figtext(0.5, 0.02, footer_text, wrap=True, horizontalalignment='center', fontsize=footer_fontsize)

    plt.show()

create_table_figure(climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d)
