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
from ambiance import Atmosphere                    
from scipy.constants import pi, g                   
from Mission_profile import Mission_profile
from Mission_profile import create_table_figure
from Drone_functions import RigidFW2
from Climb import Helicopter_properties
from BEMT_Hover import BEMT_Hover
from Drone_functions import Calculate_alpha_hover

# Baisc Input

MTOW = 40000                                                                             # Maximum Take-Off Weight (N).
W = 40000                                                                                         # Initial Weight (N).
c_r = 0.45                                                                                            # Root chord [m].
c_t = 0.35                                                                                             # Tip chord [m].
r_segn = np.linspace(0.1, 1, 19)                                                      # Dimensionless disk rotor radius.
c_tw = (c_t-c_r)/(r_segn[-1]-r_segn[0])                                                               # Twist chord [m].
Vinf_vec = np.linspace(0, 65, 10)                                                          # Asymptotic air speed [m/s].

# Input for Helicopter class.
class Helicopter:

    Omega_r_mr = 25                                                                       # Main rotor angular velocity.
    R_mr = 5.5                                                                                  # Disk rotor radius [m].
    R_hub = 0.1 * R_mr                                                                           # Hub rotor radius [m].
    c_mr = c_r + c_tw*(r_segn-r_segn[0])                                                                # Chord law [m].
    N_mr = 4                                                                                         # Number of blades.
    theta0 = 25                                                                             # Initial pitch angle [deg].
    theta_end = 15                                                                                # Tip pitch angle [deg].
    theta_tw = (theta_end-theta0)/(r_segn[-1]-r_segn[0])                                              # Pitch law [deg].
    Cla = 2 * pi                                                                            # Lift slope curve [rad^-1].
    Cd0 = 0.012                                                                            # Parasitic drag coefficient.

    # Class attributes.
    def __init__(self, Omega_r_mr, R_mr, R_hub, c_mr, N_mr, theta0, theta_tw, Cla, Cd0):
        Helicopter.Omega_r_mr = Omega_r_mr
        Helicopter.R_mr = R_mr
        Helicopter.R_hub = R_hub
        Helicopter.c_mr = c_mr
        Helicopter.N_mr = N_mr
        Helicopter.theta0 = theta0
        Helicopter.theta_tw = theta_tw
        Helicopter.Cla = Cla
        Helicopter.Cd0 = Cd0

# Input variables

R = Helicopter.R_mr
N = Helicopter.N_mr
c = Helicopter.c_mr
SFC = 0.7                                       # Specific Fuel Consumption (kg/kWh)
Pd = 500000 * 2                                 # Available power (W)
height_c = 3298.78                              # Climb altitude (m)
Voo = 20                                        # Operating speed (m/s)
Voo_d = 50                                      # Speed during descent (m/s)
height_h = 451.10                               # Hover altitude (m)
time_h = 120                                    # Hover time (s)
V_infty = 0                                     # Speed during hovering (m/s)
dist_f = 75200                                  # Forward flight distance (m)
height_f = 3298.78                              # Forward flight altitude (m)
A = math.pi * R**2                              # Rotor area
f_A = 0.00001                                   # Friction coefficient
Clalpha = 5.7                                   # Lift curve slope
Cd_mean = 0.01                                  # Mean drag
Omega = 27                                      # Rotor angular velocity (rad/s)
N_rotors = 1                                    # Number of rotors
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
    Clalpha=Clalpha,
    Cd_mean=Cd_mean,
    Omega=Omega,
    N_rotors=N_rotors,
    F_Prandtl=F_Prandtl,
    height_d=height_d,
    height_f2=height_f2,
    dist_f2=dist_f2
)

create_table_figure(climb_data, hover_data, forward_data, Descend_data, forward_data2, fuel_data, gamma_deg, gamma_deg_d)
