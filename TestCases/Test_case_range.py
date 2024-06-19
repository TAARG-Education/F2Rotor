"""
Test case for the utilization of the function evaluate_range. This test case evaluates the range of the AS 332 L1 helicopter
from AEROSPATIALE in its maximum takeoff condition.
    
    
    Authors: Beniamino Ferrara
    Date: 18/06/2024
    Version: 1.00
"""


import numpy as np
from range_fun import evaluate_range
import matplotlib.pyplot as plt

N = 100
V_inf = np.linspace(1,100,N) # cruise velocity[m/s]
W_fuel = 1620 # fuel capacity [Kg]
SFC = 0.316 # specific fuel consumption [kg/kWh]
R = 7.79 # radius [m]
B = 4    # number of blades
C = 0.6  # chord [m]
Cd_bar = 0.0121  #mean drag coefficient 
Rt = 1.525
omega = 27.85 # angular velocity [rad/s]
sigma = (B*C)/(np.pi*R) #solidity 
rho = 1.225             #density [kg/m^3]
ltv = R + Rt + 0.5 #Tail rotor distance from center of gravity
A = np.pi*R**2     #rotor Area [m^2]
W = 8600           #take off weight [kg]
mu = V_inf/(omega*R) #advance rateo
f = 0.007*A    #friction area
sigma_t = 0.209
omega_t = 204/Rt
A_t = np.pi*Rt**2
v_i =np.sqrt(-0.5*(V_inf**2) + 0.5*np.sqrt(V_inf**4+4*((W*9.81)/(2*rho*A*V_inf))**2))  #induced velocity
k = 1.15
P_i = k*v_i*W*9.81       #Induced Power
P_0 = (sigma * Cd_bar)/8 * rho * (omega**3) * (R**3) * A * (1 + 4.7*mu) #Parasite Power
P_fus = f*0.5*rho*V_inf**3         # fuselage Power

eta = 1.04                             # gearbox factor
Pr =  (P_i + P_0 + P_fus)/1000     # main rotor required Power
Q = Pr*1000/omega
Tt = Q * ltv 
v_it =np.sqrt(-0.5*(V_inf**2) + 0.5*np.sqrt(V_inf**4+4*((Tt)/(2*rho*A_t*V_inf))**2))  #induced velocity on tail rotor
Pt = k*v_i*Tt + (sigma_t * Cd_bar)/8 * rho * (omega_t**3) * (Rt**3) * A_t * (1 + 4.7*mu) 
Pt = Pt/1000
Pn = eta*(Pr + Pt) 
[range, v_max] = evaluate_range(V_inf, W_fuel, SFC, Pn)
print('the velocity for maximum range is', v_max)
v_plot = V_inf*3.6
plt.plot(v_plot, range)
plt.xlabel("V_inf(km/h)")
plt.ylabel("Range (km)")
plt.show()
plt.figure()
