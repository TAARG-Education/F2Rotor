"""
Test case for the utilization of the function evaluate_range. This test case evaluates the range of the AS 332 L1 helicopter
from AEROSPATIALE in its maximum takeoff condition.
    
    
    Authors: Beniamino Ferrara
    Date: 18/06/2024
    Version: 1.01
"""


import numpy as np
from range_fun import evaluate_range
import matplotlib.pyplot as plt

"""
Helicopter parameters. (Change parameteres accordingly to the unit of measure) 
"""
N = 100
V_inf = np.linspace(1,100,N) # cruise velocity[m/s]
W_fuel = 1620 # fuel capacity [Kg]
SFC = 0.316 # specific fuel consumption [kg/kWh]
R = 7.79 # radius [m]
B = 4    # number of blades
C = 0.6  # chord [m]
Cd_bar = 0.0121  #mean drag coefficient 
Rt = 1.525    # Tail rotor radius 
omega = 27.85 # angular velocity [rad/s]
sigma = (B*C)/(np.pi*R) #solidity 
rho = 1.225             #density [kg/m^3]
ltv = R + Rt + 0.5 #Tail rotor distance from center of gravity
A = np.pi*R**2     #rotor Area [m^2]
W = 8600           #take off weight [kg]
f = 0.007*A    #friction area
sigma_t = 0.209  # Tail rotor solidity
omega_t = 204/Rt # Tail rotor angular velocity
A_t = np.pi*Rt**2 # Tail rotor Area

"""
Calculation of the required power by the main rotor. First the induced velocity on the main rotor is evaluated to calculate the induced power.
The parasite power and fuselage power due to drag are calculated next. Finally the 3 contributions are added.
"""
mu = V_inf/(omega*R) #advance rateo
v_i =np.sqrt(-0.5*(V_inf**2) + 0.5*np.sqrt(V_inf**4+4*((W*9.81)/(2*rho*A*V_inf))**2))  #induced velocity
k = 1.15 # Correction factor for the indeced power
P_i = k*v_i*W*9.81       #Induced Power
P_0 = (sigma * Cd_bar)/8 * rho * (omega**3) * (R**3) * A * (1 + 4.7*mu) #Parasite Power
P_fus = f*0.5*rho*V_inf**3         # fuselage Power
Pr =  (P_i + P_0 + P_fus)/1000     # main rotor required Power

"""
Calculation of the required power by the tail rotor. The tail rotor required thrust to balance the main rotor momentum is calculated first,
then the required power is calculated with the same procedure as before.
"""

Q = Pr*1000/omega  #Reaction torque
Tt = Q * ltv       #Tail rotor thrust
v_it =np.sqrt(-0.5*(V_inf**2) + 0.5*np.sqrt(V_inf**4+4*((Tt)/(2*rho*A_t*V_inf))**2))  #induced velocity on tail rotor
Pt = k*v_i*Tt + (sigma_t * Cd_bar)/8 * rho * (omega_t**3) * (Rt**3) * A_t * (1 + 4.7*mu) #Tail rotor Power.
Pt = Pt/1000

"""
Finally the main rotor power and the tail rotor power are added and the "evaluate_range" function is used.
"""
eta = 1.04                             # gearbox factor
Pn = eta*(Pr + Pt)                     # Total required Power
[range, v_max] = evaluate_range(V_inf, W_fuel, SFC, Pn) 
print('the velocity for maximum range is', v_max)
v_plot = V_inf*3.6

"""
Plotting the range curve.
"""

plt.plot(v_plot, range)
plt.xlabel("V_inf(km/h)")
plt.ylabel("Range (km)")
plt.show()
plt.figure()
