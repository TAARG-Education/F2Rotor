# Organization: Universita' degli Studi di Napoli - Federico II, Dipartimento di Ingegneria Industriale, Ingegneria Aerospaziale
# Course: Aerodinamica dell'ala rotante
# Professor: Renato Tognaccini
# Supervisor: Ettore Saetta
# Academic Year: 2023-2024
# Authors: Alessio Ferrara

# Test case for the computation of the service ceilings of Agusta A109. 

# To compute the theoretical maximum altitude and the practical maximum altitude (or ceilings) of an helicopter in advance (i.e. in forward flight), is used
# the feature "ServiceCeiling.py" that takes in input a set of maximum rates of climb evaluated at different altitudes and the set of altitudes. To compute the 
# rates of climb is used the feature "Climb.py" (that can be found in F2Rotor repository on GitHub (https://github.com/TAARG-Education/F2Rotor). See that specific
# function for more details about the script). Also, the function "Climb.py" will need, among other different values that will be described when used, 
# the definition of the necessary power of the helicopter that is computed using the feature "Helicopter_Power_Advance.py", that can also be found on GitHub 
# at the same previous link (see for more details).

# Note: see also in the same repository: "TestCases/Helicopter_Power_Advance_TestCase.py" and "TestCases/Climb/Example.py" for further information.
 

# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
# Chapter 7: Il rotore rigido in volo traslato, pp.95-103.
# Giovanni Di Giorgio - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - page 89-91.




import numpy as np
import matplotlib.pyplot as plt
from Helicopter_Power_Advance import P_i_tilde_function, P_c_0_function, P_c_Fus_function
from ambiance import Atmosphere
import math
import Climb
from ServiceCeiling import ServiceCeiling_function


'''
Note: the feature "Helicopter_Power_Advance.py" computes non-dimensional forms of power relative to an helicopter in advance (i.e. in forward flight) but,
for the computation of the service ceilings, the total necessary power of the helicopter needs to be dimensional. In this first part of the test case, 
is performed the trasformation of the single power coefficients from non-dimensional form to dimensional values, where also the single power contributions
are computed for different altitudes values. Than is performed the sum of all the power dimensional contributions.

'''

# Computation of density values

max_alt = 8000                                                         # Max altitude at which the power will be computed (m)
H = 500                                                                # Number of intervals in which altitude array is divided (Keep this value high: 500 or higher)
h_values = np.linspace(0, max_alt, H)                                  # Definition of uniform altitudes array from 0 (m) to max_alt

atmosphere = Atmosphere(h_values)                                      # Computation of density values at the altitude values defined in h_values
rho_values = atmosphere.density



# Definition of helicopter characteristics

MTOW = 2600                                                            # Maximum take off weight (kg)
W = MTOW*9.81                                                          # Maximum weight (N). Will be assumumed to be equal to the Thrust because  
                                                                       # in this test case it is considered a flight condition of uniform level flight

R = 5.5                                                                # Radius of the helicopter rotor (m)
A = math.pi*R**2                                                       # Area of the disk rotor
sigma = 0.0775                                                         # Solidity of disk rotor
omega = 40.29                                                          # Angular rotor speed (1/sec)
Cd_mean = 0.012                                                         # Assumed value of drag coefficient
f_over_A = 0.009                                                       # f_over_A value


# Computation of dimensional induced power of a stiff rotor in advance.

'''
This script computes the dimensional induced power of a stiff rotor in advance, in the hypotesis of a constant thrust. The script calls the 
function P_i_tilde_function to compute the induced power in a non-dimensional form (P_i_tilde), with respect to the induced power of a stiff rotor 
in hovering. Than, the P_i_tilde function is multiplied by induced power of a stiff rotor in hovering P_hover = (W^3/(2*rho*A))^(1/2) to obtain the
induced power (P_i) in dimensional form and at different altitudes.
(Note: see "Helicopter_Power_Advance.py" for further information about the form and theory beyond the equation used)

Input:
    - V_inf_tilde: non-dimensional values of speeds defined as V_inf/w_h, where V_inf is an array of asymptotic velocity and w_h is induced velocity on 
    the disk rotor in hovering. w_h is defined as w_h = (W/(2*rho*A))^(1/2).
    - alpha_deg: Angle of Attack (AoA) of disk rotor with respect to asymptotic velocity in degrees.

'''


alpha_deg_value = 0                                                     # Assumed equal to 0. It is considered a flight condition of  level flight (W=T).
alpha_rad = math.radians(alpha_deg_value)                               # Conversion of AoA from degrees to radians

N = 300                                                                 # Number of points for V_inf interval
V_inf = np.linspace(0, 90, N)                                           # Uniform interval between 0 and 90 (m/s)

V_inf_tilde_values = np.zeros((N,H))                                    # Inizialization of the V_inf_tilde values matrix 

for j in range(H):                                                      # Loop for V_inf_tilde_values matrix computation
  
    rho = rho_values[j]                                                 # Density value at j position in the density array

    w_h = (W/(2*rho*A))**(1/2)                                          # Computation of w_h 

    for i in range(N):  
         
        V_inf_tilde_values[i,j] = V_inf[i]/w_h                          # Computation of V_inf_tilde_values



P_i_tilde = np.zeros((N,H))                                             # Inizialization of the P_i_tilde matrix

for j in range(H):                                                      # Computation loop of P_i_tilde at different altitudes  
  
    for i in range(N):                                                  

        V_inf_tilde = V_inf_tilde_values[i,j]                           # Computation of V_inf_tilde

        P_i_tilde[i,j] = P_i_tilde_function(V_inf_tilde,alpha_deg_value)          # Call to function P_i_tilde_function



# Computation of dimensional induced power

P_i = np.zeros((N,H))                                                   # Inizialization of the P_i matrix

for j in range(H):                                                      # Computation loop of P_i at different altitudes  

    rho = rho_values[j]                                                 # Density value at j position in the density array

    P_hover = ((W)**(3)/(2*rho*A))**(1/2)                               # Computation of P_hover 
  
    for i in range(N):                                                  

        P_i[i,j] = P_i_tilde[i,j]*P_hover                               # Computation of the P_i matrix



# Computation of dimensional parasite power absorbed by a stiff rotor in advance.

'''
This script computes the dimensional parasite power absorbed by a stiff rotor in advance. The hypotesis of rectangular form for the blade and 
drag coefficient of a single blade element constant along the radius and equal to a mean value, are defined.
It is called the function P_c_0_function to compute the parasite power coefficient. Than the coefficient computed is multiplied by rho*A*(omega*R)^3 
to obtain the parasite power in a dimensional form (P_o) and also at different altitudes.
(Note: see "Helicopter_Power_Advance.py" for further information about the form and theory beyond the equation used)

Input:
    - sigma: solidity of disk rotor, defined as N*c/(pi*R), where N is the number of blades, c is the mean chord of the elements and R is the
    radius of disk rotor
    - Cd_mean: mean drag coefficient of the elements of blades
    - K: coefficient which takes into account the effect of the flux's velocity along the direction of radius and other approximations. Without
    this contribution, the theory gives K = 3. Taking into account these effects, one can select a number included in the range [4,5]. 
    Stepniewski and Keys (1984) suggested K = 4.7. This is the default value in this function.
    - mu: advance ratio, defined as V_inf*cos(alpha)/(Omega*R), where V_inf is asymptotic velocity, alpha is Angle of Attack (AoA) of disk rotor
    with respect to the asymptotic velocity, Omega and R are, respectively, the angular velocity and radius of disk rotor.

'''


mu = (V_inf*math.cos(alpha_rad))/(omega*R)                               # Uniform interval for mu

P_c_0 = np.zeros(N)                                                      # Definition of the parasite power array

P_c_0 = P_c_0_function(sigma, Cd_mean, mu)                               # Call to function P_c_0_function

P_o = np.zeros((N,H))                                                    # Definition of the parasite power matrix

for j in range(H):    

    rho = rho_values[j]                                                  # Density value at j position in the density array

    for i in range(N):  

        P_o[i,j] = P_c_0[i]*rho*A*(omega*R)**(3)                         # Computation of the P_0 matrix



# Computation of dimensional parasite power absorbed by a fuselage of an helicopter in advance.

'''
This script computes the dimensional parasite power absorbed by the fuselage of an helicopter in advance. This is the dimensional form of the work of the
aerodynamic drag on the fuselage (the entire helicopter except the main rotor) in the time's unity.
It is called the function P_c_Fus_function to compute the non-dimensional parasite power coefficient. Than the coefficient computed is multiplied by
rho*A*(omega*R)^3 to obtain the dimensional parasite power (P_f) at different altitudes.
(Note: see "Helicopter_Power_Advance.py" for further information about the form and theory beyond the equation used)

Input:
    - mu: defined as V_inf/(Omega*R), where V_inf is asymptotic velocity, Omega and R are, respectively, the angular velocity and the radius
    of disk rotor
    - f/A: where f is the equivalent wetted area and A is the disk area of rotor. An optimistic value of f/A is 0.007. The default value of
    this function is 0.009, as suggested by Ing. Di Giorgio from Leonardo.

'''


P_c_Fus = np.zeros(N)                                                     # Definition of the fuselage parasite power coefficient array
                                                     
P_c_Fus = P_c_Fus_function(mu,f_over_A)                                   # Call to function P_c_Fus_function

P_f = np.zeros((N,H))                                                     # Definition of the parasite power matrix

for j in range(H):    

    rho = rho_values[j]                                                   # Density value at j position in the density array

    for i in range(N):  

        P_f[i,j] = P_c_Fus[i]*rho*A*(omega*R)**(3)                        # Computation of the P_f matrix




# Computation of the total necessary power of the helicopter

'''
It is here computed the sum of the necessary power contributions of an helicopter in advance with a condition of level flight.
For a better evaluation of the total necessary power and so a better extimation of the service ceilings, it is needed to take into account also a value 
of power associated to the auxiliary systems on board of the helicopter. This value is assumed as a constant equal to 20000 W.
The necessary power is defined as: P_n = P_i + P_0 + P_f and it is computed for different altitudes.

'''


P_n = np.zeros((N,H))                                                     # Definition of the necessary power matrix

P_aux = 20000                                                             # Auxiliary power systems value (W)

for j in range(H):    

    for i in range(N):  

        P_n[i,j] = P_i[i,j]+P_o[i,j]+P_f[i,j]+P_aux                       # Computation of the P_n matrix




# Computation of max rate of climb

'''
The following code computes maximum rates of climb at different altitudes. These values are the ones necessary for the computation of the theoretical 
and practical altitude. Because of the reduction of the air density with the increment of the altitude, it is expected an increment of necessary power
for a level flight at slow speeds and a reduction of available power (it is assumed a variation equal to the density variation with altitude).
For this reason there will be a reduction of the maximum rates of climb. The maximum rate of climb is defined as the maximum difference between the available 
power curve and the necessary power curve, divided by the weight of the helicopter. 
(Note: see "Climb.py" for further information about the form and theory beyond the equation used)

'''


P_n_h = np.zeros(N)                                                       # Definition of the total necessary power array at each altitude
P_n_min = np.zeros(H)                                                     # Definition of the array of the minimum values at each altitude


for j in range(H):                                                        # Computation of the minimum value on the necessary power curve for each altitude

    for i in range(N):  

       P_n_h[i] = P_n[i,j]                                                # Computation of the total necessary power array at each different altitude

    index_min_h = np.argmin(P_n_h)                                        # Computation of the index of the minimum value of necessary power
    
    P_n_min[j] = P_n_h[index_min_h]                                       # Computation of the array of the minimum values at each altitude
   
    

# Computation of maximum rates of climb

'''
They are now computed the maximum rates of climb using the feature "Climb.py". 

Input:
    - Pd: available power. This is the total power output provided by the helicopter's engines. It is assumed a variation equal to density variation with
    altitude.
    - P_n_min: array of the necessary power minimum values at each altitude. 
    - MTOW: maximum takeoff weight (kg). This is the maximum weight at which the helicopter is certified to take off.
    - V_inf: horizontal flight speed. This is the speed at which the helicopter is moving horizontally.

'''


Pd_sl = 544270                                                            # Available power sea level. (W)
ROC_max = np.zeros(H)                                                     # Definition of the array of the maximum rates of climb at each altitude

for j in range(H):                                                        # Computation of the maximum rate of climb at each altitude

    Pn = P_n_min[j]

    Pd = Pd_sl*rho_values[j]/rho_values[0]                                # Computation of available power at a specific altitude

    for i in range(N):

        properties, climb = Climb.Helicopter_properties(Pd, Pn, MTOW, V_inf[i])
    
        ROC, _ = climb()                                                  # Call to function climb. The ROC is expressed in ft/min. 

    ROC_max[j] = ROC*0.00508                                              # Computation of the array of the maximum rates of climb at each altitude in m/s



# Computation of theoretical and pratical altitude

'''
This final script computes the theoretical and practical altitude. It is called the function "ServiceCeiling.py".
(Note: see "ServiceCeiling.py" for further information about the form and theory beyond the equation used)

Input:
    - h_values: uniform altitudes array from 0 (m) to max_alt
    - ROC_max: the array of the maximum rates of climb at each altitude in m/s

'''
 

theo_alt, prac_alt = ServiceCeiling_function(ROC_max, h_values)


print('Theoretical altitude: ', theo_alt, 'm' )
print('-------------------')
print('Practical altitude: ', prac_alt, 'm' )
print('-------------------')





