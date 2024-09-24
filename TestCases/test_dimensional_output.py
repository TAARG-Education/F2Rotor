import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from dimensional_output_propeller import dimensional_output_propeller

# ==================================================================================================
# |Name           : test_dimensional_output_propeller.py                                           |
# |Author         : Andrea Raia, Antonino Guida.                                                   |
# |                 University of Naples Federico II.                                              |
# |Version        : 1.0                                                                            |
# |Date           : 19/06/2024                                                                     |
# |Modified       : 19/06/2024                                                                     |
# |Description    : Test case of the functions contained in 'dimensional_output.py'.               | 
# |                                                                                                |     
# |Reference      : R. Tognaccini, (a.a. 2023/2024), "Lezioni di Aerodinamica dell'ala rotante"    |                                                         
# |                                                                                                |
# |Values ​​assigned as input to the functions:                                                      |
# |                 1) Ct; thrust coefficient.		                                               |															                                                                       
# |					2) Cp; power coefficient.	                                                   |								                                                                      
# |                 3) J; advance ratio.                                                           |
# |                 4) altitude; altitudes at which to calculate the dimensional values.           |
# |                 5) D; propeller diameter.                                                      |
# |                 6) RPM; revolutions per minute.                                                |
# |                                                                                                |
# |Output quantities from functions:                                                               |
# |                 1) T; thrust.                                                                  |
# |                 2) P; power.                                                                   |
# |                 3) eta; efficiency.                                                            |
# |Note           : The plots of the calculated quantities are obtained through the appropriate    |
# |                 functions of the 'dimensional_output.py' class.                                |
## =================================================================================================



J = np.array([0.207, 0.262, 0.313, 0.365, 0.433, 0.473, 0.539])             # definition of a vector of advance ratio values
Ct = np.array([0.059, 0.0508, 0.0424,0.0338, 0.0218, 0.015, 0.004473])      # definition of a vector of thrust coefficients 
Cp = np.array([0.035, 0.0326, 0.0299, 0.02678, 0.0218, 0.0185, 0.01404])    # definition of a vector of power coefficients 


D = 0.229                                                                   # definition of a diameter                                 
R= D/2                                                                      # calculation of the radius
RPM = 3014                                                                  # definition of a RPM

# note: even if you want to perform the calculations on a single altitude, you must define it as a vector (e.g. altitude = [1000])

altitudes = [0, 1000, 2000]                                                 # definition of a vector of altitudes [m]
density = Atmosphere(altitudes).density

# note: input 5 to not save the plot in .pdf (it must be 1 if you want to save the figure as pdf)

Output1 = dimensional_output_propeller(Ct,Cp,J,D,altitudes,RPM)             # use of a 'dimensional_output_propeller' class
Output1.Power_plot(5)                                                       # use of the 'Power_plot' function to obtain the graph of power as a function of speed
Output1.Efficiency_plot(5)                                                  # use of the 'Efficiency_plot' function to obtain the graph of efficiency as a function of speed
Output1.Thrust_plot(5)                                                      # use of the 'Thrust_plot' function to obtain the graph of thrust as a function of speed
                                                                            
print(Output1.interp_Power(4))                                              # calculation of power at a speed of 4 m/s using the function 'interp_Power'
print(Output1.interp_Thrust(4))                                             # calculation of thrust at a speed of 4 m/s using the function 'interp_Thrust'
print(Output1.interp_efficiency(4))                                         # calculation of efficiency at a speed of 4 m/s using the function 'interp_efficiency'



