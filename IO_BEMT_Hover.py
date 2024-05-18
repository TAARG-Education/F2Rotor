from BEMT_hover import *
from Helicopter import Helicopter
from ambiance import Atmosphere

Altitude = float(input("Enter the helicopter altitude (m): "))
V_infty = float(input("Enter the fligth speed (m/s): "))
MTOW = float(input("Enter the Maximum Take-Off Weigth of the helicopter (kg): "))
Omega_R_mr = float(input("Enter the tip speed (m/s): "))
R_mr = float(input("Enter the main rotor radius (m): "))
R_hub = float(input("Enter the rotor hub radius (m): "))
c_mr = float(input("Enter the chord of the main rotor's blade (m): "))
N_mr = int(input("Enter the number of blades: "))
theta0 = float(input("Enter the pitch angle at the root of the rotor (°): "))
theta_tw = float(input("Enter the blade twist (tip minus root incidence) (°): "))
Cla = float(input("Enter the lift curve slope of the adopted airfoil (1/rad): "))
Cd0 = float(input("Enter the profile drag coefficient of the adopted airfoil: "))

# ATMOSPHERE PROPERTIES

atmosphere = Atmosphere(Altitude)
rho = atmosphere.density

Helicopter = Helicopter(MTOW, Omega_R_mr, R_mr, R_hub, c_mr, N_mr, theta0, theta_tw, Cla, Cd0)

[Tc, T, Qc, Q, Tc_PR, T_PR, Qc_PR, Q_PR, P, P_PR, Pc_i, P_i, Pc_0, P_0, FM, FM_PR, Disk_Loading, \
 Power_Loading, lambda_i, alpha, dTcdr_bar, dQcdr_bar, dTcdr_bar_PR, dQcdr_bar_PR, r_bar] = BEMT(Helicopter, rho, V_infty)

print("The Thrust coefficient Tc is: " + str(Tc))
print("The Thrust coefficient (Prandtl correction function) Tc_PR is:" + str(Tc_PR))
print("The Torque (Power) coefficient Qc is: " + str(Qc))
print("The Torque (Power) coefficient (Prandtl correction function) Qc_PR is: " + str(Qc_PR))
print("The induced Power (Torque) coefficient Pc_i is: " + str(Pc_i))
print("The parasite Power (Torque) coefficient Pc_0 is: " + str(Pc_0))
print("The Thrust is: " + str(T) + " N.")
print("The Thrust (Prandtl correction function) is: " + str(T_PR) + " N.")
print("The Torque is: " + str(Q) + " Nm.")
print("The Torque (Prandtl correction function) is: " + str(Q_PR) + " Nm.")
print("The Power is: " + str(P/1000) + " KW.")
print("The Power (Prandtl correction function) is: " + str(P_PR/1000) + " KW.")
print("The Induced Power is: " + str(P_i/1000) + " KW.")
print("The Parasite Power is: " + str(P_0/1000) + " KW.")
print("The Figure of Merit FM is: " + str(FM))
print("The Figure of Merit FM (Prandtl correction function) is: " + str(FM_PR))
print("The Disk Loading is: " + str(Disk_Loading) + " N/m^2.")
print("The Power Loading is: " + str(Power_Loading) + " N/W.")
