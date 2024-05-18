import json
import numpy as np
import math as ma
from scipy.constants import g
from numpy import pi
from sympy import symbols, solve

def BEMT(self, rho, V_infty):

    """BEMT is a method that executes the blade element momentum theory. It accepts as input some data of the
    helicopter's rotor, specifically:

    Altitude: the altitude at which the helicopter is hovering (m)
    V_infty: the helicopter horizontal fligth speed (m/s)
    MTOW: Maximum Take-Off Weight of the helicopter (kg);
    R_mr: main rotor radius (m);
    R_hub: rotor hub radius(m);
    c_mr: chord of the main rotor's blade (m);
    N_mr: number of blades;
    theta0: pitch angle at the root of the rotor (°);
    theta_tw: blade twist (tip minus root incidence) (°);
    Cla: lift curve slope of the adopted airfoil (1/rad);
    Cd0: profile drag coefficient of the adopted airfoil.

    It gives as output:

    Tc: Thrust coefficient;
    T: Thrust (N);
    Qc: Torque (Power) coefficient;
    Q: Torque (Nm);
    Tc_PR: Thrust coefficient (considering the Prandtl correction function);
    T_PR: Thrust (considering the Prandtl correction function) (N);
    Qc_PR: Torque (Power) coefficient (considering the Prandtl correction function);
    Q_PR: Torque (considering the Prandtl correction function) (Nm);
    P: Power  (W);
    P_PR: Power (considering the Prandtl correction function) (W);
    Pc_i: Induced Power coefficient;
    P_i: Induced Power (W);
    Pc_0: Parasite Power coefficient;
    P_0: Parasite Power (W);
    FM: Figure of Merit;
    FM_PR: Figure of Merit (considering the Prandtl correction function);
    Disk_Loading: Disk Loading T/A (N/m^2);
    Power_Loading: Power Loading T/P (N/W);
    lambda_i: Induced inflow ratio distribution along the main rotor radius;
    alpha: angle of attack distribution along the main rotor radius (rad);
    dTcdr_bar: Thrust coefficient distribution along the adimensionalized main rotor radius;
    dQcdr_bar: Torque (Power) coefficient distribution along the adimensionalized main rotor radius;
    dTcdr_bar_PR: Thrust coefficient distribution along the adimensionalized main rotor radius (considering the Prandtl
                  correction function);
    dQcdr_bar_PR: Torque (Power) coefficient distribution along the adimensionalized main rotor radius (considering the
                  Prandtl correction function);
    r_bar: dimensional stations along the main rotor radius (m).
    """

    ## INPUT

    W = self.MTOW * g  # MTOW conversion from (kg) to (N)
    Omega_R_mr = self.Omega_R_mr
    R_mr = self.R_mr
    R_hub = self.R_hub
    c_mr = self.c_mr
    N_mr = self.N_mr
    theta0 = self.theta0
    theta_tw = self.theta_tw
    Cla = self.Cla
    Cd0 = self.Cd0
    mu = V_infty/Omega_R_mr  # Rotor advance ratio


    ## GEOMETRIC PROPERTIES

    A = np.pi*R_mr**2  # Disk area (m^2)
    sol = (N_mr*c_mr)/(R_mr*ma.pi) # Rotor solidity
    r = np.linspace(R_hub, R_mr, num=19)  # Evenly spaced stations along the main rotor
    r_bar = r/R_mr  # Stations along the rotor adimensionalized with respect to the rotor radius R


    theta0 = np.deg2rad(theta0)  # Pitch angle at the root of the rotor expressed in rad
    theta_tw = np.deg2rad(theta_tw)  # Blade twist (tip minus root incidence) expressed in rad
    theta = theta0 + theta_tw*r_bar  # Pitch angle distribution along the rotor (°) (supposed linear)


    ## PHYSICAL  QUANTITY INIZIALIZATION

    F_prandtl = np.zeros(len(r_bar))  # Prandtl function value along the main rotor
    lambda_i = np.zeros(len(r_bar))  # Rotor induced inflow ratio along the main rotor
    phi = np.zeros(len(r_bar))  # Angle of inflow along the main rotor  (rad)
    alpha = np.zeros(len(r_bar))  # Angle of attach along the main rotor (rad)
    Cl = np.zeros(len(r_bar))  # Lift coefficient along the main rotor
    Cd = np.zeros(len(r_bar))  # Drag coefficient along the main rotor
    dTcdr_bar = np.zeros(len(r_bar))  # Incremental Thrust coefficient along the main rotor
    dTcdr_bar_PR = np.zeros(len(r_bar))  # Incremental Thrust coefficient along the main rotor
                                         # considering the Prandtl correction function
    dQcdr_bar = np.zeros(len(r_bar))  # Incremental Torque coefficient along the main rotor
    dQcdr_bar_PR = np.zeros(len(r_bar))  # Incremental Torque coefficient along the main rotor
                                         # considering the Prandtl correction function
    dPc_i = np.zeros(len(r_bar))  # Incremental induced Power coefficient of the rotor element
    dPc_0 = np.zeros(len(r_bar))  # Incremental parasite Power coefficient of the rotor element


## FOR CYCLE ALONG THE ROTOR

    for j in range(len(r_bar)):

     # The second degree equation for the induced inflow ratio has coefficient 1, C1 and C2

        C1 = mu + Cla * sol / 8  # First degree coefficient
        C2 = -r_bar[j] * Cla * sol / 8 * (theta[j] - mu / r_bar[j])  # Free term

        lambda_i[j] = (-C1 + np.sqrt(C1**2 - 4 * C2)) / 2  # Induced inflow ratio

        phi[j] = np.arctan(mu / r_bar[j] + lambda_i[j] / r_bar[j])  # Inflow angle (rad)
        alpha[j] = theta[j] - phi[j]  # Angle of attack (rad)

        F_prandtl[j] = 2 / np.pi * np.arccos(np.exp(N_mr / (2 * lambda_i[j]) * (r_bar[j] - 1)))  # Prandtl correction
                                                                                                 # function
        if r_bar[j] == 1:
            F_prandtl[j] = 0

        Cl[j] = Cla * alpha[j]  # Lift coefficient
        Cd[j] = Cd0  # Drag coefficient

        dTcdr_bar[j] = 0.5 * sol * Cl[j] * r_bar[j]**2  # Incremental Thrust coefficient along the main rotor

        dTcdr_bar_PR[j] = dTcdr_bar[j] * F_prandtl[j]   # Incremental Thrust coefficient along the main rotor
                                                        # considering the Prandtl correction function


        dQcdr_bar[j] = 0.5 * sol * Cd[j] * r_bar[j]**3 + lambda_i[j] * dTcdr_bar[j]  # Incremental Torque coefficient
                                                                                     # along the main rotor

        dQcdr_bar_PR[j] = 0.5 * sol * Cd[j] * r_bar[j]**3 + lambda_i[j] * dTcdr_bar_PR[j]  # Incremental Torque
                                                                                           # coefficient along the main
                                                                                           # rotor considering the
                                                                                           # Prandtl correction function

        dPc_i[j] = 0.5 * sol * Cl[j] * phi[j] * r_bar[j]**3  # Incremental induced Power coefficient of the rotor
                                                             # element
        dPc_0[j] = 0.5 * sol * Cd[j] * r_bar[j]**3  # Incremental parasite Power coefficient of the rotor element


# THRUST, TORQUE, POWER COEFFICIENTS AND FIGURE OF MERIT

    Tc = np.round(np.trapz(dTcdr_bar, r_bar), 4) # Thrust coefficient
    Tc_PR = np.round(np.trapz(dTcdr_bar_PR, r_bar),4)  # Thrust coefficient  considering the Prandtl correction
                                                               # function

    Qc = np.round(np.trapz(dQcdr_bar, r_bar),8)  # Torque (Power) coefficient
    Qc_PR = np.round(np.trapz(dQcdr_bar_PR, r_bar),8)  # Torque (Power) coefficient considering the Prandtl
                                                               # correction function

    Pc_i = np.round(np.trapz(dPc_i, r_bar),8) # Induced Power (Torque) coefficient
    Pc_0 = np.round(np.trapz(dPc_0, r_bar),8) # Parasite Power (Torque) coefficient

    FM = ((Tc**(3/2) / np.sqrt(2)))/Qc # Figure of Merit
    FM_PR = ((Tc_PR ** (3 / 2) / np.sqrt(2))) / Qc_PR  # Figure of Merit considering the Prandtl correction function

## THRUST, TORQUE, POWER

    T = Tc*rho*(Omega_R_mr**2)*A # Thrust (N)
    T_PR = Tc_PR*rho*(Omega_R_mr**2)*A # Thrust considering the Prandtl correction function (N)

    Q = Qc*rho*(Omega_R_mr**2)*A*R_mr # Torque (N m)
    Q_PR = Qc_PR*rho*(Omega_R_mr**2)*A*R_mr # Toque considering the Prandtl correction function (N m)

    P = Qc*rho*(Omega_R_mr**3)*A  # Power (kW)
    P_PR = Qc_PR*rho*(Omega_R_mr**3)*A  # Power considering the Prandtl correction function (kW)
    P_i = Pc_i*rho*(Omega_R_mr**3)*A  # Induced Power (W)
    P_0 = Pc_0*rho*(Omega_R_mr**3)*A  # Parasite Power (W)

   ## DISK LOADING AND POWER LOADING

    Disk_Loading = T/A # Disk Loading (N/m^2)
    Power_Loading = T/P # Power Loading (N/W)

    return Tc, T, Qc, Q, Tc_PR, T_PR, Qc_PR, Q_PR, P, P_PR, Pc_i, P_i, Pc_0, P_0, FM, FM_PR, Disk_Loading,\
        Power_Loading, lambda_i, alpha, dTcdr_bar, dQcdr_bar, dTcdr_bar_PR, dQcdr_bar_PR, r_bar

