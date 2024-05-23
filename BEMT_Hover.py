import numpy as np

def BEMT_Hover(Helicopter, V_infty):

    """BEMT_Hover is a method that executes the Blade Element Momentum Theory with respect for an helicopter in hover or
        (or slow climbing) condition. It accepts as input some data of the helicopter's rotor and flight conditions,
        specifically:

    V_infty: the helicopter horizontal fligth speed (m/s)
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
    Pc_i: Induced Power coefficient;
    Pc_0: Parasite Power coefficient;
    Qc: Torque (Power) coefficient;
    dTcdr_bar: Incremental Thrust coefficient along the main rotor radius considering the Prandtl correction function;
    dPcdr_bar: Incremental Power coefficient along the main rotor radius;
    FM: Figure of Merit;
    F_Prandtl: Prandtl correction function.


    The theory behind this code can be found in: Renato Tognaccini - "Lezioni di aerodinamica dell'ala rotante -
    Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04 -
    Chapter 6 - Paragraphs 6.1 and 6.2 - pages 71 to 74

    Author: Francesco Gervasio
    Date: 23/05/2024
    Version: 1.00
    """

    ## INPUT

    Omega_R_mr = Helicopter.Omega_r_mr
    R_mr = Helicopter.R_mr
    R_hub = Helicopter.R_hub
    c_mr = Helicopter.c_mr
    N_mr = Helicopter.N_mr
    theta0 = Helicopter.theta0
    theta_tw = Helicopter.theta_tw
    Cla = Helicopter.Cla
    Cd0 = Helicopter.Cd0
    mu = V_infty/Omega_R_mr  # Rotor advance ratio


    ## GEOMETRIC PROPERTIES

    sol = (N_mr*c_mr)/(R_mr*np.pi)  # Rotor solidity
    r = np.linspace(R_hub, R_mr, num=19)  # Evenly spaced stations along the main rotor
    r_bar = r/R_mr  # Stations along the rotor adimensionalized with respect to the rotor radius R

    theta0 = np.deg2rad(theta0)  # Pitch angle at the root of the rotor expressed in rad
    theta_tw = np.deg2rad(theta_tw)  # Blade twist (tip minus root incidence) expressed in rad
    theta = theta0 + theta_tw*r_bar  # Pitch angle distribution along the rotor (°) (supposed linear)


## DETERMINATION OF THE DISTRIBUTIONS OF ANGULAR AND AERODYNAMIC QUANTITY

# The second degree equation for the induced inflow ratio has coefficient 1, C1 and C2

    C1 = mu + Cla * sol / 8  # First degree coefficient
    C2 = -r_bar * Cla * sol / 8 * (theta - mu / r_bar)  # Free term

    lambda_i = (-C1 + np.sqrt(C1 ** 2 - 4 * C2)) / 2  # Induced inflow ratio

    phi = np.arctan(mu / r_bar + lambda_i / r_bar)  # Inflow angle (rad)
    alpha = theta - phi  # Angle of attack (rad)

    F_Prandtl = 2 / np.pi * np.arccos(np.exp(N_mr / (2 * np.tan(phi[-1])) * (r_bar - 1)))  # Prandtl correction
                                                                                           # function
    F_Prandtl[-1] = 0

    Cl = Cla * alpha  # Lift coefficient
    Cd = Cd0  # Drag coefficient

    dTcdr_bar = 0.5 * sol * Cl * r_bar ** 2  # Incremental Thrust coefficient along the main rotor radius

    dTcdr_bar_PR = dTcdr_bar * F_Prandtl  # Incremental Thrust coefficient along the main rotor radius
                                          # considering the Prandtl correction function

    dPc_i = 0.5 * sol * Cl * phi * r_bar ** 3  # Incremental induced Power coefficient of the rotor along the main rotor
                                               # radius
    dPc_0 = 0.5 * sol * Cd * r_bar ** 3  # Incremental parasite Power coefficient along the main rotor radius

    dPcdr_bar = dPc_i + dPc_0 # Incremental Power coefficient along the main rotor radius

# THRUST, TORQUE, POWER COEFFICIENTS AND FIGURE OF MERIT

    Tc_PR = np.round(np.trapz(dTcdr_bar_PR, r_bar),8)  # Thrust coefficient  considering the Prandtl correction
                                                                # function

    Pc_i = np.round(np.trapz(dPc_i, r_bar),8) # Induced Power (Torque) coefficient
    Pc_0 = np.round(np.trapz(dPc_0, r_bar),8) # Parasite Power (Torque) coefficient

    Qc = np.round(np.trapz(dPcdr_bar, r_bar),8)  # Torque (Power) coefficient

    FM = ((Tc_PR**(3/2)/np.sqrt(2)))/Qc  # Figure of Merit

    return Tc_PR, Pc_i, Pc_0, Qc, dTcdr_bar_PR, dPcdr_bar, FM, F_Prandtl

