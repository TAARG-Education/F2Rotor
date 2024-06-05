from BEMT_Hover import BEMT_Hover
from Helicopter_Class import Helicopter
import numpy as np

def tailrotor(Helicopter, rho, Voo):
    
    '''
    This function computes the induced power, profile power, and total power
    needed for the tail rotor operation. The tail rotor counteracts the torque
    produced by the main rotor, ensuring directional stability and control.

    Steps and parameters involved:
    1. Retrieve the main rotor torque (Qc) using Blade Element Momentum Theory (BEMT).
    2. Extract necessary tail rotor configuration parameters from the helicopter's config dictionary:
        - R_tr: Tail rotor radius
        - c_tr: Tail rotor chord
        - N_tr: Number of tail rotor blades
        - Omega_tr: Tail rotor angular velocity
        - b_tr: Arm relative to the helicopter's center of gravity for the calculation of counter-torque
    3. Calculate the tail rotor tip speed (Vtip_tr) and disk area (A_tr).
    4. Compute the tail rotor torque (Q) based on the main rotor torque, air density, and other factors.
    5. Define constants for induced power correction (k_tr) and profile drag coefficient (CD_tr).
    6. Calculate the advance ratio (mu_tr) based on forward velocity (Voo) and tip speed.
    7. Determine the tail rotor thrust (T_tr) and thrust coefficient (T_tr_c).
    8. Compute the inflow ratio (lambda_it) under two different conditions:
        - If the advance ratio is less than 0.1, use a quadratic equation.
        - Otherwise, use a simplified formula.
    9. Iterate to refine the inflow ratio (lambda_t) until the error is minimized.
    10. Calculate the induced power (P_tr_i) and profile power (P_tr_0).
    11. Sum the induced and profile power to get the total tail rotor power (P_tr).
    12. Return the induced power, profile power, and total power.

    INPUT: Qc_mr, TailRotorClass
    OTPUT: P_tr_i, P_tr_0, P_tr

    Reference: "Lezioni di aerodinamica dell'ala rotante -
    Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria", Renato Tognaccini.
    
    Authors: Fabio Beltratti, Giuseppe Russo
    Latest update: 03/06/2024
    Version: 0.9
    '''

    # INPUT
    # --------------------------------------------------
    Qc = BEMT_Hover(Helicopter, Voo)[3]

    R_tr = Helicopter.tailRotor.R_tr # tail rotor radius
    c_tr = Helicopter.tailRotor.c_tr   # tail rotor chord
    N_tr = Helicopter.tailRotor.N_tr   # tail rotor blades number
    Omega_tr = Helicopter.tailRotor.Omega_tr # tip speed
    b_tr = Helicopter.tailRotor.b_tr # distance between main rotor and tail rotor
    # ---------------------------------------------------

    Vtip_tr = Omega_tr*R_tr # tail rotor tip speed
    A_tr = np.pi*R_tr**2 # tail rotor disk area 
    Q = Qc*rho*Vtip_tr**2*R_tr*A_tr # tail rotor torque = main rotor torque
   
    k_tr = 1.15 # tip loss correction
    CD_tr = 0.01 # Cd tail rotor
    mu_tr = Voo / Vtip_tr # advance ratio tail rotor

    T_tr = Q / b_tr  # tail rotor thrust
    T_tr_c = T_tr / (rho * Vtip_tr ** 2 * A_tr) # tail rotor thrust coefficient

    # iteration for rotor induction

    if np.all(mu_tr) < 0.1:

        lambda_it = np.sqrt(-Voo ** 2 / 2 + np.sqrt(Voo ** 4 / 4
                                                     + (T_tr / (2 * rho * A_tr)) ** 2)) / Vtip_tr
    else:
        lambda_it = T_tr_c / (2 * mu_tr)

    lambda_t = lambda_it

    err = 1

    while np.all(err > 1e-3):
        lambda_it = T_tr_c / (2 * np.sqrt(mu_tr ** 2 + lambda_t ** 2))
        err = (lambda_t - lambda_it) / lambda_t
        lambda_t = lambda_it

    # OUTPUT 
    # Tail rotor power required
    P_tr_i = k_tr * lambda_it * Vtip_tr * T_tr # induced power
    P_tr_0 = 1 / 8 * rho * N_tr* c_tr * R_tr * CD_tr * Vtip_tr ** 3 * (1 + 4.7 * mu_tr ** 2) # parasite power
    P_tr = P_tr_i + P_tr_0 # required power

    return P_tr_i, P_tr_0, P_tr
