import numpy as np # type: ignore

def TailRotor(Helicopter, V_infty, Qc, k_tr = 1.20):
      
    '''
    This function computes the induced power, profile power, and total power
    coefficients needed for the tail rotor operation. The tail rotor counteracts the torque
    produced by the main rotor, ensuring directional stability and control.

    Steps and parameters involved:
    1. Accept in input:
        - Helicopter class with main and tail rotor parameters 
        - V_infty 
        - the main rotor torque coefficient (Qc) --> one can retrive it from Helicopter_power_advance function
        - k_tr which is the factor related to the induced power losses (default value: 1.20)
    2. Extract necessary tail rotor configuration parameters from the helicopter's config dictionary:
        - R_tr: Tail rotor radius
        - c_tr: Tail rotor chord
        - N_tr: Number of tail rotor blades
        - Omega_tr: Tail rotor tip speed
        - b_tr: Arm relative to the helicopter's center of gravity for the calculation of counter-torque
    3. Calculate the tail rotor disk area (A_tr).
    4. Define constants for induced power correction (k_tr) and profile drag coefficient (CD_tr).
    5. Calculate the advance ratio (mu_tr) based on forward velocity (V_infty) and tip speed.
    6. Determine the dimensionless ratio of tail rotor thrust to rho.
    7. Compute the inflow ratio (lambda_it) using a quadratic equation.
    8. Calculate the induced power (Cp_tr_i) and profile power (Cp_tr_0) coefficients.
    9. Sum the induced and profile power coefficients to get the total tail rotor power coefficient (Cp_tr).
    10. Return the induced power, profile power, and total power coefficients.

    INPUT: Qc, Helicopter class, V_infty, k_tr
    OUTPUT: Cp_tr_i, Cp_tr_0, Cp_tr

    References: 
    - "Lezioni di aerodinamica dell'ala rotante - Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria",
       Renato Tognaccini. From page 97 to 100. 
    - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024, Giovanni Di Giorgio. 

    Authors: Fabio Beltratti, Giuseppe Russo
    Latest update: 21/06/2024
    Version: 1.0
    '''
    
    # Main rotor 
    Omega_r_mr = Helicopter.Omega_r_mr                              # main rotor tip speed
    R_mr = Helicopter.R_mr                                          # main rotor radius 
    A_mr = np.pi * R_mr ** 2                                        # main rotor area

    # Tail rotor
    R_tr = Helicopter.R_tr                                          # tail rotor radius
    c_tr = Helicopter.c_tr                                          # tail rotor chord
    N_tr = Helicopter.N_tr                                          # tail rotor blades number
    Omega_r_tr = Helicopter.Omega_r_tr                              # tail rotor tip speed                 
    l_tr = Helicopter.l_tr                                          # tail rotor torque arm
    A_tr = np.pi * R_tr ** 2                                        # tail rotor area 

    CD_tr = 0.01                                                    # tail rotor profile drag coefficient                    
    mu_tr = V_infty / Omega_r_tr                                    # tail rotor advance ratio (alpha << 1)

    T_tr_over_rho = Omega_r_mr ** 2 * A_mr * R_mr * Qc / l_tr       # tail rotor thrust
    
    print(f"Main Rotor Tip Speed (Omega_R_mr): {Omega_r_mr}")
    print(f"Main Rotor Radius (R_mr): {R_mr}")
    print(f"Main Rotor Area (A_mr): {A_mr}")
    print(f"Tail Rotor Radius (R_tr): {R_tr}")
    print(f"Tail Rotor Chord (c_tr): {c_tr}")
    print(f"Tail Rotor Blades Number (N_tr): {N_tr}")
    print(f"Tail Rotor Tip Speed (Omega_r_tr): {Omega_r_tr}")
    print(f"Tail Rotor Torque Arm (l_tr): {l_tr}")
    print(f"Tail Rotor Area (A_tr): {A_tr}")
    print(f"Tail Rotor Profile Drag Coefficient (CD_tr): {CD_tr}")
    print(f"Tail Rotor Advance Ratio (mu_tr): {mu_tr}")
    print(f"Tail Rotor Thrust over Air Density (T_tr_over_rho): {T_tr_over_rho}")
    
    # induced input ratio
    lambda_it = np.sqrt(-V_infty ** 2 / 2 + 0.5 * np.sqrt(V_infty ** 4 
                     + 4 * (T_tr_over_rho / (2 * A_tr)) ** 2)) / Omega_r_tr
    

    w_i = lambda_it * Omega_r_tr

    Cp_tr_i = k_tr * lambda_it * T_tr_over_rho / (Omega_r_tr ** 2 * A_tr)
    Cp_tr_0 = 1 / 8 * N_tr * c_tr * R_tr * CD_tr * (1 + 4.7 * mu_tr ** 2) / A_tr
    Cp_tr = Cp_tr_0 + Cp_tr_i

    '''
    print(f"Induced Power Coefficient (Cp_tr_i): {Cp_tr_i}")
    print(f"Profile Power Coefficient (Cp_tr_0): {Cp_tr_0}")
    print(f"Total Power Coefficient (Cp_tr): {Cp_tr}")
    '''

    return Cp_tr_i, Cp_tr_0, Cp_tr