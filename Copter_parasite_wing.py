
#This function calculates the aerodynamic drag characteristics of a wing through a series of empirical and theoretical calculations, including factors for skin friction,
#compressibility, and interference. The function uses inputs such as wing span, chord length, wing thickness, and aspect ratio from a configuration dictionary. 
#It also takes into account freestream conditions like Mach number and Reynolds number.It computes the dimensional and non dimensional drag coefficient.

#The theory behind this code can be found in: Carlo De Nicola - "Appunti per il corso di Aerodinamica degli Aeromiboili" - a.a. 2018/2019 - Appendice E.

#Authors: Eduardo Duraccio, Samuel Filosa
#Date: 17/07/2024
#Version: 1.02



import math

def parasite_wing(helicopter):
    M_oo = helicopter.M_inf  # freestream Mach number
    Re_oo = helicopter.Re_inf  # freestream Reynolds number
    turbolentflow = helicopter.turbolentflow # turbolent flow
    
    b = helicopter.b_w  # wing span
    c = helicopter.c_w  # wing chord
    tau = helicopter.tau  # wing thickness
    AR_w = helicopter.AR  # aspect ratio of wing
    r_mr = helicopter.R_mr # radius of main rotor 

    R_LS = helicopter.R_LS # lift-induced skin friction factor 
    R_w_b = helicopter.R_w_b # interference factor wing-fusolage
    Hf = helicopter.Hf # hinge factor

    r = helicopter.r # compressibility correction factor
    gamma = helicopter.gamma # ratio of specific heats for air
    K = helicopter.k # admissible roughness
    ala_fissa = helicopter.ala_fissa # fixed or rotary wing
    profilo = helicopter.profilo # type of airfoil

    S = b**2/AR_w #wing area
    S_rotor= 3.14*r_mr**2 # reference area (rotor area)
    Swet = 2 * S  # wetted area of the wing (both sides)

    if profilo == "4-series":                       # insert "4-series" to use 4-series airfoil
        FF = 1 + 1.68 * (tau / c) + 3 * (tau / c)**2    # form factor for 4-series airfoil
    elif profilo == "biconvex":                       # insert "biconvex" to use biconvex airfoil
        FF= 1 + 1.20 * tau + 100 * tau**4               # form factor for biconvex-airfoil
    elif profilo == "6-series":                      # insert "6-series" to use 6-series airfoil
        FF = 1 + 1.44 * tau + 2 * tau**2                # form factor for 6-series airfoil
    else:
        print("Errore nell'immissione del profilo")
    
    

    if ala_fissa == False:   #select "False" if there's a rotary wing
        IF = R_LS * Hf                  # interference factor for rotary wing  
    else:                               # select "True" if there's a fixed wing
        IF = R_LS * R_w_b               # interference factor for fixed wing
    
    t = (1 + r * (gamma - 1) / 2 * M_oo**2)**(-1)  # thermodynamic temperature ratio factor
    F = 1 + 0.03916 * M_oo**2 * t  # compressibility correction factor

    Re_L1 = Re_oo  # initial Reynolds number condition
    K1 = 37.587 + 4.615 * M_oo + 2.949 * M_oo**2 + 4.132 * M_oo**3  # admissible roughness
    Re_L2 = K1 * (c / K)**(1.0489)  # second Reynolds number condition
    Re_L = min([Re_L1, Re_L2])  # select the smaller Reynolds number condition

    if turbolentflow == True:
        Cf = t * F**2 * 0.430 / (math.log10(Re_L * t**1.67 * F))**2.56  # skin friction coefficient for turbolent flow
    else:
        Cf = (1.328 / (Re_L**(1/2))) * (1 + 0.1256 * M_oo**2)**(-0.12)   #skin friction coefficient for laminar flow


    if helicopter.R_mr == 0: # select "0" if it's not an helicopter
        f_Awing = (Cf * Swet / S) * FF * IF  # drag coefficient for the wing
        f_wing = f_Awing * S  # dimensional drag coefficient for the wing
    else: 
        f_Awing = (Cf * Swet / S_rotor) * FF * IF  # drag coefficient for the wing
        f_wing = f_Awing * S_rotor  # dimensional drag coefficient for the wing

    f_Awing = round(f_Awing, 5)  # round drag coefficient for the wing
    f_wing = round(f_wing, 4)  # round dimensional drag coefficient for the wing

    return f_Awing, f_wing  # return the dimensional and  non dimensional drag coefficient and total drag
