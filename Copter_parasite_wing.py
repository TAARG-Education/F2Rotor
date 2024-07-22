"""
This function calculates the aerodynamic drag characteristics of a wing through a series of empirical and theoretical calculations, including factors for skin friction,
compressibility, and interference. The function uses inputs such as wing span, chord length, wing thickness, and aspect ratio from a configuration dictionary. 
It also takes into account freestream conditions like Mach number and Reynolds number.It computes the dimensional and non dimensional drag coefficient.

The theory behind this code can be found in: Carlo De Nicola - "Appunti per il corso di Aerodinamica degli Aeromiboili" - a.a. 2018/2019 - Appendice E.

    Authors: Eduardo Duraccio, Samuel Filosa
    Date: 17/07/2024
    Version: 1.01

"""

def wing_comp(self):
    b = self.config["Lift-Compound"][0]["b_w"]  # wing span
    c = self.config["Lift-Compound"][0]["c_w"]  # rectangular wing chord
    tau = self.config["Lift-Compound"][0]["tau"]  # wing thickness
    AR_w = self.config["Lift-Compound"][0]["AR"]  # aspect ratio of rectangular wing
    r_mr = self.config["Lift-Compound"][0]["r_mr"] # radius of main rotor 
    R_LS = self.config["Lift-Compound"][0]["R_LS"] # lift-induced skin friction factor 
    Hf = self.config["Lift-Compound"][0]["Hf"] # hinge factor
    r = self.config["Lift-Compound"][0]["r"] # compressibility correction factor
    gamma = self.config["Lift-Compound"][0]["gamma"] # ratio of specific heats for air
    K = self.config["Lift-Compound"][0]["K"] # admissible roughness
    M_oo = self.M_inf  # freestream Mach number
    Re_oo = self.Re_inf  # freestream Reynolds number

    S = b**2/AR_w #wing area
    S_rotor= 3.14*r_mr**2 # reference area (rotor area)
    Swet = 2 * S  # wetted area of the wing (both sides)
    FF = 1 + 1.68 * (tau / c) + 3 * (tau / c)**2  # form factor for 4-series airfoil
    IF = R_LS * Hf # interference factor for horizontal stabilizer  
    
    t = (1 + r * (gamma - 1) / 2 * M_oo**2)**(-1)  # thermodynamic temperature ratio factor
    F = 1 + 0.03916 * M_oo**2 * t  # compressibility correction factor

    Re_L1 = Re_oo  # initial Reynolds number condition
    K1 = 37.587 + 4.615 * M_oo + 2.949 * M_oo**2 + 4.132 * M_oo**3  # # admissible roughness
    Re_L2 = K1 * (c / K)**(1.0489)  # second Reynolds number condition
    Re_L = min([Re_L1, Re_L2])  # select the smaller Reynolds number condition

    Cf = t * F**2 * 0.430 / (math.log10(Re_L * t**1.67 * F))**2.56  # skin friction coefficient for turbolent flow
    f_Awing = (Cf * Swet / S_rotor) * FF * IF  # drag coefficient for the wing
    f_wing = f_Awing * S_rotor  # dimensional drag coefficient for the wing

    f_Awing = round(f_Awing, 5)  # round drag coefficient for the wing
    f_wing = round(f_wing, 4)  # round dimensional drag coefficient for the wing

    return f_Awing, f_wing  # return the dimensional and  non dimensional drag coefficient and total drag