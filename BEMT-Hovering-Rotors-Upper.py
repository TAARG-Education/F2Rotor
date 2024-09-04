# BEMT COUNTER-ROTATING ROTORS IN HOVERING
#
# Description:
# This function is part of the set of functions that computes the performances
# of a counter-rotating configuration of two rotors in hovering.
# Specifically, this function compute the lambda_i_u of the upper rotor, where: lambda_i_u =  w/(Omega*R); in which:
# - 'w' is the axial induction.
# - 'Omega' is the Rotor angular speed.
# - 'R' is the radius of the blade.
# 
# Furthermore the following quantities are computed:
# - Thrust Coefficient per unit length, with Prandtl correction applied and without.
# - Torque Coefficient per unit length, with Prandtl correction applied and without.
# - Elemental Power per unit length. 
# - Overall Thrust coefficient, with Prandtl correction applied and without.
# - Overall Torque coefficient, with Prandtl correction applied and without.
# - Overall Power coefficient, with Prandtl correction applied and without.
#
# The lambda_i_u output will then be used to compute the performances of the lower rotor.
#
# The Hypothesis assumed are:
# - Inviscid Flow.
# - Incompressible Flow.
# - It is assumed to be working in the linear trend of the lift curve, whose slope is assumed to be 2*pi.
# - a' = 0. 
# - phi << 1.
# - The lower rotor does not influence by any mean the upper rotor.  
#
# Authors of the Function: Giovanni Perez, Catello Meglio.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic year 2023/2024.
#
# References: 'Lezioni di Aerodinamica dell'ala rotante', Renato Tognaccini: Chapter 6 - Paragraph 6.10.2. 

import numpy as np
from scipy.integrate import trapezoid

# The input parameters are:
# - r_segn vector: adimentionalized distance of the blade element in respect with the blade radius.
# - N: number of blades
# - solidity vector
# - twist angle vector
# - Cl_alpha (assumed 2*pi)
# - Cd0 

def TeoriaElemPalaRotoreUpper(r_segn, N, solidity, theta, Clalpha, Cd0):
    
    # The vectors that are going to be used are here initialized.
    # This function compute the Prandtl Correction factor iteratively. Two vectors are needed to store the old and the new values.
    F_prandtl_u = np.zeros(len(r_segn))
    F_prandtl_new_u = np.zeros(len(r_segn))

    lambda_i_u = np.zeros(len(r_segn))

    # 'phi' is the inflow angle.
    # 'alpha' is the angle of attack
    phi_u = np.zeros(len(r_segn))
    alpha_u = np.zeros(len(r_segn))

    # lift coeffficient and drag coefficient
    Cl_u = np.zeros(len(r_segn))
    Cd_u = np.zeros(len(r_segn))

    dTcdr_segn_u = np.zeros(len(r_segn))
    dTcdr_segn_PR_u = np.zeros(len(r_segn))
    dQcdr_segn_u = np.zeros(len(r_segn))
    dQcdr_segn_PR_u = np.zeros(len(r_segn))
    dPc_i_u = np.zeros(len(r_segn))
    dPc_0_u = np.zeros(len(r_segn))

    # These vectors store the coefficients of the second order linear equation in the unknown 'lambda_i_u'.
    # 'C1' is the coefficient of the linear term.
    # 'C2' is the costant term.  
    C1 = np.zeros(len(r_segn))
    C2 = np.zeros(len(r_segn))

    for j in range(len(r_segn)):

        # The error and Prandtl factor values are initially set to 1.    
        err = 1
        F_prandtl_u[j] = 1

        # Control condition on the error.
        while err > 1e-3:

            # 'C1' and 'C2' are computed. 
            # The discriminant of the 2nd order equation is computed. Only real positive values of 'lambda_i_u' shall be taken.
            C1[j] = Clalpha*solidity[j]/(8*F_prandtl_u[j])
            C2[j] = -theta[j]*r_segn[j]*Clalpha*solidity[j]/(8*F_prandtl_u[j])
            Discriminant = np.pow(C1[j],2)-4*C2[j]

            if Discriminant < 0:
                lambda_i_u[j] = 0
            elif Discriminant > 0 and (-C1[j]+np.sqrt(Discriminant))/2 < 0:
                lambda_i_u[j] = 0
            elif Discriminant > 0 and (-C1[j]+np.sqrt(Discriminant))/2 > 0:
                lambda_i_u[j] = (-C1[j]+np.sqrt(Discriminant)) / 2

            # At the end of the blade 'lambda_i_u' value must be null.
            lambda_i_u[-1] = 0

            # Inflow angle and attack angle are computed.
            phi_u[j]   = lambda_i_u[j]/r_segn[j]
            alpha_u[j] = theta[j]-phi_u[j]

            # Prandtl Correction Factor in computed through the use o the new computed values.
            # The error values is updated.
            F_prandtl_new_u[j] = (2/np.pi)*np.arccos(np.exp(N/(2*lambda_i_u[j])*(r_segn[j]-1)))
            err = abs((F_prandtl_new_u[j] - F_prandtl_u[j]) / F_prandtl_u[j])
            F_prandtl_u[j] = F_prandtl_new_u[j]

            # At the end of the blade 'F_Prandtl' value must be null
            F_prandtl_u[-1] = 0 

        # Drag coefficient is assumed to be constant to Cd0
        Cl_u[j] = Clalpha*alpha_u[j]
        Cd_u[j] = Cd0

        dTcdr_segn_u[j] = 0.5*solidity[j]*Cl_u[j]*np.pow(r_segn[j],2) 
        dTcdr_segn_PR_u[j] = dTcdr_segn_u[j]*F_prandtl_u[j]

        dQcdr_segn_u[j] = 0.5*solidity[j]*Cd_u[j]*np.pow(r_segn[j],3) + lambda_i_u[j]*dTcdr_segn_u[j]
        dQcdr_segn_PR_u[j] = 0.5*solidity[j]*Cd_u[j]*np.pow(r_segn[j],3) + lambda_i_u[j]*dTcdr_segn_PR_u[j]

        dPc_i_u[j] = 0.5*solidity[j]*Cl_u[j]*phi_u[j]*np.pow(r_segn[j],3)
        dPc_0_u[j] = 0.5*solidity[j]*Cd_u[j]*np.pow(r_segn[j],3)

    # Overall coefficients are computed via trapezoidal integration
    Tc_PR_u   = trapezoid(dTcdr_segn_PR_u, r_segn)
    Tc_u      = trapezoid(dTcdr_segn_u, r_segn)

    Qc_u      = trapezoid(dQcdr_segn_u, r_segn)
    Qc_PR_u   = trapezoid(dQcdr_segn_PR_u, r_segn)

    Pc_i_u    = trapezoid(dPc_i_u, r_segn)
    Pc_0_u    = trapezoid(dPc_0_u, r_segn)

    FM_u = (np.pow(Tc_PR_u,3/2)/np.sqrt(2))/Qc_PR_u        

    return [dTcdr_segn_u, dTcdr_segn_PR_u, dQcdr_segn_u, dQcdr_segn_PR_u, dPc_i_u, dPc_0_u, Tc_PR_u, Tc_u, Qc_u, Qc_PR_u, Pc_i_u, Pc_0_u, FM_u, F_prandtl_u, lambda_i_u]
