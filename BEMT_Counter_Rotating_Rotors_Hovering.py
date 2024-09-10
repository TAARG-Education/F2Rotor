# BEMT COUNTER-ROTATING ROTORS IN HOVERING
#
# Description:
# This function computes the performances of a counter-rotating configuration of two rotors in , displayed on the same axis.
# Specifically, lambda_i is computed, where: lambda_i_l =  w/(Omega*R); in which:
# - 'w' is the axial induction.
# - 'Omega' is the Rotor angular speed.
# - 'R' is the radius of the blade.
# by resolving a 2nd order linear equation. This parameter will then be used to evaluate the performances of both rotors.
#
# The Hypothesis assumed are:
# - Inviscid Flow.
# - Incompressible Flow.
# - It is assumed to be working in the linear trend of the lift curve, whose slope is assumed to be 2*pi.
# - a' = 0. 
# - phi << 1.
# - The lower rotor does not influence by any mean the upper rotor. 
# - The Influence of the upper rotor is consistent. An interference factor is assumed: r_segn_interf = 0.71. 
#   The interference factor accounts the elements of the lower rotor's blades that are affected by the wake of the upper rotor.
#   This values may be changed by the user if needed. 
#
# Authors of the Function: Giovanni Perez, Catello Meglio.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic year 2023/2024.
#
# References: 'Lezioni di Aerodinamica dell'ala rotante', Renato Tognaccini: Chapter 6 - Paragraph 6.10.2. 

import numpy as np
from scipy.integrate import trapezoid

# The notation 'u' will denote the upper rotor.
# The notation 'l' will denote the upper rotor.

# Number of points with which discretize the domain
N_points = 100

# Example variables for the upper rotor
r_segn_u = np.linspace(0.1, 1, N_points) # Adimensionalized element blase distance
N_u      = 2 # Number of blades
solidity_u = 0.027*np.ones(N_points) # Solidity vector
theta_u = np.deg2rad(np.linspace(0, 10, N_points)) # Twist angle vector.
Clalpha_u = 2*np.pi 
Cd0_u = 0.011 

# Example variables for the lower rotor
r_segn_l = np.linspace(0.1, 1, N_points) # Adimensionalized element blase distance
N_l      = 2 # Number of blades
solidity_l = 0.027*np.ones(N_points) # Solidity vector
theta_l = np.deg2rad(np.linspace(0, 10, N_points)) # Twist angle vector.
Clalpha_l = 2*np.pi 
Cd0_l = 0.011 

r_segn_interf = 0.81 # Distance of the element of the blades that are affected by the wake of the upper rotor.

# The inputs of the function are:
# - (r_segn_u, r_segn_l): arrays that define the normalised blade element distance from the axis of the rotor or the hub
# - (N_u, N_l): Number of blades
# - (solidity_u, solidity_l): solidity arrays 
# - (theta_u, theta_l): twist angle arrays
# - (Clalpha_u, Clalpha_l): it is assumed to work in the linear part of the Cl-alpha trend
# - (Cd0_u, Cd0_l)
# - r_segn_interf: represents the last blade element of the lower rotor which is affected by the wake of the upper rotor
# - N_points: Number of points with which discretize the domain

def BEMT_Counter_Rotating_Rotors_Hovering(r_segn_u, r_segn_l, N_u, N_l, solidity_u, solidity_l, theta_u, theta_l, Clalpha_u, Clalpha_l, Cd0_u, Cd0_l, r_segn_interf, N_points):

    # The script works with matrices, in which the first row (i = 0) contains the values of the variable relative to the upper rotor.
    # Whereas the second row of the matrix (i = 1) contains the values of the vaiables of the lower rotor.

    # np.vstack is a function imported from the library 'numpy' that stacks the arrays to form a matrix.
    r_segn = np.vstack((r_segn_u, r_segn_l))
    N = np.vstack((N_u, N_l))
    solidity = np.vstack((solidity_u, solidity_l))
    theta = np.vstack((theta_u, theta_l))
    Clalpha = np.vstack((Clalpha_u, Clalpha_l))
    Cd0 = np.vstack((Cd0_u, Cd0_l))
    
    # The following variables are initialised as 'zeros'.
    lambda_i = np.zeros((2,N_points))

    # 'phi' is the inflow angle.
    # 'alpha' is the angle of attack
    phi = np.zeros((2,N_points))
    alpha = np.zeros((2,N_points))

    Cl = np.zeros((2,N_points))
    Cd = np.zeros((2,N_points))
    
    # 'dTcdr_segn' and 'dTcdr_segn_PR' accomodates the values of the Thrust coefficient over r_segn, 
    # without Prandtl correction and with Prandtl correction respectively
    dTcdr_segn = np.zeros((2,N_points))
    dTcdr_segn_PR = np.zeros((2,N_points))
    Tc = [0,0]
    Tc_PR = [[0],[0]]

    # 'dQcdr_segn' and 'dQcdr_segn_PR' accomodates the values of the Thrust coefficient over r_segn, 
    # without Prandtl correction and with Prandtl correction respectively
    dQcdr_segn = np.zeros((2,N_points))
    dQcdr_segn_PR = np.zeros((2,N_points))
    Qc = [[0],[0]]
    Qc_PR = [[0],[0]]


    # 'dPc_i' and 'dP_0' accomodates the values of the inducted power and parasite power 
    dPc_i = np.zeros((2,N_points))
    dPc_0 = np.zeros((2,N_points))
    Pc_i = [[0],[0]]
    Pc_0 = [[0],[0]]

    FM = [[0],[0]]

    # 'C1' and 'C2' matrices accomodates the coefficients of the linear term and constant term of the 2nd order equation, respectively.
    # 'Discriminant' matrix accomodates the values of the discriminant of the equation.
    C1 = np.zeros((2,N_points))
    C2 = np.zeros((2,N_points))
    Discriminant = np.zeros((2,N_points))

    # 'C1_int' and 'C2_int' matrices accomodates the coefficients of the lements that are affected by the wake.
    # 'Discriminant_int' matrix has the same purpose of 'Discriminant'.
    C1_int = np.zeros((2,N_points))
    C2_int = np.zeros((2,N_points))
    Discriminant_int = np.zeros((2,N_points))

    # error matrix and F_prandtl matrix are initialised as 'ones' due to the control condition.
    # F_prandtl_new accomodates the newer values of F_prandtl, so it is initialised as zeros.
    err = np.ones((2,N_points))
    F_prandtl = np.ones((2,N_points))
    F_prandtl_new = np.zeros((2,N_points))

    # The first cycle discriminates the upper from the lower, so that 'i == 0' is related to the upper, 'i == 1' to the lower
    for i in range(0,2):
        for j in range(0,N_points-1):

            # Control condition on the error.
            while err[i,j] > 1e-3:        

                # The 2nd order linear equation in the unknown lambda_i is solved for both rotors. Note that the solution must 
                # have only real positive values, whence the control on the discriminant.
                if i == 0:
                    C1[i,j] = Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j])
                    C2[i,j] = -theta[i,j]*r_segn[i,j]*Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j])
                    Discriminant[i,j] = np.pow(C1[i,j],2)-4*C2[i,j]

                    if Discriminant[i,j] < 0: lambda_i[i,j] = 0
                    elif Discriminant[i,j] > 0 and (-C1[i,j]+np.sqrt(Discriminant[i,j]))/2 < 0: lambda_i[i,j] = 0
                    elif Discriminant[i,j] > 0 and (-C1[i,j]+np.sqrt(Discriminant[i,j]))/2 > 0: 
                        lambda_i[i,j] = (-C1[i,j]+np.sqrt(Discriminant[i,j])) / 2

                elif i == 1:

                    C1[i,j] = Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j])
                    C2[i,j] = -r_segn[i,j]*Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j])*theta[i,j]     
                    
                    C1_int[i,j] = Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j]) - lambda_i[0,j]/(np.pow(r_segn_interf, 2))
                    C2_int[i,j] = -r_segn[i,j]*Clalpha[i]*solidity[i,j]/(8*F_prandtl[i,j])*theta[i,j]

                    Discriminant[i,j] = np.pow(C1[i,j],2)-4*C2[i,j]
                    Discriminant_int[i,j] = np.pow(C1_int[i,j],2)-4*C2_int[i,j]

                    if r_segn[i,j] < r_segn_interf:
                        if Discriminant_int[i,j] < 0:
                            lambda_i[i,j] = 0
                        elif Discriminant_int[i,j] > 0 and (-C1_int[i,j]+np.sqrt(np.pow(C1_int[i,j],2)-4*C2_int[i,j]))/2 < 0:
                            lambda_i[i,j] = 0
                        elif Discriminant_int[i,j] > 0 and (-C1_int[i,j]+np.sqrt(np.pow(C1_int[i,j],2)-4*C2_int[i,j]))/2 > 0:
                            lambda_i[i,j] = (-C1_int[i,j]+np.sqrt(Discriminant_int[i,j]))/2
                    else:
                        if Discriminant[i,j] < 0:
                            lambda_i[i,j] = 0
                        elif Discriminant[i,j] > 0 and (-C1[i,j]+np.sqrt(np.pow(C1[i,j],2)-4*C2[i,j]))/2 < 0:
                            lambda_i[i,j] = 0
                        elif Discriminant[i,j] > 0 and (-C1[i,j]+np.sqrt(np.pow(C1[i,j],2)-4*C2[i,j]))/2 > 0:
                            lambda_i[i,j] = (-C1[i,j]+np.sqrt(Discriminant[i,j]))/2

                # At the end of the blade 'lambda_i_u' value must be null.
                lambda_i[i,-1] = 0

                # Inflow angle and attack angle are computed.
                phi[i,j]   = lambda_i[i,j]/r_segn[i,j]
                alpha[i,j] = theta[i,j]-phi[i,j]

                # Prandtl Correction Factor in computed through the use o the new computed values.
                # The error values is updated.
                F_prandtl_new[i,j] = (2/np.pi)*np.arccos(np.exp(N[i]/(2*lambda_i[i,j])*(r_segn[i,j]-1)))
                err[i,j] = abs((F_prandtl_new[i,j] - F_prandtl[i,j]) / F_prandtl[i,j])
                F_prandtl[i,j] = F_prandtl_new[i,j]

                # At the end of the blade 'F_Prandtl' value must be null
                F_prandtl[i,-1] = 0 

            # Drag coefficient is assumed to be constant to Cd0
            Cl[i,j] = Clalpha[i]*alpha[i,j]
            Cd[i,j] = Cd0[i]

            dTcdr_segn[i,j] = 0.5*solidity[i,j]*Cl[i,j]*np.pow(r_segn[i,j],2)
            dTcdr_segn_PR[i,j] = dTcdr_segn[i,j]*F_prandtl[i,j]

            dQcdr_segn[i,j] = 0.5*solidity[i,j]*Cd[i,j]*np.pow(r_segn[i,j],3) + lambda_i[i,j]*dTcdr_segn[i,j]
            dQcdr_segn_PR[i,j] = 0.5*solidity[i,j]*Cd[i,j]*np.pow(r_segn[i,j],3) + lambda_i[i,j]*dTcdr_segn_PR[i,j]

            dPc_i[i,j] = 0.5*solidity[i,j]*Cl[i,j]*phi[i,j]*np.pow(r_segn[i,j],3)
            dPc_0[i,j] = 0.5*solidity[i,j]*Cd[i,j]*np.pow(r_segn[i,j],3)

        # Overall coefficients are computed via trapezoidal integration
        Tc_PR[i]   = trapezoid(dTcdr_segn_PR[i], r_segn[i,:])
        Tc[i]      = trapezoid(dTcdr_segn[i], r_segn[i,:])

        Qc[i]      = trapezoid(dQcdr_segn[i], r_segn[i,:])
        Qc_PR[i]   = trapezoid(dQcdr_segn_PR[i], r_segn[i,:])

        Pc_i[i]    = trapezoid(dPc_i[i], r_segn[i,:])
        Pc_0[i]    = trapezoid(dPc_0[i], r_segn[i,:])

        FM[i] = (np.pow(Tc_PR[i],3/2)/np.sqrt(2))/Qc_PR[i]

    return [r_segn, dTcdr_segn, dTcdr_segn_PR, dQcdr_segn, dQcdr_segn_PR, dPc_i, dPc_0, Tc, Tc_PR, Qc, Qc_PR, Pc_i, Pc_0, FM, F_prandtl]


   