# BEMT COUNTER-ROTATING ROTORS IN HOVERING
#
# Description:
# This function is part of the set of functions that computes the performances of 
# a counter-rotating configuration of two rotors in hovering.
# Specifically, this function compute the lambda_i_l of the lower rotor, where: lambda_i_l =  w/(Omega*R); in which:
# - 'w' is the axial induction.
# - 'Omega' is the Rotor angular speed.
# - 'R' is the radius of the blade.
#
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

# The input parameters are:
# - r_segn_vector: adimentionalized distanceof the blade element in respect with the blade radius.
# - r_segn_interf.
# - N: number of blades.
# - solidity vector.
# - twist angle vector.
# - Cl_alpha (assumed 2*pi).
# - Cd0. 
# - lambda_i_u vector of the upper rotor.

def TeoriaElemPalaRotoreLower(r_segn, r_segn_interf, N, solidity, theta, Clalpha, Cd0, lambda_i_u):

   # The vectors that are going to be used are here initialized.
   # This function compute the Prandtl Correction factor iteratively. Two vectors are needed to store the old and the new values.
   F_Prandtl_l = np.zeros(len(r_segn))
   F_Prandtl_new_l = np.zeros(len(r_segn))

   lambda_i_l = np.zeros(len(r_segn))

   # 'phi' is the inflow angle.
   # 'alpha' is the angle of attack
   phi_l = np.zeros(len(r_segn))
   alpha_l = np.zeros(len(r_segn))

   # lift coeffficient and drag coefficient
   Cl_l = np.zeros(len(r_segn))
   Cd_l = np.zeros(len(r_segn))

   dTcdr_segn_l = np.zeros(len(r_segn))
   dTcdr_segn_PR_l = np.zeros(len(r_segn))
   dQcdr_segn_l = np.zeros(len(r_segn))
   dQcdr_segn_PR_l = np.zeros(len(r_segn))
   dPc_i_l = np.zeros(len(r_segn))
   dPc_0_l = np.zeros(len(r_segn))

   # These vectors store the coefficients of the second order linear equation in the unknown 'lambda_i_l'. 
   # Specifically, based on the interference distance along the lower rotor's blade, two different eqaution are solved for both
   # the affected and unaffected elemens of the blade. 
   # The discriminants of the two 2nd order equations are computed. Only real positive values of 'lambda_i_l' shall be taken.
   # 'C1' is the coefficient of the linear term for the unaffected blade element.
   # 'C2' is tthe costant term for the unaffected blade element. 
   # 'C1_int' is the coefficient of the linear term for the affected blade element.
   # 'C2_int' is tthe costant term for the affected blade element.  
   C1 = np.zeros(len(r_segn))
   C2 = np.zeros(len(r_segn))
   C1_int = np.zeros(len(r_segn))
   C2_int = np.zeros(len(r_segn))
   Discriminant = np.zeros(len(r_segn))
   Discriminant_int = np.zeros(len(r_segn))

   for j in range(len(r_segn)):

      # The error and Prandtl factor values are initially set to 1.
      err = 1
      F_Prandtl_l[j] = 1
      
      # Control condition on the error.
      while err > 1e-3:
         
         # 'C1', 'C2', 'C1_int', 'C2_int'  and the discriminants are computed. 
         C1[j] = Clalpha*solidity[j]/(8*F_Prandtl_l[j])
         C2[j] = -r_segn[j]*Clalpha*solidity[j]/(8*F_Prandtl_l[j])*theta[j]     
         C1_int[j] = Clalpha*solidity[j]/(8*F_Prandtl_l[j]) - lambda_i_u[j]/(np.pow(r_segn_interf, 2))
         C2_int[j] = -r_segn[j]*Clalpha*solidity[j]/(8*F_Prandtl_l[j])*theta[j]
         Discriminant[j] = np.pow(C1[j],2)-4*C2[j]
         Discriminant_int[j] = np.pow(C1_int[j],2)-4*C2_int[j]
      
         if r_segn[j] < r_segn_interf:
            if Discriminant_int[j] < 0:
               lambda_i_l[j] = 0
            elif Discriminant_int[j] > 0 and (-C1_int[j]+np.sqrt(np.pow(C1_int[j],2)-4*C2_int[j]))/2 < 0:
               lambda_i_l[j] = 0
            elif Discriminant_int[j] > 0 and (-C1_int[j]+np.sqrt(np.pow(C1_int[j],2)-4*C2_int[j]))/2 > 0:
               lambda_i_l[j] = (-C1_int[j]+np.sqrt(np.pow(C1_int[j],2)-4*C2_int[j]))/2
            
         else:
            if Discriminant[j] < 0:
               lambda_i_l[j] = 0
            elif Discriminant[j] > 0 and (-C1[j]+np.sqrt(np.pow(C1[j],2)-4*C2[j]))/2 < 0:
               lambda_i_l[j] = 0
            elif Discriminant[j] > 0 and (-C1[j]+np.sqrt(np.pow(C1[j],2)-4*C2[j]))/2 > 0:
               lambda_i_l[j] = (-C1[j]+np.sqrt(np.pow(C1[j],2)-4*C2[j]))/2
         
         # At the end of the blade 'lambda_i_l' value must be null.
         lambda_i_l[-1] = 0    

         # Inflow angle and attack angle are computed.
         phi_l[j]   = lambda_i_l[j]/r_segn[j]
         alpha_l[j] = theta[j] - phi_l[j]
            
         # Prandtl Correction Factor in computed through the use o the new computed values.
         # The error values is updated.
         F_Prandtl_new_l[j] = 2/np.pi * np.arccos(np.exp(N / (2*lambda_i_l[j]) * (r_segn[j] - 1)))
         err = np.abs((F_Prandtl_new_l[j] - F_Prandtl_l[j]) / F_Prandtl_l[j])
         F_Prandtl_l[j]     = F_Prandtl_new_l[j]
         
         # At the end of the blade 'F_Prandtl' value must be null.
         F_Prandtl_l[-1] = 0

      # Drag coefficient is assumed to be constant to Cd0
      Cl_l[j] = Clalpha*alpha_l[j]
      Cd_l[j] = Cd0

      dTcdr_segn_l[j] = 0.5*solidity[j]*Cl_l[j]*np.pow(r_segn[j],2) 
      dTcdr_segn_PR_l[j] = dTcdr_segn_l[j]*F_Prandtl_l[j]

      dQcdr_segn_l[j] = 0.5*solidity[j]*Cd_l[j]*np.pow(r_segn[j],3) + lambda_i_l[j]*dTcdr_segn_l[j]
      dQcdr_segn_PR_l[j] = 0.5*solidity[j]*Cd_l[j]*np.pow(r_segn[j],3) + lambda_i_l[j]*dTcdr_segn_PR_l[j]

      dPc_i_l[j] = 0.5*solidity[j]*Cl_l[j]*phi_l[j]*np.pow(r_segn[j],3)
      dPc_0_l[j] = 0.5*solidity[j]*Cd_l[j]*np.pow(r_segn[j],3)

   # Overall coefficients are computed via trapezoidal integration
   Tc_PR_l   = trapezoid(dTcdr_segn_PR_l, r_segn)
   Tc_l      = trapezoid(dTcdr_segn_l, r_segn)

   Qc_l      = trapezoid(dQcdr_segn_l, r_segn)
   Qc_PR_l   = trapezoid(dQcdr_segn_PR_l, r_segn)

   Pc_i_l    = trapezoid(dPc_i_l, r_segn)
   Pc_0_l    = trapezoid(dPc_0_l, r_segn)

   FM_l = (np.pow(Tc_l,3/2)/np.sqrt(2))/Qc_l

   return [dTcdr_segn_l, dTcdr_segn_PR_l, dQcdr_segn_l, dQcdr_segn_PR_l, dPc_i_l, dPc_0_l, Tc_PR_l, Tc_l, Qc_l, Qc_PR_l, Pc_i_l, Pc_0_l, FM_l, F_Prandtl_l]
