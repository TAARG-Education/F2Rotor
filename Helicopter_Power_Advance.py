# Organization: Universita' degli Studi di Napoli - Federico II, Dipartimento di Ingegneria Industriale, Ingegneria Aerospaziale
# Course: Aerodinamica dell'ala rotante
# Professor: Renato Tognaccini
# Supervisor: Ing. Ettore Saetta
# Academic Year: 2023-2024

# Description: this file defines four functions which calculate different non-dimensional form of power relative to an helicopter in advance (i.e. in forward flight).
# In particular:
# - w_tilde_function: solves the fourth grade's equation that describes the operation's curve of a rotor in advance in the hypotesis fo constant thrust. It gives in
#                     output w_tilde, defined as w/w_h, where w is induced velocity for a stiff rotor in advance and w_h is the induced velocity for a stiff rotor in
#                     hovering.
# - P_i_tilde_function: compute P_i_tilde, which is defined as P_i/P_ih, where P_i is the induced power of a stiff rotor in advance and P_ih is the induced power of a
#                       a stiff rotor in hovering.
# - P_c_0_function: computes the parasite power's coefficient of a stiff rotor in advance with respect to a reference's power in some simplyfing hyposthesis.
# - P_c_Fus_function: computes the non-dimensional form of power absorbed by a fuselage of an helicopter in advance with respect to a reference's power.
#                     The power absorbed by the fuselage in forward flight is the work fulfilled by aerdynamic drag on the fuselage in the time's unity.
# More details are in the definitions of the single functions.

#  Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
# Chapter 7: Il rotore rigido in volo traslato, pp.97-105.

# Authors: Afflitto Fernando, Bruno Fabio
# Last change: 27/07/2024


import numpy as np
import math

def w_tilde_function(V_inf_tilde, alpha_deg):

    '''
    Description: this function solves the fourth grade's equation which describes the operation's curve of a stiff rotor in advance. The equation
    is obtained in the hypotesis of rotor's operation with a constant thrust (a dual theory in the hypotesis of constant power exists but it is not
    widespread). The equation is:

    w_tilde^4 + 2*V_inf_tilde*sin(alpha)*w_tilde^3 + V_inf_tilde^2*w_tilde^2 - 1 = 0

    The unknown is w_tilde. In general there are four roots for every couple (V_inf_tilde, alpha), but this function selects only real and positive
    solutions. These solutions are ordered in ascending order and then the function chooses the smallest one.

    Input:
    - V_inf_tilde: defined as V_inf/w_h, where V_inf is asymptotic velocity and w_h is induced velocity on disk rotor in hovering
    - alpha_deg: Angle of Attack (AoA) of disk rotor with respect to asymptotic velocity in degrees.

    Output:
    - w_tilde: defined as w/w_h, where w is induced velocity on disk rotor and w_h is induced velocity on disk rotor in hovering. The function choose
      the smallest real positive solution of the equation.


    Documentation: 
    - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
      Chapter 7: Il rotore rigido in volo traslato, Subsection 7.1.1.: Funzionamento a spinta costante, pp. 98-99.

    Authors: Afflitto Fernando, Bruno Fabio
    Last change: 27/07/2024

    '''

    alpha_rad = math.radians(alpha_deg)      # Conversion of AoA from degrees to radians
    Coefficients = [1, 2*V_inf_tilde*math.sin(alpha_rad), V_inf_tilde**2, 0, -1] # List of coefficients of fourth's grade equation
    Roots = np.roots(Coefficients)           # Roots of equation for the selected couple (V_inf_tilde, alpha_deg)
    
    w_tilde = []                             # Definition of an empty list for w_tilde

    for Root in Roots:
        if np.isreal(Root) and Root > 0:     # Selection of real positive roots 
            w_tilde.append(Root)             # Addition to the list

            w_tilde = np.real(w_tilde)       # Selected real roots have a zero imaginary part (0*j). Function takes into account only real part
            w_tilde.sort()                   # Ordering real positive roots in ascending order
            w_tilde = w_tilde[0]             # Choice of the smallest real positive root
            
            return w_tilde


def P_i_tilde_function(V_inf_tilde, alpha_deg):

    '''
    Description: this function computes P_i_tilde, induced power of a stiff rotor in advance in a non-dimensional form with respect to
    induced power of a stiff rotor in hovering. The relation is:

    P_i_tilde = V_inf_tilde*sin(alpha) + w_tilde

    where w_tilde is the smallest real positive root of the fourth grade's equation which describes the operation's curve of a stiff rotor in
    advance, in the hypotesis of constant thrust. w_tilde is calculated by function "w_tilde_function".

    Input:
    - V_inf_tilde: defined as V_inf/w_h, where V_inf is asymptotic velocity and w_h is induced velocity on disk rotor in hovering
    - alpha_deg: Angle of Attack (AoA) of disk rotor with respect to asymptotic velocity in degrees.

    Output:
    - P_i_tilde: defined as P_i/P_i_h where P_i is induced power of a stiff rotor in advance and P_i_h is induced power of a stiff rotor
      in hovering

    Documenation:
    - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
      Chapter 7: Il rotore rigido in volo traslato, Subsection 7.1.1: Funzionamento a spinta costante, pp. 98-99.
    
    Authors: Afflitto Fernando, Bruno Fabio
    Last change: 27/07/2024

    '''

    alpha_rad = math.radians(alpha_deg)                   # Conversion of AoA from degrees to radians
    w_tilde = w_tilde_function(V_inf_tilde,alpha_deg)     # Resolution of fourth grade's equation for w_tilde and choice of the smallest real 
                                                          # positive roots
    P_i_tilde = V_inf_tilde*math.sin(alpha_rad) + w_tilde # Calculation of P_i_tilde

    return P_i_tilde

def P_c_0_function(sigma, Cd_mean, mu, K = 4.7):

    '''
    Description: this function computes the parasite power's coefficient P_c_0 of a stiff rotor in advance. The parasite power absorbed
    by a generic blade element at the station r along the radius R, with the blade at the generic angular position psi (psi: azimuth angle
    of the blade), is integrated with respect r (from 0 to R) and with respect to psi (from 0 to 2pi). Then this quantity is multiplied for
    /2pi (with N number of blades). In this way the mean parasite power absorbed by the entire stiff rotor in advance, during a single rotation,
    is obtained. Dividing by a reference's power, in the hypotesis of rectangular form in plant for the blade and introducing a mean drag
    coefficient, one can obtain the following non-dimensional expression for the parasite power absorbed by a stiff rotor in advance:

    P_c_0 = (sigma*Cd/8)*(1 + K*mu^2)

    Input:
    - sigma: solidity of disk rotor, defined as N*c/(pi*R), where N is the number of blades, c is the mean chord of the elements and R is the
             radius of disk rotor
    - Cd_mean: mean drag coefficient of the elements of blades
    - K: coefficient which takes into account the effect of the flux's velocity along the direction of radius and other approximations. Without
         this contribution, the theory gives K = 3. Taking into account these effects, one can select a number included in the range [4,5]. 
         Stepniewski and Keys (1984) suggested K = 4.7. This is the default value in this function.
    - mu: advance ratio, defined as V_inf*cos(alpha)/(Omega*R), where V_inf is asymptotic velocity, alpha is Angle of Attack (AoA) of disk rotor
          with respect to the asymptotic velocity, Omega and R are, respectively, the angular velocity and radius of disk rotor.

    Output:
    - P_c_0: Parasite power's coefficient of a stiff rotor in advance.

    Documentation: 
    - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato TOgnaccini, a.a. 2023-2024, vers. 2.05,
      Chapter 7: Il rotore rigido in volo traslato, Section 7.3: Potenza parassita in volo traslato, pp. 100-102.

    Authors: Afflitto Fernando, Bruno Fabio
    Last change: 27/07/2024
                   
    '''

    P_c_0  = (sigma*Cd_mean/8)*(1 + K*mu**2) # Compute of P_c_0

    return P_c_0

def P_c_Fus_function(mu, f_over_A = 0.009):

    '''
    Description: this function computes parasite power's coefficient absorbed by the fuselage of an helicopter in advance. This is the
    non-dimensional form of the work of the aerodynamic drag on the fuselage (the entire helicopter except the main rotor) in the time's unity.
    Starting from the dimensional relation of a power:

    P_fus = f*(1/2)*rho*V_inf^3

    where:
    - f is the equivalent wetted area, has a dimension of a surface and is defined as the sum of C_d_n*S_n. C_d_n and S_n are, respectively,
    the drag coefficient and the reference surface of the n-element of fuselage (generic element of the fuselage).
    - rho: density of the air
    - V_inf: asymptotic velocity
    
    And considering a reference's power, defined as rho*(Omega*R)^3*A, 
    where:
    - Omega is the angular velocity of disk rotor
    - R is the radius of disk rotor
    - A is the area of disk rotor (pi*R^2).

    One can obtain the correpsonding non-dimensional form:
     
    P_c_Fus = (1/2)*(f/A)*mu^3

    where mu is defined as V_inf/(Omega R), with V_inf asymptotic velocity, Omega and R, respctively, the angular velocity and the radius of
    the disk rotor.

    Input:
    - mu: defined as V_inf/(Omega*R), where V_inf is asymptotic velocity, Omega and R are, respectively, the angular velocity and the radius
     of disk rotor
    - f/A: where f is the equivalent wetted area and A is the disk area of rotor. An optimistic value of f/A is 0.007. The default value of
    this function is 0.009, as suggested by Ing. Di Giorgio from Leonardo.

    Output:
    - P_c_Fus: Parasite power's coefficient absorbed by the fuselage of the helicopter in advance (the entire helicopter except the main rotor)

    Documentation: 
    - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
      Chapter 7: Il rotore rigido in volo traslato, Section 7.5: Potenza parassita della fusoliera, pp. 104-105.

    Authors: Afflitto Fernando, Bruno Fabio
    Last change: 27/07/2024

    '''

    P_c_Fus = (f_over_A/2)*mu**3 # Calculation of P_c_Fus
    
    return P_c_Fus

    





            



