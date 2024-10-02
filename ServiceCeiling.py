# Organization: Universita' degli Studi di Napoli - Federico II, Dipartimento di Ingegneria Industriale, Ingegneria Aerospaziale
# Course: Aerodinamica dell'ala rotante
# Professor: Renato Tognaccini
# Supervisor: Ing. Ettore Saetta, Ing. Michele Massa
# Academic Year: 2023-2024

# This function computes the theoretical maximum altitude and the practical maximum altitude (or ceilings) of an helicopter in advance 
# (i.e. in forward flight). 
# The function takes in input a set of maximum rates of climb (here defined as Vcmax (maximum climb speeds)) evaluated at different altitudes 
# and the set of altitudes.  
# The function output is:
#   -SC_theo = the theoretical maximum altitude, considering fixed weight, at which the maximum climb speed (Vcmax_theo) is zero.
#    In other words is the maximum altitude that the helicopter can reach in theory.
#   -SC_prac = the practical maximum altitude, considering fixed weight, at which the maximum climb speed (Vcmax_prac) is 100 ft/min. 
#
# 
# Documentation: "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - page 90-91.
# 
# Author: Alessio Ferrara
# Last change: 02/10/2024


import numpy as np


def ServiceCeiling(Vcmax,h):

    Vcmax_theo = 0                                              # Climb speed at theoretical maximum altitude (ft/min)
    Vcmax_prac = 100                                            # Climb speed at practical maximum altitude (ft/min)


    index_theo = np.argmin(np.abs(Vcmax - Vcmax_theo))          # Computation of the position in the climb speed array of the minimum absolute value of the
                                                                # difference between the maximum climb speed and the theoretical maximum climb speed

    index_prac = np.argmin(np.abs(Vcmax - Vcmax_prac))          # Computation of the position in the climb speed array of the minimum absolute value of the
                                                                # difference between the maximum climb speed and the practical maximum climb speed


    SC_theo = h[index_theo]                                     # Extraction of the theoretical maximum altitude from altitudes array
    SC_prac = h[index_prac]                                     # Extraction of the practical maximum altitude from altitudes array

    
    
    return SC_theo, SC_prac                                     # Output: SC_theo = (theoretical maximum altitude), SC_prac = (practical maximum altitude)

    