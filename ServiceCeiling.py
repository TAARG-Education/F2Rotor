# Organization: Universita' degli Studi di Napoli - Federico II, Dipartimento di Ingegneria Industriale, Ingegneria Aerospaziale
# Course: Aerodinamica dell'ala rotante
# Professor: Renato Tognaccini
# Supervisor: Ing. Ettore Saetta, Ing. Michele Massa
# Academic Year: 2023-2024
#
# This function computes the theoretical maximum altitude and the practical maximum altitude (or ceilings) of an helicopter in advance 
# (i.e. in forward flight). It computes the absolute value of the minimum difference between a set of maximum rates of climb, computed 
# at different altitudes, and the rate of climb defined at the theoretical or the practical altitude. When this condition is found, 
# the function computes the altitude corresponding to that rate of climb value. 
# The maximum rate of climb is defined as the maximum difference between the available power curve and the necessary power curve, divided by 
# the weight of the helicopter. 
# So basically, the theoretical maximum altitude or the practical maximum altitude are the altitudes in which the rates of climb are 0ft/min 
# or 100ft/min (0m/s or 0,508m/s), as definition, or the altitudes at which the available power curve is tangent to the minimum value of the 
# necessary power curve. 
#
# The function takes in input a set of maximum rates of climb (here defined as Vcmax (maximum climb speeds)) evaluated at different altitudes 
# and the set of altitudes.  
# The function output is:
#   -SC_theo = the theoretical maximum altitude, considering fixed weight, at which the maximum climb speed (Vcmax_theo) is 0 ft/min.
#    In other words is the maximum altitude that the helicopter can reach in theory.
#   -SC_prac = the practical maximum altitude, considering fixed weight, at which the maximum climb speed (Vcmax_prac) is 100 ft/min. 
#
# 
# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
# Chapter 7: Il rotore rigido in volo traslato, pp.95-103.
# Giovanni Di Giorgio - "Lezioni integrative dell'insegnamento di Aerodinamica dell'Ala Rotante" - a.a. 2023/2024 - page 89-91.
#
#  
# Author: Alessio Ferrara
# Last change: 07/10/2024


import numpy as np


def ServiceCeiling_function(Vcmax,h):

    Vcmax_theo = 0                                              # Climb speed at theoretical maximum altitude (m/s)
    Vcmax_prac = 0.508                                          # Climb speed at practical maximum altitude (m/s)


    index_theo = np.argmin(np.abs(Vcmax - Vcmax_theo))          # Computation of the position in the climb speed array of the minimum absolute value of the
                                                                # difference between the maximum climb speed and the theoretical maximum climb speed

    index_prac = np.argmin(np.abs(Vcmax - Vcmax_prac))          # Computation of the position in the climb speed array of the minimum absolute value of the
                                                                # difference between the maximum climb speed and the practical maximum climb speed


    SC_theo = h[index_theo]                                     # Extraction of the theoretical maximum altitude from altitudes array
    SC_prac = h[index_prac]                                     # Extraction of the practical maximum altitude from altitudes array

    
    
    return SC_theo, SC_prac                                     # Output: SC_theo = (theoretical maximum altitude), SC_prac = (practical maximum altitude)

    