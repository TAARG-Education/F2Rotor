# Organization: Universita' degli Studi di Napoli - Federico II, Dipartimento di Ingegneria Industriale, Ingegneria Aerospaziale
# Course: Aerodinamica dell'ala rotante
# Professor: Renato Tognaccini
# Supervisor: Ettore Saetta
# Academic Year: 2023-2024
# Authors: Afflitto Fernando, Bruno Fabio


# Test Cases for the following functions, defined in the file "Helicopter_Power_Advance.py" :

# - w_tilde_function
# - P_i_tilde_function
# - P_c_0_function
# - P_c_Fus_function



import numpy as np
import matplotlib.pyplot as plt
from Helicopter_Power_Advance import w_tilde_function, P_i_tilde_function, P_c_0_function, P_c_Fus_function


# w_tilde_function

# Brief description of the function: this function solves the fourth grade's equation that describes the operating's curve for a stiff rotor
# in advance in the hypotesis of constant thrust. The function selects the smallest real and positive root of the equation for every couple
# (V_inf, alpha).

# Test case: the function w_tilde_function is evaluated in function of a vector of values of asymptotic velocity V_inf_tilde (in a particular 
# non-dimensional form, with respect to the induced velocity on the disk rotor in hovering) and some values of AoA, given in input. This test case
# aims to reproduce the graphic of Figure 7.1: Curve di funzionamento w_tilde(V_inf_tilde) a spinta costante per eliche in flusso non assiale,
# at pp.99 in "Lezioni di Aerodinamica dell'ala rotante".

# Input:
# - V_inf_tilde: uniform interval between 0 to 10
# - alpha_deg: four values of AoA between 0 deg and 90 deg (0, 30, 60, 90 deg)

# Output:
# - Curves which describe w_tilde in function of V_inf_tilde for some values of alpha_deg

# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
#  Chapter 7: Il rotore in volo traslato, Subsection 7.1.1: Funzionamento a spinta costante, pp. 98-99.
# - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
#  Figure 7.1: Curve di funzionamento w_tilde(V_inf_tilde) a spinta costante per eliche in flusso non assiale, pp.99.

# Authors: Afflitto Fernando, Bruno Fabio
# Last change: 27/07/2024


N = 300                                                    # Number of points in the interval for V_inf_tilde
M = 4                                                      # Number of points in the interval for alpha_deg

V_inf_tilde_values = np.linspace(0,10,N)                   # Uniform interval for V_inf_tilde
alpha_deg_values = np.linspace(0,90,M)                     # Uniform interval for alpha_deg

w_tilde = np.zeros((N,M))                                  # Inizialization for the matrix of w_tilde

for j in range(M):                                         # Loop for alpha_deg

    alpha = alpha_deg_values[j]

    for i in range(N):                                     # Loop for V_inf_tilde

        V_inf_tilde = V_inf_tilde_values[i]

        w_tilde[i,j] = w_tilde_function(V_inf_tilde,alpha)  # Call to function

plt.figure(1) # Graphic
plt.plot(V_inf_tilde_values, w_tilde[:,0], label=r'$\alpha = 0\ deg$')
plt.plot(V_inf_tilde_values, w_tilde[:,1], label =r'$\alpha = 30\ deg$')
plt.plot(V_inf_tilde_values,w_tilde[:,2], label=r'$\alpha = 60\ deg$')
plt.plot(V_inf_tilde_values,w_tilde[:,3], label=r'$\alpha = 90\ deg$')
plt.xlabel(r'$\tilde{V}_{\infty}$')
plt.ylabel(r'$\tilde{w}$')
plt.legend()
plt.grid()


# P_i_tilde_function

# Brief description of the function: this function computes the induced power of a stiff rotor in advance in a non-dimensional form with respect
# to the induced power rotor in hovering, in the hypotesis of a constant thrust. This function calls the function w_tilde_function.

# Test case: the function P_i_tilde_function is evaluated in function of a vector of values for asymptotic velocity V_inf_tilde (in a particular
# non-dimensional form, with respect to the induced velocity on the disk rotor in hovering) and ten values of AoA (from 0 to 90 degrees), given
# in input. This test case aims to reproduce the graphic in Figure 7.2: Curve P_i_tilde(V_inf_tilde) a spinta costante per eliche in flusso non
# assiale, at pp. 100 in "Lezioni Aerodinamica dell'ala rotante".

# Input:
# - V_inf_tilde: uniform interval between 0 and 10
# - alpha_deg: uniform interval between 0 and 90 (ten values)

# Output:
# - Curves which describe P_i_tilde in function of V_inf_tilde for some values of alpha_deg

# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
#  Chapter 7: Il rotore rigido in volo traslato, Subsection 7.1.1: Funzionamento a spinta costante, pp. 98-99.
# - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
#  Figure 7.2: Curve P_i_tilde(V_inf_tilde) a spinta costante per eliche in flusso non assiale, pp.100.

# Authors: Afflitto Fernando, Bruno Fabio
# Last change: 27/07/2024


N = 300                                                          # Number of points of interval for V_inf_tilde
M = 10                                                           # Number of points of interval for alpha_deg

alpha_deg_values = np.linspace(0,90,M)                           # Uniform interval for alpha_deg
V_inf_tilde_values = np.linspace(0,10,N)                         # Uniform interval for V_inf_tilde

P_i_tilde = np.zeros((N,M))                                      # Inizialization for the matrix of P_i_tilde

for j in range(M):                                               # Loop for alpha_deg

    alpha = alpha_deg_values[j]

    for i in range(N):                                           # Loop for V_inf_tilde

        V_inf_tilde = V_inf_tilde_values[i]

        P_i_tilde[i,j] = P_i_tilde_function(V_inf_tilde,alpha)  # Call to function P_i_tilde_function

plt.figure(2) # Graphic
plt.plot(V_inf_tilde_values,P_i_tilde[:,0],label=r'$\alpha = 0\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,1],label=r'$\alpha = 10\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,2],label=r'$\alpha = 20\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,3],label=r'$\alpha = 30\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,4],label=r'$\alpha = 40\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,5],label=r'$\alpha = 50\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,6],label=r'$\alpha = 60\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,7],label=r'$\alpha = 70\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,8],label=r'$\alpha = 80\ deg$')
plt.plot(V_inf_tilde_values,P_i_tilde[:,9],label=r'$\alpha = 90\ deg$')
plt.xlabel(r'$\tilde{V}_{\infty}$')
plt.ylabel(r'$\tilde{P}_{i}$')
plt.legend()
plt.grid()


# P_c_0_function

# Brief description of the function: this function computes the parasite power's coefficient absorbed by a stiff rotor in advance. In the
# hypotesis of rectangular form for the blade and drag coefficient of a single blade element constant along the radius and equal to a mean value,
# one can obtain the non-dimensional form of the parasite power absorbed by a stiff rotor in advance, which is a parabolic formula in function
# of mu. mu is the advance ratio, defined as V_inf*cos(alpha)/(Omega*R). By a factor K > 3 and 4 < K < 5, one can consider the effect on parasite
# power's coefficient of the velocity's component of flux along the direction of radius (K = 4.7 is the default value in this function).

# Test case: the function P_c_Fus_function is evaluated in function of a vector of values for mu, for some different values 
# of sigma and Cd_mean. A parabolic function of mu is obtained, as excepted.

# Input:
# - sigma: 0.10, 0.15, 0.20
# - Cd_mean: 0.02, 0.04, 0.06
# - K: 4.7
# - mu: Uniform interval between 0 and 2

# Output:
# - Curves which describe P_c_0 in function of mu for some values of Cd_mean and sigma = 0.15
# - Curves which describe P_c_0 in function of mu for some values of sigma and Cd_mean = 0.04

# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof.Renato Tognaccini, a.a. 2023-2024, vers. 2.05,
# Chapter 7: Il rotore rigido in volo traslato, Section 7.3: Potenza parassita in volo traslato, pp.100-102.

# Authors: Afflitto Fernando, Bruno Fabio
# Last change: 27/07/2024


N = 50                                      # Number of points of interval for mu
sigma_values = [0.10, 0.15, 0.20]           # Values of Solidity
M = len(sigma_values)                       # Length of solidity's vector
Cd_mean_values = [0.02, 0.04, 0.06]         # Values of Cd
P = len(Cd_mean_values)                     # Lenght of mean drag coefficients' vector
mu = np.linspace(0,2,N)                     # Uniform interval for mu

P_c_0 = np.zeros((N,3,3))                   # Inizialization for matrix of parasite power's coefficients

for k in range(P):                          # Loop for Cd_mean
    Cd_mean = Cd_mean_values[k]

    for j in range(M):                      # Loop for sigma
     sigma = sigma_values[j]

     P_c_0[:,j,k] = P_c_0_function(sigma,Cd_mean,mu)        # Call to function


plt.figure(3)                                                # Graphic
plt.plot(mu,P_c_0[:,1,0],label=r'$\bar{C_{d}}\ =\ 0.02$')
plt.plot(mu,P_c_0[:,1,1],label=r'$\bar{C_{d}}\ =\ 0.04$')
plt.plot(mu,P_c_0[:,1,2],label=r'$\bar{C_{d}}\ =\ 0.06$')
plt.xlabel(r'$\mu = V_{\infty}\cos(\alpha)/(\Omega R)$')
plt.ylabel(r'$P_{c,0}$')
plt.legend()
plt.title(r'$\sigma\ =\ 0.15$')
plt.grid()

plt.figure(4)                                                # Graphic
plt.plot(mu,P_c_0[:,0,1],label=r'$\sigma\ =\ 0.10$')
plt.plot(mu,P_c_0[:,1,1],label=r'$\sigma\ =\ 0.15$')
plt.plot(mu,P_c_0[:,2,1],label=r'$\sigma\ =\ 0.20$')
plt.xlabel(r'$\mu = V_{\infty}\cos(\alpha)/(\Omega R)$')
plt.ylabel(r'$P_{c,0}$')
plt.legend()
plt.title(r'$\bar{C_{d}}\ =\ 0.04$')
plt.grid()


# P_c_Fus_function

# Brief description of the function: this function computes the parasite power's coeffcient absorbed by a fuselage of an helicopter in advance.
# Starting from the classic relation for the power, function of density rho, cube of asymptotic velocity V_inf and equivalent wetted area f,
# one can obtain a non-dimensional form with respect to a reference power. This non dimensional form is function of f/A and the cube of mu,
# where f is the equivalent wetted area, A is the disk area and mu is defined as V_inf/(Omega*R). In this function the default value of
# f/A is 0.009.

# Test case: the function P_c_Fus_function is evaluated in function of a vector of values of mu, for some values of f/A
# (the first is the default value 0.009). A cubic function of mu is obtained, as excepted.

# Input:
# - mu: Uniform interval between 0 and 2
# - f_over_A: 0.009, 0.012, 0.015

# Output:
# - Curves which describe P_c_Fus in function of mu, for some values of f_over_A


# Documentation:
# - Lezioni di Aerodinamica dell'ala rotante, Prof. Renato Tognaccini, a.a. 2023/2024, vers.2.05,
#  Chapter 7: Il rotore rigido in volo traslato, Section 7.5: Potenza parassita della fusoliera, pp. 104-105.

# Authors: Afflitto Fernando, Bruno Fabio
# Last change: 27/07/2024


N = 50                                                  # Number of points of interval for mu
mu = np.linspace(0,2,N)                                 # Uniform interval for mu
f_over_A_values = [0.009, 0.012, 0.015]                 # Values for f_over_A
M = len(f_over_A_values)                                # Length of vector of values for f_over_A

P_c_Fus = np.zeros((N,M))                               # Inizialization for the matrix of parasite power's coefficients of fuselage

for i in range(M):                                      # Loop for f_over_A
    f_over_A = f_over_A_values[i]

    P_c_Fus[:,i] = P_c_Fus_function(mu,f_over_A)        # Call to function

plt.figure(5)                                           # Graphic
plt.plot(mu,P_c_Fus[:,0],label=r'$f/A\ =\ 0.009$')
plt.plot(mu,P_c_Fus[:,1],label=r'$f/A\ =\ 0.012$')
plt.plot(mu,P_c_Fus[:,2],label=r'$f/A\ =\ 0.015$')
plt.xlabel(r'$\mu = V_{\infty}/(\Omega R)$')
plt.ylabel(r'$P_{c,Fus}$')
plt.legend()
plt.grid()
plt.show()








         

        