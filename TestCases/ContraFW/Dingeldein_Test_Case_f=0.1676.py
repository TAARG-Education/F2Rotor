from ambiance import Atmosphere
import numpy as np
from sympy import symbols, Eq, solve, atan, sin, cos
from scipy.constants import g, pi
 
def ContraFW1(V_inf, N_contra, R, M, w_h, height, f, Omega, Pci_h, Pc0_h):
    """
    OUTPUT:
        - P_old: Lista delle potenze totali richieste per ciascuna velocità di volo
        - P_i_old: Lista delle potenze indotte
        - P_0_old: Lista delle potenze di profilo
        - P_fus_old: Lista delle potenze di fusoliera
    INPUT:
        - V_inf: Array o lista delle velocità di volo
        - N_contra: Numero di rotori
        - R: Raggio del rotore
        - M: Massa del veicolo
        - w_h: Velocità indotta all'hover
        - height: Altitudine
        - f: Area equivalente
        - Omega: Velocità angolare del rotore
        - Pci_h: Coefficiente di potenza indotta all'hover
        - Pc0_h: Coefficiente di potenza di profilo all'hover
    """
   
    # Calcolo dell'area del disco del rotore
    A = pi * R**2
 
  # Function to obtain air density based on altitude
    def air_density_at_altitude(h):
        """
        Calculate air density at a given altitude using the ambiance library
        (based on the International Standard Atmosphere - ISA model).
        
        """
        if h < 0:
            raise ValueError("Altitude cannot be negative.")
        # Use the ambiance library to calculate air density.
        atmosphere = Atmosphere(h)
        return atmosphere.density[0]  # Estrae la densità in kg/m^3
 
    rho_inf = air_density_at_altitude(height)
 
    # Calcolo del peso per rotore
    W = M * g / N_contra
 
    # Inizializzazione delle liste per i risultati
    P_old = []
    P_i_old = []
    P_0_old = []
    P_fus_old = []
 
    # Variabile simbolica per risolvere l'equazione per w
    w = symbols('w', real=True, positive=True)
 
    # Ciclo su tutte le velocità di volo in V_inf
    for Vinf in V_inf:
        # Calcolo della forza di trascinamento della fusoliera
        D_fus = 0.5 * rho_inf * Vinf**2 * f
        # Calcolo dell'angolo d'attacco del rotore
        alpha = atan((D_fus / N_contra) / W)
 
        # Definizione dell'equazione simbolica per la velocità indotta w
        eqn1 = Eq((Vinf / w_h * w / w_h * sin(alpha) + (w / w_h)**2)**2 +
                  (Vinf / w_h)**2 * (w / w_h)**2 * (cos(alpha))**2 - 1, 0)
       
        # Risoluzione dell'equazione per w e selezione della soluzione reale e positiva
        sol1 = solve(eqn1, w)
        w_value = None
        for sol in sol1:
            if sol.is_real and sol > 0:
                w_value = float(sol)
                break
       
        # Calcolo delle potenze richieste
        if w_value is not None:
            P_i_old_value = (w_value / w_h) * Pci_h * rho_inf * Omega**3 * R**3 * A
            P_i_old_value += 0.16*P_i_old_value
            P_0_old_value = Pc0_h * rho_inf * A * Omega * R * (Omega**2 * R**2 + 4.7 * Vinf**2)
            P_fus_old_value = 0.5 * rho_inf * Vinf**3 * f
 
            # Calcolo della potenza totale e salvataggio dei risultati
            P_total = P_i_old_value + P_0_old_value + P_fus_old_value / N_contra
            P_old.append(P_total)
            P_i_old.append(P_i_old_value)
            P_0_old.append(P_0_old_value)
            P_fus_old.append(P_fus_old_value)
        else:
            # Se non c'è soluzione valida per w
            P_old.append(None)
            P_i_old.append(None)
            P_0_old.append(None)
            P_fus_old.append(None)
 
    return P_old, P_i_old, P_0_old, P_fus_old
 
# Input Data for Dingeldein Coaxial-rotor:
V_inf = [0, 18, 21, 24, 27.4, 32]     # Flight Speed (m/s)
N_contra = 2                          # Rotors Number
R = 3.81                              # Rotor Radius (m)
M = 650                               # Aircraft Weight (kg)
w_h = 5.34                            # Hovering Induced Speed (m/s)
height = 0                            # Height (m)
f = 0.1676                            # Equivalent Wet Area (m^2)
Omega = 40.0                          # Rotor Angular Speed (rad/s)
Pci_h = 0.00017226                    # Hovering Induced Power Coefficient        
Pc0_h = 0.00011475                    # Hovering Profile Power Coefficient
 
# Executing the function with the Kamov-Ka226 input data
P_old, P_i_old, P_0_old, P_fus_old = ContraFW1(V_inf, N_contra, R, M, w_h, height, f, Omega, Pci_h, Pc0_h)
 
# Printing of the results for verification
print("Total Required Power (P_old):", P_old)
print("Induced Power (P_i_old):", P_i_old)
print("Profile Power (P_0_old):", P_0_old)
print("Fuselage Power (P_fus_old):", P_fus_old) 