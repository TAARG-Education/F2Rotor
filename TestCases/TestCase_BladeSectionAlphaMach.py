"""

    Test case for the BladeSection_alpha_Mach function.

    This test case verifies the function's ability to compute aerodynamic properties,
    including angle of attack, velocity components, Mach number, and stall regions,
    for different radial and azimuthal positions along a rotor blade.

    Author: Andrea Malafronte
    Rotary Wing Aerodynamics course, prof. Renato Tognaccini
    University of Naples Federico II
    Academic Year 2023-2024

"""

import numpy as np
from bladeSection_alpha_mach import BladeSection_alpha_Mach

# Define test parameters for the function
Omega = 30.0  # Rotor rotational speed in rad/s
R = 10.0  # Rotor radius in meters
lambda_ = 0.1  # Inflow ratio (non-dimensional)
r_segn = np.linspace(0.2, 1.0, 5)  # Radial positions along the blade (non-zero range)
beta = np.radians([10.0, 12.0, 11.0, 9.0, 10.0, 11.0, 9.0, 12.0])  # Array with the same length as psi
dbeta = np.radians([0.5, 0.4, 0.3, 0.2, 0.1, 0.3, 0.2, 0.4])  # Array with the same length as psi
mu = 0.2  # Advance ratio (non-dimensional)
psi = np.linspace(0, 2 * np.pi, 8)  # Azimuthal angles in radians
theta = np.radians([15.0, 14.5, 15.5, 14.0, 15.0])  # Pitch angles along the blade in radians
alpha_stall_up = 15.0  # Upper stall angle in degrees
alpha_stall_lo = -15.0  # Lower stall angle in degrees
h = 1000.0  # Altitude in meters

# Execute BladeSection_alpha_Mach function with the test parameters
alpha_e, phi, u_P, u_T, u_R, V_eff, non_stall, M = BladeSection_alpha_Mach(
    Omega, R, lambda_, r_segn, beta, dbeta, mu, psi, theta, alpha_stall_up, alpha_stall_lo, h
)

# Print some test results to verify function output
print("Effective angle of attack (alpha_e):")
print(alpha_e)

print("\nIncidence angle (phi):")
print(phi)

print("\nNormal velocity component to the disk plane (u_P):")
print(u_P)

print("\nTangential velocity component to the disk plane (u_T):")
print(u_T)

print("\nRadial velocity component (u_R):")
print(u_R)

print("\nEffective velocity (V_eff):")
print(V_eff)

print("\nNon-stalled regions (non_stall):")
print(non_stall)

print("\nMach number (M):")
print(M)
