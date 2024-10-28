"""
Rotor Blade Aerodynamic Analysis Module

This module provides functions to compute atmospheric properties and aerodynamic
characteristics of a rotor blade at different radial and azimuthal positions. It includes:

1. `atmosisa(h)`: Computes atmospheric properties at a given altitude.
2. `BladeSection_alpha_Mach(...)`: Computes the angle of attack, effective velocity,
   Mach number, and other parameters of a blade section at various positions.

Dependencies:
- numpy
- scipy.constants
"""

import numpy as np
from scipy.constants import convert_temperature
from scipy.constants import atmosphere
def atmosisa(h):
    """
    Computes atmospheric properties based on altitude, assuming International Standard Atmosphere (ISA) conditions.

    Parameters:
        h (float): Altitude in meters.

    Returns:
        T (float): Temperature at altitude in Kelvin.
        a_inf (float): Speed of sound at altitude in m/s.
        rho (float): Air density at altitude in kg/m^3.
        p (float): Atmospheric pressure at altitude in kPa.
    """
    # Calculate temperature and pressure based on altitude
    if h < 11000:
        T = 15.04 - 0.00649 * h
        p = 101.29 * ((T + 273.1) / 288.08) ** 5.256
    else:
        T = -56.46
        p = 22.65 * np.exp(1.73 - 0.000157 * h)

    T += 273.15  # Convert from degrees Celsius to Kelvin
    rho = p / (0.2869 * T)  # Calculate air density
    a_inf = np.sqrt(1.4 * 287.05 * T)  # Calculate speed of sound in air
    return T, a_inf, rho, p


def BladeSection_alpha_Mach(Omega, R, lambda_, r_segn, beta, dbeta, mu, psi, theta, alpha_stall_up, alpha_stall_lo, h):
    """
    Computes the aerodynamic characteristics of rotor blade sections at different radial and azimuthal positions.

    Parameters:
        Omega (float): Rotational speed of the rotor in rad/s.
        R (float): Rotor radius in meters.
        lambda_ (float): Inflow ratio (non-dimensional).
        r_segn (array): Array of radial positions along the blade.
        beta (array): Array of collective pitch angles at different azimuth angles in radians.
        dbeta (array): Variation in pitch angle with respect to azimuth.
        mu (float): Advance ratio (non-dimensional).
        psi (array): Array of azimuthal angles in radians.
        theta (array): Array of pitch angles at different radial positions in radians.
        alpha_stall_up (float): Upper stall angle in degrees.
        alpha_stall_lo (float): Lower stall angle in degrees.
        h (float): Altitude in meters.

    Returns:
        alpha_e (2D array): Effective angle of attack at each blade section.
        phi (2D array): Incidence angle at each blade section.
        u_P (2D array): Normal velocity component to the disk plane at each section.
        u_T (2D array): Tangential velocity component to the disk plane at each section.
        u_R (2D array): Radial velocity component at each section.
        V_eff (2D array): Effective velocity at each blade section.
        non_stall (2D array): Blade regions not in stall (NaN for stalled sections).
        M (2D array): Mach number at each blade section.
        M_cr (2D array): Placeholder array for critical Mach numbers (not calculated).
    """
    # Mesh grid creation for radial and azimuthal coordinates
    r_2d, psi_2d = np.meshgrid(np.linspace(0, r_segn[0], len(r_segn)), psi)
    x_hub = r_2d * np.cos(psi_2d - np.pi / 2)
    y_hub = r_2d * np.sin(psi_2d - np.pi / 2)

    r_2d, psi_2d = np.meshgrid(r_segn, psi)
    x = r_2d * np.cos(psi_2d - np.pi / 2)
    y = r_2d * np.sin(psi_2d - np.pi / 2)

    # Initialize variables for velocity components, incidence angle, and angle of attack
    u_P = np.zeros((len(r_segn), len(psi)))  # Normal velocity component to the disk plane
    u_T = np.zeros((len(r_segn), len(psi)))  # Tangential velocity component to the disk plane
    u_R = np.zeros((len(r_segn), len(psi)))  # Radial velocity component
    phi = np.zeros((len(r_segn), len(psi)))  # Incidence angle
    alpha_e = np.zeros((len(r_segn), len(psi)))  # Effective angle of attack
    V_eff = np.zeros((len(r_segn), len(psi)))  # Effective velocity

    _, a_inf, _, _ = atmosisa(h)

    # Compute velocity components and angles for each blade section
    for i in range(len(psi)):
        for j in range(len(r_segn)):
            u_P[j, i] = lambda_ + r_segn[j] * dbeta[i] + beta[i] * mu * np.cos(psi[i])
            u_T[j, i] = r_segn[j] + mu * np.sin(psi[i])
            u_R[j, i] = mu * np.cos(psi[i])
            phi[j, i] = np.arctan(u_P[j, i] / u_T[j, i])
            alpha_e[j, i] = theta[j] - phi[j, i]
            V_eff[j, i] = Omega * R * np.sqrt(u_P[j, i] ** 2 + u_T[j, i] ** 2)

    # Initialize arrays for stall region, non-stall region, and Mach number
    stall = np.zeros((len(r_segn), len(psi)))  # Rotor blade's stalled region
    non_stall = np.full((len(r_segn), len(psi)), np.nan)  # Non-stalled region of the blade
    M = np.full((len(r_segn), len(psi)), np.nan)
    M_cr = np.full((len(r_segn), len(psi)), np.nan)

    # Determine stalled and non-stalled sections, calculate Mach numbers
    for iii in range(len(r_segn)):
        for jjj in range(len(psi)):
            if alpha_e[iii, jjj] - np.radians(alpha_stall_up) > 0 or alpha_e[iii, jjj] - np.radians(alpha_stall_lo) < 0:
                stall[iii, jjj] = alpha_e[iii, jjj]
            else:
                non_stall[iii, jjj] = alpha_e[iii, jjj]
            M[iii, jjj] = V_eff[iii, jjj] / a_inf

    return alpha_e, phi, u_P, u_T, u_R, V_eff, non_stall, M, M_cr
