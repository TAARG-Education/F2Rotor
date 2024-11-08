"""

    Estimation of aerodynamic parameters for rotor blade sections in forward flight based on radial and azimuthal positions.
    This function calculates the effective angle of attack, velocity components, Mach number, and stall regions across the blade, considering:
    - Rotor rotational speed Omega and radius R
    - Inflow ratio lambda_ and advance ratio mu
    - Arrays of azimuthal angles psi, radial positions r_segn, and pitch angles theta
    - Upper and lower stall angles (alpha_stall_up, alpha_stall_lo)
    - Collective pitch angles beta and variations in pitch angle dbeta with respect to azimuth
    - Atmospheric properties at a given altitude h

    Procedure:
        1. Computes velocity components (normal, tangential, radial) and incidence angles across blade sections.
        2. Calculates the effective angle of attack and the resulting effective velocity at each section.
        3. Identifies stalled and non-stalled blade regions based on stall angles.
        4. Calculates Mach numbers at each blade section, factoring in atmospheric speed of sound based on altitude h.
    Outputs:
    - Effective angle of attack (alpha_e) and incidence angle (phi) for each section
    - Velocity components u_P (normal), u_T (tangential), u_R (radial), and effective velocity V_eff
    - Mach number M and non-stalled regions (non_stall), with stalled sections marked as NaN.

    References:
    - Carlo de Nicola, (2018-2019), Notes on Aircraft Aerodynamics, Appendix E: http://wpage.unina.it/denicola/AdA/DOWNLOAD/Appunti_AdA_2018_2019.pdf
    - Renato Tognaccini, (2023-2024), Rotary Wing Aerodynamics Lectures, section 7.5

    Author: Andrea Malafronte
    Latest Update: 08/11/2024
    Version: 1.0

"""

import numpy as np
from ambiance import Atmosphere  # Importing ambiance package

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
    # Use ambiance to calculate atmospheric properties at the given altitude
    atmosphere = Atmosphere(h)
    T = atmosphere.temperature[0]  # Temperature in Kelvin
    a_inf = atmosphere.speed_of_sound[0]  # Speed of sound in m/s
    rho = atmosphere.density[0]  # Air density in kg/m^3
    p = atmosphere.pressure[0]  # Pressure in Pa

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
