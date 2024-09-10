# FLAPPING COEFFICIENTS FOR ROTOR BLADES
#
# Description: This function computes the flapping coefficients (beta_0, beta_1c, beta_1s, and beta) 
#              for rotor blades based on the given input parameters including angular velocity (omega), 
#              Lock number (gamma), collective pitch angle (theta_0), blade twist angle (theta_tw), advance 
#              ratio (mu) and inflow ratio (lamb).
#              The flapping equations used are from standard rotor dynamics.
#              The study was conducted with the following assumptions: the rotor angular velocity (omega) is 
#              constant, the eccentricity of the flapping hinge (e = 0) is neglected, it is assumed that there 
#              is no feathering hinge, and that the asymptotic velocity is constant and only slightly inclined 
#              with respect to the control plane. 
#              The blade motion is assumed to be described with respect to the control plane, which is the 
#              plane perpendicular to the direction of thrust.
#
# References: "Lezioni di Aerodinamica dell'Ala Rotante", Renato Tognaccini, 
#             Chapter 8  Section 8, Calculations of Flapping Coefficients.
#
# Authors: El Mehdi Rahmaoui
# Rotary Wing Aerodynamics Course
# University of Naples Federico II
# Academic Year 2023-2024
#
# 
#
# Version: 1.0.0

import numpy as np
import matplotlib.pyplot as plt


def compute_beta_coefficients(gamma, theta_0, theta_tw, mu, lamb, omega, t):
    """
    Compute the flapping coefficients for the rotor blade based on time (t).
    The azimuth angle psi = omega * t is used for calculating the overall flapping coefficient beta.
    
    Parameters:
    - gamma: Lock number (measure of the relative importance between the aerodynamic and inertial forces acting on the blade): 
      For an articulated rotor, the Lock number  typically ranges from 8 to 10. For a rotor without hinges, gamma typically ranges from 5 to 7.
    - theta_0: Collective pitch angle (radians)
    - theta_tw: Blade twist angle (radians)
    - mu: Advance ratio
    - lamb: Inflow ratio
    - omega: Angular velocity (rad/s)
    - t: Time (seconds)
    
    Returns:
    - beta_0: Conicity coefficient
    - beta_1c: Longitudinal flapping coefficient
    - beta_1s: Lateral flapping coefficient
    - beta: Overall flapping coefficient
    """
    psi = omega * t  # Azimuth angle psi (radians)

    # Calculate beta_0: Conicity coefficient
    beta_0 = gamma * ((theta_0 / 8) * (1 + mu ** 2) + 
                      (theta_tw / 10) * (1 + (5 / 6) * mu ** 2) - 
                      (lamb / 6))

    # Calculate beta_1c: Longitudinal flapping coefficient
    beta_1c = -2 * mu * ((4 / 3) * theta_0 + theta_tw - lamb) / (1 - mu ** 2 / 2)

    # Calculate beta_1s: Lateral flapping coefficient
    beta_1s = -(4 / 3) * mu * beta_0 / (1 + mu ** 2 / 2)

    # Calculate overall beta as a function of time (psi = omega * t)
    beta = beta_0 + beta_1s * np.cos(psi) + beta_1c * np.sin(psi)

    return beta_0, beta_1c, beta_1s, beta
