import numpy as np


def BEMT_Hover(N, r_segn, solidity, theta, Clalpha, Cd0, mu):

    """
       Calculation of the coefficients Tc, Qc, Pc_i, Pc_0 and the merit figure FM.

       Args:
           N: Number of blades
           r_segn: Normalized vector of radial positions (r/R)
           solidity: Rotor solidity
           theta: Pitch angle distribution  [rad]
           Clalpha: Slope of the load-angle curve  [rad^-1]
           Cd0: Parasitic resistance coefficient
           mu: Progress report

       Returns:
           Tc: Total thrust coefficient
           Qc: Torque coefficient
           Pc_i: Induced power (coefficient)
           Pc_0: Parasitic power (coefficient)
           FM: Figure of merit
       """

    epsilon = 1e-6  # Avoid /0
    C1 = mu + Clalpha * solidity / 8
    C2 = -r_segn * Clalpha * solidity / 8 * (theta - mu / (r_segn + epsilon))

    # Ensure discriminant > 0
    discriminant = C1 ** 2 - 4 * C2
    discriminant[discriminant < 0] = 0  # Avoid complex values
    lambda_i = (-C1 + np.sqrt(discriminant)) / 2

    # Calculate phi e alpha
    phi = np.arctan(np.clip(mu / (r_segn + epsilon) + lambda_i / (r_segn + epsilon), -1e3, 1e3))
    alpha = np.clip(theta - phi, -np.radians(90), np.radians(90))

    # Prandtl factor
    F_Prandtl = np.where(lambda_i > 0, 2 / np.pi * np.arccos(np.exp(N / (2 * lambda_i) * (r_segn - 1))), 0)

    # Calculate Cl and Cd
    Cl = Clalpha * alpha
    Cd = Cd0

    # Thrust and power distributions
    dTcdr_segn = 0.5 * solidity * Cl * r_segn ** 2
    dTcdr_segn_PR = dTcdr_segn * F_Prandtl

    dPc_i = 0.5 * solidity * Cl * phi * r_segn ** 3
    dPc_0 = 0.5 * solidity * Cd * r_segn ** 3
    dPcdr_segn = dPc_i + dPc_0

    # Integrate with trapezoidal method
    Tc = float(np.trapz(dTcdr_segn_PR, r_segn))
    Pc_i = float(np.trapz(dPc_i, r_segn))
    Pc_0 = float(np.trapz(dPc_0, r_segn))
    Qc = float(np.trapz(dPcdr_segn, r_segn))

    # Calculate figure of merit (FM)
    FM = (Tc ** (3 / 2) / np.sqrt(2)) / Qc if abs(Qc) > epsilon else np.nan

    return Tc, Qc, Pc_i, Pc_0, FM
