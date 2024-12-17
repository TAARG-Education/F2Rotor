# DRONE IN FORWARD FLIGHT F2ROTOR
#
# Description: this functions have the aim of determine the drone performance, in terms of power contributions, thrust
#              and torque coefficients. Specifically, the BEMT and the simple momentum theory are applied, and a
#              comparison of the results is carried out.
#
# References:  Renato Tognaccini - Lezioni per il corso di Aerodinamica dell'Ala Rotante - Eliche rotori ed aeromotori
#              con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04.
#
# Authors:     Pasquale De Riso, Samuel Filosa.
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini.
# University of Naples Federico II.
# Academic Year 2023-2024.
#
# Date: 11/28/2024.
#
# Version: 1.0.0.

import numpy as np
from math import atan, sin, pi
from ambiance import Atmosphere
from sympy import symbols, Eq, solve

def RotorSolving(r_segn, psi, Vinf, Omega_new, R, A, N, N_rotors, w, alpha, theta, Clalpha, Cd_mean, c, rho_inf, f, W,
                 F_Prandtl):
    """
            Calculation of Tc, Qc and power contributions, determining the dimensionless error
            in order to reach horizontal flight condition (T = W).

            Args:
                r_segn: Normalized vector of radial positions (r/R).
                psi: Azimuth angle [rad].
                Vinf: airspeed [m/s].
                Omega_new: Angular velocity calculated with bisection method [rad/s].
                R: Rotor radius [m].
                A: Rotor area [m^2].
                N: number of blades.
                N_rotors: number of rotors.
                w: induced velocity [m/s].
                alpha: Angle of attack [rad].
                theta: pitch angle [rad].
                Clalpha: Slope of the load-angle curve [rad^-1].
                Cd_mean: Parasitic resistance coefficient.
                c: Mean geometric chord [m].
                rho_inf: air density [Kg/m^3].
                f: equivalent wet area (% chord).
                W: Drone weight [N].
                F_Prandtl: Prandtl function.

            Returns:
                err: Dimensionless error.
                T: Thrust [N].
                Q: Torque [N].
                H: Rotor drag [N].
                Y: Cross force [N].
                P: Total power by impulsive theory [W].
                P_i: Induced power by impulsive theory [W].
                P_0: Parasitic power by impulsive theory [W].
                P_fus: Power absorbed by "fuselage" by impulsive theory [W].
                alpha_e: Effective angle of attack [rad].
                V_e: Effective velocity [m/s].
                dT_dr_dpsi: Derivative of thrust versus r_segn and psi.
                U_P: Perpendicular velocity component [m/s].
                dH_dpsi_0: Derivative of H versus psi, valued at 0.
                dQ_dpsi_0: Derivative of Q versus psi, valued at 0.
                dT_dr: Derivative of T versus r (thrust distribution).
                dQ_dr: Derivative of Q versus r (torque distribution).
    """

    dT_dr = np.zeros(len(r_segn))                                                       # Derivative of T versus radius.
    dQ_dr = np.zeros(len(r_segn))                                                       # Derivative of Q versus radius.

    psi_grid, R_grid = np.meshgrid(psi, r_segn)                                        # Make psi and r_segn arrays.

    mu = Vinf / (Omega_new * R)                                                             # Advance ratio, eq. (7.14).

    U_T = (Omega_new * R_grid * R + Vinf * np.cos(alpha) * np.sin(psi_grid))   # Parallel velocity component, eq. (8.5).
    U_P = Vinf * sin(alpha) + w                        # Perpendicular velocity component, eq. (8.5) with no flap [m/s].

    V_e = np.sqrt(U_T ** 2 + U_P ** 2)                                    # Effective velocity in the blade plane [m/s].
    phi = np.arctan(U_P / U_T)                                                  # Inflow angle in the blade plane [rad].

    alpha_e = theta[:, np.newaxis] - phi                                              # Effective angle of attack [rad].

    # If alpha_e is greater or minor than 15 deg, set it respectively to 15 deg and -15 deg,
    # and turn it in rad (stall control).
    alpha_e = np.clip(alpha_e, np.deg2rad(-15), np.deg2rad(15))

    Cl = Clalpha * alpha_e                                       # Lift coefficient calculation according linear theory.
    Cd = Cd_mean                                                                                     # Drag coefficient.

    # If chord is variable with r_segn use c[ii], else use c.
    chord = c[:] if np.ndim(c) > 0 else c
    dL = 0.5 * rho_inf * V_e ** 2 * chord[:, np.newaxis] * Cl                                       # Differential lift.
    dD = 0.5 * rho_inf * V_e ** 2 * chord[:, np.newaxis] * Cd                                       # Differential drag.

    Fx = dL * np.sin(phi) + dD * np.cos(phi)                                     # x component of the thrust, eq. (8.8).
    Fz = dL * np.cos(phi) - dD * np.sin(phi)                                     # z component of the thrust, eq. (8.8).
    Fr = 0                                                                                         # Flapping neglected.

    dT_dr_dpsi = Fz * F_Prandtl[:, np.newaxis]                      # Derivative of T versus r_segn and psi, eq. (8.11).
    dH_dr_dpsi = Fx * np.sin(psi_grid) + Fr * np.cos(psi_grid)      # Derivative of H versus r_segn and psi, eq. (8.11).
    dH_dr_dpsi_0 = dD * np.cos(phi) * np.sin(psi_grid)    # Derivative of H parasitic versus r_segn and psi, eq. (8.11).
    dY_dr_dpsi = -Fx * np.cos(psi_grid) + Fr * np.sin(psi_grid)     # Derivative of Y versus r_segn and psi, eq. (8.11).
    dQ_dr_dpsi = R_grid * R * Fx                                    # Derivative of Q versus r_segn and psi, eq. (8.11).
    dQ_dr_dpsi_0 = R_grid * R * dD * np.cos(phi)          # Derivative of Q parasitic versus r_segn and psi, eq. (8.11).

    # Integration along the dimensionless radius.
    dT_dpsi = np.trapz(dT_dr_dpsi, r_segn * R, axis = 0)                          # Thrust integration along the radius.
    dH_dpsi = np.trapz(dH_dr_dpsi, r_segn * R, axis = 0)                      # Rotor drag integration along the radius.
    dH_dpsi_0 = np.trapz(dH_dr_dpsi_0, r_segn * R, axis = 0)        # Rotor parasitic drag integration along the radius.
    dY_dpsi = np.trapz(dY_dr_dpsi, r_segn * R, axis = 0)                     # Cross force integration along the radius.
    dQ_dpsi = np.trapz(dQ_dr_dpsi, r_segn * R, axis = 0)                          # Torque integration along the radius.
    dQ_dpsi_0 = np.trapz(dQ_dr_dpsi_0, r_segn * R, axis = 0)            # Parasitic torque integration along the radius.

# Final computation looking at (8.11) formulas.
    T = N / (2 * pi) * np.trapz(dT_dpsi, psi)                                                              # Thrust [N].
    H = N / (2 * pi) * np.trapz(dH_dpsi, psi)                                                          # Rotor drag [N].
    H_0 = N / (2 * pi) * np.trapz(dH_dpsi_0, psi)                                            # Parasitic rotor drag [N].
    Y = N / (2 * pi) * np.trapz(dY_dpsi, psi)                                                         # Cross force [N].
    Q = N / (2 * pi) * np.trapz(dQ_dpsi, psi)                                                             # Torque [Nm].
    Q_0 = N / (2 * pi) * np.trapz(dQ_dpsi_0, psi)                                               # Parasitic torque [Nm].

    D_fus = f * 0.5 * rho_inf * Vinf ** 2                                                                    # Drag [N].

    Tc = T / (rho_inf * Omega_new ** 2 * R ** 2 * A)                                    # Thrust coefficient, eq. (2.5).
    Hc_0 = H_0 / (rho_inf * Omega_new ** 2 * R ** 2 * A)                             # Parasitic rotor drag coefficient.

    Qc = Q / (rho_inf * Omega_new ** 2 * R ** 3 * A)
    Qc_0 = Q_0 / (rho_inf * Omega_new ** 2 * R ** 3 * A)                                 # Parasitic torque coefficient.

    # Power coefficient contributions according to (8.35) equation.
    Pc_0 = Qc_0 + mu * Hc_0                                                                 # Profile power coefficient.
    Pc_fus = mu * D_fus / W * Tc                                                          # Parasitic power coefficient.
    Pc = Qc                                                                   # Single rotor required power coefficient.
    Pc_i = Pc - Pc_0 - Pc_fus/N_rotors

    # Power contributions.
    P_i = Pc_i * rho_inf * Omega_new ** 3 * R ** 3 * A                                              # Induced power [W].
    P_0 = Pc_0 * rho_inf * Omega_new ** 3 * R ** 3 * A                                              # Profile power [W].
    P_fus = Pc_fus * rho_inf * Omega_new ** 3 * R ** 3 * A                                        # Parasitic power [W].

    P = Q*Omega_new                                                                                   # Total power [W].

    err = (T - W)/W                                                                      # Defining dimensionless error.

    return (err, T, Q, H, Y, P, P_i, P_0, P_fus, alpha_e, V_e, dT_dr_dpsi, dQ_dr_dpsi, U_P, dH_dpsi_0, dQ_dpsi_0,
            dT_dr, dQ_dr
            )



def RigidFW2(
        Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A, P_ih, f, Pc0,
        Clalpha, Cd_mean, Omega, N_rotors, F_Prandtl):
    """
        Calculation of Tc, Qc and power contributions through a bisection method.
        Note that the comments refers to the equation contained in the following book:
        Renato Tognaccini. Lezioni per il corso di Aerodinamica dellAla Rotante. 2023.

        Args:
            Vinf_vec: Airspeed [m/s].
            N: Number of blades.
            r_segn: Normalized vector of radial positions (r/R).
            c: Mean geometric chord [m].
            theta: Pitch angle distribution [rad].
            W: Drone weight [N].
            w_h: Induced velocity in hover conditions [m/s].
            height: Altitude [m].
            R: Rotor radius [m].
            A: Rotor area [m^2].
            P_ih: Induced power in hovering [W].
            f: Equivalent wet area (% chord).
            Pc0: Power parasitic coefficient.
            Clalpha: Slope of the load-angle curve [rad^-1].
            Cd_mean: Parasitic resistance coefficient.
            Omega: Angular velocity [rad/s].
            N_rotors: Number of rotors.
            F_Prandtl: Prandtl function.

        Returns:
            P_i: Induced power by impulsive theory [W].
            P_0: Parasitic power by impulsive theory [W].
            P_fus: Power absorbed by "fuselage" by impulsive theory[W].
            P: Total power by impulsive theory [W].
            P_new: Total power by BEMT [W].
            P_new_i: Induced power by BEMT [W].
            P_new_0: Profile power by BEMT [W].
            P_new_fus: Parasitic power absorbed by "fuselage" by BEMT [W].
            alpha_vec: Angle of attack [rad].
            Omega_new: Angular velocity calculated with bisection method [rad/s].
            alpha_e_vec: Effective angle of attack [rad].
            Mach_vec: Mach number.
            dTc: Differential thrust coefficient.
            U_P_vec: U_P [m/s].
            dH_dpsi_0: Derivative of H versus psi, valued at 0.
            dQ_dpsi_0: Derivative of Q versus psi, valued at 0.
            T_vec: Thrust [N].
            Q_vec: Torque [N].
            Y_vec: Cross force [N].
            dT_dr_dpsi: Derivative of T versus r (thrust distribution).
            dQ_dr_dpsi: Derivative of Q versus r (torque distribution).
            psi: Azimuth angle [rad].
        """

    atmosphere = Atmosphere(height)
    rho_inf = atmosphere.density[0]                                                              # Air density [kg/m^3].

    psi = np.linspace(0, 2 * np.pi, 19)                                                      # Azimuth angle.

    # Initialization.
    Omega_new = np.zeros(len(Vinf_vec))                                                              # Angular velocity.
    alpha_vec = np.zeros(len(Vinf_vec))                                                               # Angle of attack.
    P_new_i = np.zeros(len(Vinf_vec))                                              # Induced power calculated with BEMT.
    P_new_0 = np.zeros(len(Vinf_vec))                                              # Profile power calculated with BEMT.
    P_new_fus = np.zeros(len(Vinf_vec))                                          # Parasitic power calculated with BEMT.
    P_new = np.zeros(len(Vinf_vec))                                                  # Total power calculated with BEMT.
    alpha_e_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                          # Effective angle of attack.
    dTc = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                            # Differential thrust coefficient.
    dT_dr_dpsi = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                     # Differential thrust coefficient.
    dQ_dr_dpsi = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                     # Differential thrust coefficient.
    U_P_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                       # Perpendicular velocity component.
    dT_dr = np.zeros((len(Vinf_vec), len(r_segn)))                                      # Derivative of T versus radius.
    dQ_dr = np.zeros((len(Vinf_vec), len(r_segn)))                                      # Derivative of Q versus radius.

    T_vec = np.zeros(len(Vinf_vec))                                                                # Thrust coefficient.
    Q_vec = np.zeros(len(Vinf_vec))                                                                # Torque coefficient.
    Y_vec = np.zeros(len(Vinf_vec))                                                                       # Cross force.

    P_i = np.zeros(len(Vinf_vec))                                                    # Induced power calculated with IT.
    P_0 = np.zeros(len(Vinf_vec))                                                    # Profile power calculated with IT.
    P_fus = np.zeros(len(Vinf_vec))                                                # Parasitic power calculated with IT.
    P = np.zeros(len(Vinf_vec))                                  # Single rotor total power required calculated with IT.
    P_tot = np.zeros(len(Vinf_vec))                                           # Total power required calculated with IT.


    for i, Vinf in enumerate(Vinf_vec):
        alpha = 0                                                                                 # Initial alpha value.
        err = 1                                                                                     # Initial err value.

        while err > 5e-6:

            w = symbols('w', real=True, positive=True)                                  # Symbolic solution in w.
            eqn = (
                    (Vinf / w_h * w / w_h * np.sin(alpha) + (w / w_h) ** 2) ** 2
                    + (Vinf / w_h) ** 2 * (w / w_h) ** 2 * (np.cos(alpha)) ** 2
                    - 1
            )                                                                             # Fourth grade equation (7.9).
            sol = solve(Eq(eqn, 0), w)                         # Fourth grade equation solution with RHS 0.

            for s in sol:
                if s.is_real and s > 0:                            # Take into account only real and positive solutions.
                    w = float(s)
                    break

            if alpha == 0:
                w_0 = w

            "Starting bisection method"
            Omega_a, Omega_b = 1, 10000                                                          # Initial Omega values.

            # Errors calculation with Omega_a and Omega_b.
            err_a, *_ = RotorSolving(r_segn, psi, Vinf, Omega_a, R, A, N,N_rotors, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, F_Prandtl)

            err_b, *_ = RotorSolving(r_segn, psi, Vinf, Omega_b, R, A, N,N_rotors, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, F_Prandtl)

            # Check the errors signs.
            if err_a * err_b > 0:
                raise ValueError(f"Redefine the search range: "                                          # Error signal.
                                 f"err_a: {err_a}, err_b: {err_b}, range: [{Omega_a}, {Omega_b}]")

            # Updating variables.
            it = 0                                                                                   # Starting counter.
            err_m = 1                                                                                      # Mean error.
            tol = 5e-6                                                                                      # Tolerance.

            while abs(err_m) > tol:
                it += 1                                                                           # Improve the counter.

                Omega_New = (Omega_a + Omega_b) / 2                                                  # Update Omega_New.
                result = RotorSolving(                                                   # RotorSolving function output.
                    r_segn, psi, Vinf, Omega_New, R, A, N,N_rotors, w, alpha, theta,
                    Clalpha, Cd_mean, c, rho_inf, f, W, F_Prandtl)

                (
                    err_m, T_v, Q_v, H, Y_v, P_New, P_New_i, P_New_0, P_New_fus, alpha_e, V_e, dT_dr_dpsi_temp,
                    dQ_dr_dpsi_temp, U_P, dH_dpsi_0, dQ_dpsi_0, dT, dQ) = result  # RotorSolving function output saving.

                # Bisection method instructions.
                if err_m * err_a < 0:
                    Omega_b = Omega_New                                                              # Updating Omega_b.
                    err_m = err_b                                                                      # Updating err_m.
                else:
                    Omega_a = Omega_New                                                              # Updating Omega_a.
                    err_a = err_m                                                                      # Updating err_a.

            D_fus = 0.5 * rho_inf * Vinf ** 2 * f                                                   # Fuselage drag [N].
            alpha_new = atan((D_fus / N_rotors + H) / W)                # Angle of attack calculation, eq. (8.33) [rad].

            err = abs((alpha_new - alpha) / (alpha_new+1e-6))                                          # Updating error.
            alpha = alpha_new                                                                    # Updating alpha [rad].

        # Power contributions calculation with BEMT.
        Omega_new[i] = Omega_New                                                                 # Omega update [rad/s].
        dT_dr[i,:] = dT                                                                           # Thrust distribution.
        dQ_dr[i,:] = dQ                                                                           # Torque distribution.
        P_new[i] = P_New                                                                              # Total power [W].
        P_new_i[i] = P_New_i                                                                        # Induced power [W].
        P_new_0[i] = P_New_0                                                                        # Profile power [W].
        P_new_fus[i] = P_New_fus                                                                  # Parasitic power [W].

        alpha_e_vec[i, :, :] = alpha_e                                                # Effective angle of attack [rad].
        dTc[i, :, :] = dT_dr_dpsi_temp / (rho_inf * Omega_new[i] ** 2 * R ** 2 * A)   # Differential thrust coefficient.
        dT_dr_dpsi[i, :, :] = dT_dr_dpsi_temp
        dQ_dr_dpsi[i, :, :] = dQ_dr_dpsi_temp

        U_P_vec[i, :, :] = U_P                                                 # Perpendicular component velocity [m/s].
        alpha_vec[i] = alpha                                                                    # Angle of attack [rad].

        # Force calculations.
        T_vec[i] = T_v                                                                                     # Thrust [N].
        Q_vec[i] = Q_v                                                                                     # Torque [N].
        Y_vec[i] = Y_v                                                                                # Cross force [N].

        # Power contributions calculation with IT.
        P_i[i] = w_0 / w_h * P_ih                                                                   # Induced power [W].
        P_0[i] = Pc0 * rho_inf * A * Omega * R * (Omega ** 2 * R ** 2 + 4.7 * Vinf ** 2)            # Profile power [W].
        P_fus[i] = 0.5 * rho_inf * Vinf ** 3 * f                                      # Parasitic power, eq. (7.25) [W].
        P[i] = P_i[i] + P_0[i] + P_fus[i] / N_rotors                            # Single rotor total power required [W].
        P_tot[i] = N_rotors * P[i]                                                           # Total power required [W].

    return (
        P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
        Omega_new, alpha_e_vec, dTc, U_P_vec, dH_dpsi_0, dQ_dpsi_0,
        T_vec, Y_vec, Q_vec, dT_dr_dpsi, dQ_dr_dpsi, psi
    )

def Calculate_alpha_hover(mu, Cla, sol, r_bar, theta):
    """
            Calculation of alpha in hover conditions, based on BEMT_Hover function.

            Args:
                mu: Advance ratio.
                Cla: Lift curve slope [rad^-1].
                sol: Rotor solidity.
                r_bar: Dimensionless rotor radius.
                theta: Pitch angle distribution [rad].

            Returns:
                alpha_hover: angle of attack in hovering [rad].
    """

    C1 = mu + Cla * sol / 8                                                                  # First degree coefficient.
    C2 = -r_bar * Cla * sol / 8 * (theta - mu / r_bar)                                                      # Free term.

    lambda_i = (-C1 + np.sqrt(C1 ** 2 - 4 * C2)) / 2                                             # Induced inflow ratio.

    phi = np.arctan(mu / r_bar + lambda_i / r_bar)                                                 # Inflow angle [rad].
    alpha_hover = theta - phi                                                                   # Angle of attack [rad].

    return alpha_hover