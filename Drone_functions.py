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
from math import atan, sin, cos, pi, exp, acos
from ambiance import Atmosphere
from sympy import symbols, Eq, solve

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

    epsilon = 1e-6                                                                                           # Avoid /0.

    ## DETERMINATION OF THE DISTRIBUTIONS OF ANGULAR AND AERODYNAMIC QUANTITY.

    # The second degree equation for the induced inflow ratio has coefficient 1, C1 and C2.

    C1 = mu + Clalpha * solidity / 8                                                         # First degree coefficient.
    C2 = -r_segn * Clalpha * solidity / 8 * (theta - mu / (r_segn + epsilon))                               # Free term.

    # Ensure discriminant > 0.
    discriminant = C1 ** 2 - 4 * C2
    discriminant[discriminant < 0] = 0                                                           # Avoid complex values.
    lambda_i = (-C1 + np.sqrt(discriminant)) / 2                                                 # Induced inflow ratio.

    # Calculate phi and alpha.
    phi = np.arctan(np.clip(mu / (r_segn + epsilon) + lambda_i / (r_segn + epsilon), -1e3, 1e3)) # Inflow angle
    # [rad].
    alpha = np.clip(theta - phi, -np.radians(90), np.radians(90))                               # Angle of attack [rad].

    # Prandtl correction function.
    F_Prandtl = np.where(lambda_i > 0, 2 / np.pi * np.arccos(np.exp(N / (2 * lambda_i) * (r_segn - 1))), 0)

    # Calculate Cl and Cd.
    Cl = Clalpha * alpha                                                                             # Lift coefficient.
    Cd = Cd0                                                                                         # Drag coefficient.

    # Thrust and power distributions.
    dTcdr_segn = 0.5 * solidity * Cl * r_segn ** 2         # Incremental thrust coefficient along the main rotor radius.
    dTcdr_segn_PR = dTcdr_segn * F_Prandtl # Incremental thrust coefficient along the main rotor radius considering the
    # Prandtl correction function.

    dPc_i = 0.5 * solidity * Cl * phi * r_segn ** 3 # Incremental induced power coefficient of the rotor along the main
    # rotor radius.
    dPc_0 = 0.5 * solidity * Cd * r_segn ** 3      # Incremental parasite power coefficient along the main rotor radius.
    dPcdr_segn = dPc_i + dPc_0                              # Incremental power coefficient along the main rotor radius.

    # Integrate with trapezoidal method
    Tc = float(np.trapz(dTcdr_segn_PR, r_segn))       # Thrust coefficient  considering the Prandtl correction function.
    Pc_i = float(np.trapz(dPc_i, r_segn))                                                   # Induced power coefficient.
    Pc_0 = float(np.trapz(dPc_0, r_segn))                                                  # Parasite power coefficient.
    Qc = float(np.trapz(dPcdr_segn, r_segn))                                                       # Torque coefficient.

    # Calculate figure of merit (FM).
    FM = (Tc ** (3 / 2) / np.sqrt(2)) / Qc if abs(Qc) > epsilon else np.nan

    return Tc, Qc, Pc_i, Pc_0, FM



# Convang function.
def convang(angle, from_unit, to_unit):
    if from_unit == 'deg' and to_unit == 'rad':                                       # If [angle] = [deg], use deg2rad.
        return np.deg2rad(angle)
    elif from_unit == 'rad' and to_unit == 'deg':                                     # If [angle] = [rad], use rad2deg.
        return np.rad2deg(angle)
    else:
        raise ValueError("Non supported units")

def RotorSolving(r_segn, psi, Vinf, Omega_new, R, A, N, w, alpha, theta, Clalpha, Cd_mean, c, rho_inf, f, W,
                      linear_inflow):
    """
            Calculation of Tc, Qc and power contributions, determining the dimensionless error
            in order to reach horizontal flight condition (T = W).

            Args:
                r_segn: Normalized vector of radial positions (r/R)
                psi: Azimuth angle [rad]
                Vinf: airspeed [m/s]
                Omega_new: Angular velocity calculated with bisection method [rad/s]
                R: Rotor radius [m]
                A: Rotor area [m^2]
                N: number of blades
                w: induced velocity [m/s]
                alpha: Angle of attack [rad]
                theta: pitch angle [rad]
                Clalpha: Slope of the load-angle curve [rad^-1]
                Cd_mean: Parasitic resistance coefficient
                c: Mean geometric chord [m]
                rho_inf: air density [Kg/m^3]
                f: equivalent wet area (% chord)
                W: Drone weight [N]
                linear_inflow: Non identified function

            Returns:
                err: Dimensionless error
                T: Thrust [N]
                H: Rotor drag [N]
                P: Total power by impulsive theory [W]
                P_i: Induced power by impulsive theory [W]
                P_0: Parasitic power by impulsive theory [W]
                P_fus: Power absorbed by "fuselage" by impulsive theory [W]
                alpha_e: Effective angle of attack [rad]
                V_e: Effective velocity [m/s]
                dT_dr_dpsi: Derivative of thrust
                U_P: U_P velocity [m/s]
                dT_dr_dpsi: Derivative of thrust versus r_segn and psi
                dH_dpsi_0: Derivative of H versus psi, valued at 0
                dQ_dpsi_0: Derivative of Q versus psi, valued at 0
                Tc: Thrust coefficient
                Y: Cross force [N]
                Qc: Torque coefficient
    """

    # Initialization
    F_prandtl = np.zeros(len(r_segn))                                                                # Prandtl Function.
    U_R = np.zeros((len(r_segn), len(psi)))                                                 # Radial velocity component.
    U_T = np.zeros((len(r_segn), len(psi)))                                             # Tangential velocity component.
    U_P = np.zeros((len(r_segn), len(psi)))                                          # Perpendicular velocity component.
    V_e = np.zeros((len(r_segn), len(psi)))                                                        # Effective velocity.
    phi = np.zeros((len(r_segn), len(psi)))                                                             # Azimuth angle.
    alpha_e = np.zeros((len(r_segn), len(psi)))                                              # Effective angle of attack.
    Cl = np.zeros((len(r_segn), len(psi)))                                                           # Lift coefficient.
    dL = np.zeros((len(r_segn), len(psi)))                                                          # Differential lift.

    dT_dr_dpsi = np.zeros((len(r_segn), len(psi)))                              # Derivative of T versus r_segn and psi.
    dH_dr_dpsi = np.zeros(len(r_segn))                                          # Derivative of H versus r_segn and psi.
    dH_dr_dpsi_0 = np.zeros(len(r_segn))                              # Derivative of H parasitic versus r_segn and psi.
    dY_dr_dpsi = np.zeros(len(r_segn))                                          # Derivative of H versus r_segn and psi.
    dQ_dr_dpsi = np.zeros(len(r_segn))                                          # Derivative of Q versus r_segn and psi.
    dQ_dr_dpsi_i = np.zeros(len(r_segn))                                # Derivative of Q induced versus r_segn and psi.
    dQ_dr_dpsi_0 = np.zeros(len(r_segn))                              # Derivative of Q parasitic versus r_segn and psi.

    dT_dpsi = np.zeros(len(psi))                                                           # Derivative of T versus psi.
    dH_dpsi = np.zeros(len(psi))                                                           # Derivative of H versus psi.
    dH_dpsi_0 = np.zeros(len(psi))                                               # Derivative of H parasitic versus psi.
    dY_dpsi = np.zeros(len(psi))                                                           # Derivative of Y versus psi.
    dQ_dpsi = np.zeros(len(psi))                                                           # Derivative of Q versus psi.
    dQ_dpsi_i = np.zeros(len(psi))                                                 # Derivative of Q induced versus psi.
    dQ_dpsi_0 = np.zeros(len(psi))                                               # Derivative of Q parasitic versus psi.

    for jj in range(len(psi)):
        for ii in range(len(r_segn)):
            mu = Vinf / (Omega_new * R)                                                     # Advance ratio, eq. (7.14).
            lambda_i = w / (Omega_new * R)                                                  # Axial interference factor.

            # Ignore the following instructions.
            ##########
            if linear_inflow == 1:
                Tc = W / (rho_inf * Omega_new ** 2 * R ** 2 * A)
                lambda_i = linear_inflow(alpha, mu, 0, Tc, r_segn[ii] * R, psi[jj], lambda_i)
            ##########

            F_prandtl[ii] = 2 / pi * acos(exp(N / (2 * lambda_i) * (r_segn[ii] - 1)))    # Prandtl Function, eq. (6.79).

            U_R[ii, jj] = Vinf * cos(alpha) * cos(psi[jj])                # Radial velocity component, eq. (7.17) [m/s].
            U_T[ii, jj] = (Omega_new * r_segn[ii] *
                           R + Vinf * cos(alpha) * sin(psi[jj]))      # Tangential velocity component, eq. (7.17) [m/s].

            if linear_inflow == 1:                                                            # Ignore this instruction.
                U_P[ii, jj] = lambda_i * Omega_new * R
            else:
                U_P[ii, jj] = Vinf * sin(alpha) + w    # Perpendicular velocity component, eq. (8.5) with no flap [m/s].

            V_e[ii, jj] = np.sqrt(U_T[ii, jj] ** 2 + U_P[ii, jj] ** 2)    # Effective velocity in the blade plane [m/s].
            phi[ii, jj] = atan(U_P[ii, jj] / U_T[ii, jj])                      # Azimuth angle in the blade plane [rad].
            alpha_e[ii, jj] = theta[ii] - phi[ii, jj]                                           # Angle of attack [rad].

            # If alpha_e is greater or minor than 15 deg, set it respectively to 15 deg and -15 deg,
            # and turn it in rad.
            if alpha_e[ii, jj] < convang(-15, 'deg', 'rad'):
                alpha_e[ii, jj] = convang(-15, 'deg', 'rad')
            if alpha_e[ii, jj] > convang(15, 'deg', 'rad'):
                alpha_e[ii, jj] = convang(15, 'deg', 'rad')

            Cl[ii, jj] = Clalpha * alpha_e[ii, jj]               # Lift coefficient calculation according linear theory.
            Cd = Cd_mean                                                                             # Drag coefficient.

            # If chord is variable with r_segn use c[ii], else use c.
            chord = c[ii] if np.ndim(c) > 0 else c
            dL[ii, jj] = 0.5 * rho_inf * V_e[ii, jj] ** 2 * chord * Cl[ii, jj]                      # Differential lift.
            dD = 0.5 * rho_inf * V_e[ii, jj] ** 2 * chord * Cd                                      # Differential drag.

            Fx = dL[ii, jj] * sin(phi[ii, jj]) + dD * cos(phi[ii, jj])           # x component of the thrust, eq. (8.8).
            Fz = dL[ii, jj] * cos(phi[ii, jj]) - dD * sin(phi[ii, jj])           # z component of the thrust, eq. (8.8).
            Fr = 0                                                                                 # Flapping neglected.

            dT_dr_dpsi[ii, jj] = Fz * F_prandtl[ii]                 # Derivative of T versus r_segn and psi, eq. (8.11).
            dH_dr_dpsi[ii] = Fx * sin(psi[jj]) + Fr * cos(psi[jj])  # Derivative of H versus r_segn and psi, eq. (8.11).
            dH_dr_dpsi_0[ii] = dD * cos(phi[ii, jj]) * sin(psi[jj]) # Derivative of H parasitic versus r_segn and psi,
            # eq. (8.11).
            dY_dr_dpsi[ii] = -Fx * cos(psi[jj]) + Fr * sin(psi[jj]) # Derivative of Y versus r_segn and psi, eq. (8.11).
            dQ_dr_dpsi[ii] = r_segn[ii] * R * Fx                    # Derivative of Q versus r_segn and psi, eq. (8.11).
            dQ_dr_dpsi_i[ii] = r_segn[ii] * R * dL[ii, jj] * sin(phi[ii, jj]) # Derivative of Q induced versus r_segn
            # and psi, eq. (8.11).
            dQ_dr_dpsi_0[ii] = r_segn[ii] * R * dD * cos(phi[ii, jj]) # Derivative of Q parasitic versus r_segn and psi,
            # eq. (8.11).

        # Integration
        dT_dpsi[jj] = np.trapz(dT_dr_dpsi[:, jj], r_segn * R)                     # Thrust integration along the radius.
        dH_dpsi[jj] = np.trapz(dH_dr_dpsi, r_segn * R)                        # Rotor drag integration along the radius.
        dH_dpsi_0[jj] = np.trapz(dH_dr_dpsi_0, r_segn * R)          # Rotor parasitic drag integration along the radius.
        dY_dpsi[jj] = np.trapz(dY_dr_dpsi, r_segn * R)                       # Cross force integration along the radius.
        dQ_dpsi[jj] = np.trapz(dQ_dr_dpsi, r_segn * R)                            # Torque integration along the radius.
        dQ_dpsi_i[jj] = np.trapz(dQ_dr_dpsi_i, r_segn * R)                # Induced torque integration along the radius.
        dQ_dpsi_0[jj] = np.trapz(dQ_dr_dpsi_0, r_segn * R)              # Parasitic torque integration along the radius.

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
    Hc_0 = 1.7 / 2 * Hc_0
    Qc = Q / (rho_inf * Omega_new ** 2 * R ** 3 * A)                                    # Torque coefficient, eq. (2.5).
    Qc_0 = Q_0 / (rho_inf * Omega_new ** 2 * R ** 3 * A)                                 # Parasitic torque coefficient.

    # Power coefficient contributions according to (8.35) equation.
    Pc_i = 1.15 * lambda_i * Tc                                                             # Induced power coefficient.
    Pc_0 = Qc_0 + mu * Hc_0                                                                 # Profile power coefficient.
    Pc_fus = mu * D_fus / W * Tc                                                          # Parasitic power coefficient.
    Pc = Pc_i + Pc_0 + Pc_fus / 4                                             # Single rotor required power coefficient.

    # Power contributions.
    P_i = Pc_i * rho_inf * Omega_new ** 3 * R ** 3 * A                                              # Induced power [W].
    P_0 = Pc_0 * rho_inf * Omega_new ** 3 * R ** 3 * A                                              # Profile power [W].
    P_fus = Pc_fus * rho_inf * Omega_new ** 3 * R ** 3 * A                                        # Parasitic power [W].
    P = Pc * rho_inf * Omega_new ** 3 * R ** 3 * A                                    # Single rotor power required [W].

    err = (T - W) / W
    return err, T, H, P, P_i, P_0, P_fus, alpha_e, V_e, dT_dr_dpsi, U_P, dH_dpsi_0, dQ_dpsi_0, Tc, Y, Qc



def RigidFW2(
        Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A, P_ih, f, Pc0,
        Clalpha, Cd_mean, Omega, N_rotors, linear_inflow
):
    """
        Calculation of Tc, Qc and power contributions through a bisection method.
        Note that the comments refers to the equation contained in the following book:
        Renato Tognaccini. Lezioni per il corso di Aerodinamica dellAla Rotante. 2023.

        Args:
            Vinf_vec: Airspeed [m/s]
            N: Number of blades
            r_segn: Normalized vector of radial positions (r/R)
            c: Mean geometric chord [m]
            theta: Pitch angle distribution [rad]
            W: Drone weight [N]
            w_h: Induced velocity in hover conditions [m/s]
            height: Altitude [m]
            R: Rotor radius [m]
            A: Rotor area [m^2]
            P_ih: Induced power in hovering [W]
            f: Equivalent wet area (% chord)
            Pc0: Power parasitic coefficient
            Clalpha: Slope of the load-angle curve [rad^-1]
            Cd_mean: Parasitic resistance coefficient
            Omega: Angular velocity [rad/s]
            N_rotors: Number of rotors
            linear_inflow: Non identified function

        Returns:
            P_i: Induced power by impulsive theory [W]
            P_0: Parasitic power by impulsive theory [W]
            P_fus: Power absorbed by "fuselage" by impulsive theory[W]
            P: Total power by impulsive theory [W]
            P_new: Total power by BEMT [W]
            P_new_i: Induced power by BEMT [W]
            P_new_0: Profile power by BEMT [W]
            P_new_fus: Parasitic power absorbed by "fuselage" by BEMT [W]
            alpha_vec: Angle of attack [rad]
            Omega_new: Angular velocity calculated with bisection method [rad/s]
            alpha_e_vec: Effective angle of attack [rad]
            Mach_vec: Mach number
            dTc: Differential thrust coefficient
            U_P_vec: U_P [m/s]
            dH_dpsi_0: Derivative of H versus psi, valued at 0
            dQ_dpsi_0: Derivative of Q versus psi, valued at 0
            Tc_vec: Thrust coefficient
            Y_vec: Cross force [N]
            Qc_vec: Torque coefficient

        """

    atmosphere = Atmosphere(height)
    rho_inf = atmosphere.density[0]                                                              # Air density [kg/m^3].
    a_inf = atmosphere.speed_of_sound[0]                                                         # Speed of sound [m/s].

    psi = np.linspace(0, 2 * np.pi, 60)                                                       # Azimuth angle.

    # Initialization.
    Omega_new = np.zeros(len(Vinf_vec))                                                              # Angular velocity.
    alpha_vec = np.zeros(len(Vinf_vec))                                                               # Angle of attack.
    P_new_i = np.zeros(len(Vinf_vec))                                              # Induced power calculated with BEMT.
    P_new_0 = np.zeros(len(Vinf_vec))                                              # Profile power calculated with BEMT.
    P_new_fus = np.zeros(len(Vinf_vec))                                          # Parasitic power calculated with BEMT.
    P_new = np.zeros(len(Vinf_vec))                                                  # Total power calculated with BEMT.
    alpha_e_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                          # Effective angle of attack.
    Mach_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                                           # Mach number.
    dTc = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                            # Differential thrust coefficient.
    U_P_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))                       # Perpendicular velocity component.

    Tc_vec = np.zeros(len(Vinf_vec))                                                               # Thrust coefficient.
    Qc_vec = np.zeros(len(Vinf_vec))                                                               # Torque coefficient.
    Y_vec = np.zeros(len(Vinf_vec))                                                                       # Cross force.

    P_i = np.zeros(len(Vinf_vec))                                                    # Induced power calculated with IT.
    P_0 = np.zeros(len(Vinf_vec))                                                    # Profile power calculated with IT.
    P_fus = np.zeros(len(Vinf_vec))                                                # Parasitic power calculated with IT.
    P = np.zeros(len(Vinf_vec))                                  # Single rotor total power required calculated with IT.
    P_tot = np.zeros(len(Vinf_vec))                                            # Total power required calculated with IT.


    for i, Vinf in enumerate(Vinf_vec):
        alpha = 0                                                                                 # Initial alpha value.
        err = 1                                                                                     # Initial err value.

        while err > 1e-6:

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
            Omega_a, Omega_b = 600, 1500                                                         # Initial Omega values.

            # Errors calculation with Omega_a and Omega_b.
            err_a, *_ = RotorSolving(r_segn, psi, Vinf, Omega_a, R, A, N, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow)

            err_b, *_ = RotorSolving(r_segn, psi, Vinf, Omega_b, R, A, N, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow)

            # Check the errors signs.
            if err_a * err_b > 0:
                raise ValueError(f"Redefine the search range: "                                          # Error signal.
                                 f"err_a: {err_a}, err_b: {err_b}, range: [{Omega_a}, {Omega_b}]")

            # Updating variables.
            it = 0                                                                                   # Starting counter.
            err_m = 1                                                                                      # Mean error.
            tol = 5e-3                                                                                      # Tolerance.

            while abs(err_m) > tol:
                it += 1                                                                           # Improve the counter.

                Omega_New = (Omega_a + Omega_b) / 2                                                  # Update Omega_New.

                result = RotorSolving(                                                   # RotorSolving function output.
                    r_segn, psi, Vinf, Omega_New, R, A, N, w, alpha, theta,
                    Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow
                )

                (
                    err_m, _, H, P_New, P_New_i, P_New_0, P_New_fus, alpha_e, V_e, dT_dr_dpsi,
                    U_P, dH_dpsi_0, dQ_dpsi_0, Tc_vec[i], Y_vec[i], Qc_vec[i]
                ) = result                                                        # RotorSolving function output saving.

                # Bisection method instructions.
                if err_m * err_a < 0:
                    Omega_b = Omega_New                                                              # Updating Omega_b.
                    err_m = err_b                                                                      # Updating err_m.
                else:
                    Omega_a = Omega_New                                                              # Updating Omega_a.
                    err_a = err_m                                                                      # Updating err_a.

            D_fus = 0.5 * rho_inf * Vinf ** 2 * f                                                   # Fuselage drag [N].
            alpha_new = np.arctan((D_fus / N_rotors + H) / W)           # Angle of attack calculation, eq. (8.33) [rad].

            err = abs((alpha_new - alpha) / alpha_new)                                                 # Updating error.
            alpha = alpha_new                                                                    # Updating alpha [rad].

        # Power contributions calculation with BEMT.
        Omega_new[i] = Omega_New
        P_new[i] = P_New                                                                              # Total power [W].
        P_new_i[i] = P_New_i                                                                        # Induced power [W].
        P_new_0[i] = P_New_0                                                                        # Profile power [W].
        P_new_fus[i] = P_New_fus                                                                  # Parasitic power [W].

        alpha_e_vec[i, :, :] = alpha_e                                                # Effective angle of attack [rad].
        Mach_vec[i, :, :] = V_e / a_inf                                                                   # Mach number.
        dTc[i, :, :] = dT_dr_dpsi / (rho_inf * Omega_new[i] ** 2 * R ** 2 * A)        # Differential thrust coefficient.
        U_P_vec[i, :, :] = U_P                                                 # Perpendicular component velocity [m/s].
        alpha_vec[i] = alpha                                                                    # Angle of attack [rad].

        # Power contributions calculation with IT.
        P_i[i] = w_0 / w_h * P_ih                                                                   # Induced power [W].
        P_0[i] = Pc0 * rho_inf * A * Omega * R * (Omega ** 2 * R ** 2 + 4.7 * Vinf ** 2)            # Profile power [W].
        P_fus[i] = 0.5 * rho_inf * Vinf ** 3 * f                                      # Parasitic power, eq. (7.25) [W].
        P[i] = P_i[i] + P_0[i] + P_fus[i] / N_rotors                            # Single rotor total power required [W].
        P_tot[i] = N_rotors * P[i]                                                           # Total power required [W].

    return (
        P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
        Omega_new, alpha_e_vec, Mach_vec, dTc, U_P_vec, dH_dpsi_0, dQ_dpsi_0,
        Tc_vec, Y_vec, Qc_vec
    )