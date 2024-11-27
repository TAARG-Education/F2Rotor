from math import atan, sqrt, sin, cos, pi, exp, acos
import numpy as np

def convang(angle, from_unit, to_unit):
    if from_unit == 'deg' and to_unit == 'rad':
        return np.deg2rad(angle)
    elif from_unit == 'rad' and to_unit == 'deg':
        return np.rad2deg(angle)
    else:
        raise ValueError("Non supported units")

def RotorSolving(r_segn, psi, Vinf, Omega_new, R, A, N, w, alpha, theta, Clalpha, Cd_mean, c, rho_inf, f, W,
                      linear_inflow):
    """
            Calculation of Tc, Qc and power contributions, determining the dimensionless error
            in order to reach horizontal flight condition (T = W)

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
                H:
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

    def atmosisa(height):
        """
       International Standard Atmosphere model (ISA)
        Calculation of density, pressure and temperature at certain altitudes (excluding stratosphere).

        Parameters:
            height (float): Altitude [m].

        Returns:
            rho (float): Air density [kg/m^3]
            P (float): Atmospheric pressure [Pa]
            T (float): Temperature [K]
            a_inf (float): sound speed [m/s]
        """

    F_prandtl = np.zeros(len(r_segn))
    U_R = np.zeros((len(r_segn), len(psi)))
    U_T = np.zeros((len(r_segn), len(psi)))
    U_P = np.zeros((len(r_segn), len(psi)))
    V_e = np.zeros((len(r_segn), len(psi)))
    phi = np.zeros((len(r_segn), len(psi)))
    alpha_e = np.zeros((len(r_segn), len(psi)))
    Cl = np.zeros((len(r_segn), len(psi)))
    dL = np.zeros((len(r_segn), len(psi)))

    dT_dr_dpsi = np.zeros((len(r_segn), len(psi)))
    dH_dr_dpsi = np.zeros(len(r_segn))
    dH_dr_dpsi_0 = np.zeros(len(r_segn))
    dY_dr_dpsi = np.zeros(len(r_segn))
    dQ_dr_dpsi = np.zeros(len(r_segn))
    dQ_dr_dpsi_i = np.zeros(len(r_segn))
    dQ_dr_dpsi_0 = np.zeros(len(r_segn))

    dT_dpsi = np.zeros(len(psi))
    dH_dpsi = np.zeros(len(psi))
    dH_dpsi_0 = np.zeros(len(psi))
    dY_dpsi = np.zeros(len(psi))
    dQ_dpsi = np.zeros(len(psi))
    dQ_dpsi_i = np.zeros(len(psi))
    dQ_dpsi_0 = np.zeros(len(psi))

    for jj in range(len(psi)):
        for ii in range(len(r_segn)):
            mu = Vinf / (Omega_new * R)
            lambda_i = w / (Omega_new * R)

            if linear_inflow == 1:
                Tc = W / (rho_inf * Omega_new ** 2 * R ** 2 * A)
                lambda_i = linear_inflow(alpha, mu, 0, Tc, r_segn[ii] * R, psi[jj], lambda_i)

            F_prandtl[ii] = 2 / pi * acos(exp(N / (2 * lambda_i) * (r_segn[ii] - 1)))

            U_R[ii, jj] = Vinf * cos(alpha) * cos(psi[jj])
            U_T[ii, jj] = Omega_new * r_segn[ii] * R + Vinf * cos(alpha) * sin(psi[jj])

            if linear_inflow == 1:
                U_P[ii, jj] = lambda_i * Omega_new * R
            else:
                U_P[ii, jj] = Vinf * sin(alpha) + w

            V_e[ii, jj] = np.sqrt(U_T[ii, jj] ** 2 + U_P[ii, jj] ** 2)
            phi[ii, jj] = atan(U_P[ii, jj] / U_T[ii, jj])
            alpha_e[ii, jj] = theta[ii] - phi[ii, jj]

            if alpha_e[ii, jj] < convang(-15, 'deg', 'rad'):
                alpha_e[ii, jj] = convang(-15, 'deg', 'rad')
            if alpha_e[ii, jj] > convang(15, 'deg', 'rad'):
                alpha_e[ii, jj] = convang(15, 'deg', 'rad')

            Cl[ii, jj] = Clalpha * alpha_e[ii, jj]
            Cd = Cd_mean

            chord = c[ii] if np.ndim(c) > 0 else c
            dL[ii, jj] = 0.5 * rho_inf * V_e[ii, jj] ** 2 * chord * Cl[ii, jj]
            dD = 0.5 * rho_inf * V_e[ii, jj] ** 2 * chord * Cd

            Fx = dL[ii, jj] * sin(phi[ii, jj]) + dD * cos(phi[ii, jj])
            Fz = dL[ii, jj] * cos(phi[ii, jj]) - dD * sin(phi[ii, jj])
            Fr = 0  # Flapping neglected

            dT_dr_dpsi[ii, jj] = Fz * F_prandtl[ii]
            dH_dr_dpsi[ii] = Fx * sin(psi[jj]) + Fr * cos(psi[jj])
            dH_dr_dpsi_0[ii] = dD * cos(phi[ii, jj]) * sin(psi[jj])
            dY_dr_dpsi[ii] = -Fx * cos(psi[jj]) + Fr * sin(psi[jj])
            dQ_dr_dpsi[ii] = r_segn[ii] * R * Fx
            dQ_dr_dpsi_i[ii] = r_segn[ii] * R * dL[ii, jj] * sin(phi[ii, jj])
            dQ_dr_dpsi_0[ii] = r_segn[ii] * R * dD * cos(phi[ii, jj])

        # Ensure dT_dpsi is correctly calculated for integration
        dT_dpsi[jj] = np.trapz(dT_dr_dpsi[:, jj], r_segn * R)

        dH_dpsi[jj] = np.trapz(dH_dr_dpsi, r_segn * R)
        dH_dpsi_0[jj] = np.trapz(dH_dr_dpsi_0, r_segn * R)
        dY_dpsi[jj] = np.trapz(dY_dr_dpsi, r_segn * R)
        dQ_dpsi[jj] = np.trapz(dQ_dr_dpsi, r_segn * R)
        dQ_dpsi_i[jj] = np.trapz(dQ_dr_dpsi_i, r_segn * R)
        dQ_dpsi_0[jj] = np.trapz(dQ_dr_dpsi_0, r_segn * R)

    # Final computation of T, H, Y, Q, etc.
    T = N / (2 * pi) * np.trapz(dT_dpsi, psi)
    H = N / (2 * pi) * np.trapz(dH_dpsi, psi)
    H_0 = N / (2 * pi) * np.trapz(dH_dpsi_0, psi)
    Y = N / (2 * pi) * np.trapz(dY_dpsi, psi)
    Q = N / (2 * pi) * np.trapz(dQ_dpsi, psi)
    Q_0 = N / (2 * pi) * np.trapz(dQ_dpsi_0, psi)

    D_fus = f * 0.5 * rho_inf * Vinf ** 2

    Tc = T / (rho_inf * Omega_new ** 2 * R ** 2 * A)
    Hc_0 = H_0 / (rho_inf * Omega_new ** 2 * R ** 2 * A)
    Hc_0 = 1.7 / 2 * Hc_0
    Qc = Q / (rho_inf * Omega_new ** 2 * R ** 3 * A)
    Qc_0 = Q_0 / (rho_inf * Omega_new ** 2 * R ** 3 * A)

    Pc_i = 1.15 * lambda_i * Tc
    Pc_0 = Qc_0 + mu * Hc_0
    Pc_fus = mu * D_fus / W * Tc
    Pc = Pc_i + Pc_0 + Pc_fus / 4

    P = Pc * rho_inf * Omega_new ** 3 * R ** 3 * A
    P_i = Pc_i * rho_inf * Omega_new ** 3 * R ** 3 * A
    P_0 = Pc_0 * rho_inf * Omega_new ** 3 * R ** 3 * A
    P_fus = Pc_fus * rho_inf * Omega_new ** 3 * R ** 3 * A

    err = (T - W) / W
    return err, T, H, P, P_i, P_0, P_fus, alpha_e, V_e, dT_dr_dpsi, U_P, dH_dpsi_0, dQ_dpsi_0, Tc, Y, Qc




