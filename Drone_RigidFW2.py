import numpy as np
from scipy.optimize import brentq
from sympy import symbols, Eq, solve
from scipy.constants import g, R
from Drone_RotorSolving import RotorSolving


def RigidFW2(
        Vinf_vec, N, r_segn, c, theta, W, w_h, height, R, A, P_ih, f, Pc0,
        Clalpha, Cd_mean, Omega, N_rotors, linear_inflow
):
    """
        Calculation of Tc, Qc and power contributions through a bisection method.

        Args:
            Vinf_vec: airspeed [m/s]
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
            f: equivalent wet area (% chord)
            Pc0:
            Clalpha: Slope of the load-angle curve [rad^-1]
            Cd_mean: Parasitic resistance coefficient
            Omega: angular velocity [rad/s]
            N_rotors: Number of rotors
            linear_inflow: Non identified function

        Returns:
            P_i: Induced power by impulsive theory [W]
            P_0: Parasitic power by impulsive theory [W]
            P_fus: Power absorbed by "fuselage" by impulsive theory[W]
            P: Total power by impulsive theory [W]
            P_new: Total power by BEMT [W]
            P_new_i: Induced power by BEMT [W]
            P_new_0: Parasitic power by "fuselage" [W]
            P_new_fus: Power absorbed by "fuselage" by BEMT [W]
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
        # Standard constants
        T0 = 288.15     # Temperature at sea level [K]
        P0 = 101325     # Pressure at sea level [Pa]
        L = 0.0065      # Temperature gradient (K/m)
        R_air = 287.05  # Specific constant for dry air (J/(kg·K))
        g0 = 9.80665    # Gravity acceleration (m/s^2)
        gamma = 1.4     # Adiabatic expansion constant

        if height <= 11000:  # Troposphere (altitude <= 11 km)
            Temp = T0 - L * height
            Pressure = P0 * (Temp / T0) ** (g0 / (R_air * L))
        else:
            raise ValueError("Altitude beyond the limits of the model implemented.")

        # Air density and sound speed (perfect gas)
        rho = Pressure / (R_air * Temp)
        a_inf = (gamma * R_air * Temp) ** 0.5
        return rho, Pressure, Temp, a_inf

    rho_inf, pressure_inf, temperature_inf, a_inf = atmosisa(height)

    psi = np.linspace(0, 2 * np.pi, 60)

    # Initialization
    Omega_new = np.zeros(len(Vinf_vec))
    alpha_vec = np.zeros(len(Vinf_vec))
    P_new_i = np.zeros(len(Vinf_vec))
    P_new_0 = np.zeros(len(Vinf_vec))
    P_new_fus = np.zeros(len(Vinf_vec))
    P_new = np.zeros(len(Vinf_vec))
    alpha_e_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))
    Mach_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))
    dTc = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))
    U_P_vec = np.zeros((len(Vinf_vec), len(r_segn), len(psi)))

    Tc_vec = np.zeros(len(Vinf_vec))
    Qc_vec = np.zeros(len(Vinf_vec))
    Y_vec = np.zeros(len(Vinf_vec))

    P_i = np.zeros(len(Vinf_vec))
    P_0 = np.zeros(len(Vinf_vec))
    P_fus = np.zeros(len(Vinf_vec))
    P = np.zeros(len(Vinf_vec))
    P_tot = np.zeros(len(Vinf_vec))


    for i, Vinf in enumerate(Vinf_vec):
        alpha = 0
        err = 1

        while err > 1e-6:

            w = symbols('w', real=True, positive=True)
            eqn = (
                    (Vinf / w_h * w / w_h * np.sin(alpha) + (w / w_h) ** 2) ** 2
                    + (Vinf / w_h) ** 2 * (w / w_h) ** 2 * (np.cos(alpha)) ** 2
                    - 1
            )
            sol = solve(Eq(eqn, 0), w)

            for s in sol:
                if s.is_real and s > 0:
                    w = float(s)
                    break

            if alpha == 0:
                w_0 = w

            # Starting bisection method
            Omega_a, Omega_b = 600, 1500

            err_a, *_ = RotorSolving(r_segn, psi, Vinf, Omega_a, R, A, N, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow)

            err_b, *_ = RotorSolving(r_segn, psi, Vinf, Omega_b, R, A, N, w, alpha, theta,
                                          Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow)

            # Debug output
            print(f"Omega_a: {Omega_a}, err_a: {err_a}")
            print(f"Omega_b: {Omega_b}, err_b: {err_b}")

            if err_a * err_b > 0:
                raise ValueError(f"Redefine the search range: "
                                 f"err_a: {err_a}, err_b: {err_b}, range: [{Omega_a}, {Omega_b}]")

            it = 0
            err_m = 1
            tol = 5e-3

            while abs(err_m) > tol:
                it += 1

                Omega_New = (Omega_a + Omega_b) / 2

                result = RotorSolving(
                    r_segn, psi, Vinf, Omega_New, R, A, N, w, alpha, theta,
                    Clalpha, Cd_mean, c, rho_inf, f, W, linear_inflow
                )

                (
                    err_m, _, H, P_New, P_New_i, P_New_0, P_New_fus, alpha_e, V_e, dT_dr_dpsi,
                    U_P, dH_dpsi_0, dQ_dpsi_0, Tc_vec[i], Y_vec[i], Qc_vec[i]
                ) = result

                if err_m * err_a < 0:
                    Omega_b = Omega_New
                    err_m = err_b
                else:
                    Omega_a = Omega_New
                    err_a = err_m

            D_fus = 0.5 * rho_inf * Vinf ** 2 * f
            alpha_new = np.arctan((D_fus / N_rotors + H) / W)

            err = abs((alpha_new - alpha) / alpha_new)
            alpha = alpha_new

        # Saving results
        Omega_new[i] = Omega_New

        P_new[i] = P_New
        P_new_i[i] = P_New_i
        P_new_0[i] = P_New_0
        P_new_fus[i] = P_New_fus
        alpha_e_vec[i, :, :] = alpha_e
        Mach_vec[i, :, :] = V_e / a_inf
        dTc[i, :, :] = dT_dr_dpsi / (rho_inf * Omega_new[i] ** 2 * R ** 2 * A)
        U_P_vec[i, :, :] = U_P
        alpha_vec[i] = alpha

        P_i[i] = w_0 / w_h * P_ih
        P_0[i] = Pc0 * rho_inf * A * Omega * R * (Omega ** 2 * R ** 2 + 4.7 * Vinf ** 2)
        P_fus[i] = 0.5 * rho_inf * Vinf ** 3 * f

        P[i] = P_i[i] + P_0[i] + P_fus[i] / N_rotors
        P_tot[i] = N_rotors * P[i]

    return (
        P_i, P_0, P_fus, P, P_new, P_new_i, P_new_0, P_new_fus, alpha_vec,
        Omega_new, alpha_e_vec, Mach_vec, dTc, U_P_vec, dH_dpsi_0, dQ_dpsi_0,
        Tc_vec, Y_vec, Qc_vec
    )
