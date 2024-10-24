import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import convert_temperature
from scipy.constants import atmosphere


def atmosisa(h):
    # Calcola la velocità del suono e altre proprietà atmosferiche a una certa altezza
    # Assumiamo h in metri (altitudine)
    # Restituisce (T, a_inf, rho, p)
    if h < 11000:
        T = 15.04 - 0.00649 * h
        p = 101.29 * ((T + 273.1) / 288.08) ** 5.256
    else:
        T = -56.46
        p = 22.65 * np.exp(1.73 - 0.000157 * h)

    T += 273.15  # Converti da gradi Celsius a Kelvin
    rho = p / (0.2869 * T)
    a_inf = np.sqrt(1.4 * 287.05 * T)  # Velocità del suono in aria
    return T, a_inf, rho, p


def BladeSection_alpha_Mach(Omega, R, lambda_, r_segn, beta, dbeta, mu, psi, theta, alpha_stall_up, alpha_stall_lo, h):
    # Creazione della griglia mesh.
    r_2d, psi_2d = np.meshgrid(np.linspace(0, r_segn[0], len(r_segn)), psi)
    x_hub = r_2d * np.cos(psi_2d - np.pi / 2)
    y_hub = r_2d * np.sin(psi_2d - np.pi / 2)

    r_2d, psi_2d = np.meshgrid(r_segn, psi)
    x = r_2d * np.cos(psi_2d - np.pi / 2)
    y = r_2d * np.sin(psi_2d - np.pi / 2)

    # Inizializzazione delle variabili.
    u_P = np.zeros(
        (len(r_segn), len(psi)))  # velocità dell'aria della sezione della pala, perpendicolare al piano del disco.
    u_T = np.zeros((len(r_segn), len(psi)))  # velocità dell'aria della sezione della pala, tangente al piano del disco.
    u_R = np.zeros((len(r_segn), len(psi)))  # velocità radiale dell'aria della sezione della pala.
    phi = np.zeros((len(r_segn), len(psi)))  # angolo d'incidenza della sezione.
    alpha_e = np.zeros((len(r_segn), len(psi)))  # angolo d'attacco della sezione della pala.
    V_eff = np.zeros((len(r_segn), len(psi)))

    _, a_inf, _, _ = atmosisa(h)

    for i in range(len(psi)):
        for j in range(len(r_segn)):
            u_P[j, i] = lambda_ + r_segn[j] * dbeta[i] + beta[i] * mu * np.cos(psi[i])
            u_T[j, i] = r_segn[j] + mu * np.sin(psi[i])
            u_R[j, i] = mu * np.cos(psi[i])
            phi[j, i] = np.arctan(u_P[j, i] / u_T[j, i])
            alpha_e[j, i] = theta[j] - phi[j, i]
            V_eff[j, i] = Omega * R * np.sqrt(u_P[j, i] ** 2 + u_T[j, i] ** 2)

    stall = np.zeros((len(r_segn), len(psi)))  # regione stallata della pala del rotore.
    non_stall = np.full((len(r_segn), len(psi)), np.nan)  # regione non stallata della pala del rotore.
    M = np.full((len(r_segn), len(psi)), np.nan)
    M_cr = np.full((len(r_segn), len(psi)), np.nan)

    for iii in range(len(r_segn)):
        for jjj in range(len(psi)):
            if alpha_e[iii, jjj] - np.radians(alpha_stall_up) > 0 or alpha_e[iii, jjj] - np.radians(alpha_stall_lo) < 0:
                stall[iii, jjj] = alpha_e[iii, jjj]
            else:
                non_stall[iii, jjj] = alpha_e[iii, jjj]
            M[iii, jjj] = V_eff[iii, jjj] / a_inf

    # Plot
    plt.figure(1)
    plt.fill(x_hub, y_hub, 'w')
    plt.axis('equal')
    plt.axis('off')
    plt.grid(False)
    plt.text(-r_segn[-10], r_segn[-1], r'$\mu = $' + str(mu), color='r', fontsize=12)
    contourf = plt.contourf(x, y, np.degrees(non_stall.T), 20, cmap='gray')
    plt.colorbar(contourf)
    plt.show()

    return alpha_e, phi, u_P, u_T, u_R, V_eff, non_stall, M, M_cr

