
    
#This test case analyzes the flapping behavior of helicopter rotor blades by computing and plotting
#the flapping angle (denoted as β) for different collective pitch angles (θ_0). The steps involved are:

#1. Parameters Setup:
#   - Key flight and rotor parameters such as the Lock number (gamma), blade twist angle (θ_tw), advance ratio (mu),
#     inflow ratio (lambda), and rotor angular velocity (Omega) are defined.
#   - The azimuth angle (psi) is computed over a full rotor cycle.


#2. Flapping Coefficients Calculation:   
#    -For three different collective pitch angles (θ_0 = 5°, 10°, 15°), the flapping coefficients β_0, β_1c, β_1s,
#     and the flapping angle beta are computed as functions of time using the `compute_beta_coefficients` function.
#
#3. Data Visualization:
#    - The computed flapping angles are converted from radians to degrees and plotted as a function of the azimuth
#      angle for each θ_0 value.
    
#4. Results:
#
#     - The plot illustrates how the collective pitch angle affects the flapping behavior of the rotor blades, comparing
#       beta vs. azimuth angle psi for different values of θ_0.



# Authors: El Mehdi Rahmaoui
# Rotary Wing Aerodynamics Course
# University of Naples Federico II
# Academic Year 2023-2024

import numpy as np
import matplotlib.pyplot as plt

from copter_flapping_coeff import compute_beta_coefficients

# Test parameters
gamma = 8        # Lock number
theta_0_values = [np.radians(5), np.radians(10), np.radians(15)]  # Peak collective pitch angle in radians
theta_tw = np.radians(2)  # Blade twist angle in radians
mu = 0.1            # Advance ratio
lamb = 0.0008       # Inflow ratio
omega = 30.0        # Angular velocity in rad/s
t = np.linspace(0, 2 * np.pi / omega, 200)  # Time interval for a full cycle
psi = np.degrees(omega * t)  # Azimuth angle converted to degrees

# Plot beta for different values of theta_0
plt.figure(figsize=(12, 6))

for theta_0 in theta_0_values:
    # Compute flapping coefficients beta_0, beta_1c, beta_1s, and beta over time
    beta_0, beta_1c, beta_1s, beta = compute_beta_coefficients(gamma, theta_0, theta_tw, mu, lamb, omega, t)

    # Convert beta to degrees for plotting
    beta_deg = np.degrees(beta)

    # Print the computed coefficients in degrees
    print(f"Theta_0: {np.degrees(theta_0):.1f}°")
    print(f"  beta_0: {np.degrees(beta_0):.1f}°")
    print(f"  beta_1c: {np.degrees(beta_1c):.1f}°")
    print(f"  beta_1s: {np.degrees(beta_1s):.1f}°")
    
    # Plot beta as a function of azimuth angle psi (in degrees)
    plt.plot(psi, beta_deg, label=f'$\\theta_0$ = {np.degrees(theta_0):.1f}°')

# Labeling the axes and setting plot details
plt.xlabel('$\psi$ [deg]')
plt.ylabel(' $\\beta$ [deg]')
plt.title('Flapping Coefficient $\\beta$ vs Azimuth Angle for different $\\theta_0$ values')
plt.legend()
plt.grid(True)
plt.show()