class Helicopter:
    def __init__(self, Omega_R_mr, R_mr, R_hub, c_mr, N_mr, theta0, theta_tw, Cla, Cd0):
        self.Omega_r_mr = Omega_R_mr # Main rotor tip speed (m/s)
        self.R_mr = R_mr  # Main rotor radius (m)
        self.R_hub = R_hub  # Rotor hub radius (m)
        self.c_mr = c_mr  # Rotor's chord (m)
        self.N_mr = N_mr  # Number of blades of the helicopter
        self.theta0 = theta0  # Pitch angle at the root (°)
        self.theta_tw = theta_tw  # Blade twist (tip minus root incidence) (°)
        self.Cla = Cla  # Lift curve slope of the airfoil (1/rad)
        self.Cd0 = Cd0  # Profile drag coefficient of the airfoil

        # Developing the possibility of integration within the feature_geometry
