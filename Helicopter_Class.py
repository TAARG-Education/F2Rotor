
class Helicopter:
    def __init__(self, Omega_r_mr, R_mr, R_hub, c_mr, N_mr, theta0, theta_tw, Cla, Cd0, Omega_r_tr, R_tr, R_hub_tr, c_tr, N_tr, l_tr):

        # Main rotor
        self.Omega_r_mr = Omega_r_mr # Main rotor tip speed (m/s)

        self.R_mr = R_mr  # Main rotor radius (m)
        self.R_hub = R_hub  # Rotor hub radius (m)
        self.c_mr = c_mr  # Rotor's chord (m)
        self.N_mr = N_mr  # Number of blades of the helicopter
        self.theta0 = theta0  # Pitch angle at the root (°)
        self.theta_tw = theta_tw  # Blade twist (tip minus root incidence) (°)
        self.Cla = Cla  # Lift curve slope of the airfoil (1/rad)
        self.Cd0 = Cd0  # Profile drag coefficient of the airfoil

         # Tail rotor
        self.Omega_r_tr = Omega_r_tr         # Tail rotor tip speed (m/s)
        self.R_tr = R_tr                     # Tail rotor radius (m)
        self.R_hub_tr = R_hub_tr             # Tail Rotor hub radius (m)
        self.c_tr = c_tr                     # Tail rotor chord (m)
        self.N_tr = N_tr                     # Number of blades of tail rotor
        self.l_tr = l_tr                     # Distance between tail rotor and main rotor

        # Developing the possibility of integration within the feature_geometry
