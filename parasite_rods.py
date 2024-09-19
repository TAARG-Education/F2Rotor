def landing_rods(self):
    '''Methods that evaluates the parasite drag area of the landing rods'''

    # Extracting the parameters
    R_b = self.config["ParasiteArea"][4]["R"]
    L_b = self.config["ParasiteArea"][4]["L"]
    alpha = self.config["ParasiteArea"][4]["alpha"]
    S = self.S

    # Calculating the wet Area
    # the wet area corresponds to the lateral area of a cylinder
    S_wet = ma.pi*R_b*L_b
    # Calculating the drag coefficient
    # this is an empirical formula. 1.1 is the basic drag coefficient of the cylinder
    # at alpha=90 degrees
    CD = 1.1*(ma.sin(alpha))**3+0.02
    # Calculating the adimensional parasite drag area of the hub
    f_A_lr_b = CD*S_wet/S
    # Calculating the parasite drag area of the hub
    # alternatively we can see this formula as:
    # f_cs= CD*S_wet
    f_lr_b = f_A_lr_b*S
    # Rounding of f_Acs and f_cs
    f_A_lr_b = round(f_A_lr_b,5)
    f_lr_b = round(f_lr_b,4)

    R_f = self.config["ParasiteArea"][5]["R"]
    L_f = self.config["ParasiteArea"][5]["L"]
    alpha = self.config["ParasiteArea"][5]["alpha"]
    S = self.config["ParasiteArea"][5]["S"]

    S_wet = ma.pi*R_f*L_f
    CD = 1.1*(ma.sin(alpha))**3+0.02

    f_A_lr_f = CD*S_wet/S
    f_lr_f = f_A_lr_f*S

    f_A_lr_f = 2*round(f_A_lr_f,5)
    f_lr_f = 2*round(f_lr_f,4)

    f_A_lr = f_A_lr_b + f_A_lr_f
    f_lr = f_lr_b + 2*f_lr_f
    # Results
    return f_A_lr, f_lr