#PARASITE RODS
#Description: This function calculates the parasite drag area
#             of the landing gear rods, assuming it is a circular section cylinder
#
#   WARNING: This code operates on a helicopter configuration that uses TWO rods
#            for the FRONT landing gear and ONE for the REAR landing gear.
#            If the helicopter you are working with does not follow this configuration,
#            it is necessary to modify the code internally by multiplying the values
#            of f_A and f_A_lr by the correct number of rods.
#
#References:-"Lezioni Aerodinamica dell'ala rotante", v2.05, Renato Tognaccini,
#             section 7.5, "Potenza parassita della fusoliera"
#           -"Fluid-Dynamic Drag", Hoerner, chapter 3, section 5, "Drag of round bodies"
#
#Author: Catello Del Sorbo, Marcello Zammarrelli
#Rotary Wing Aerodynamics course, prof. Renato Tognaccini
#University of Naples Federico II
#Academic Year 2023-2024
#Date:17/09/2024
# 
#Version:1.0.0


import math as ma

def landing_rods(self):
    '''Methods that evaluates the parasite drag area of the landing rods.

    Input:
        -R:     rod radius
        -L:     rod length
        -alpha: angle of attack [rad]
        -S:     reference area

    Output:
        f_A_lr: Combined parasite drag area of the landing rods.
        f_lr: Total parasite drag of the landing rods.
    '''

    def calculate_drag_for_section(section_index):
        ''' Helper function to calculate drag for a given section of the landing rods '''
        R = self.config["ParasiteArea"][section_index]["R"]
        L = self.config["ParasiteArea"][section_index]["L"]
        alpha = self.config["ParasiteArea"][section_index]["alpha"]
        S = self.config["ParasiteArea"][section_index]["S"]

        # Calculating the wet area
        # The wet area corresponds to the lateral area of a cylinder
        S_wet = ma.pi * R * L

        # Calculating the drag coefficient
        # This is an empirical formula from "Fluid-Dynamic Drag", Hoerner,
        # chapter 3, section 5, "Drag of round bodies"
        CD = 1.1 * (ma.sin(alpha)) ** 3 + 0.02

        # Calculating the adimensional parasite drag area of the hub
        f_A_lr = CD * S_wet / S
        f_lr = f_A_lr * S

        # Rounding of f_Acs and f_cs
        return round(f_A_lr, 5), round(f_lr, 4)


    # Calculate for the back section (index 4)
    f_A_lr_b, f_lr_b = calculate_drag_for_section(4)

    # Calculate for the front section (index 5)
    f_A_lr_f, f_lr_f = calculate_drag_for_section(5)

    # Adjust front section for symmetry (x2)
    f_A_lr_f *= 2
    f_lr_f *= 2

    # Combine both sections
    f_A_lr = f_A_lr_b + f_A_lr_f
    f_lr = f_lr_b + f_lr_f

    return f_A_lr, f_lr
