#PARASITE HUB
#Description: This function calculates the parasite drag area
#             of the hub, assuming it is a circular section cylinder
#
#References:-"Lezioni Aerodinamica dell'ala rotante", v2.05, Renato Tognaccini,
#             section 7.5, "Potenza parassita della fusoliera"
#           -"Fluid-Dynamic Drag", Hoerner, chapter 3, section 5, "Drag of round bodies"
#
#Author: Arianna Palumbo
#Rotary Wing Aerodynamics course, prof. Renato Tognaccini
#University of Naples Federico II
#Academic Year 2023-2024
#Date:14/07/2024
#
#Version:1.0.0



import math as ma

def hub(config, S):
        '''
        This function calculates the parasite drag area of the hub, assuming it as a circular section
        cylinder.
        
        Input:
        -config: dictionary containing the parameters:
                 -R:hub radius
                 -Lcs: hub height
                 -alpha: angle of attack
        -S: reference area
        
        Output:
        -f_Acs:adimensional parasite drag area of the hub
        -f_cs:parasite drag area of the hub
        '''

        # Extracting the parameters
        R = config["ParasiteArea"][3]["R"]
        Lcs = config["ParasiteArea"][3]["Lcs"]
        alpha = config["ParasiteArea"][3]["alpha"]

        # Calculating the wet Area
        #the wet area corresponds to the lateral area of a cylinder
        S_wet = ma.pi * R * Lcs

        # Calculating the drag coefficient
        #this is an empirical formula. 1.1 is the basic drag coefficient of the cylinder
        # at alpha=90 degrees
        CD = 1.1 * (ma.sin(alpha)) ** 3 + 0.02

        # Calculating the adimensional parasite drag area of the hub
        f_Acs = CD * S_wet / S

        # Calculating the parasite drag area of the hub
        #alternatively we can see this formula as:
        # f_cs= CD*S_wet
        f_cs = f_Acs * S

        # Rounding of f_Acs and f_cs
        f_Acs = round(f_Acs, 5)
        f_cs = round(f_cs, 4)

        # Results
        return f_Acs, f_cs
