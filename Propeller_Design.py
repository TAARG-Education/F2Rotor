import numpy as np
from scipy.integrate import simps
from scipy.interpolate import pchip_interpolate

class Propeller_Design():
    def __init__(self, C_t:float, J:float, N_b:int, r_hub_R:float, N_stations:int = 100, C_l = 1, r_bar_Cl = 1):
        self.C_t = C_t;
        self.J = J;
        self.N_b = N_b;
        self.r_hub_R = r_hub_R;
        self.C_l = np.array(C_l);
        self.N_stations = N_stations;
        self.r_bar_Cl = r_bar_Cl;
        if np.size(self.r_bar_Cl) == 1:
            self.r_bar_Cl = np.linspace(self.r_hub_R, 1, self.N_stations);
    
    def computation(self, a):
        phi = np.arctan(self.J/(np.pi*self.r_bar)*(1 + a));
        a_vec = a*(np.cos(phi))**2;
        a_prime = self.J/(np.pi*self.r_bar)*a*np.cos(phi)*np.sin(phi);

        F = 2/np.pi*np.arccos(np.exp(- self.N_b*(1 - self.r_bar)/(2*phi[self.N_stations - 1])));

        NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F;
        dCt_dr_bar = 0.25*NGF_OmR2*np.pi**2*(1 - a_prime)*self.r_bar;
        C_t= simps(dCt_dr_bar, self.r_bar);
    
        return C_t, a_vec, a_prime, F, NGF_OmR2, phi

    
    def iteration(self):
        a_0 = - 0.5 + np.sqrt(0.25 + 2*self.C_t/(self.J**2 * np.pi));
        self.r_bar = np.linspace(self.r_hub_R, 1, self.N_stations);

        if np.size(self.C_l) == 1:
            self.C_l = self.C_l*np.ones(self.N_stations);
        elif np.size(self.C_l) > 1:
            r_interp = np.linspace(self.r_hub_R, 1, np.size(self.C_l));
            print(np.size(r_interp), np.size(self.C_l));
            self.C_l = pchip_interpolate(r_interp, self.C_l, self.r_bar);

        C_t_0, _, _, _, _, _ = Propeller_Design.computation(self, a_0);

        err_0 = (C_t_0 - self.C_t)/self.C_t;

        a_1 = 1.1*a_0;
        C_t_1, _, _, _, _, _ = Propeller_Design.computation(self, a_1);

        err_1 = (C_t_1 - self.C_t)/self.C_t;
        err = err_1

        while np.absolute(err) > 1e-12:
            a_n = (a_0*err - a_1*err_0)/(err - err_0);

            C_t_n, a_vec, a_prime, F, NGF_OmR2, phi = Propeller_Design.computation(self, a_n);
            err_0 = err;
            err = (C_t_n - self.C_t)/self.C_t;

            a_0 = a_1;
            a_1 = a_n;

        sigma_Cl = NGF_OmR2*1/(np.pi*self.r_bar**2*np.sqrt((((self.J)/(np.pi*self.r_bar))*(1 + a_vec))**2 + \
                                                      (1 - a_prime)**2));

        sigma = sigma_Cl/self.C_l;
        theta = self.C_l/(2*np.pi) + phi;

        return sigma, theta, a_vec, a_prime, F, NGF_OmR2, self.C_l
    
    def sigma(self):
        sig, _, _, _, _, _, _ = Propeller_Design.iteration(self);
        return sig
    
    def pitch(self):
        _, theta, _, _, _, _, _ = Propeller_Design.iteration(self);
        return theta
    
    def Prandtl_correction(self):
        _, _, _, _, F, _, _ = Propeller_Design.iteration(self);
        return F
    
    def axial_induction(self):
        _, _, a_vec, _, _, _, _ = Propeller_Design.iteration(self);
        return a_vec
    
    def rotational_induction(self):
        _, _, _, a_prime, _, _, _ = Propeller_Design.iteration(self);
        return a_prime
    
    def lift_coefficient(self):
        _, _, _, _, _, _, self.C_l = Propeller_Design.iteration(self);
        return self.C_l


    









    