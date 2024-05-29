############ Propeller Design F2Rotor ############
# Authors: Antonio Mazzara, Antonio D'Onofrio
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini
# University of Naples Federico II
# Academic Year 2023-2024
# Date: 29/05/2024
#
# Description: This class computes the solidity, pitch and load distributions 
#              for an optimal propeller. The procedure is non-dimensional.
#              Three different tip-loss corrections can be applied
#              (Prandtl, Goldstein or no correction).
#
# References: "Lezioni di Aerodinamica dell'Ala Rotante", vsn 2.01, Renato Tognaccini,
#              section 3.8, "Progetto dell'elica"
#             "Wind Turbine Aerodynamics and vorticity-based methods", Volume 7, Emmanuel Branlard,
#              Chapter 14, Springer, 2017


import numpy as np
from scipy.integrate import simps
from scipy.interpolate import pchip_interpolate

class Propeller_Design():

    ######## Initialization method ########
        # The input parameters are: thrust coefficient "C_t" (Renard definition); 
        #                           advance ratio "J" (= V/(n*D));
        #                           number of blades "N_b";
        #                           hub radius-tip radius ratio "r_hub_R", number of radial stations "N_stations" (default values is 100);
        #                           lift coefficient distribution along radius "C_l" (default value is 1 for each station);
        #                           radial stations in which the lif coefficient distribution is imposed "r_bar_Cl" (default value is 1);
        #                           tip-loss correction: "none" ---> no correction is applied;
        #                                                "p"    ---> Prandtl correction is applied;
        #                                                "g"    ---> Goldstein correction is applied.


    def __init__(self, C_t:float, J:float, N_b:int, r_hub_R:float, N_stations:int = 100, C_l = 1, r_bar_Cl = 1, correction = 'none'):

        self.C_t = C_t;

        self.J = J;

        self.N_b = N_b;

        self.r_hub_R = r_hub_R;

        self.C_l = np.array(C_l);

        self.N_stations = N_stations;

        self.r_bar_Cl = r_bar_Cl;

        if np.size(self.r_bar_Cl) == 1:
            self.r_bar_Cl = np.linspace(self.r_hub_R, 1, self.N_stations);
        
        self.correction = correction;
    

    ################## Goldstein correction ##################
    # The procedure to compute the Goldstein correction consists of 3 methods:
    # - fGoldsteinFactor;
    # - fGoldsteinCirculation;
    # - fUi_HelixNTheory.
    # The Goldstein factor "G" is computed using a superposition of helic filaments.
    
    def fGoldsteinFactor(self):

        # Computes the Goldstein factor "G".
        # l_bar : dimepnsionless torsional parameter h/(2*pi*R)
        
        w = 1;

        l_bar = self.J/np.pi;

        G = Propeller_Design.fGoldsteinCirculation(self, l_bar, w);
        G = G*self.N_b/(2*np.pi*l_bar*w);
    
        return G
    
    def fGoldsteinCirculation(self, l_bar, w):

        # Computes Goldstein circulation using superposition of helix.

        # Definition of the radial stations.
        vr = np.linspace(0,1,self.N_stations);
        vr = vr[1:self.N_stations];
        N = self.N_stations - 1;

        # Definition of the control points (CPs) in between vortices.
        vrCP = np.arange(3/(2*N),1,1/N);

        # Calculation of matrices of influence at CPs for a helical filament.
        A = np.zeros((N,N));

        # Circulation is unitary.
        Gamma_h = 1; 

        #  Azimutal position is 0 (on the Lifting line).
        psi_h = 0;

        # Loop on helical vortex radii.
        for j in range(1,N): 

            # Loop on Control points radii.
            for i in range(1,N-1):
                
                A[i,j] = Propeller_Design.fUi_HelixNTheory(Gamma_h, vrCP[i], vr[j], l_bar, psi_h, self.N_b);
    
        #  Boundary conditions values on the vortex sheet CP.
        U= np.zeros(N);

        U[1:(N -1)] = w*1./(1+ (l_bar**2)/(vrCP[1:( N -1)])**2);

        # Condition : sum of gamma = 0
        A[-1,:] = 1;
        U[-1] = 0;

        # Solving for individual trailed vorticity and cumulative one.
        A_pinv = np.linalg.pinv(A);
        Gamma_t = np.dot(A_pinv, U);
        Gamma_Goldstein = np.cumsum(Gamma_t);
        Gamma_Goldstein = np.append(0, Gamma_Goldstein);
    
        return Gamma_Goldstein

    def fUi_HelixNTheory(Gamma, r , r0 , l ,psih , B):

        # Computes induction from N - helices.

        C0z = (l**2+ r0**2)**(1/4)/(l**2+ r**2)**(1/4) ;

        C1z = l/24*((3*r**2 - 2*l**2)/(l**2 + r**2)**(3/2) + (2*l**2 + 9*r0**2)/( l**2 + r0**2)**(3/2));

        pexi = r/r0*(l + np.sqrt(l**2+ r0**2))/(l + np.sqrt (l**2+ r**2) )*np.exp(np.sqrt(l**2 + r**2)/l)/np.exp(np.sqrt(l**2 + r0**2)/l);

        mexi = 1/pexi; 

        t = psih;

        if abs(r)<r0:

            tmp = 1/(( mexi*np.exp(-1j* t))**B - 1);
            vz = 1/(2*np.pi*l) + 1/(2*np.pi* l)*C0z*np.real (tmp + C1z/B*np.log(1 + tmp ));  
         
        elif abs(r)>r0:

            tmp = 1/(( pexi*np.exp(-1j* t))**B - 1 );
            vz = 0 + 1/( 2*np.pi*l )*C0z*np.real( - tmp + C1z /B*np.log(1 + tmp ) );
        else:
            vz = 0;
    
        uz =  - B*Gamma*vz;
    
        return uz
    

    ######## Computation method #########
    # This method is used in the iteration one. The input parameters are self and the 
    # axial induction factor in the current iteration. The output parameters are:
    # C_t;
    # the axial induction factor along the blade "a_vec";
    # the rotational induction factor along the blade "a_prime";
    # the correction factor "F";
    # the non-dimensional load distribution "NGF_OmR2";
    # the inflow angle distribution along the blade "phi";
    # the thrust coefficient distribution "dCt_dr_bar".
    
    def computation(self, a):

        # Computation of the inflow angle.
        phi = np.arctan(self.J/(np.pi*self.r_bar)*(1 + a));

        # Computation of the axial induction factor along the blade.
        a_vec = a*(np.cos(phi))**2;

        # Computation of the axial induction factor along the blade.
        a_prime = self.J/(np.pi*self.r_bar)*a*np.cos(phi)*np.sin(phi);

        # Computation of the load distribution depending on the correction in input.
        if self.correction.lower() == 'p':
            F = 2/np.pi*np.arccos(np.exp(- self.N_b*(1 - self.r_bar)/(2*phi[self.N_stations - 1])));
            NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F;
        elif self.correction.lower() == 'g':
            G = Propeller_Design.fGoldsteinFactor(self);
            vr = np.linspace(0,1,self.N_stations);
            G = pchip_interpolate(vr, G, self.r_bar);
            NGF_OmR2 = (2*G*a_vec*self.J**2)*(1 + a_vec)/(np.pi*(1 - a_prime));
            F = G;
        else:
            F = 1;
            NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F;

        # Computation of the thrust coefficient distribution along the blade.
        dCt_dr_bar = 0.25*NGF_OmR2*np.pi**2*(1 - a_prime)*self.r_bar;

        # Computation of the thrust coefficient through Cavalieri-Simpson integration.
        C_t = simps(dCt_dr_bar, self.r_bar);
    
        return C_t, a_vec, a_prime, F, NGF_OmR2, phi, dCt_dr_bar
    

    ######## Iteration method #########
    #  Iteration using the false position method.
    #  Based on the error from the thrust coefficient given in input.
    
    def iteration(self):
        
        # Computation of the first attempt axial induction factor from the Momentum Theory.
        a_0 = - 0.5 + np.sqrt(0.25 + 2*self.C_t/(self.J**2 * np.pi));
        
        # Definition of the radial stations vector.
        self.r_bar = np.linspace(self.r_hub_R, 1, self.N_stations);

        # Defintion of the lift coefficient distribution vector.
        if np.size(self.C_l) == 1:
            self.C_l = self.C_l*np.ones(self.N_stations);
        elif np.size(self.C_l) > 1:
            r_interp = np.linspace(self.r_hub_R, 1, np.size(self.C_l));
            print(np.size(r_interp), np.size(self.C_l));
            self.C_l = pchip_interpolate(r_interp, self.C_l, self.r_bar);

        # Computation of the first attempt thrust coefficient.
        C_t_0, _, _, _, _, _, _ = Propeller_Design.computation(self, a_0);

        # Evaluation of the error based on thrust coefficient given in input.
        err_0 = (C_t_0 - self.C_t)/self.C_t;

        # Second attempt: the first axial induction factor is increased of 10%.
        a_1 = 1.1*a_0;

        # Computation of the second attempt thrust coefficient.
        C_t_1, _, _, _, _, _, _ = Propeller_Design.computation(self, a_1);

        # Evaluation of the second error.
        err_1 = (C_t_1 - self.C_t)/self.C_t;
        err = err_1;

        # Loop in which the false position method is implemented.
        while np.absolute(err) > 1e-12:
            a_n = (a_0*err - a_1*err_0)/(err - err_0);

            C_t_n, a_vec, a_prime, F, NGF_OmR2, phi, dCt_dr_bar = Propeller_Design.computation(self, a_n);
            err_0 = err;
            err = (C_t_n - self.C_t)/self.C_t;

            a_0 = a_1;
            a_1 = a_n;

        # Computation of the solidity.
        sigma_Cl = NGF_OmR2*1/(np.pi*self.r_bar**2*np.sqrt((((self.J)/(np.pi*self.r_bar))*(1 + a_vec))**2 + \
                                                      (1 - a_prime)**2));
        sigma = sigma_Cl/self.C_l;

        # Computation of the pitch.
        theta = self.C_l/(2*np.pi) + phi;

        return sigma, theta, a_vec, a_prime, F, NGF_OmR2, self.C_l, dCt_dr_bar
    
    ##################################
    ######### Output methods #########
    ##################################
       
    # This method returns the solidity.
    def sigma(self):
        sig, _, _, _, _, _, _, _ = Propeller_Design.iteration(self);
        return sig
    
    # This method returns the pitch.
    def pitch(self):
        _, theta, _, _, _, _, _, _= Propeller_Design.iteration(self);
        return theta
    
    # This method returns the correction factor.
    def correction_factor(self):
        _, _, _, _, F, _, _, _ = Propeller_Design.iteration(self);
        return F
    
    # This method returns the axial induction factor.
    def axial_induction(self):
        _, _, a_vec, _, _, _, _, _ = Propeller_Design.iteration(self);
        return a_vec
    
    # This method returns the rotational induction factor.
    def rotational_induction(self):
        _, _, _, a_prime, _, _, _, _ = Propeller_Design.iteration(self);
        return a_prime
    
    # This method returns the lift coefficient distribution.
    def lift_coefficient(self):
        _, _, _, _, _, _, self.C_l, _ = Propeller_Design.iteration(self);
        return self.C_l
    
    # This method returns the load distribution.
    def load_distribution(self):
        _, _, _, _, _, NGF_OmR2, _, _ = Propeller_Design.iteration(self);
        return NGF_OmR2
    
    # This method returns the thrust coefficient distribution.
    def thrust_coeff_distribution(self):
        _, _, _, _, _, _, _, dCt_dr_bar = Propeller_Design.iteration(self);
        return dCt_dr_bar


    









    