# PROPELLER DESIGN F2ROTOR
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
#             "Aerodynamics of V/STOL Flight", McCormick, B. W., 1967
#
# Authors: Antonio Mazzara, Antonio D'Onofrio
# Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini
# University of Naples Federico II
# Academic Year 2023-2024
#
# Date: 02/06/2024
#
# Version: 1.0.0

import numpy as np
from scipy.integrate import simps
from scipy.interpolate import pchip_interpolate

class Propeller_Design():

    def __init__(self, C_t:float, J:float, N_b:int, r_hub_R:float, N_stations:int = 100, C_l = 1, r_bar_Cl = 1, correction = 'none'):
        
        '''
        The input parameters are: 
        - thrust coefficient "C_t" (Renard definition); 
        - advance ratio "J" (= V/(n*D));
        - number of blades "N_b";
        - hub radius-tip radius ratio "r_hub_R", number of radial stations "N_stations" (default values is 100);
        - lift coefficient distribution along radius "C_l" (default value is 1 for each station);
        - radial stations in which the lif coefficient distribution is imposed "r_bar_Cl" (default value is 1);
        - tip-loss correction: "none" ---> no correction is applied;
                               "p"    ---> Prandtl correction is applied;
                               "g"    ---> Goldstein correction is applied.

        The output parameters, obtained are
        '''

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
    
        self.sigma, self.theta, self.a_vec, self.a_prime, self.F, self.NGF_OmR2, self.C_l, self.dCt_dr_bar = Propeller_Design.iteration(self);
    
    def fGoldsteinFactor(self):

        '''
        This function computes the tip-loss factor according to the Goldstein theory, using a superposition of
        helical vortex filaments of identical pitch (see "Wind Turbine Aerodynamics and vorticity-based methods", 
        Volume 7, Emmanuel Branlard, Chapter 14, Springer, 2017).
        The procedure to compute the Goldstein correction is divided in 3 further functions in a hierarchical way:

        fGoldsteinFactor 
                 |
                 V
        fGoldsteinCirculation
                 |
                 V
        fUi_HelixNTheory

        The latter computes the induction from N-helices, the second one evaluates the Goldstein circulation and
        the main one determines the Goldstein factor "G" defined as:
        
        G = Gamma/(h * w)

        where Gamma is the circulation, h is the pitch of the helix and w is the axial induced velocity.

        The input parameter is a struct variable with the following fields:
        - the advance ratio "J" (= V/(n * D));
        - the number of blades "N_b".

        The output variable is "G", the Goldstein tip-loss correction factor.

        Authors: Antonio Mazzara, Antonio D'Onofrio

        Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini
        University of Naples Federico II
        Academic Year 2023-2024

        Date: 02/06/2024

        Version: 1.0.0

        '''

        def fGoldsteinCirculation(self, l_bar, w):

            # Computes Goldstein circulation using superposition of helix.

            def fUi_HelixNTheory(Gamma, r , r0 , l , psih, B):

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

                    A[i,j] = fUi_HelixNTheory(Gamma_h, vrCP[i], vr[j], l_bar, psi_h, self.N_b);
        
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
        
        w = 1;

        # l_bar : dimepnsionless torsional parameter h/(2 * pi * R)
        l_bar = self.J/np.pi;

        G = fGoldsteinCirculation(self, l_bar, w);

        G = G*self.N_b/(2*np.pi*l_bar*w);
    
        return G
    
    def iteration(self):

        '''
        This function evaluates iteratively the solidity, pitch and load distributions for an optimal propeller.
        The iterative method used here is the false-position one, where the error is defined as:

        err = C_t_it - C_t  

        where "C_t_it" is the thrust coefficient computed at the current iteration and "C_t" is the thrust coefficient given in input.

        The input parameter is a struct variable, whose fields are:
        - thrust coefficient "C_t" (Renard definition); 
        - advance ratio "J" (= V/(n*D));
        - number of blades "N_b";
        - hub radius-tip radius ratio "r_hub_R", number of radial stations "N_stations" (default values is 100);
        - lift coefficient distribution along radius "C_l" (default value is 1 for each station);
        - radial stations in which the lif coefficient distribution is imposed "r_bar_Cl" (default value is 1);
        - tip-loss correction: "none" ---> no correction is applied;
                               "p"    ---> Prandtl correction is applied;
                               "g"    ---> Goldstein correction is applied. 

        The output parameters are:
        - the solidity distribution "sigma";
        - the pitch angle distribution "theta";
        - the axial induction factor distribution "a_vec";
        - the rotational induction factor distribution "a_prime";
        - the tip-loss correction factor "F";
        - the non-dimensional load distribution "NGF_OmR2";
        - the interpolated lift coefficient distribution "C_l";
        - the thrust coefficient per unit length "dCt_dr_bar".

        Note: the load is non-dimensionalized with "Omega * R^2", where "Omega" is the angular velocity of the propeller and "R" is the 
        tip radius.

        Authors: Antonio Mazzara, Antonio D'Onofrio

        Rotary Wing Aerodynamics Course, Prof. Renato Tognaccini
        University of Naples Federico II
        Academic Year 2023-2024

        Date: 02/06/2024

        Version: 1.0.0

        '''

        def computation(self, a):

            # Computation of the inflow angle.
            phi = np.arctan(self.J/(np.pi*self.r_bar)*(1 + a));

            # Computation of the axial induction factor along the blade.
            a_vec = a*(np.cos(phi))**2;

            # Computation of the axial induction factor along the blade.
            a_prime = self.J/(np.pi*self.r_bar)*a*np.cos(phi)*np.sin(phi);

            # Computation of the load distribution depending on the correction in input.
            if self.correction.lower() == 'p':

                # Prandtl tip-loss correction factor.
                # F = 2/np.pi*np.arccos(np.exp(- self.N_b*(1 - self.r_bar)/(2*phi[self.N_stations - 1])));
                F = 2/np.pi*np.arccos(np.exp(- self.N_b*(1 - self.r_bar)/(2*np.arctan(self.J/np.pi))));

                # Tip-loss correction of the load.
                NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F;
            
            elif self.correction.lower() == 'g':

                # Goldstein factor.
                G = Propeller_Design.fGoldsteinFactor(self);

                vr = np.linspace(0,1,self.N_stations);

                # Interpolation of the Goldstein factor on the real radial stations vector.
                G = pchip_interpolate(vr, G, self.r_bar);

                # Local speed ratio (Omega * r / V). 
                lambda_r = np.pi*self.r_bar/self.J;

                # Goldstein tip-loss correction factor.
                F_Go = G*(1 + lambda_r**2)/(lambda_r**2);
                # F_Go = G*((1 + 2*a_vec)/(1 - 2*a_prime)*(1 + lambda_r**2)/(lambda_r**2);

                # Tip-loss correction of the load.
                NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F_Go;

                F = F_Go;
            
            else:

                F = 1;

                # Load without tip-loss correction.
                NGF_OmR2 = 4*np.pi*(self.r_bar**2)*a_prime*F;

            # Computation of the thrust coefficient distribution along the blade.
            dCt_dr_bar = 0.25*NGF_OmR2*np.pi**2*(1 - a_prime)*self.r_bar;

            # Computation of the thrust coefficient through Cavalieri-Simpson integration.
            C_t = simps(dCt_dr_bar, self.r_bar);
        
            return C_t, a_vec, a_prime, F, NGF_OmR2, phi, dCt_dr_bar
        
        # Computation of the first attempt axial induction factor from the Momentum Theory.
        a_0 = - 0.5 + np.sqrt(0.25 + 2*self.C_t/(self.J**2 * np.pi));
        
        # Definition of the radial stations vector.
        self.r_bar = np.linspace(self.r_hub_R, 1, self.N_stations);

        # Defintion of the lift coefficient distribution vector.
        if np.size(self.C_l) == 1:

            self.C_l = self.C_l*np.ones(self.N_stations);
        
        elif np.size(self.C_l) > 1:

            r_interp = np.linspace(self.r_hub_R, 1, np.size(self.C_l));

            # Interpolation of the lift coefficient distribution on the real radial stations vector.
            self.C_l = pchip_interpolate(r_interp, self.C_l, self.r_bar);

        # Computation of the first attempt thrust coefficient.
        C_t_0, _, _, _, _, _, _ = computation(self, a_0);

        # Evaluation of the error based on thrust coefficient given in input.
        err_0 = (C_t_0 - self.C_t)/self.C_t;

        # Second attempt: the first axial induction factor is increased of 10%.
        a_1 = 1.1*a_0;

        # Computation of the second attempt thrust coefficient.
        C_t_1, _, _, _, _, _, _ = computation(self, a_1);

        # Evaluation of the second error.
        err_1 = (C_t_1 - self.C_t)/self.C_t;

        err = err_1;

        # Loop in which the false position method is implemented.
        while np.absolute(err) > 1e-12:

            # Computation of the new Lagrange moltiplicator value based on the false position method.
            a_n = (a_0*err - a_1*err_0)/(err - err_0);

            # Computation of the new variables.
            C_t_n, a_vec, a_prime, F, NGF_OmR2, phi, dCt_dr_bar = computation(self, a_n);

            # Updating the stored values for the next iteration.
            err_0 = err;

            a_0 = a_1;

            a_1 = a_n;

            # Computation of the total thrust coefficient error with respect to the input value.
            err = (C_t_n - self.C_t)/self.C_t;

            
        # Computation of the solidity.
        sigma_Cl = NGF_OmR2*1/(np.pi*self.r_bar**2*np.sqrt((((self.J)/(np.pi*self.r_bar))*(1 + a_vec))**2 + (1 - a_prime)**2));

        sigma = sigma_Cl/self.C_l;

        # Computation of the pitch.
        theta = self.C_l/(2*np.pi) + phi;

        return sigma, theta, a_vec, a_prime, F, NGF_OmR2, self.C_l, dCt_dr_bar
