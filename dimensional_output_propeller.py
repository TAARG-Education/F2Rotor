import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

class dimensional_output_propeller:
    # This class is created to produce the dimensional output of a propeller given
    # non-dimensional parameters as Ct, Cp and J
    # For the dimensional output, the user has to define:
    # 1) altitudes
    # 2) propeller diameter
    # 3) RPM

    def __init__(self,Ct,Cp,J,D,altitude,RPM) :
        self.Ct = Ct #Thrust coefficient: Ct = T/(rho*n^2*D^4)
        self.Cp = Cp #Power coefficient: Cp = P/(rho*n^3*D^5)
        self.J = J # advance ratio J = V_inf/(n*D)
        self.D = D # propeller diameter [m]
        self.altitude = altitude # [m]
        self.RPM = RPM    # revolutions per minute 
        self.n =self.RPM/60 #revolutions per second
        self.V_infty = self.J*self.n*self.D  # air speed [m/s]
        self.density = Atmosphere(self.altitude).density #air density [kg/m^3]

        


    def Thrust(self):
        # This function returns thrust [N] for given altitudes and velocities
        T=np.zeros((len(self.altitude),len(self.V_infty)))
        for i in range(len(self.altitude)):
            T[i,:] = self.Ct*self.density[i]*self.n**2*self.D**4
        return T

    
    def Power(self):
         # This function returns power [W] for given altitudes and velocities
         P=np.zeros((len(self.altitude),len(self.V_infty)))
         for i in range(len(self.altitude)):
            P[i,:] = self.Cp*self.density[i]*self.n**3*self.D**5
         return P
    
    def Efficiency(self):
        # This function returns efficiency for given velocities
        eta = self.J*self.Ct/self.Cp
        return eta


    def Thrust_plot(self,a):
        T=self.Thrust()
        #This function generates plots of thrust as function of velocity and altitude
        for i in range(len(self.altitude)):
            plt.plot(self.V_infty,T[i,:], linewidth=2, linestyle='solid',marker = 'o',label="quota: " + str(self.altitude[i])+ "m")
            plt.xlabel('$V_\infty [m/s]$')
            plt.ylabel ('Thrust [N]')
            plt.grid(True)  
            
        
        plt.xlim([0,max(self.V_infty)*1.1])
        plt.ylim([0,np.max(np.max(T))*1.1])
        plt.legend()
        # Set a=1 to save the figure as pdf
        if a==1:
            plt.savefig('Thrust.pdf')

        plt.show()
        


    def Power_plot(self,a):
        #This function generates plots of power as function of velocity and altitude
        P = self.Power()
        for i in range(len(self.altitude)):
            plt.plot(self.V_infty,P[i,:], linewidth=2, linestyle='solid',marker = 'o',label="quota: " + str(self.altitude[i])+"m")
        plt.xlabel('$V_\infty [m/s]$')
        plt.ylabel ('Power [W]')
        
        plt.grid(True)
        plt.xlim([0,max(self.V_infty)*1.1])
        plt.ylim([0,np.max(np.max(P))*1.1])
        plt.legend()
        # Set a=1 to save the figure as pdf
        if a==1:
            plt.savefig('Power.pdf')
        plt.show()

    

    def Efficiency_plot(self,a):
        #This function generates plot of efficieny as function of velocity 
        eta = self.Efficiency()
        plt.plot(self.V_infty,eta,color='black', linewidth=2, linestyle='solid',marker = 'o')
        plt.xlabel('$V_\infty [m/s]$')
        plt.ylabel ('Efficiency')
        plt.xlim([0,max(self.V_infty)*1.1])
        plt.ylim([0,1.2])
        plt.grid(True)
        # Set a=1 to save the figure as pdf
        if a==1:
            plt.savefig('Efficiency.pdf')
        plt.show()
        

    def interp_Thrust(self,V_input):
        # This function returns thrust for a specific velocity given by the user
        T=self.Thrust()
        T_output = np.zeros(len(self.altitude))
        for i in range(len(self.altitude)):
            T_output[i] = np.interp(V_input, self.V_infty, T[i,:], left=None, right=None, period=None)
        return T_output
    
    def interp_Power(self,V_input):
        # This function returns power for a specific velocity given by the user
        P=self.Power()
        P_output = np.zeros(len(self.altitude))
        for i in range(len(self.altitude)):
            P_output[i] = np.interp(V_input, self.V_infty, P[i,:], left=None, right=None, period=None)
        return P_output
    
    def interp_efficiency(self,V_input):
        # This function returns efficiency for a specific velocity given by the user
        eta=self.Efficiency()
        eta_output = np.interp(V_input, self.V_infty, eta, left=None, right=None, period=None)
        return eta_output

   

    
        
