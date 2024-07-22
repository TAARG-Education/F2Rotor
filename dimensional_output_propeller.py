
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere



# ==================================================================================================
# |Name           : dimensional_output_propeller.py                                                |
# |Author         : Andrea Raia, Antonino Guida.                                                   |
# |                 University of Naples Federico II.                                              |
# |Version        : 1.0                                                                            |
# |Date           : 31/05/2024                                                                     |
# |Modified       : 31/05/2024                                                                     |
# |Description    : Functions that produce the dimensional output of a propeller.                  | 
# |                                                                                                |     
# |Reference      : R. Tognaccini, (a.a. 2023/2024), "Lezioni di Aerodinamica dell'ala rotante"    |                                                         
# |                                                                                                |
# |Non dimensional                                                                                 |
# |input variables:	1) Ct; thrust coefficient.		                                                |															                                                                       
# |					   2) Cp; power coefficient.	                                                   |								                                                                      
# |                  3) J; advance ratio.                                                          |
# |                                                                                                |
# |Dimensional                                                                                     |
# |input variables: 1) altitude; altitudes at which to calculate the dimensional values.           |
# |                 2) D; propeller diameter.                                                      |
# |                 3) RPM; revolutions per minute.                                                |
# |                                                                                                |
# |Output         : 1) T; thrust.                                                                  |
# |                 2) P; power.                                                                   |
# |                 3) eta; efficiency.                                                            |
# |Note           :                                                                                |
## =================================================================================================




class dimensional_output_propeller:

    def __init__(self,Ct,Cp,J,D,altitude,RPM) :
        self.Ct = Ct                                                     # Thrust coefficient: Ct = T/(rho*n^2*D^4)
        self.Cp = Cp                                                     # Power coefficient: Cp = P/(rho*n^3*D^5)
        self.J = J                                                       # Advance ratio J = V_inf/(n*D)
        self.D = D                                                       # Propeller diameter [m]
        self.altitude = altitude                                         # Altitude [m]
        self.RPM = RPM                                                   # Revolutions per minute 
        self.n =self.RPM/60                                              # Revolutions per second
        self.V_infty = self.J*self.n*self.D                              # Air speed [m/s]
        self.density = Atmosphere(self.altitude).density                 # Air density [kg/m^3]

        


    def Thrust(self):
         """Thrust is a function that calculates the thrust of a propeller at given altitudes and at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.                                        

         It gives as output the dimensional thrust T [N].
          
         The theory behind this code can be found in: Renato Tognaccini - "Lezioni di aerodinamica dell'ala rotante -
         Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04 -
         Chapter 2 - Paragraphs 2.1 - page 21
          
         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         T=np.zeros((len(self.altitude),len(self.V_infty)))         # Initialize a matrix so that each row can contain the thrust vector 
                                                                    # at a given altitude and at the given velocities
         for i in range(len(self.altitude)):                        # For cycle to repeat the calculation at each altitude
            T[i,:] = self.Ct*self.density[i]*self.n**2*self.D**4    # Calculate the thrust   
         return T                                                   

    
    def Power(self):
         """Power is a function that calculates the power of a propeller at given altitudes and at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.                                        

         It gives as output the dimensional power P [W].
          
         The theory behind this code can be found in: Renato Tognaccini - "Lezioni di aerodinamica dell'ala rotante -
         Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04 -
         Chapter 2 - Paragraphs 2.1 - page 21
          
         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         P=np.zeros((len(self.altitude),len(self.V_infty)))         # Initialize a matrix so that each row can contain the power vector 
                                                                    # at a given altitude and for the given velocities
         for i in range(len(self.altitude)):                        # For cycle to repeat the calculation at each altitude
            P[i,:] = self.Cp*self.density[i]*self.n**3*self.D**5    # Calculate the power 
         return P                                                  
    
    def Efficiency(self):
         """Efficiency is a function that calculates the efficiency of a propeller at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]
            (the altitude doesn't affect the efficiency)        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.                                        

         It gives as output the  efficency eta.
          
         The theory behind this code can be found in: Renato Tognaccini - "Lezioni di aerodinamica dell'ala rotante -
         Eliche rotori ed aeromotori con un'introduzione all'aerodinamica instazionaria" - a.a. 2023/2024 - vsn 2.04 -
         Chapter 2 - Paragraphs 2.1 - page 21
          
         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         eta = self.J*self.Ct/self.Cp                               # Calculate the efficiency
         return eta


    def Thrust_plot(self,a):
         """Thrust:plot is a function that generates plots of the thrust of a propeller at given altitudes and at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) a; (it must be 1 if the user wants to save the figure as pdf)                                    

         It gives the plot of the thrust curves as a function of speed and at different altitudes. 

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         T=self.Thrust()                                            # Calculate thrust whith 'Thrust' function           
         for i in range(len(self.altitude)):                        # Plot the thrust curves as a function of speed at different altitudes 
            plt.plot(self.V_infty,T[i,:], linewidth=2,         
                     linestyle='solid',marker = 'o',
                     label="altitude: " +str(self.altitude[i])+"m")
            plt.xlabel('$V_\infty [m/s]$')
            plt.ylabel ('Thrust [N]')
            plt.grid(True)  
            
        
         plt.xlim([0,max(self.V_infty)*1.1])
         plt.ylim([0,np.max(np.max(T))*1.1])
         plt.legend()
         if a==1:                                                   # Save the figure in pdf format in the case a = 1
            plt.savefig('Thrust.pdf')

         plt.show()
        


    def Power_plot(self,a):
         """Power_plot is a function that generates plots of the power of a propeller at given altitudes and at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) a; (it must be 1 if the user wants to save the figure as pdf)                                    

         It gives the plot of the power curves as a function of speed and at different altitudes. 

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         P = self.Power()                                           # Calculate power whith 'Power' function      
         for i in range(len(self.altitude)):
            plt.plot(self.V_infty,P[i,:], linewidth=2,              # Plot the power curves as a function of speed and at different altitudes      
                      linestyle='solid',marker = 'o',
                      label="altitude: "
                      +str(self.altitude[i])+"m")
         plt.xlabel('$V_\infty [m/s]$')
         plt.ylabel ('Power [W]')
        
         plt.grid(True)
         plt.xlim([0,max(self.V_infty)*1.1])
         plt.ylim([0,np.max(np.max(P))*1.1])
         plt.legend()
         if a==1:                                                   # Save the figure in pdf format in the case a = 1
            plt.savefig('Power.pdf')
         plt.show()

    

    def Efficiency_plot(self,a):
         """Efficiency_plot is a function that generates plots of the efficiency of a propeller at given velocities. 
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m] 
            (the altitude doesn't affect the efficiency)       
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) a; (it must be 1 if the user wants to save the figure as pdf)                                    

         It gives the plot of the efficiency curve as a function of speed. 

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         eta = self.Efficiency()                                     # Calculate the efficiency with 'Efficiency' function
         plt.plot(self.V_infty,eta,color='black',                    # Plot the efficiency curve as a function of speed
                  linewidth=2, linestyle='solid',
                  marker = 'o')
         plt.xlabel('$V_\infty [m/s]$')
         plt.ylabel ('Efficiency')
         plt.xlim([0,max(self.V_infty)*1.1])
         plt.ylim([0,1.2])
         plt.grid(True)                                              # Save the fiugre in pdf format in the case a = 1
         if a==1:
            plt.savefig('Efficiency.pdf')
         plt.show()
        

    def interp_Thrust(self,V_input):
         """interp_Thrust interpolates the thrust values ​​at the various speeds so as to 
         be able to calculate the thrust value at any speed in the range of those studied.
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) V_input; speed at which the user wants to know the thrust.          [m/s]
                                         

         It gives the values of the thrust at the given speed and at the different altitudes. 

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         T=self.Thrust()                                           # Calculate thrust whith 'Thrust' function                                                     
         T_output = np.zeros(len(self.altitude))                   # Initialize a array so that it can contain the thrust values 
                                                                   # at given altitudes and at the given speed
         for i in range(len(self.altitude)):                       # For cycle to repeat the calculation at each altitude      
            T_output[i] = np.interp(V_input, self.V_infty,        
                 T[i,:], left=None, right=None, period=None)       # Through linear interpolation, it gives the thrust value at 
                                                                   # the desired speed 
         return T_output
    
    def interp_Power(self,V_input):
         """interp_Power interpolates the power values at the various speeds so as to 
         be able to calculate the power value at any speed in the range of those studied.
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]        
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) V_input; speed at which the user wants to know the power.          [m/s]
                                         

         It gives the values of the power at the given speed and at the different altitudes. 

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         P=self.Power()                                            # Calculate power whith 'Power' function
         P_output = np.zeros(len(self.altitude))                   # Initialize a array so that it can contain the power values 
                                                                   # at given altitudes and at the given speed
         for i in range(len(self.altitude)):                       # For cycle to repeat the calculation at each altitude
            P_output[i] = np.interp(V_input, self.V_infty, 
                    P[i,:], left=None, right=None, period=None)    # Through linear interpolation, it gives the power value at 
                                                                   # the desired speed 
         return P_output
    
    def interp_efficiency(self,V_input):
         """interp_efficiency interpolates the efficiency values at the various speeds so as to 
         be able to calculate the efficiency value at any speed in the range of those studied.
         It accepts as input:
         1) Ct; thrust coefficient.		                                               															                                                                       
         2) Cp; power coefficient.	                                                   								                                                                      
         3) J; advance ratio.                                                          
         4) altitude; altitudes at which to calculate the dimensional values.   [m]    
            (the altitude doesn't affect the efficiency)
         5) D; propeller diameter.                                              [m]                          
         6) RPM; revolutions per minute.    
         7) V_input; speed at which the user wants to know the power.          [m/s]
                                         

         It gives the value of the efficiency at the given speed.

         Author: Andrea Raia, Antonino Guida
         Date: 31/05/2024
         Version: 1.00
         """
         eta=self.Efficiency()                                     # Calculate efficiency whith 'Efficiency' function
         eta_output = np.interp(V_input, self.V_infty,
                 eta, left=None, right=None, period=None)          # Through linear interpolation, it gives the efficiency value at 
                                                                   # the desired speed     
         return eta_output
   

    
        
