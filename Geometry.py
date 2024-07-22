#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#x This code includes the geometrical parameters of the designed blade and some operational parameters within a python class.                                   x
#x In addition, there are functions capable of manipulating the geometry and providing an example of the blade visualization.                                   x
#x The generation of the geometry is based on the input of an arbitrary number of blade sections: chord adim, span adim,                                        x
#x twist and airfoil must be provided for the same section. Regarding airfoils, NACA 4 and 5 digit are directly implemented, but custom airfoil                 x
#x can be added manually. All input airfoil are reproduced in output as xFoil txt to be used in other functions.                                                x
#x The rest of the sections are obtained by interpolation.																	                                    x
#x		                                                                                                                                                        x
#x  NOTE: REFERENCE FRAME: X (chordwise direction), Y (spanwise direction), Z (thickness direction), origin in the centre of the hub.                           x
#x  NOTE: lengths are in meters, twists are in deg                                                                                                              x
#x									                                                                                                                            x
#x Input variables:	                                                                                                                                            x
#x																	                                                                                            x
#x R                =   global radius (float)                                                                                                                   x
#x r_hub            =   hub radius (float)                                                                                                                      x
#x N                =   number of blades (int)                                                                                                                  x
#x RPM              =   round per minute (int)                                                                                                                  x
#x r_R_known        =   known adimensional RADIUS in ascending order (array)                                                                                    x
#x c_R_known        =   known adimensional CHORDS in ascending order (array)                                                                                    x
#x beta_known       =   known TWIST in ascending order WITH RESPECT TO THE PROPELLER PLANE (array)                                                              x
#x beta75           =   assigned nominal twist at 75% span (float)                                                                                              x
#x pitch            =   blade aerodymanic pitch (float)															                                                x
#x airfoil_known    =   array containing airfoil info. It can include:                                                                                          x
#x                      1) 4 digit NACA. ex: '0012'                                                                                                             x
#x                      2) 5 digit NACA. ex: '23012'                                                                                                            x
#x                      3) arbitrary airfoil. In this case in mandatory to type 'custom'                                                                        x
#x interp_kind      =   string containing the kind of the interpolator. You must use 'linear' if you have less than 4 stations besause                          x
#                       'cubic' needs at least 4 station to work, otherwise an error is dispayed 																x                  
#x																			                                                                                    x
#x  WARNING 1:  hub(r_hub/R) and tip(1) sections MUST be set. Extrapolation doesn't work well.                                                                  x
#x                                                                                                                                                              x
#x  WARNING 2:  Only NACA airfoils are built-in generated. If you want to include CUSTOM airfoils you MUST  girst create                                        x
#x              a direcotory called 'airfoil_input' in the same path you launch the script and manually insert the txt file in this directory.                  x
#x              You also MUST rename the txt file as 'airfoil_i' where i indicstes the position in airfoil_known.                                               x
#x              ex. if 4 station are assigned: airfoil_known = ['0012', '23012', 'custom', '2412'] so the custom airfoil must be airfoil_3                      x
#x              Start the geometry class only after this procedure .                                                                                            x
#x                                                                                                                                                              x
#x	WARNING 3:  Open the working folder in your code editor for proper functioning.																			    x
#x 																			                                                                                    x
#x Author: Antonio Brunaccini.														                                                                            x
#x																				                                                                                x
#x Version: 1.1.0	FIXED STL GENERATION SINGULARITY, ADDED CHORD COMPUTATION FUNCTION																	        x
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, splprep, splev
import pyvista as pv
import os
import pandas as pd
from math import cos, sin, atan, pi, pow, sqrt, dist


class Geometry():

    def __init__(self, R, r_hub, N, RPM, r_R_known, c_R_known, beta_known, beta75, airfoil_known, pitch, kind):

        self.R = R                                                             
        self.r_hub = r_hub                      
        self.N = N                              
        self.RPM = RPM                          
        self.r_R_known = r_R_known              
        self.c_R_known = c_R_known              
        self.beta_known = beta_known            
        self.beta75 = beta75                     
        self.airfoil_known =  airfoil_known     
        self.pitch = pitch                      
        self.airfoil_dir()
        self.x, self.z = self.init_blade()
        self.interp_kind = kind

        '''
        Utility functions:
        fc,fb,fx,fz are all 1D interpolator generated when the Geometry class is defined. 
        They are used as tools for BEMT module (fc,fb) and blade geometry generation (fx,fz)
        '''
        self.fc = interp1d(self.r_R_known,self.c_R_known, kind=self.interp_kind, fill_value='none') 
        self.fb = interp1d(self.r_R_known,self.beta_known, kind=self.interp_kind, fill_value='none')
        self.fx = interp1d(self.R*self.r_R_known,self.x, kind=self.interp_kind, fill_value='none')
        self.fz = interp1d(self.R*self.r_R_known,self.z, kind=self.interp_kind, fill_value='none')


    def airfoil_dir(self):

        '''
        This function creates a directory called 'airfoil input' where airfoil (AF) coordinates text file are stored (Xfoil type).
        
        output:
        - path_dir: it is the path of the directory, it can be used in another function to generate and store txt 

        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        name = 'airfoil_input'              # name of the directory where txt airfoil files are stored
        dir = os.getcwd()                   # get the current working directory 
        path_dir = os.path.join(dir, name)  # concatenate airfoil_input name to the working directory to make a complete path
        if not os.path.exists(path_dir):    # create the directory if it doesn't already exist
            os.makedirs(path_dir)           

        return path_dir                     # return the path for further actions
    

    def init_blade(self):
        '''
        This function initialize the blade generation by creating airfoils (AF) imposed in airfoil_known.
        NACA airfoils generation is built-in while custom airfoil must be manually insert
        
        output:
        - x,z : multidimensional array containing x and z coordinates, for all assigned station (every column is a station)
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        n = 200                                                                                     # number of points

        x = []                                                                                      # chordwise coordinates matrix initialization
        z = []                                                                                      # thickness coordinates matrix initialization

        if type(self.r_R_known) != 'array': self.r_R_known = np.array(self.r_R_known)               # type control: make r_R_known a np.array if it not already 

        for i, airfoil in enumerate(self.airfoil_known):                                            # evaluate what is written in each element of the vector 
            
            if airfoil.isdigit() and len(airfoil) == 4:                                             # if there are 4 digits

                AF = self.naca4(airfoil)                                                            # then generate that profile with naca4 function
                AF = self.adapt_AF_points(AF, n)                                                    # and adapt it to n points (for generating the blade the code needs the same number of elements for each airfoil)

            elif airfoil.isdigit() and len(airfoil) == 5:                                           # if there are 4 digits

                AF = self.naca5(airfoil)                                                            # then generate that profile with naca4 function
                AF = self.adapt_AF_points(AF, n)                                                    # and adapt it to n points (for generating the blade the code needs the same number of elements for each airfoil)


            elif airfoil.lower() == 'custom':                                                       # if it is custom

                path_dir = self.airfoil_dir()                                                       # look at the path of "airfoil_input"
                name = f"airfoil_{i+1}.txt"                                                         # get the txt file name corresponding to the position of the custom airfoil
                path_file = os.path.join(path_dir, name)                                            # concatenate the txt file name to the working directory to make a complete path
                if os.path.exists(path_file):                                                       # if this path exists
                    data = pd.read_csv(path_file, skiprows=1, delim_whitespace=True, header=None)   # read the corresponding txt file 
                    AF = data.to_numpy()                                                            # convert coordinates to numpy elements
                    AF = np.array(AF)                                                               # convert coordinates to array
                    AF = self.adapt_AF_points(AF, n)                                                # and adapt it to n points (for generating the blade the code needs the same number of elements for each airfoil)


            AF = self.AF_scale(AF,self.R*self.c_R_known[i])                                         # scale the airfoil 
            AF = self.AF_trasl(AF)                                                                  # translate the airfoil 
            AF = self.AF_rot(AF,self.beta_known[i])                                                 # rotate the airfoil 

            x.append(AF[:,0])                                                                       # append chordwise elements for each airfoil 
            z.append(AF[:,1])                                                                       # append spanwise elements for each airfoil 

        x = np.column_stack(x)                                                                      # store in previous initialized matrix
        z = np.column_stack(z)                                                                      # store in previous initialized matrix

        return x, z                                                                                 # return both matrices


    def chord(self,AF):
        '''
        This function provides the chord computation without the need to provide the span station but based on the Euclidean distance 
        
        input variables: 
        - AF = two column array containing airfoil coordinates
        
        output:
        Computed chord length
        
        Author: Antonio Brunaccini
        Date: 27/06/2024
        Version: 1
        '''
        x1 = min(AF[:,0])                   # find the minimum x
        i1 = np.argmin(AF[:,0])             # find the corresponding index
        y1 = AF[i1,1]                       # use the index to find corresponding y

        x2 = max(AF[:,0])                   # find the maximum x
        i2 = np.argmax(AF[:,0])             # find the corresponding index
        y2 = AF[i2,1]                       # use the index to find corresponding y

        chord = dist([x1,y1], [x2,y2])      # compute the euclidean distance

        return chord
    

    def twist_rspct_75(self):
        '''
        This function returns the twist distribution for a given beta75 for an arbitrary beta distribution.
        

        input variables: 
        - beta75: internally declared
        
        output:
        - beta_known_sign: twist distribution with respect to adim radius and given blade nominal twist 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        beta75_star = self.fb(0.75)                                     # interpolate the 75% value of beta according to arbitrary beta_known array input
        beta_known_sign_zero = self.beta_known - beta75_star            # translate beta_known to obtain a twist angle of the 75% section equal to zero
        beta_known_sign = beta_known_sign_zero + self.beta75            # translate again to obtain the twist distribution according to nominal twist input beta75

        return beta_known_sign                                          # return the required twist distribution
    
    
    def gen_AF_txt(self):
        '''
        This function generates NACA airfoil txt coordinate files and stores them in 'airfoil_input'. 
        There is a warning for custom airfoil which have to be manually insert
        
        output:
        - AF txt files in Xfoil format: one row header and tab space x and z coordinated for all remaining rows.
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        path_dir = self.airfoil_dir()                                       # find the path according to self.airfoil_dir

        for i, airfoil in enumerate(self.airfoil_known):                    # evaluate what is written in each element of the vector 

            if airfoil.isdigit() and len(airfoil) == 4:                     # if there are 4 digits

                AF = self.naca4(airfoil, n=64)                              # then generate that profile with naca4 function
                df = pd.DataFrame(AF, columns=['x', 'y'])                   # create a two columns dataframe to store chordwise and thickness coordinates
                name = f"airfoil_{i+1}.txt"                                 # name the i_th airfoil
                path_file = os.path.join(path_dir, name)                    # concatenate the txt file name to the working directory to make a complete path
                with open(path_file, 'w') as file:                          # creates the file as specified in path_file
                    file.write(f"NACA {airfoil}\n")                         # writes the airfoil name in the header
                    for index, row in df.iterrows():                        # for each row in the dataframe df
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")      # write in the txt file the corresponding df coordinates
                        

            elif airfoil.isdigit() and len(airfoil) == 5:                   # if there are 5 digits

                AF = self.naca5(airfoil, n=64)                              # then generate that profile with naca4 function
                df = pd.DataFrame(AF, columns=['x', 'y'])                   # create a two columns dataframe to store chordwise and thickness coordinates
                name = f"airfoil_{i+1}.txt"                                 # name the i_th airfoil
                path_file = os.path.join(path_dir, name)                    # concatenate the txt file name to the working directory to make a complete path
                with open(path_file, 'w') as file:                          # creates the file as specified in path_file
                    file.write(f"NACA {airfoil}\n")                         # writes the airfoil name in the header
                    for index, row in df.iterrows():                        # for each row in the dataframe df
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")      # write in the txt file the corresponding df coordinates
                           
            elif airfoil.lower() == 'custom':                               # if it is custom
                print(f"you have to insert airfoil_{i+1} manually")         # you have to insert it manually



    def naca4(self,number, n=240, finite_TE = False, half_cosine_spacing=True):
        '''
        This function generates two columns array containing x anz z coordinates of the input NACA 4-digit airfoil
    
        input variables: 
        - number: 4 digits in a string. ex: '0012'
        - n: number of points to be generated
        - finite_TE: True if you want a closed trailing edge
        - half_cosine_spacing:True for fitting law to improve leading and trailing edge resolution
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''
        
        m = float(number[0])/100.0                                                          # max curvature in chord percentage
        p = float(number[1])/10.0                                                           # distance of the point of maximum curvature from the leading edge expressed as a chord percentage and in multiples of 10
        t = float(number[2:])/100.0                                                         # max thickness in chord percentage

        a0 = +0.2969                                                                        # coeffcients from NACA4 equation
        a1 = -0.1260                                                                        
        a2 = -0.3516                                                                        
        a3 = +0.2843                                                                        

        if finite_TE:                                                                       # last coefficient depending on the choice on the trailing edge
            a4 = -0.1015 
        else:
            a4 = -0.1036 

        if half_cosine_spacing:                                                             # chordwise direction spacing depending on the chosen spacing law
            beta = np.linspace(0,pi,n+1)
            x = [(0.5*(1.0-cos(xx))) for xx in beta]                                        # cosine spacing
        else:
            x = np.linspace(0,1.0,n+1)                                                      # uniform spacing

        yt = [5*t*(a0*sqrt(xx)+a1*xx+a2*pow(xx,2)+a3*pow(xx,3)+a4*pow(xx,4)) for xx in x]   #NACA4 generation equation 

        xc1 = [xx for xx in x if xx <= p]                                                   # split x with respect to p in two halves
        xc2 = [xx for xx in x if xx > p]

        if p == 0:                                                                          # if the airfoil is symmetric
            xu = x                                                                          # x upper is equal to x
            yu = yt                                                                         # y upper is equal to yt

            xl = x                                                                          # x lower is equal to x
            yl = [-xx for xx in yt]                                                         # y lower is equal to the inverse of yt

            xc = xc1 + xc2                                                                  # concatenate the two halves
            zc = [0]*len(xc)                                                                # array of zeros


        else:                                                                               # if the airfoil is NOT symmetric
            yc1 = [m/pow(p,2)*xx*(2*p-xx) for xx in xc1]                                    # represent the ordinates of the mean line of the NACA4, calculated separately for two chord segments: 
            yc2 = [m/pow(1-p,2)*(1-2*p+xx)*(1-xx) for xx in xc2]                            # one up to the point of maximum curvature and the other beyond this point
            zc = yc1 + yc2                                                                  # concatenate the two halves

            dyc1_dx = [m/pow(p,2)*(2*p-2*xx) for xx in xc1]                                 # curvature derivative first chord segment
            dyc2_dx = [m/pow(1-p,2)*(2*p-2*xx) for xx in xc2]                               # curvature derivative second chord segment
            dyc_dx = dyc1_dx + dyc2_dx                                                      # concatenate curvature derivatives

            theta = [atan(xx) for xx in dyc_dx]                                             # angle between x axis and tangents

            xu = [xx - yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]                        # upper x coordinates
            yu = [xx + yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]                       # upper thickness coordinates

            xl = [xx + yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]                        # lower x coordinates
            yl = [xx - yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]                       # lower thickness coordinates

        X = xu[::-1] + xl[1:]                                                               # concatenate the coordinates of the upper and lower points into a single list 
        Z = yu[::-1] + yl[1:]

        AF_xz = np.column_stack((X, Z))                                                     # make a 2d array: first column (chordwise), second column (thickness)

        return AF_xz                                                                        # return the 2d array 
    

    
    def naca5(self, number, n=240, finite_TE = False, half_cosine_spacing = True):
        '''
        This function generates two columns array containing x anz z coordinates of the input NACA 5-digit airfoil
    
        input variables: 
        - number: 5 digits in a string. ex: '23012'
        - n: number of points to be generated
        - finite_TE: True if you want a closed trailing edge
        - half_cosine_spacing:True for fitting law to improve leading and trailing edge resolution
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''
        
        naca1 = int(number[0])                                                                  # first digit, multiplied by 0.15 gives the theoretical optimum lift coefficient of design with ideal angle of attack
        naca23 = int(number[1:3])                                                               # second digit, multiplied by 5, gives the relative position, in percentage, of the point of maximum curvature along the chord with respect to the leading edge;
                                                                                                # third digit indicates whether the curvature is simple (0) or reflected (1);
        naca45 = int(number[3:])                                                                # fourth and fifth digits provide the maximum profile thickness as a percentage of the chord, as in NACA4.

        cld = naca1*(3.0/2.0)/10.0                                                              # conversion to apply procedure
        p = 0.5*naca23/100.0
        t = naca45/100.0

        a0 = +0.2969                                                                            # coeffcients from NACA4 equation
        a1 = -0.1260
        a2 = -0.3516
        a3 = +0.2843

        if finite_TE:                                                                           # last coefficient depending on the choice on the trailing edge
            a4 = -0.1015
        else:
            a4 = -0.1036  

        if half_cosine_spacing:                                                                 # chordwise direction spacing depending on the chosen spacing law
            beta = np.linspace(0.0,pi,n+1)
            x = [(0.5*(1.0-cos(x))) for x in beta]                                              # cosine spacing
        else:
            x = np.allcloselinspace(0.0,1.0,n+1)                                                # uniform spacing

        yt = [5*t*(a0*sqrt(xx)+a1*xx+a2*pow(xx,2)+a3*pow(xx,3)+a4*pow(xx,4)) for xx in x]       # NACA4 generation equation

        P = [0.05,0.1,0.15,0.2,0.25]                                                            # camber-line profile coefficients
        M = [0.0580,0.1260,0.2025,0.2900,0.3910]
        K = [361.4,51.64,15.957,6.643,3.230]

        f1=interp1d(P, M, kind = 'cubic', fill_value='none')                                    # 1D P and M interpolator
        m = f1([p])[0]                                                                          # interpolated m value for a giben p (1 element array)

        f2=interp1d(P, K, kind = 'cubic', fill_value='none')                                    # 1D P and K interpolator
        k1 = f2([m])[0]                                                                         # interpolated k1 value for a giben m (1 element array)

        xc1 = [xx for xx in x if xx <= p]                                                       # split x with respect to p in two halves
        xc2 = [xx for xx in x if xx > p]
        xc = xc1 + xc2

        if p == 0:                                                                              # if the airfoil is symmetric
            xu = x                                                                              # x upper is equal to x
            yu = yt                                                                             # y upper is equal to yt

            xl = x                                                                              # x lower is equal to x
            yl = [-x for x in yt]                                                               # y lower is equal to the inverse of yt

            zc = [0]*len(xc)                                                                    # array of zeros

        else:                                                                                   # if the airfoil is NOT symmetric

            yc1 = [k1/6.0*(pow(xx,3)-3*m*pow(xx,2)+ pow(m,2)*(3-m)*xx) for xx in xc1]           # represent the ordinates of the mean line of the NACA4, calculated separately for two chord segments:
            yc2 = [k1/6.0*pow(m,3)*(1-xx) for xx in xc2]                                        # one up to the point of maximum curvature and the other beyond this point
            zc  = [cld/0.3 * xx for xx in yc1 + yc2]                                            # concatenate the two halves

            dyc1_dx = [cld/0.3*(1.0/6.0)*k1*(3*pow(xx,2)-6*m*xx+pow(m,2)*(3-m)) for xx in xc1]  # curvature derivative first chord segment
            dyc2_dx = [cld/0.3*-(1.0/6.0)*k1*pow(m,3)]*len(xc2)                                 # curvature derivative second chord segment

            dyc_dx = dyc1_dx + dyc2_dx
            theta = [atan(xx) for xx in dyc_dx]                                                 # angle between x axis and tangents

            xu = [xx - yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]                            # upper x coordinates
            yu = [xx + yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]                           # upper thickness coordinates

            xl = [xx + yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]                            # lower thickness coordinates
            yl = [xx - yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]                           # lower thickness coordinates


        X = xu[::-1] + xl[1:]                                                                   # concatenate the coordinates of the upper and lower points into a single list 
        Z = yu[::-1] + yl[1:]

        AF_xz = np.column_stack((X, Z))                                                         # make a 2d array: first column (chordwise), second column (thickness)

        return AF_xz                                                                            # return the 2d array 


    def AF_scale(self,AF_xz,c):
        '''
        This function scale the airfoil according to the reference chord C (chord at which the airfoil is generated)
    
        input variables: 
        - c = dimensional chord of the airfoil
        - AF_xz = two column array containing airfoil coordinates
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        C = 1                           # reference chird to which the airfoil was generated    
        AF_xz_scl = AF_xz * (c/C)       # scaled airfoil with respect to the ratio of the chords

        return AF_xz_scl                # return the 2d array 
        

    def AF_trasl(self,AF):
        '''
        This function translate the airfoil by 25% of the chord (Centre Of Rotation, COR)
        It is needed to allineate airfoils' focuses and reduce hinge forces
    
        input variables: 
        - AF = two column array containing airfoil coordinates
        -x = set to 25% of chord airfoil
        - z = height coordinate as DoF
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 27/06/2024
        Version: 2.00
        '''

        chord = self.chord(AF)  # compute che chord length

        x = 0.25*chord
        z = 0               

        COR = (x,z)             # amplitude of how long you want the translation, originally a quarter of the chord
        AF_trasl = AF - COR     # translated airfoil with respect to COR

        return AF_trasl         # return the 2d array 
    
    
    def AF_rot(self,AF_xz,beta_query):
        '''
        This function rotates the airfoil by applying a rotation matrix.
    
        input variables: 
        - AF_xz = two column array containing airfoil coordinates
        - beta_query = twist angle of which the profile is rotated (deg)
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        beta_query = -np.radians(beta_query)                                        # angle of which you want to rotate the airfoil

        rot_matrix = np.array([[np.cos(beta_query), -np.sin(beta_query)],           # define a 2d rotation matrix
                                    [np.sin(beta_query), np.cos(beta_query)]])

        AF_xz_rot = np.dot(AF_xz,rot_matrix.T)                                      # matrix multiplication between the two-column airfoil vector and the rotation matrix to obtain the rotated airfoil
        
        return AF_xz_rot                                                            # return the 2d array
        
    
    
    def adapt_AF_points(self, AF, n):
        '''
        This function change the number of points of the input airfoil according to the input number n.
        According to this result, it is possible to import AF of arbitrary number of point and have in output Af described by the same nunmber of points, 
        necessary for blade generation (interpolation based on fixed number of points)
        It is mostly useful for custom airfoil beacuse, for example using AirfoilTools, AF coordinates are provided in a highly variable way.
        
        The function applies e double cosine spacing in order to have a fitted LE and TE, comuting a Bi-spline interpolation
        (in fact it is impossible ti order points taking into account that the x vector start from 1, comes to 0 and returns to 1)

        input variables: 
        - AF_xz = two column array containing airfoil coordinates
        - n = number of points 
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        N=n//2                                                      # define half of input number to equally divide upper and lower 

        x1 = 0.25 * (1 - np.cos(np.linspace(0, np.pi, N)))          # cosine spacing for the first half
        x2 = 0.5 + 0.25 * (1 - np.cos(np.linspace(0, np.pi, N)))    # cosine spacing for the second half
        x = np.hstack((x1,x2))                                      # concatenate to have fitted LE and TE for both halves
                                                    
        tck, u = splprep(AF.T, u=None, s=0.0, per=1)                # tck: a tuple containing the B-spline coefficients, u: calculated curve parameters
                                                                    # s:  specifies that the spline must pass exactly through the data points
                                                                    # per: specifies that the B-spline is periodic (suitable for closed forms such as an airfoil)
        x_new, y_new = splev(x, tck, der=0)                         # compute the interpolated coordinates with respect to B-spline coefficients
        AF = np.column_stack((x_new, y_new))                        # make a 2d array: first column (chordwise), second column (thickness)


        return AF                                                   # return the 2d array

    
    def AF_span_interp(self, y_query):
        '''
        This function provides interpolated coordinates for a generic span station in a two columns array
        performing 1D iterpolators fx and fz. 
        
        input variables: 
        - y_query = dimensional span station 
        
        output:
        - AF_in: 2D column array containing coordinates. Every row describes a point. Product of interpolation.
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        xint=self.fx(y_query)                   # interpolate chordwise points according to the input span
        xint = np.array(xint.reshape(-1, 1))    # ensures xint is an array with a single column 

        zint=self.fz(y_query)                   # interpolate thickness points according to the input span
        zint = np.array(zint.reshape(-1, 1))    # ensures zint is an array with a single column 

        AFint = np.hstack([xint,zint])          # concatenate in 2d array 

        return AFint                            # return the 2d array
    

    
    def AF_max_tk(self,y_query):
        '''
        This function computes the max thickness for a generic span station. It first performes the interpolation 
        of the AF. Then, thanking to the adapt points function (acring in the blade initialization), computes the 
        thickness for every station x (for each x station there is an upper and lower z, LOWER<--x-->UPPER) and 
        picks the maximum

        input variables: 
        - y_query = dimensional span station 
        
        output:
        - tk_max: max thickness
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        AF = self.AF_span_interp(y_query)                       # interpolate airfoil coordinates according to span input

        half_len = len(AF) // 2                                 # round to half the number of points of the airfoil 

        p1 = AF[:half_len]                                      # define first half
        p2 = AF[-1:-half_len-1:-1]                              # define second half

        tk_max = max(np.sqrt(np.sum((p2 - p1) ** 2, axis=1)))   # find max thickness comparing coupled upper and lower points (works properly even if the number of points is odd)

        return tk_max                                           # return the max thickness value


    def gen_blade(self,nspan=100, wireframe = True):
        '''
        This function provides the interpolation for an arbitrary number of station 
        and the visualization of the blade in two ways.
        It also provides the STL file as output for further applications 
        
        input variables: 
        - nspan = number of span points 
        - wireframe = True for a second visualization type 
        
        output:
        Blade visualizatin and STL file
        
        Author: Antonio Brunaccini
        Date: 27/06/2024
        Version: 1.10
        '''

        yint_v = np.linspace(self.r_hub, self.R, nspan)                     # set the number of span interpolation stations (unfiform spacing)

        x = []                                                              # init chordwise coordinates matrix
        z = []                                                              # init thickness coordinates matrix
        y = []                                                              # init spanwise coordinates matrix

        for i in range(len(yint_v)):                                        # for each span station 

            y_int = yint_v[i]                                               # set the span statin from the array
            AFint = self.AF_span_interp(y_int)                              # interpolate the airfoil
            AFint = np.delete(AFint,len(AFint)//2,axis=0)                   # delete the repeated point (from previous actions, repetition neeeded) 

            x.append(AFint[:,0])                                            # append coordinates in each matrix
            z.append(AFint[:,1])
            y.append(y_int*np.ones(len(AFint)))


        X = np.column_stack(x)                                              # define numpy matrices
        Z = np.column_stack(z)
        Y = np.column_stack(y)

        wing = pv.StructuredGrid(X, Y, Z)                                   # creates a structured grid using interpolated data X, Y, Z     
        wing = wing.triangulate()                                           # convert the grid from quad to tria elements as preparation to STL
        pv.save_meshio('blade.stl', wing)                                   # save the trias grid as STL file

        p1 = pv.Plotter()                                                   # creates a visualization object for Pyvista
        p1.add_mesh(wing, color='grey', show_edges=False)                    # adds the grid (wing) to the plotter with the gray color. True to show trias 
        p1.show()                                                           # show the grid

        if wireframe:                                                       # another kind of visualization 

            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})         # creates a figure and a 3D axis using Matplotlib
            ax.plot_wireframe(X,Y,Z, color='k', rstride=5, cstride=5)   # draws a wireframe of the surface defined by the coordinates X, Y, Z. rstride and cstride control the density of the wireframe lines.
            ax.axis('equal'); ax.axis('off')                                # sets the axes so that the units are equal in all directions, maintaining the correct proportions of the geometry
            plt.show()                                                      # show the wireframe figure
        



    

