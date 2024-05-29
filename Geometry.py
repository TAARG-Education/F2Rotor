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
#x pitch            =   blade pitch (float)															                                                            x
#x airfoil_known    =   array containing airfoil info. It can include:                                                                                          x
#x                      1) 4 digit NACA. ex: '0012'                                                                                                             x
#x                      2) 5 digit NACA. ex: '23012'                                                                                                            x
#x                      3) arbitrary airfoil. In this case in mandatory to type 'custom'                                                                        x
#x																				                                                                                x
#x																			                                                                                    x
#x  WARNING 1:  hub(r_hub/R) and tip(1) sections MUST be set. Extrapolation doesn't work well.                                                                  x
#x                                                                                                                                                              x
#x  WARNING 2:  Only NACA airfoils are built-in generated. If you want to include CUSTOM airfoils you MUST  create                                              x
#x              a direcotory called 'airfoil_input' in the same path you launch the script and manually insert the txt file in this directory.                  x
#x              You also MUST rename the txt file as 'airfoil_i' where i indicstes the position in airfoil_known.                                               x
#x              ex. if 4 station are assigned: airfoil_known = ['0012', '23012', 'custom', '2412'] so the custom airfoil must be airfoil_3                      x
#x                                                                                                                                                              x
#x                                                                                                                                                              x
#x																				                                                                                x
#x 																			                                                                                    x
#x Author: Antonio Brunaccini.														                                                                            x
#x																				                                                                                x
#x Version: 1.2.0																		                                                                        x
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, splprep, splev
import pyvista as pv
import os
import pandas as pd
from math import cos, sin, atan, pi, pow, sqrt


class Geometry():

    def __init__(self, R, r_hub, N, RPM, r_R_known, c_R_known, beta_known, beta75, airfoil_known, pitch):
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

        # utility functions
        # fc,fb,fx,fz are all 1D interpolator generated when the Geometry class is defined. 
        # They are used as tools for BEMT module (fc,fb) and blade geeneration (fx,fz)
        self.fc = interp1d(self.r_R_known,self.c_R_known, kind='cubic', fill_value='none') 
        self.fb = interp1d(self.r_R_known,self.beta_known, kind='cubic', fill_value='none')
        self.fx = interp1d(self.R*self.r_R_known,self.x, kind='cubic', fill_value='none')
        self.fz = interp1d(self.R*self.r_R_known,self.z, kind='cubic', fill_value='none')


    def airfoil_dir(self):

        '''
        This function creates a directory called 'airfoil input' where airfoil (AF) coordinates text file are stored (Xfoil type).
        
        output:
        - path_dir: it is the path of the directory, it can be used in another function to generate and store txt 

        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        name = 'airfoil_input'  
        dir = os.getcwd()
        path_dir = os.path.join(dir, name)
        if not os.path.exists(path_dir):
            os.makedirs(path_dir)

        return path_dir

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

        n = 200

        x = []
        z = []

        if type(self.r_R_known) != 'array': self.r_R_known = np.array(self.r_R_known)

        for i, airfoil in enumerate(self.airfoil_known):
            
            if airfoil.isdigit() and len(airfoil) == 4:

                AF = self.naca4(airfoil)
                AF = self.adapt_AF_points(AF, n)

            elif airfoil.isdigit() and len(airfoil) == 5:

                AF = self.naca5(airfoil)
                AF = self.adapt_AF_points(AF, n)

            elif airfoil.lower() == 'custom':

                path_dir = self.airfoil_dir()
                name = f"airfoil_{i+1}.txt"
                path_file = os.path.join(path_dir, name)
                if os.path.exists(path_file):
                    data = pd.read_csv(path_file, skiprows=1, delim_whitespace=True, header=None)
                    AF = data.to_numpy()
                    AF = np.array(AF)
                    AF = self.adapt_AF_points(AF, n)


            AF = self.AF_scale(AF,self.R*self.c_R_known[i])
            AF = self.AF_trasl(AF)
            AF = self.AF_rot(AF,self.beta_known[i])

            x.append(AF[:,0])
            z.append(AF[:,1])

        x = np.column_stack(x)
        z = np.column_stack(z)

        return x, z


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

        beta75_star = self.fb(0.75)
        beta_known_sign_zero = self.beta_known - beta75_star
        beta_known_sign = beta_known_sign_zero + self.beta75

        return beta_known_sign
    
    
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

        path_dir = self.airfoil_dir()

        for i, airfoil in enumerate(self.airfoil_known):

            if airfoil.isdigit() and len(airfoil) == 4:

                AF = self.naca4(airfoil, n=64)
                df = pd.DataFrame(AF, columns=['x', 'y'])
                name = f"airfoil_{i+1}.txt"
                path_file = os.path.join(path_dir, name)
                with open(path_file, 'w') as file:
                    file.write(f"NACA {airfoil}\n")
                    for index, row in df.iterrows():
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")
                        

            elif airfoil.isdigit() and len(airfoil) == 5:

                AF = self.naca5(airfoil, n=64)
                df = pd.DataFrame(AF, columns=['x', 'y'])
                name = f"airfoil_{i+1}.txt"
                path_file = os.path.join(path_dir, name)
                with open(path_file, 'w') as file:
                    file.write(f"NACA {airfoil}\n")
                    for index, row in df.iterrows():
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")
                           
            elif airfoil.lower() == 'custom':
                print(f"you have to insert airfoil_{i+1} manually")



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
        
        m = float(number[0])/100.0
        p = float(number[1])/10.0
        t = float(number[2:])/100.0

        a0 = +0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = +0.2843

        if finite_TE:
            a4 = -0.1015 
        else:
            a4 = -0.1036 

        if half_cosine_spacing:
            beta = np.linspace(0,pi,n+1)
            x = [(0.5*(1.0-cos(xx))) for xx in beta]  
        else:
            x = np.linspace(0,1.0,n+1)

        yt = [5*t*(a0*sqrt(xx)+a1*xx+a2*pow(xx,2)+a3*pow(xx,3)+a4*pow(xx,4)) for xx in x]

        xc1 = [xx for xx in x if xx <= p]
        xc2 = [xx for xx in x if xx > p]

        if p == 0:
            xu = x
            yu = yt

            xl = x
            yl = [-xx for xx in yt]

            xc = xc1 + xc2
            zc = [0]*len(xc)


        else:
            yc1 = [m/pow(p,2)*xx*(2*p-xx) for xx in xc1]
            yc2 = [m/pow(1-p,2)*(1-2*p+xx)*(1-xx) for xx in xc2]
            zc = yc1 + yc2

            dyc1_dx = [m/pow(p,2)*(2*p-2*xx) for xx in xc1]
            dyc2_dx = [m/pow(1-p,2)*(2*p-2*xx) for xx in xc2]
            dyc_dx = dyc1_dx + dyc2_dx

            theta = [atan(xx) for xx in dyc_dx]

            xu = [xx - yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
            yu = [xx + yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

            xl = [xx + yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
            yl = [xx - yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

        X = xu[::-1] + xl[1:]
        Z = yu[::-1] + yl[1:]

        AF_xz = np.column_stack((X, Z))

        return AF_xz
    

    
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
        
        naca1 = int(number[0])
        naca23 = int(number[1:3])
        naca45 = int(number[3:])

        cld = naca1*(3.0/2.0)/10.0
        p = 0.5*naca23/100.0
        t = naca45/100.0

        a0 = +0.2969
        a1 = -0.1260
        a2 = -0.3516
        a3 = +0.2843

        if finite_TE:
            a4 = -0.1015
        else:
            a4 = -0.1036  

        if half_cosine_spacing:
            beta = np.linspace(0.0,pi,n+1)
            x = [(0.5*(1.0-cos(x))) for x in beta]  
        else:
            x = np.allcloselinspace(0.0,1.0,n+1)

        yt = [5*t*(a0*sqrt(xx)+a1*xx+a2*pow(xx,2)+a3*pow(xx,3)+a4*pow(xx,4)) for xx in x]

        P = [0.05,0.1,0.15,0.2,0.25]
        M = [0.0580,0.1260,0.2025,0.2900,0.3910]
        K = [361.4,51.64,15.957,6.643,3.230]

        f1=interp1d(P, M, kind = 'cubic', fill_value='none')
        m = f1([p])[0]

        f2=interp1d(P, K, kind = 'cubic', fill_value='none')
        k1 = f2([m])[0]

        xc1 = [xx for xx in x if xx <= p]
        xc2 = [xx for xx in x if xx > p]
        xc = xc1 + xc2

        if p == 0:
            xu = x
            yu = yt

            xl = x
            yl = [-x for x in yt]

            zc = [0]*len(xc)
        else:
            yc1 = [k1/6.0*(pow(xx,3)-3*m*pow(xx,2)+ pow(m,2)*(3-m)*xx) for xx in xc1]
            yc2 = [k1/6.0*pow(m,3)*(1-xx) for xx in xc2]
            zc  = [cld/0.3 * xx for xx in yc1 + yc2]

            dyc1_dx = [cld/0.3*(1.0/6.0)*k1*(3*pow(xx,2)-6*m*xx+pow(m,2)*(3-m)) for xx in xc1]
            dyc2_dx = [cld/0.3*-(1.0/6.0)*k1*pow(m,3)]*len(xc2)

            dyc_dx = dyc1_dx + dyc2_dx
            theta = [atan(xx) for xx in dyc_dx]

            xu = [xx - yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
            yu = [xx + yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

            xl = [xx + yy * sin(zz) for xx,yy,zz in zip(x,yt,theta)]
            yl = [xx - yy * cos(zz) for xx,yy,zz in zip(zc,yt,theta)]


        X = xu[::-1] + xl[1:]
        Z = yu[::-1] + yl[1:]

        AF_xz = np.column_stack((X, Z))

        return AF_xz


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

        C = 1                       
        AF_xz_scl = AF_xz * (c/C)

        return AF_xz_scl
        

    def AF_trasl(self,AF_xz, z_COR=0):
        '''
        This function translate the airfoil by applying a rotation matrix with respect to a Centre Of Rotation (COR)
        It is needed to allineate airfoils' focuses and reduce hinge forces
    
        input variables: 
        - AF_xz = two column array containing airfoil coordinates
        - z_COR = height coordinate as DoF
        
        output:
        - AF_xz: 2D column array containing coordinates. Every row describes a point 
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''
         
        chord = max(AF_xz[:,0]) - min(AF_xz[:,0])
        COR = (0.25*chord,z_COR)
        AF_xz_trasl = AF_xz - COR

        return AF_xz_trasl
    
    
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

        beta_query = -np.radians(beta_query)

        rot_matrix = np.array([[np.cos(beta_query), -np.sin(beta_query)],
                                    [np.sin(beta_query), np.cos(beta_query)]])

        AF_xz_rot = np.dot(AF_xz,rot_matrix.T)
        
        return AF_xz_rot
        
    
    # returns the coordinates of each airfoil on a common n number of points so as to give freedom on the input txt
    # needed for interpolation.
    
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

        N=n//2

        x1 = 0.25 * (1 - np.cos(np.linspace(0, np.pi, N)))
        x2 = 0.5 + 0.25 * (1 - np.cos(np.linspace(0, np.pi, N)))
        x = np.hstack((x1,x2))

        pts = AF
        tck, u = splprep(pts.T, u=None, s=0.0, per=1) 
        u_new = x
        x_new, y_new = splev(u_new, tck, der=0)
        AF = np.column_stack((x_new, y_new))

        return AF

    
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

        #npoints = len(self.x[:,0])
        #AFint = np.full((npoints, 2), np.nan)

        xint=self.fx(y_query)
        xint = np.array(xint.reshape(-1, 1))

        zint=self.fz(y_query)
        zint = np.array(zint.reshape(-1, 1))

        AFint = np.hstack([xint,zint])

        return AFint
    

    
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

        AF = self.AF_span_interp(y_query)

        half_len = len(AF) // 2

        p1 = AF[:half_len]  
        p2 = AF[-1:-half_len-1:-1]  

        tk_max = max(np.sqrt(np.sum((p2 - p1) ** 2, axis=1)))

        return tk_max


    def gen_blade(self,n=50, wireframe = True):
        '''
        This function provides the interpolation for an arbitrary number of station 
        and the visualization of the blade in two ways.
        It also provides the STL file as output for further applications 
        
        input variables: 
        - n = number of span points 
        - wireframe = True for a second visualization type 
        
        output:
        Blade visualizatin and STL file
        
        Author: Antonio Brunaccini
        Date: 19/05/2024
        Version: 1.00
        '''

        yint_v = np.linspace(self.r_hub, self.R, n)

        x = []
        z = []
        y = []

        for i in range(len(yint_v)):

            y_int = yint_v[i]
            AFint = self.AF_span_interp(y_int)

            x.append(AFint[:,0])
            z.append(AFint[:,1])
            y.append(y_int*np.ones(len(AFint)))


        X = np.column_stack(x)
        Z = np.column_stack(z)
        Y = np.column_stack(y)

        wing = pv.StructuredGrid(X, Y, Z)
        wing = wing.triangulate()
        pv.save_meshio('blade.stl', wing)

        p1 = pv.Plotter()
        p1.add_mesh(wing, color='grey', show_edges=True)
        p1.show()

        if wireframe:

            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            ax.plot_wireframe(X,Y,Z, color='k', rstride=10, cstride=10)
            ax.axis('equal'); ax.axis('off')
            plt.show()
        