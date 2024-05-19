""" 
File: Geometry.py
Author: Antonio Brunaccini
Creation Date: May 19, 2024

This code includes the geometrical parameters of the designed blade and some operational parameters within a python class.
In addition, there are functions capable of manipulating the geometry and providing an example of the blade visualization.
The generation of the geometry is based on the input of an arbitrary number of blade sections: chord adim, span adim, 
twist and airfoil must be provided for the same section. Regarding airfoils, NACA 4 and 5 digit are directly implemented, but custom airfoil
can be added manually. All input airfoil are reproduced in output as xFoil txt to be used in other functions. 
The rest of the sections are obtained by interpolation.
"""

#NOTE: REFERENCE FRAME: X (chordwise direction), Y (spanwise direction), Z (thickness direction)

#NOTE: lengths are in meters, twists are in deg

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pyvista as pv
import os
import pandas as pd

from math import cos, sin
from math import atan
from math import pi
from math import pow
from math import sqrt

# INPUTS 
class Geometry():

    def __init__(self, R, r_hub, N, RPM, r_R_known, c_R_known, beta_known, beta75, airfoil_known, pitch):
        self.R = R                              # global radius (float)
        self.r_hub = r_hub                      # hub radius (float)
        self.N = N                              # number of blades (int)
        self.RPM = RPM                          # round per minute (float)
        self.r_R_known = r_R_known              # known adimensional RADIUS in ascending order (array)
        self.c_R_known = c_R_known              # known adimensional CHORDS in ascending order (array)
        self.beta_known = beta_known            # known TWIST in ascending order WITH RESPECT TO THE PROPELLER PLANE (array)
        self.beta75 = beta75                    # assigned nominal twist at 75% span 
        self.airfoil_known =  airfoil_known     # array containing NACA 4 digit airfoil as strings
        self.pitch = pitch                      # blade pitch (float)
        self.airfoil_dir()

        # cubic interp adim radius-adim chord/twist (used in BEMT module)
        self.fc = interp1d(self.r_R_known,self.c_R_known, kind='cubic', fill_value='extrapolate') 
        self.fb = interp1d(self.r_R_known,self.beta_known, kind='cubic', fill_value='extrapolate')


    # this function return the twist distribution for a fixed twist at 75% span and
    # preserving the twist distribution function shape
    def twist_rspct_75(self):

        beta75_star = self.fb(0.75)
        beta_known_sign_zero = self.beta_known - beta75_star
        beta_known_sign = beta_known_sign_zero + self.beta75

        return beta_known_sign


    # creation of the directory where txt input are stored 
    # (created in the same path where the code is located)
    def airfoil_dir(self):
        name = 'airfoil_input'  
        dir = os.getcwd()
        path_dir = os.path.join(dir, name)
        if not os.path.exists(path_dir):
            os.makedirs(path_dir)

        return path_dir


    # read the array "airfoil known" and generates the corresponding .txt coordinate file.
    # NACA 4 and 5 digits are built-in. 
    # custom airfoil must be manually added
    # it doesn't matter in the number of points of each txt file ins't the same
    # all airfoil will be named "airfoil_{i}" where i stands for the index position in "airfoil_known" (from 0). ex: airfoil_0, airfoil_1 ... 
    # delete .txt file manually to restart
    def gen_AF_txt(self):

        path_dir = self.airfoil_dir()

        for i, airfoil in enumerate(self.airfoil_known):

            if airfoil.isdigit() and len(airfoil) == 4:

                AF = self.naca4(airfoil, n=64)
                df = pd.DataFrame(AF, columns=['x', 'y'])
                name = f"airfoil_{i}.txt"
                path_file = os.path.join(path_dir, name)
                with open(path_file, 'w') as file:
                    file.write(f"NACA {airfoil}\n")
                    for index, row in df.iterrows():
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")
                        

            elif airfoil.isdigit() and len(airfoil) == 5:

                AF = self.naca5(airfoil, n=64)
                df = pd.DataFrame(AF, columns=['x', 'y'])
                name = f"airfoil_{i}.txt"
                path_file = os.path.join(path_dir, name)
                with open(path_file, 'w') as file:
                    file.write(f"NACA {airfoil}\n")
                    for index, row in df.iterrows():
                        file.write(f"{row['x']:.6f} {row['y']:.6f}\n")
                           
            elif airfoil.lower() == 'custom':
                print(f"you have to insert airfoil_{i} manually")



    # NACA 4 DIGIT COORDINATE GENERATOR. 
    #Two columns(X and Z) array as output
    # number = 4 digit string. ex: '0012'
    # n = number of points
    # finite_TE -> for finite(true) or zero(false) thickness trailing edge
    # half_cosine_spacing -> X spacing based on half cosine law (true) or uniform (false) 
    def naca4(self,number, n=240, finite_TE = False, half_cosine_spacing=True):
        
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
    

    # NACA 5 DIGIT COORDINATE GENERATOR. 
    #Two columns(X and Z) array as output
    # number = 4 digit string. ex: '0012'
    # n = number of points
    # finite_TE -> for finite(true) or zero(false) thickness trailing edge
    # half_cosine_spacing -> X spacing based on half cosine law (true) or uniform (false) 
    def naca5(self, number, n=240, finite_TE = False, half_cosine_spacing = True):
        
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


    # SCALE AIRFOIL FUNCTION
    #AF_xz = two column array (X and Z)
    #C = reference chord to which the airfoil is generated 
    def AF_scale(self,AF_xz,c):

        C = 1                       
        AF_xz_scl = AF_xz * (c/C)

        return AF_xz_scl
        

    # TRASLATION AIRFOIL FUNCTION
    #AF_xz = two column array (X and Z)
    #COR = centre of rotation. Airfoil is traslated by COR to allineate airfoils' aerodynamic centres
    #z_COR = COR height, default is zero (not true for non symmetric airfoils)
    def AF_trasl(self,AF_xz, z_COR=0):
            
        chord = max(AF_xz[:,0]) - min(AF_xz[:,0])
        COR = (0.25*chord,z_COR)
        AF_xz_trasl = AF_xz - COR

        return AF_xz_trasl
    
    
    # ROTATION AIRFOIL FUNCTION
    #AF_xz = two column array (X and Z)
    #beta_query = twist angle of which the profile is rotated 

    def AF_rot(self,AF_xz,beta_query):

        beta_query = -np.radians(beta_query)

        rot_matrix = np.array([[np.cos(beta_query), -np.sin(beta_query)],
                                    [np.sin(beta_query), np.cos(beta_query)]])

        AF_xz_rot = np.dot(AF_xz,rot_matrix.T)
        
        return AF_xz_rot
        
    
    # returns the coordinates of each airfoil on a common n number of points so as to give freedom on the input txt
    # needed for interpolation.
    def adapt_AF_points(self, AF, n):

        a = int(len(AF[:,0])/2 +1)

        xa = np.linspace(0,pi,n+1)
        x = [(0.5*(1.0-cos(xx))) for xx in xa]  
        xx = x[::-1]
        f_interpolate = interp1d(AF[0:a,0], AF[0:a,1], kind='cubic')
        zz = f_interpolate(xx)
        AF1 = np.vstack((xx, zz)).T

        xx = x
        f_interpolate = interp1d(AF[a-1:,0], AF[a-1:,1], kind='cubic')
        zz = f_interpolate(xx)
        AF2 = np.vstack((xx, zz)).T

        AF = np.vstack((AF1, AF2[1:]))

        return AF


    # SPAN AIRFOIL FUNCTION
    # returns the X and Z co-ordinates at an arbitrary y station of the span using 1D cubic interpolation along the two directions.
    # 1) reads the airfoil_known array
    # 2) generates internally if the profiles are NACA on the common number of points n
    #    and/or adjusts the number of points of the custom profiles and finally interpolates 
    # 3) applies geometric transformations 
    # 4) interpolates
    # AF_xz = two column array (X and Z)
    # y_query = station to which interpolation is requested (dimensional)
    def AF_span_interp(self,y_query):

        n = 100

        x = []
        z = []
        if type(self.r_R_known) != 'array': self.r_R_known = np.array(self.r_R_known)
        y = self.R*self.r_R_known

        for i, airfoil in enumerate(self.airfoil_known):
            
            if airfoil.isdigit() and len(airfoil) == 4:

                AF = self.naca4(airfoil)
                AF = self.adapt_AF_points(AF, n)

            elif airfoil.isdigit() and len(airfoil) == 5:

                AF = self.naca5(airfoil)
                AF = self.adapt_AF_points(AF, n)

            elif airfoil.lower() == 'custom':

                path_dir = self.airfoil_dir()
                name = f"airfoil_{i}.txt"
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

        npoints = len(x[:,0])
        AFint = np.full((npoints, 2), np.nan)

        for i in range(npoints):

            xx = x[i,:]
            f=interp1d(y, xx, kind = 'cubic', fill_value='extrapolate')
            xint=f(y_query)
            AFint[i,0]=xint

            zz = z[i,:]
            f=interp1d(y, zz, kind = 'cubic', fill_value='extrapolate')
            zint=f(y_query)
            AFint[i,1]=zint
            
        return AFint

    
    # BLADE GENERATION FUNCTION
    # provides the visualization of the blade interpolating a 3D grid where X and Z are provided by generative functions and 
    # Y is a vector of unifrom spacing from hub radius to global radius.
    # It also provide an STL file of the blade for third-party software visualization
    #n = number of station for interpolation 
    #wireframe -> alternative blade visualization method
    def gen_blade(self,n=50, wireframe = True):

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
        