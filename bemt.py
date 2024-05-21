# ########################## INPUT VARIABLES #####################################
# input variables: 
# CONDITION: height, Vinf, rpm, design pitch
#
# GEOMETRY: diameter, blades number
#           hub radius, pitch
#           calettamento distribution (as function)
#           chord distribution (as function)
#
# METHOD: 1. momentum theory
#         2. vortical theory and small disturbance
#         3. vortical theory 
#
# AERO:  cl, cd for each section of the blade   
#

import numpy as np
from scipy.optimize import newton
from ambiance import Atmosphere
# import aero
from modules import aero
import os 


# utility functions
def initialize_vars(z,geom,dx,J):
    # start function
    atmosphere = Atmosphere([z*1000])
    T, rho,mu,a_sound = atmosphere.temperature[0], atmosphere.density[0], atmosphere.dynamic_viscosity[0], atmosphere.speed_of_sound[0]
    ni = mu/rho
   
    D = 2*geom.R
    omega = geom.RPM*2*np.pi/60    # rotational velocity
    Vt = omega*geom.R        # tangential velocity

    Vinf = 2*geom.R*geom.RPM*J/60

    lam = Vinf/Vt       # lambda

    x_hub = geom.r_hub/geom.R      # hub fraction in percentage of total blade radius

    x_vec = np.arange(x_hub, 1+ dx/geom.R, dx/geom.R)   #vector station
    if x_vec[-1] >= 1: x_vec[-1] = 0.99
    
    return T,rho,ni,a_sound,D,omega,Vt,Vinf,lam,x_vec,x_hub

def section_characteristic(x,geom,Vinf,omega,ni,a_sound,aero):
    beta = np.deg2rad(geom.fb(x))
    beta +=np.deg2rad(geom.pitch)
        
    # evaluate chord at r/R station (multiply for R because the interpolation returns c/R)
    chord = geom.fc(x)
    Vr = np.sqrt(Vinf**2 + (omega*geom.R*x)**2)  # effective velocity
    phi = np.arctan(Vinf/(omega*geom.R*x))
    
    # calculate solidity
    sigma = (geom.N*chord)/(np.pi*geom.R)
    # evaluate aerodynamic parameters
    cla, bo = aero.eval_lift_properties(Re_ref =1e+6 ,Re =Vr*chord/ni,M = Vr/a_sound, alpha_0_lift = np.deg2rad(0), Cl_alpha = 6.28)
    
    # # twist correction
    # beta += bo
    
    return sigma, chord, cla, Vr, phi, beta

def section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, phi, Vinf, Vr, wa, wt, omega, ni, a_sound, curvature=True, thickness=True):
    # curvature effect
    if curvature:
        dalpha_c = 0.25 * (np.arctan((Vinf + wa)/(omega*x*geom.R-2*wt)) - np.arctan((Vinf + wa)/(omega*x*geom.R)))
    else:
        dalpha_c = 0
        
    # thickness effect
    tmax_c = float(geom.airfoil_known[0][2:])/100
    if thickness:
        dalpha_t = (4/15) * (lam*sigma/(lam**2 + x**2)) * (tmax_c)
    else:
        dalpha_t = 0
   
    # calculate effective AoA
    alpha = beta - phi - alpha_i - dalpha_c - dalpha_t
    print("alpha:",np.rad2deg(alpha))
    
    ###############################################################################
    ### AERODYNAMICS (CIRO AND DANIELE)
    # calculate airfoil coefficients 
    cl, cd = aero.eval_coeffs(Re_ref =1e+6 ,Re =Vr*chord/ni ,M = Vr/a_sound,AoA=alpha,
                              alpha_0_lift = np.deg2rad(0), Cl_alpha = 6.28, 
                              Cl_alpha_stall = 0.100, Cl_max = 1.5, Cl_min = -0.5, Cl_incr_to_stall = 0.100,
                              Cd_min = 0.0130, Cl_at_cd_min = 0.500, dCd_dCl2 = 0.004, Mach_crit = 0.8, Re_scaling_exp = -0.15,
                              polar_filename = 'airfoil_polar.dat')
 
    ###############################################################################
    
    ## PERFORMANCE ##
    delct = (np.pi/8)*((J**2) + ((np.pi*x)**2))*sigma*((cl*np.cos(phi+alpha_i+dalpha_c+dalpha_t))
                                                    - (cd*np.sin(phi+alpha_i+dalpha_c+dalpha_t)))*dx
        
    delcp = (np.pi/8)*(np.pi*x)*((J**2)+((np.pi*x)**2))*sigma*((cl*np.sin(phi+alpha_i+dalpha_c+dalpha_t))
                                                                + (cd*np.cos(phi+alpha_i+dalpha_c+dalpha_t)))*dx
        
    
    return delct, delcp, cl, cd
    
def blade_performance(ct,cp,rho,Vinf,J,D,geom, Hub_Loss=True):
    
    # hub loss
    if Hub_Loss: ct -= 0.5*np.pi*((geom.r_hub*J)**3)/(D**2)

    thrust = ct*rho*((geom.RPM/60)**2)*D**4
    pin = cp*rho*((geom.RPM/60)**3)*D**5
    puse = thrust*Vinf
    eta = J*ct/cp
    
    return ct,cp,thrust,pin,puse,eta


##########################################################################################################################################################

def BEMT_timp(z,J,dx,geom,aero, curvature=True, thickness=True):
    
    # initialize variables
    T,rho,ni,a_sound,D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    for x in x_vec:
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta = section_characteristic(x,geom,Vinf,omega,ni,a_sound,aero)
        
        # calculate alpha
        alpha_i = (0.5*np.sqrt((((lam/x)+((sigma*cla*Vr)/(8*Vt*x**2)))**2)+
               ((sigma*cla*Vr*(beta-phi))/(2*Vt*x**2)))) - 0.5*((lam/x)+((sigma*cla*Vr)/(8*Vt*x**2)))
    
        wt = Vr*alpha_i*np.sin(phi+alpha_i)
        wa = Vr*alpha_i*np.cos(phi+alpha_i)
        
        # evaluate section performance 
        delct, delcp, _, _ = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, phi, Vinf, Vr, wa, wt, omega, ni, a_sound, curvature, thickness)
        ct += delct
        cp += delcp
        
        # print("x: {:.2f}, beta: {:.2f}, AoA: {:.2f}, chord: {:.2f}, ct: {:.2f}, cp: {:.2f}".format(x, geom.fb(x), np.rad2deg(alpha), chord, ct, cp))
    
    # evaluate blade performance
    ct,cp,thrust,pin,puse,eta = blade_performance(ct,cp,rho,Vinf,J,D,geom) 
    
    return ct,cp,thrust,pin,puse,eta,rho

##########################################################################################################################################################

def BEMT_tvorpd(z,J,dx,geom,aero,curvature=True, thickness=True):
    
    # initialize variables
    T,rho,ni,a_sound,D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    for x in x_vec:
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta = section_characteristic(x,geom,Vinf,omega,ni,a_sound,aero)
        
        # calculate alpha
        # calculate inflow angle at blade tip
        phi_t = np.arctan(lam)
        
        # calculate Prandtl's correction function
        if phi_t == 0:
            f = 1
        else:
            f = (2/np.pi)*np.arccos(np.exp(((x-1)*geom.N)/(2*np.sin(phi_t))))
        
        # calculate induced AoA 
        alpha_i = 0.5*((np.sqrt((((lam/x)+((cla*sigma)/(8*x*f*np.cos(phi))))**2)
                            + ((sigma*cla*(beta-phi))/(2*x*f*np.cos(phi)))))-
                                ((lam/x)+((cla*sigma)/(8*x*f*np.cos(phi)))))
        
        # calculate induced velocities
        wt = Vr*alpha_i*np.sin(phi+alpha_i)
        wa = Vr*alpha_i*np.cos(phi+alpha_i)
        
        # evaluate section performance 
        delct, delcp, _, _ = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, phi, Vinf, Vr, wa, wt, omega, ni, a_sound, curvature, thickness)
        ct += delct
        cp += delcp
        
        # print("x: {:.2f}, beta: {:.2f}, AoA: {:.2f}, chord: {:.2f}, ct: {:.2f}, cp: {:.2f}".format(x, geom.fb(x), np.rad2deg(alpha), chord, ct, cp))
    
    # evaluate blade performance
    ct,cp,thrust,pin,puse,eta = blade_performance(ct,cp,rho,Vinf,J,D,geom) 
    
    return ct,cp,thrust,pin,puse,eta,rho


##########################################################################################################################################################

def BEMT_tvor(z,J,dx,geom,aero,curvature=True, thickness=True):
    
    # initialize variables
    T,rho,ni,a_sound,D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    cl_values = []
    cd_values = []
    for x in x_vec:
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta = section_characteristic(x,geom,Vinf,omega,ni,a_sound,aero)
        
        # calculate alpha
        # calculate inflow angle at blade tip
        phi_t = np.arctan(lam)
        
        # calculate Prandtl's correction function
        if phi_t == 0:
            f = 1
        else:
            f = (2/np.pi)*(np.arccos(np.exp(((x-1)*geom.N)/(2*np.sin(phi_t)))))
        
        # calculate induced AoA 
        alpha_i = 0.5*(np.sqrt((((lam/x)+((cla*sigma)/(8*x*f*np.cos(phi))))**2)
                            + (4*((sigma*cla*(beta-phi))/(8*x*f*np.cos(phi)))))-
                                ((lam/x)+((cla*sigma)/(8*x*f*np.cos(phi)))))
        if alpha_i<0:alpha_i = 0.001
        # calculate induced velocities
        wt = Vr*alpha_i*np.sin(phi+alpha_i)
        wa = Vr*alpha_i*np.cos(phi+alpha_i)
        alfaisa = alpha_i
        wtsa = wt
        wasa = wa

        i = 1
        while i <= 100:
            test1 = Vinf**2 + 4 * wt * ((omega * geom.R * x) - wt)
            test2 = lam**2 +(4 * (wt / Vt) * (x - (wt / Vt)))
            
            if test1 < 0 or test2 < 0:
                wt = wtsa
                wa = wasa
                break
            
            term1 = (np.sqrt((Vinf**2) + (4 * wt * ((omega * geom.R * x) - wt))) - Vinf)
            term2 = lam + np.sqrt((lam**2) + (4 * (wt / Vt) * (x - (wt / Vt))))
            
            aa = beta - (np.arctan((2 * wt) / term1))
            bb = np.sqrt((0.25 * (term2**2)) + ((x - (wt / Vt))**2))
            cc = 8 * x * f * wt / Vt
            
            # Additional calculations for da, db, dc
            da = ((-2) / ((term1**2) + (4 * wt**2))) * (term1 - (wt * (((2 * omega * geom.R * x) - (4 * wt)) / (term1 + Vinf))))
            db = (((term2 * ((x / Vt) - (2 * wt / (Vt**2)))) / (2 * (term2 - lam))) + (wt / Vt**2) - (x / Vt))/bb
            dc = 8 * x * f / Vt
            
            y = (sigma * cla * aa * bb) - cc
            dy = (sigma * cla * ((aa * db) + (bb * da))) - dc
            
            if abs(y / dy) <= 1e-6:
                wt -= y / dy
                break
            else:
                i += 1
                if i > 100:
                    print("Non ci sono soluzioni per wt")
                    break
                wt -= y / dy
        
        # Calcolo del valore esatto di wa
        wa = 0.5 * (np.sqrt(Vinf**2 + 4 * wt * (omega * x * geom.R - wt)) - Vinf)
        
        # Calcolo del valore esatto dell’alfa indotto
        # Calcolo del valore esatto dell’alfa indotto
        alpha_i = np.arctan(wt / wa) - phi
        if (alpha_i<0):
            wt = wtsa
            wa = wasa
            alpha_i = alfaisa
        
        
        # evaluate section performance 
        delct, delcp, cl,cd = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, phi, Vinf, Vr, wa, wt, omega, ni, a_sound, curvature, thickness)
        cl_values.append(cl)
        cd_values.append(cd)
        ct += delct
        cp += delcp
    
        # print("x: {:.2f}, beta: {:.2f}, AoA: {:.2f}, chord: {:.2f}, ct: {:.2f}, cp: {:.2f}".format(x, geom.fb(x), np.rad2deg(alpha), chord, ct, cp))
    
    # evaluate blade performance
    ct,cp,thrust,pin,puse,eta = blade_performance(ct,cp,rho,Vinf,J,D,geom) 
    
    return ct,cp,thrust,pin,puse,eta,rho,cl_values,cd_values