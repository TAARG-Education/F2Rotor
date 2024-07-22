#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#x Calculate the performance of a propeller using various methods:												                                                x
#x         1. momentum theory																	                                                                x
#x         2. vortical theory and small disturbance														                                                        x
#x         3. vortical theory 																	                                                                x
#x																				                                                                                x
#x Reference: McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight"												                                            x
#x																				                                                                                x
#x Input variables:																		                                                                        x
#x z = height																			                                                                        x
#x dx = Pitch to move along the blade (Alert: The pitch is dimensionless with dx/R where R is the radius at the tip.)						                    x
#x J = advance ratio																		                                                                    x
#x Geom: objects of the "geom" library 																                                                            x
#x Aero: objects of the "Aero" library 																                                                            x
#x																				                                                                                x
#x Output variables:																		                                                                    x
#x ct = thrust coefficient 																	                                                                    x
#x cp = power coefficient																	                                                                    x
#x																				                                                                                x
#x Authors: Daniele Di Somma, Emanuele Viglietti.														                                                        x
#x																				                                                                                x
#x Version: 1.4.0																		                                                                        x
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

import numpy as np
import os 
from scipy.integrate import trapezoid as trapz

# utility functions
def initialize_vars(z,geom,dx,J):
    '''
    This function evaluates the blade characteristics.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    
    input variables: 
    - z: altitude
    - J: advance ratio 
    - D: blade diameter
    - geom class: further details in geometry class
    
    output:
    - D: blade diameter
    - omega: angular velocity
    - Vt: tangential velocity
    - Vinf: asymptotic velocity
    - lam
    - x_vec: blade discretization vector
    - x_hub: hub radius ratio
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    
    '''
 
    # diameter
    D = 2*geom.R
    x_hub = geom.r_hub/geom.R      			# hub fraction in percentage of total blade radius
    x_vec = np.arange(x_hub, 1+ dx/geom.R, dx/geom.R)   # blade sections
    if x_vec[-1] >= 1: x_vec[-1] = 0.99
    
    # velocities
    omega = geom.RPM*2*np.pi/60    # rotational velocity
    Vt = omega*geom.R              # tangential velocity
    Vinf = 2*geom.R*geom.RPM*J/60  # asymptotic velocity

    lam = Vinf/Vt       # lambda
    
    return D,omega,Vt,Vinf,lam,x_vec,x_hub


def section_characteristic(x,geom,Vinf,omega,aero):
    '''
    This function calculates the main characteristics of each blade section.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    input variables: 
    - x: radius ratio
    - omega: blade angular velocity
    - Vinf: asymptotic velocity 
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    
    output:
    - sigma: solidity
    - chord
    - cla: lift curve slope
    - Vr: effective velocity
    - phi: inflow angle
    - beta: twist angle
    - bo: zero-lift angle
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    '''
    
    beta = np.deg2rad(geom.fb(x))
    beta +=np.deg2rad(geom.pitch)
        
    # evaluate chord at r/R station (multiply for R because the interpolation returns c/R)
    chord = geom.fc(x)
    Vr = np.sqrt(Vinf**2 + (omega*geom.R*x)**2)  # effective velocity
    phi = np.arctan(Vinf/(omega*geom.R*x))
    
    # calculate solidity
    sigma = (geom.N*chord)/(np.pi*geom.R)

    # evaluate Lift curve slope and Alpha zero Lift using aero class object.
    if aero.aero_method == 1 or aero.aero_method == 2:
        cla = aero.aero_params['Cl_alpha']
        # if aero.M_corr:
        #     cla/=(1-(Vr/a_sound)**2)**0.5
        bo = aero.aero_params['alpha_0_lift']
    else: 
        cla, bo = aero.eval_lift_properties(x = x, Vr = Vr)
    
    # define beta with respect to zero-lift line
    beta -= bo

    return sigma, chord, cla, Vr, phi, beta, bo


def section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, bo, phi, Vinf, Vr, wa, wt, omega, curvature=True, thickness=True):
    '''
    This function allows the calculation of thrust and power coefficients for a given blade section adopting the momentum theory,
    with or without taking into account curvature and thickness effects.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    input variables: 
    - z: altitude
    - J: advance ratio 
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    
    output:
    - x: radius ratio 
    - dx: radius increment
    - J: advance ratio
    - lam:
    - sigma: solidity
    - chord
    - alpha_i: induced angle of attack
    - beta: twist angle
    - bo: zero-lift angle
    - phi: inflow angle
    - Vinf: asymptotic velocity
    - Vr: effective velocity
    - wa: axial induction
    - wt: tangential induction
    - omega: angular velocity
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    '''
    
    # curvature effect
    if curvature:
        dalpha_c = 0.25 * (np.arctan((Vinf + wa)/(omega*x*geom.R-2*wt)) - np.arctan((Vinf + wa)/(omega*x*geom.R)))
    else:
        dalpha_c = 0
        
    # thickness effect
    if thickness:
        tmax_c = geom.AF_max_tk(x*geom.R)/chord
        dalpha_t = (4/15) * (lam*sigma/(lam**2 + x**2)) * (tmax_c)
    else:
        dalpha_t = 0
    
    # calculate effective AoA with respect to zero lift line
    alpha = beta - phi - alpha_i - dalpha_c - dalpha_t
    # correct AoA for aligning with the chord
    alpha += bo
    
    # calculate airfoil coefficients using aero class object
    if aero.aero_method == 1:
        cl, cd = aero.clcd1(Re_ref =1e+6 , Vr = Vr ,chord = chord, AoA=alpha)
    elif aero.aero_method == 2:
        cl,cd = aero.clcd2(Re_ref = 1e+6, Vr = Vr ,chord = chord, AoA = alpha)
    else: 
        cl, cd = aero.clcd3(AoA = np.rad2deg(alpha), Re_ref = 1e+6, Vr = Vr ,chord = chord, x = x)
    
    # calculate section contribute of ct and cp
    delct = (np.pi/8)*((J**2) + ((np.pi*x)**2))*sigma*((cl*np.cos(phi+alpha_i+dalpha_c+dalpha_t))
                                                    - (cd*np.sin(phi+alpha_i+dalpha_c+dalpha_t)))
        
    delcp = (np.pi/8)*(np.pi*x)*((J**2)+((np.pi*x)**2))*sigma*((cl*np.sin(phi+alpha_i+dalpha_c+dalpha_t))
                                                                + (cd*np.cos(phi+alpha_i+dalpha_c+dalpha_t)))
    return delct, delcp

    
def hub_loss(ct,J,D,geom):
    '''
    This function corrects the thrust coefficient introducing hub effect.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    
    input variables: 
    - ct: altitude
    - J: advance ratio 
    - D: blade diameter
    - geom class: further details in geometry class
    
    output:
    - ct: thrust coefficient
    - cp: power coefficient
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    
    '''
    
    # hub loss
    ct -= 0.5*np.pi*((geom.r_hub*J)**3)/(D**2)
    
    return ct


##########################################################################################################################################################

def BEMT_timp(z,J,dx,geom,aero, curvature=True, thickness=True, hub_corr=True):
    
    '''
    This function allows the calculation of thrust and power coefficients adopting the momentum theory.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    input variables: 
    - z: altitude
    - J: advance ratio 
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    - hub_corr:  boolean var, if 'True' applies hub correction
    
    output:
    - ct: thrust coefficient
    - cp: power coefficient
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    '''
    
    # initialize variables
    D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    delct = np.zeros(x_vec.shape)
    delcp = np.zeros(x_vec.shape)
    for int_x,x in enumerate(x_vec):
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta, bo = section_characteristic(x,geom,Vinf,omega,aero)
        
        # calculate alpha
        alpha_i = (0.5*np.sqrt((((lam/x)+((sigma*cla*Vr)/(8*Vt*x**2)))**2)+
               ((sigma*cla*Vr*(beta-phi))/(2*Vt*x**2)))) - 0.5*((lam/x)+((sigma*cla*Vr)/(8*Vt*x**2)))
    
        wt = Vr*alpha_i*np.sin(phi+alpha_i)
        wa = Vr*alpha_i*np.cos(phi+alpha_i)
        
        # evaluate section performance 
        delct[int_x], delcp[int_x] = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, bo, phi, Vinf, Vr, wa, wt, omega, curvature, thickness)
        
    # trapezoidal rule for integration
    ct = trapz(delct, x=x_vec)
    cp = trapz(delcp, x=x_vec)
    
    # hub correction
    if hub_corr:
        ct = hub_loss(ct,J,D,geom)
    
    return ct,cp

##########################################################################################################################################################

def BEMT_tvorpd(z,J,dx,geom,aero,curvature=True, thickness=True, hub_corr=True):

    '''
    This function allows the calculation of thrust and power coefficients adopting the vortex small disturbances theory.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    input variables: 
    - z: altitude
    - J: advance ratio 
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    - hub_corr:  boolean var, if 'True' applies hub correction

    output:
    - ct: thrust coefficient
    - cp: power coefficient
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    '''
    
    # initialize variables
    D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    delct = np.zeros(x_vec.shape)
    delcp = np.zeros(x_vec.shape)
    for int_x,x in enumerate(x_vec):
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta, bo = section_characteristic(x,geom,Vinf,omega,aero)
        
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
        delct[int_x], delcp[int_x] = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, bo, phi, Vinf, Vr, wa, wt, omega, curvature, thickness)
        
    # trapezoidal rule for integration
    ct = trapz(delct, x=x_vec)
    cp = trapz(delcp, x=x_vec)
    
    # hub correction
    if hub_corr:
        ct = hub_loss(ct,J,D,geom)
    
    return ct,cp

##########################################################################################################################################################

def BEMT_tvor(z,J,dx,geom,aero,curvature=True, thickness=False, hub_corr=True):
    '''
    This function allows the calculation of thrust and power coefficients adopting the vortex theory.
    Theory in McCormick, B. W., (1967), "Aerodynamics of V/STOL Flight", p. ...
    input variables: 
    - z: altitude
    - J: advance ratio 
    - geom class: further details in geometry class
    - aero class: further details in aero class
    - curvature: boolean var, if 'True' applies curvature correction
    - thickness: boolean var, if 'True' applies thickness correction
    - hub_corr:  boolean var, if 'True' applies hub correction
    
    output:
    - ct: thrust coefficient
    - cp: power coefficient
    
    Authors: Daniele Di Somma, Emanuele Viglietti
    Date: 24/05/2024
    Version: 1.00
    '''
    
    # initialize variables
    D,omega,Vt,Vinf,lam,x_vec,x_hub = initialize_vars(z,geom,dx,J)

    # start loop over blade sections
    cp, ct = 0.,0.
    delct = np.zeros(x_vec.shape)
    delcp = np.zeros(x_vec.shape)
    for int_x,x in enumerate(x_vec):
        # evaluate section characteristics
        sigma, chord, cla, Vr, phi, beta, bo = section_characteristic(x,geom,Vinf,omega,aero)
        
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
        
        # calculate of exact values for wa and induced angle of attack (alpha_i)
        wa = 0.5 * (np.sqrt(Vinf**2 + 4 * wt * (omega * x * geom.R - wt)) - Vinf)
    
        alpha_i = np.arctan(wt / wa) - phi
        if (alpha_i<0):
            wt = wtsa
            wa = wasa
            alpha_i = alfaisa
        
        # evaluate section performance 
        delct[int_x], delcp[int_x] = section_performance(x, dx, J, lam, sigma, chord, geom, aero, alpha_i, beta, bo, phi, Vinf, Vr, wa, wt, omega, curvature, thickness)
        
    # trapezoidal rule for integration
    ct = trapz(delct, x=x_vec)
    cp = trapz(delcp, x=x_vec)
    
    # hub correction
    if hub_corr:
        ct = hub_loss(ct,J,D,geom)
    
    return ct,cp