import Geometry 
import matplotlib.pyplot as plt
import numpy as np


# Test case according to data from paper "WindPACT Reference Wind Turbines"    https://www.nrel.gov/docs/fy18osti/67667.pdf



# INPUT DATA
R = 25   
R_hub = 2.04      
N = 3          
RPM = 100
r_R_known = [0.0816,   0.1448,    0.2080,    0.2716,    0.3348,    0.3980,    0.4616,    0.5248,    0.5880,    0.6516,    0.7148,    0.7780,   0.8416,    0.9048,    0.9680,    1.0000]
c_R_known = [ 0.0556,    0.0648,    0.0736,    0.0784,    0.0736,    0.0688,    0.0640,    0.0592,    0.0548,    0.0500,    0.0448,    0.0408,    0.0364,    0.0323,    0.0280,         0.025]
beta_known = [11.1,11.1,11.1,10.41,8.38,6.35,4.33,2.85,2.22,1.58,0.95,0.53,0.38,0.23,0.08,0] 
beta75 = 20
airfoil_known = ['0030', 'custom','custom','custom','custom','custom','custom','custom','custom','custom','custom','custom','custom','custom','custom','23012']
pitch = 0. 
interp_kind = 'cubic'

# Class iniziatization
g = Geometry.Geometry( R, R_hub, N, RPM, r_R_known, c_R_known, beta_known, beta75, airfoil_known, pitch,interp_kind)


# DEMO of all functions

##################### CHORD AND TWIST INTERPOLANT DISTRIBUTIONS,  TWIST_RSPCT_75 #####################
r_R_interp = np.linspace(R_hub/R, max(r_R_known), 100)
c_R_interp = g.fc(r_R_interp)
beta_interp = g.fb(r_R_interp)
beta_new = g.twist_rspct_75()


plt.figure(figsize=(10, 6))
plt.plot(r_R_known, c_R_known, 'o', label='Original data')
plt.plot(r_R_interp, c_R_interp, '-', label='Interpolated data')
plt.xlabel('adimensional RADIUS')
plt.ylabel('adimensional CHORD')
plt.legend()
plt.title('Chord distribution')
#plt.show()

plt.figure(figsize=(10, 6))
plt.plot(r_R_known, beta_known, 'o', label='Original data')
plt.plot(r_R_interp, beta_interp, '-', label='Interpolated data')
plt.plot(r_R_known, beta_new, 'o', label='Shifted twist with respect to c_R = 0.75')

plt.xlabel('adimensional RADIUS')
plt.ylabel('TWIST')
plt.legend()
plt.title('Twist distribution')
#plt.show()




##################### AIRFOIL GENERATION #####################

#NACA0030 (built-in)
naca0030 = g.naca4('0030')
plt.figure(figsize=(10, 6))
plt.plot(naca0030[:,0],naca0030[:,1],'-')
plt.xlabel('x')
plt.ylabel('z')
plt.title('NACA0030 coordinates')
plt.axis('equal')
#plt.show()

#NACA23012 (built-in)
naca23012 = g.naca5('23012')
plt.figure(figsize=(10, 6))
plt.plot(naca23012[:,0],naca23012[:,1],'-')
plt.xlabel('x')
plt.ylabel('z')
plt.title('NACA23012 coordinates')
plt.axis('equal')
#plt.show()

# Custom airfoil are inserted as input files




##################### TXT GENERATION FOR NACA AIRFOILS  #####################

g.gen_AF_txt()

# just check the directory "airfoil_input"



#####################   ADAPT AIRFOIL NUMBER OF POINTS  #####################

naca23012_adapted = g.adapt_AF_points(naca23012, 40)

plt.figure(figsize=(10, 6))
plt.plot(naca23012[:,0],naca23012[:,1],'o', label='Generated', markersize=2)
plt.plot(naca23012_adapted[:,0],naca23012_adapted[:,1],'o', label='Adapted', markersize=4)
plt.xlabel('x')
plt.ylabel('z')
plt.title('NACA23012 points comparison')
plt.legend()
plt.axis('equal')
#plt.show()



#####################   BODY TRANSFORMATIONS  #####################

naca_23012_scaled = g.AF_scale(naca23012,0.4)
naca_23012_trasl = g.AF_trasl(naca_23012_scaled)
naca_23012_rot = g.AF_rot(naca_23012_trasl,20)

plt.figure(figsize=(10, 6))
plt.plot(naca23012[:,0],naca23012[:,1],'-', label='Scaled')
plt.plot(naca_23012_scaled[:,0],naca_23012_scaled[:,1],'-', label='Scaled')
plt.plot(naca_23012_trasl[:,0],naca_23012_trasl[:,1],'-', label='Translated')
plt.plot(naca_23012_rot[:,0],naca_23012_rot[:,1],'-', label='Rotated')
plt.plot(0, 0, 'ko',markersize= 4)  
plt.xlabel('x')
plt.ylabel('z')
plt.title('Body transformations')
plt.legend()
plt.axis('equal')
#plt.show()




#####################   SPAN INTERPOLATION  #####################

y_span = np.linspace(R,R_hub,10)

plt.figure(figsize=(10, 6))
for i in range(len(y_span)):

    y = y_span[i]
    airfoil = g.AF_span_interp(y)

    c = plt.cm.Reds(i / len(y_span))  # i è normalizzato tra 0 e 1


    plt.plot(airfoil[:,0],airfoil[:,1],'-', color=c ,label=f'y: {y:.2f}')
    plt.axis('equal')


plt.xlabel('x')
plt.ylabel('z')
plt.title('Interpolated airfoil over 10 stations')
plt.plot(0, 0, 'ko',markersize= 4)  
plt.legend()
#plt.show()



#####################   CHORD AND MAX THICKNESS COMPUTATION  #####################

y_span = np.linspace(R_hub,R,100)
chords = []
tks = []

for i in range(len(y_span)):

    y = y_span[i]
    airfoil = g.AF_span_interp(y)

    chord = g.chord(airfoil)
    chords.append(chord)

    tk = g.AF_max_tk(y)
    tks.append(tk)


chords = np.array(chords)
tks = np.array(tks)

c_R_interp2 = g.fc(y_span/R)
c_interp = R*c_R_interp2



plt.figure(figsize=(10, 6))
plt.plot(y_span, chords,'-',label='computed')
plt.plot(y_span, c_interp,'.',label='interpolated')

plt.xlabel('y')
plt.ylabel('c')
plt.title('Computated chords as distance')
#plt.show()

plt.figure(figsize=(10, 6))
plt.plot(y_span, tks,'-')
plt.xlabel('y')
plt.ylabel('Max thickness')
plt.title('Max thickness')


plt.show()


#####################   BLADE GENERATION  #####################

# close all previous plots to make the blade appear


g.gen_blade()

# for the .stl file just check the worjing folder. It is named as "blade.stl"


