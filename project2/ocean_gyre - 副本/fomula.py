# -*- coding: utf-8 -*-
"""
student number: 25806676
"""
from __future__ import division
from parameter import *
import matplotlib.pyplot as plt


def Arakawa(day, question):    # day is the number of the running days
    
    nt = day*int(86400/dt)
    # Initially the flow is at rest and velocity and elevation are zero
    enta0 = np.zeros((nd,nd))# enta is surface elevation
    u0 = np.zeros((nd,nd+1)) # u is zonal velocity
    v0 = np.zeros((nd+1,nd)) # v is meridional velocity
    
    # Create the arrays for Arakawa-C grid in forward-backward time scheme
    enta1 = enta0.copy()
    u1 = u0.copy()
    v1 = v0.copy()
    
    En = np.zeros(nt)  # En is total energy
    ee = np.zeros((nd,nd)) # ee is a 2-D array for storing the energy of each grid
    
    # Boundary Conditions
    u0[:,0] = 0.; u0[:,nd] = 0.
    v0[0,:] = 0.; v0[nd,:] = 0.    
    
    for t in range(0,nt,2):
        enta1 = enta_cal(enta0, u0, v0, nd) # calculate eta first       
        u1 = u_cal(enta1, u0, v0, nd) # use new result of eta to solve u
        u1[:,0] = 0.; u1[:,nd] = 0. # Reset the boundary conditions
        v1 = v_cal(enta1, u1, v0, nd)# use new result of eta and u to solve v
        v1[0,:] = 0.; v1[nd,:] = 0.
        
       
        En[t] ,ee = TaskC(enta1, u1, v1, nd) 
        
        # Re-calculate the fomula again using the previou results but v before u
        enta0 = enta_cal(enta1, u1, v1, nd)
        v0 = v_cal(enta0, u1, v1, nd) # at this time step calculate v before u
        v0[0,:] = 0.; v0[nd,:] = 0.
        u0 = u_cal(enta0, u1, v0, nd)
        u0[:,0] = 0.; u0[:,nd] = 0.
        
   
        En[t+1] ,ee = TaskC(enta0, u0, v0, nd)
        
    if (question == 'B' or question == 'D1' or question == 'D2'):
        return enta0, u0, v0
    if (question == 'C'):
        return En
    



def enta_cal(enta, u, v, nd):
# small function for calculating the value of eta
    d = dict['d']
    for i in range(nd):
        for j in range(nd):
            enta_next = enta[i,j] - H*dt*((u[i,j+1]-u[i,j])/d+(v[i+1,j]-v[i,j])/d)
            enta[i,j] = enta_next
            
    return enta
    
def u_cal(enta, u, v, nd):
# small function for calculating the value of u
    d = dict['d']
    for i in range(0,nd):
        for j in range(1,nd):
            y = (i+0.5)*d
            u_next = u[i,j] + f(y)*dt*(v[i,j]+v[i+1,j-1]+v[i,j-1]+v[i+1,j])/4.\
                     - g*dt*(enta[i,j] - enta[i,j-1])/d - gama*dt*u[i,j]\
                     + ta(y)*dt/(ro*H)
            u[i,j] = u_next
    return u
    
def v_cal(enta, u, v, nd):
# small function for calculating the value of v
    d = dict['d']
    for i in range(1,nd):
        for j in range(0,nd):
            y = (i-0.5)*d
            v_next = v[i,j] - f(y)*dt*(u[i,j]+u[i-1,j+1]+u[i-1,j]+u[i,j+1])/4.\
                     - g*dt*(enta[i,j] - enta[i-1,j])/d - gama*dt*v[i,j] 
            v[i,j] = v_next
            
    return v 
    
def TaskC(enta, u, v, nd):
# In this fuction just use numerical way to finish the integral work for calculate 
# the total energy of the perturbation E, and the energy 2-D distribution ee
    d = dict['d']
    ee = np.zeros((nd,nd))
    a = 0.
    for i in range(nd):
        for j in range(nd):
            ee[i,j] = 0.5*ro*(d**2)*(H*((0.5*(u[i,j+1]+u[i,j]))**2 + \
                     (0.5*(v[i+1,j]+v[i,j]))**2) + g*enta[i,j]**2)
            a = a + ee[i,j]
     
    E = a 
    
    return E, ee
    
def TaskD(day,question):
# This function aims to calculate the exact solution of ocean gyre, calculate 
# the difference between it and numerical results, and get the energy calculated 
# from the difference fields ee, and the total energy error E.
    d = dict['d']
    enta, u, v = Arakawa(day,question)
    entast = enta.copy(); ust = u.copy(); vst = v.copy()
    entap = enta.copy(); up = u.copy(); vp = v.copy()    
    
    C = 0.2/(np.pi*gama*ro*H)
    
    for i in range(nd):
        for j in range(nd):
            x = (j+0.5)*d
            y = (i+0.5)*d
            f1, f2 = f12(x/L)
            
            a = 0.5*L/d
            if (int(a) - a == 0.):
                eta = 0.5*(enta[int(a),0] + enta[int(a)-1,0])
            else:
                eta = enta[int(a),0]

            entast[i,j] = eta + C*(f0*L/g)*((gama/(f0*np.pi))*\
            f2*np.cos(np.pi*y/L) + f1/np.pi*(np.sin(np.pi*y/L)*(1+beta*y/f0)+\
            beta*L/(f0*np.pi)*np.cos(np.pi*y/L)))
            
            entap[i,j] = enta[i,j] - entast[i,j]
            
    for i in range(0,nd):
        for j in range(0,nd+1):
            y = (i+0.5)*d
            x = j*d
            f1, f2 = f12(x/L)
            ust[i,j] = -C*f1*np.cos(np.pi*y/L)
            up[i,j] = u[i,j] - ust[i,j]
    
    for i in range(0,nd+1):
        for j in range(0,nd):
            y = i*d
            x = (j+0.5)*d
            f1, f2 = f12(x/L)
            vst[i,j] = C*f2*np.sin(np.pi*y/L)
            vp[i,j] = v[i,j] - vst[i,j]
            
    E, ee = TaskC(entap, up, vp, nd)
    
    if (question == 'D1'):
        return enta, entast
    if (question == 'D2'):    
        return enta, entast, E, ee

    
def f12(x):
# small function to provide the parameter for task D    
    e = gama/(L*beta)
    
    a = (-1 - np.sqrt(1 + (2*np.pi*e)**2))/(2*e)
    b = (-1 + np.sqrt(1 + (2*np.pi*e)**2))/(2*e)
    
    f1 = np.pi*(1 + ((np.exp(a)-1)*np.exp(b*x)+(1-np.exp(b))*np.exp(a*x)) \
         /(np.exp(b) - np.exp(a)))
        
    f2 = ((np.exp(a)-1)*np.exp(b*x)+(1-np.exp(b))*np.exp(a*x)) \
         /(np.exp(b)-np.exp(a))
    
    return f1, f2




   