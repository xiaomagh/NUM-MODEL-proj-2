# -*- coding: utf-8 -*-
"""
Student ID: 25806676
"""

from __future__ import division
from parameter import *
import matplotlib.pyplot as plt
import fomula


def interpolation(enta, u, v, uold, vold, r):
# Function for interpolation of eta, u and v on the depature points    
    uu = 1.5*u - 0.5*uold
    vv = 1.5*v -0.5*vold # simple linear ways to get u(*)
    if (r == 'enta'):
        s_b = s = enta.copy()
        m = len(s[0,:]); n = len(s[:,0])
        for i in range(n):
            for j in range(m):
            # the distance here means relative to the arrival points, not the grids!!
                x = j*d - 0.5*(uu[i,j+1]+uu[i,j])*dt # Calculate the depature
                y = i*d - 0.5*(vv[i+1,j]+vv[i,j])*dt # points' positions
                if (x <= 0):
                    if (y <= 0):
                        s_b[i,j] = s[0,0] # inside eta should keep the same with outside eta
                    elif (y >= n*d-d):
                        s_b[i,j] = s[nd-1,0] # inside eta should keep the same with outside eta
                    else:
                        a = y/d - int(y/d) # interpolation
                        s_b[i,j] = a*s[int(y/d)+1,0] + (1-a)*s[int(y/d),0]
                elif (x >= m*d-d):
                    if (y <= 0):
                        s_b[i.j] = s[0,nd-1] # inside eta should keep the same with outside eta
                    elif (y >= n*d-d):
                        s_b[i,j] = s[nd-1,nd-1] # inside eta should keep the same with outside eta
                    else:
                        a = y/d - int(y/d) # interpolation
                        s_b[i,j] = a*s[int(y/d)+1,nd-1] + (1-a)*s[int(y/d),nd-1]
                else:
                    if (y <= 0):
                        b = x/d - int(x/d) # interpolation
                        s_b[i,j] = b*s[0,int(x/d)+1] + (1-b)*s[0,int(x/d)]
                    elif (y >= n*d-d):
                        b = x/d - int(x/d) # interpolation
                        s_b[i,j] = b*s[nd-1,int(x/d)+1] + (1-b)*s[nd-1,int(x/d)]
                    else:
                        b = x/d - int(x/d)
                        a = y/d - int(y/d) # interpolation
                        s_b[i,j] = a*b*s[int(y/d)+1,int(x/d)+1] +\
                                   a*(1-b)*s[int(y/d)+1,int(x/d)] +\
                                   (1-a)*b*s[int(y/d),int(x/d)+1] +\
                                   (1-a)*(1-b)*s[int(y/d),int(x/d)]
    if (r == 'u'):
        s_b = s = u.copy() # the similar process to eta's
        m = len(s[0,:]); n = len(s[:,0])
        tax_b = np.zeros((n,m))
        for i in range(n):
            for j in range(m):
                x = j*d - uu[i,j]*dt
                if (j == 0):
                    y = i*d - 0.25*(vv[i,j]+0.+0.+vv[i+1,j])*dt
                elif (j == m-1):
                    y = i*d - 0.25*(vv[i,j-1]+0.+0.+vv[i+1,j-1])*dt
                else:
                    y = i*d - 0.25*(vv[i,j]+vv[i+1,j-1]+vv[i,j-1]+vv[i+1,j])*dt
                if (y <= -0.5*d or y >= n*d - 0.5*d or x <= 0 or x >= L):
                    s_b[i,j] = 0.
                elif (y > -0.5*d and y <= 0):
                    b = x/d - int(x/d);a = (y+0.5*d)/(0.5*d)
                    s_b[i,j] = a*(b*s[0,int(x/d)+1] + (1-b)*s[0,int(x/d)])
                elif (y < n*d - 0.5*d and y >= n*d - d):
                    b = x/d - int(x/d);a = (y-(n-1)*d)/(0.5*d)
                    s_b[i,j] = (1-a)*(b*s[n-1,int(x/d)+1] + (1-b)*s[n-1,int(x/d)])
                else:
                    b = x/d - int(x/d)
                    a = y/d - int(y/d)
                    s_b[i,j] = a*b*s[int(y/d)+1,int(x/d)+1] +\
                               a*(1-b)*s[int(y/d)+1,int(x/d)] +\
                               (1-a)*b*s[int(y/d),int(x/d)+1] +\
                               +(1-a)*(1-b)*s[int(y/d),int(x/d)]
                tax_b[i,j] = ta(y)
        
    if (r == 'v'):
        s_b = s = v.copy() # the similar process to eta's
        m = len(s[0,:]); n = len(s[:,0])
        for i in range(n):
            for j in range(m):
                if (i == 0):
                    x = j*d - 0.25*(0.+uu[i,j+1]+uu[i,j]+0.)*dt
                elif (i == n-1):
                    x = j*d - 0.25*(0.+uu[i-1,j+1]+uu[i-1,j]+0.)*dt 
                else:
                    x = j*d - 0.25*(uu[i,j]+uu[i-1,j+1]+uu[i-1,j]+uu[i,j+1])*dt
                y = i*d - vv[i,j]*dt
                if (x <= -0.5*d or x >= m*d - 0.5*d or y <= 0 or y >= L):
                    s_b[i,j] = 0.
                elif (x > -0.5*d and x <= 0):
                    a = y/d - int(y/d);b = (x+0.5*d)/(0.5*d)
                    s_b[i,j] = b*(a*s[0,int(x/d)+1] + (1-a)*s[0,int(x/d)])
                elif (x < m*d - 0.5*d and x >= m*d - d):
                    a = y/d - int(y/d);b = (x-(m-1)*d)/(0.5*d)
                    s_b[i,j] = (1-b)*(a*s[n-1,int(y/d)+1] + (1-a)*s[n-1,int(y/d)])
                else:
                    b = x/d - int(x/d)
                    a = y/d - int(y/d)
                    s_b[i,j] = a*b*s[int(y/d)+1,int(x/d)+1] +\
                               a*(1-b)*s[int(y/d)+1,int(x/d)] +\
                               (1-a)*b*s[int(y/d),int(x/d)+1] +\
                               (1-a)*(1-b)*s[int(y/d),int(x/d)]
        
    
    return s_b
            
            
#def calculate(day):
#    
#    nt = day*int(86400/dt)
#    nd = int(L/d)
#    # Initially the flow is at rest and velocity and elevation are zero
#    enta0 = np.zeros((nd,nd))
#    u0 = np.zeros((nd,nd+1))
#    v0 = np.zeros((nd+1,nd))
#    
#    # Create the arrays for Arakawa-C grid in forward-backward time scheme
#    enta1 = entab = enta0.copy()
#    u1 = ub = u0.copy()
#    v1 = vb = v0.copy()
#    
#    E = np.zeros(nt)
#    ee = np.zeros((nd,nd))
#    
#    # Boundary Conditions
#    u0[:,0] = 0.; u0[:,nd] = 0.
#    v0[0,:] = 0.; v0[nd,:] = 0.    
#    
#        
#    
#    
#    for t in range(2,nt,2):
#        enta1 = enta_cal(enta0, u0, v0, nd)        
#        u1 = u_cal(enta1, u0, v0, nd)
#        u1[:,0] = 0.; u1[:,nd] = 0.
#        v1 = v_cal(enta1, u1, v0, nd)
#        v1[0,:] = 0.; v1[nd,:] = 0.
#        
#        E[t] ,ee = TaskC(enta1, u1, v1, nd)
#        
#        # Re-calculate the fomula again using the previou results but v before u
#        enta0 = enta_cal(enta1, u1, v1, nd)
#        v0 = v_cal(enta0, u1, v1, nd)
#        v0[0,:] = 0.; v0[nd,:] = 0.
#        u0 = u_cal(enta0, u1, v0, nd)
#        u0[:,0] = 0.; u0[:,nd] = 0.
#        
#        E[t+1] ,ee = TaskC(enta0, u0, v0, nd)
#        
#
#    plt.figure(1)
#    plt.contourf(enta0,15)
#    plt.colorbar()
#    
#    plt.figure(2)
#    plt.plot(E, linewidth=2, color='black')
#    plt.show()
#    
#    plt.figure(3)
#    plt.contourf(ee,20)
#    plt.colorbar()
    
    
    
    
    
    
    
    