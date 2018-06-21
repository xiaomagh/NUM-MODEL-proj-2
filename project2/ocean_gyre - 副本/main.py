# -*- coding: utf-8 -*-
"""
Student ID: 25806676
The main code for plotting
"""
import matplotlib.pyplot as plt
from fomula import *
from parameter import *
import numpy as np
def main(day, question):  
    if (question == 'B'):
        enta, u, v = Arakawa(day, question)
        plt.figure(figsize = (14,5))
        plt.clf()
        plt.ion()
        p1 = plt.subplot(121)
        p1.plot(u[1,:],linewidth=2,color='blue')
        p1.set_xlabel('East-West grid index')
        p1.set_ylabel('Zonal velocity(m/s)')
        p2 = plt.subplot(122)
        p2.plot(v[:,1],linewidth=2,color='black')
        p2.set_xlabel('South-North grid index')
        p2.set_ylabel('Meridional velocity(m/s)')
        plt.show()
        plt.figure(figsize = (14,5))
        plt.clf()
        plt.ion()
        x = np.zeros(nd)
        for i in range(nd):
            x[i] = i+0.5
        p1 = plt.subplot(121)
        p1.plot(x, 0.5*(enta[nd/2,:]+enta[nd/2+1,:]),linewidth=2,color='r')
        p1.set_xlabel('East-West grid index')
        p1.set_ylabel('Surface elevation(m)')
        p2 = plt.subplot(122)
        et = p2.contourf(enta, cmap = plt.cm.get_cmap('RdYlBu_r'))
        p2.set_xlabel('East-West grid index')
        p2.set_ylabel('South-North grid index')
        plt.colorbar(et)
        plt.show()

    if (question == 'C'):
        En = Arakawa(day, question)# plot the time series of total energy
        plt.figure(1)
        plt.clf()
        plt.ion()
        plt.plot(En,linewidth = 2)
        plt.xlabel('Runnning Time-steps(dimensionless)');plt.ylabel('Total energy($J$)')
        plt.show()

    if (question == 'D1'):
        enta,entast = TaskD(day,question)
        plt.figure(figsize = (14,5))# plot the surface elevation of exact result
        plt.clf()                   # and numerical model
        plt.ion()
        p1 = plt.subplot(121)
        et = p1.contourf(enta, cmap = plt.cm.get_cmap('RdYlBu_r'))
        p1.set_xlabel('East-West grid index')
        p1.set_ylabel('South-North grid index')
        plt.colorbar(et)
        p2 = plt.subplot(122)
        etst = p2.contourf(entast, cmap = plt.cm.get_cmap('RdYlBu_r'))
        p2.set_xlabel('East-West grid index')
        p2.set_ylabel('South-North grid index')
        plt.colorbar(etst)
        plt.show()
    if (question == 'D2'):
        enta,entast,E,ee = TaskD(day,question)
        plt.figure(figsize = (14,5))# plot the total energy error and the difference
        plt.clf()                  # field between exact result and modeling
        plt.ion() 
        p1 = plt.subplot(121)
        et = p1.contourf(enta-entast, cmap = plt.cm.get_cmap('RdYlBu_r'))
        p1.set_xlabel('East-West grid index')
        p1.set_ylabel('South-North grid index')
        plt.colorbar(et)
        p2 = plt.subplot(122)
        etst = p2.contourf(ee, cmap = plt.cm.get_cmap('RdYlBu_r'))
        p2.set_xlabel('East-West grid index')
        p2.set_ylabel('South-North grid index')
        plt.colorbar(etst)
        plt.show()
        
        