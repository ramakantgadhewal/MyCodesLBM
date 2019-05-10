#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:02:35 2019

@author: maruthinh
"""
import numpy as np

def Grid1D(xl, xu, N):
    '''generates 1d mesh (N points) by taking xl=lower limit adnd xu=upper limi
    t of domain'''
    return np.linspace(xl, xu, N)


def Init1DRiemann( x, x0, rhol, ul, pl, rhor, ur, pr, rho, u, p):
    '''Initilalizes Riemann problem by taking values to the left and right states
    of the diaphragm'''
    for i in range(np.size(x)):
        if(x[i]<x0):
            rho[i] = rhol
            u[i] = ul
            p[i] = pl
        else:
            rho[i] = rhor
            u[i] = ur
            p[i] = pr

        
 
xll = 0.0
xup = 1.0    
N = 100
x0 = 0.5

rhol = 1.0; ul = 0.0; pl = 1.0
rhor = 0.125; ur = 0.0; pr = 0.1

x=Grid1D(xll, xup, N)

rho=np.zeros(np.size(x))
u=np.zeros(N)
p=np.zeros(N)

Init1DRiemann(x, x0, rhol, ul, pl, rhor, ur, pr, rho, u, p)

U = np.zeros((N,3))
G = np.zeros((N,3))
f1 = np.zeros((N,3))
f1eq = np.zeros((N,3))
f2 = np.zeros((N,3))
f2eq = np.zeros((N,3))
E = np.zeros(N)


a = np.sqrt((1.4*p)/rho)
E = p/(rho*0.4) + 0.5*np.multiply(u, u)

Lambda = np.zeros(N)

for i in range(N):
    Lambda[i] = max(abs(u[i]-a[i]), u[i], abs(u[i]+a[i]))
    
U[:,0] = rho
U[:,1] = np.multiply(rho,u)
U[:,2] = np.multiply(rho,E)

G[:,0] = np.multiply(rho,u)
G[:,1] = p + np.multiply(rho, u, u)
G[:,2] = np.multiply(p, u) + np.multiply(rho,u, E)

f1eq[:,0] = U[:,0]/2.0 - G[:,0]/(2.0*Lambda)
f1eq[:,1] = U[:,1]/2.0 - G[:,1]/(2.0*Lambda)
f1eq[:,2] = U[:,2]/2.0 - G[:,2]/(2.0*Lambda)

f2eq[:,0] = U[:,0]/2.0 + G[:,0]/(2.0*Lambda)
f2eq[:,1] = U[:,1]/2.0 + G[:,1]/(2.0*Lambda)
f2eq[:,2] = U[:,2]/2.0 + G[:,2]/(2.0*Lambda)


print("lambda=", Lambda)    

#f1[:,0] = rho/2.0 - (rho*u)/(2.0*)





    
dt=x[2]-x[1]
print(dt)

t=0




#f1eq[0] = 
while(t<0.2):
    t=t+dt


        

    
    

