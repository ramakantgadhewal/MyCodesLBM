#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 13:30:36 2019
@author: maruthinh
This is to simulate Lid driven cavity using D2Q9 LBM method

Some important formulas:
    1. formula for kinematic viscosity nu = ((2*tau - 1)*(dx^2))/(6*dt)
    2. lattice speed: lambda = dx/dt
    3. equilibrium distribution function is: 
        feq_i = w_i*{3*[(e_i.u)/Lambda] + (9/2)*[(e_i.u)^2/Lambda^2] - (3/2)*[(e_i.u)/Lambda^2]}        
"""

def WallBcOnGridBounceBack(f):
    #bottom wall
    f[5,0] = f[7,0]
    f[2,0] = f[4,0]
    f[6,0] = f[8,0]
    
    #left wall 
    f[8,:,0] = f[6,:,0]
    f[1,:,0] = f[3,:,0]
    f[5,:,0] = f[7,:,0]
    
    #right wall
    f[6,:,-1] = f[8,:,-1]
    f[3,:,-1] = f[1,:,-1]
    f[7,:,-1] = f[7,:,-1]
    

def MovingWallBcOnGridBounceBack(f, Uwall):    
    '''implemented from Q. Zou, X.He: Physc Fluids, 1997'''
    rhoN = f[0,-1,:]+f[2,-1,:]+f[3,-1,:]+2*(f[2,-1,:]+f[6,-1,:]+f[5,-1,:])
    f[4,-1,:] = f[2,-1,:]
    f[7,-1,:] = f[5,-1,:] + 0.5*(f[1,-1,:]-f[3,-1,:])-0.5*rhoN*Uwall
    f[8,-1,:] = f[6,-1,:] + 0.5*(f[3,-1,:]-f[1,-1,:])+0.5*rhoN*Uwall
    
       
import numpy as np
import D2Q9

Nx = 5; Ny = 5;
D = 2
Q = 9
Re = 100
nu = 0.03
dx = 1.0
dt = 1.0
Lambda = 1.0;
tau = (6.0*nu+1.0)/2.0

rho_ini = 1.0;
u_ini = 0.0
v_ini = 0.0

#initialization of arrays required 
rho = np.zeros((Nx, Ny))
x = np.zeros((D, Nx, Ny))
vel = np.zeros((D, Nx, Ny))
feq = np.zeros((Q, Nx, Ny))
f = np.zeros((Q, Nx, Ny))
#initialization
rho.fill(rho_ini)
#setting u velocity to 1 and v velocity to zero
vel[0,-1,:] = 1.0
vel[1,:,:] = 0.0


######First iteration#####
D2Q9.CalEqbDistrFun(rho, vel, Q, Lambda, feq)
f.fill(1.0)

#CalEqbDistrFun(vel, Q, Lambda, feq)
#
#f = fnew
#WallBcOnGridBounceBack(f)
#MovingWallBcOnGridBounceBack(f,1.0)
#fnew = f - (1/tau)*(f-feq)
#GetDenFromFeq(feq)
#GetVelFromFeq(feq, vel)
#CalEqbDistrFun(vel, Q, Lambda, feq)

iter=0
while(iter<10):
    iter = iter+1
    print(iter)
    
    fnew = f - (1/tau)*(f-feq)
    f = fnew
    WallBcOnGridBounceBack(f)
    MovingWallBcOnGridBounceBack(f,0.5)
    D2Q9.GetDenFromFeq(feq)
    D2Q9.GetVelFromFeq(feq, vel)
    
    
