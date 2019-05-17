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


    
       
import numpy as np
import D2Q9
import BCs

Nx = 5; Ny = 5;
D = 2
Q = 9
Ma = 0.1
a = 1/3**0.5
Re = 100
kn = Ma/Re
RefLen = 1.0
UtopWall = Ma*a
nu = (RefLen*UtopWall)/Re
dx = 1.0
dt = RefLen/Nx
tau = (6.0*nu+1.0)/2.0
tauNdim = tau/dt
beta = 1.0/(2.0*tauNdim+1.0)
Lambda = 1.0;

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
while(iter<1):
    iter = iter+1
    print(iter)
    
    fnew = f - (1/tau)*(f-feq)
    f = fnew
    BCs.WallBcOnGridBounceBack(f)
    BCs.MovingWallBcOnGridBounceBack(f,0.5)
    D2Q9.GetDenFromFeq(feq)
    D2Q9.GetVelFromFeq(feq, vel)
    
    
