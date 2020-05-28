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
import matplotlib.pyplot as plt
import D2Q9 
import BCs
import variable_prandtl as vp
# import data_output_functions as dop

c = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0], [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
w = np.array([4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0])

    
Nx = 50; Ny = 50;
GhostX = 2
GhostY = 2
TotNx = Nx + GhostX
TotNy = Ny + GhostY
D = 2
Q = 9
Ma = 0.1
T0=1/3.0
a = 1/3**0.5
Re = 100
kn = Ma/Re
RefLen = 1.0
UtopWall = Ma*a
nu = (RefLen*UtopWall)/Re
dx = 1.0
dt = RefLen/Nx
tau = 3.0*nu
Pr=0.710
tau2 = 4.0*tau/Pr
# tau2=tau
tauNdim = tau/dt
omega1 = (2.0*dt)/(2.0*tau+dt)
beta = 1.0/(2.0*tauNdim+1.0)


SimulTime = 100000#(40*Nx)/UtopWall

rho_ini = 1.0;
p=rho_ini*T0
u_ini = 0.0
v_ini = 0.0


#initialization of arrays required 
Grid = np.zeros((Q, TotNx, TotNx))
PrimVar = np.zeros((4, TotNx, TotNx))
feq = np.zeros((Q, TotNx, TotNx))
fquasi_eq = np.zeros((Q, TotNx, TotNx))
#initialization

#setting u velocity to 1 and v velocity to zero

PrimVar[0] = rho_ini
PrimVar[1] = u_ini
PrimVar[2] = v_ini
PrimVar[3] = p

u0 = np.zeros((TotNx, TotNx))
v0 = np.zeros((TotNx, TotNx))
######First iteration#####
# D2Q9.CalEqbDistrFun(c, w, PrimVar, feq)
vp.feq_iso(c, w, PrimVar, feq)
# vp.quasi_feq(dt, tau2, c, w, PrimVar, Grid, feq, fquasi_eq)
Grid=feq
iter=0; err=0.0

while(iter<int(SimulTime)):
# while(err>1e-14):
    
#while(iter<int(100)):
    iter = iter+1
    # print(iter)
    
    # D2Q9.Collision(beta, c, w, tau, Grid, PrimVar, feq)
    vp.collision(dt, omega1, tau, tau2, c, w, Grid, PrimVar, feq, fquasi_eq)
    
    BCs.PrepareWall(Grid)
    D2Q9.Advection(Grid)
    # BCs.WallBcOnGridBounceBack(Grid)
    BCs.bc_periodic_x(Grid)
    BCs.bc_bounceBack_bottomWall(Grid)
    BCs.MovingWallBcOnGridBounceBack(UtopWall, w, c, PrimVar, Grid)
    vp.GetMoments(c, Grid, PrimVar)
    
    if(iter%100==0):        
        e1=np.sum(np.sqrt((PrimVar[1]-u0)**2+(PrimVar[2]-v0)**2))
        e2=np.sum(np.sqrt(PrimVar[1]**2+PrimVar[2]**2))
        err=e1/e2
        print('error is = ', err)
    u0=PrimVar[1]
    v0=PrimVar[2]
    
#to plot contours of normal velocity
x = np.linspace(1, Nx, Nx)    
y = np.linspace(1, Ny, Ny)    
X, Y = np.meshgrid(x, y)
vp.GetMoments(c, Grid, PrimVar)
VelNorm = np.sqrt(PrimVar[1, 1:-1, 1:-1]*PrimVar[1, 1:-1, 1:-1] + PrimVar[2, 1:-1, 1:-1]*PrimVar[2, 1:-1, 1:-1])

cp=plt.contourf(X, Y, np.transpose(VelNorm), 30)
plt.colorbar(cp)
plt.xlabel('x')
plt.ylabel('y')

plt.figure()
plt.plot(PrimVar[1,int(Nx/2),1:-1])
