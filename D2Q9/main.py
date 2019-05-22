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

c = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0], [0.0, 0.0, 
                  1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
w = np.array([4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 
                  1.0/36.0, 1.0/36.0, 1.0/36.0])

    
Nx = 46; Ny = 46;
GhostX = 2
GhostY = 2
TotNx = Nx + GhostX
TotNy = Ny + GhostY
D = 2
Q = 9
Ma = 0.1
a = 1/3**0.5
Re = 2000
kn = Ma/Re
RefLen = 1.0
UtopWall = Ma*a
nu = (RefLen*UtopWall)/Re
dx = 1.0
dt = RefLen/Nx
tau = 3.0*nu
tauNdim = tau/dt
beta = 1.0/(2.0*tauNdim+1.0)
Lambda = 1.0;

SimulTime = (60*Nx)/UtopWall

rho_ini = 1.0;
u_ini = 0.0
v_ini = 0.0


#initialization of arrays required 

Grid = np.zeros((Q, TotNx, TotNx))
PrimVar = np.zeros((3, TotNx, TotNx))
feq = np.zeros((Q, TotNx, TotNx))

#initialization

#setting u velocity to 1 and v velocity to zero
PrimVar[0, 1:-1, 1:-1] = rho_ini
PrimVar[1, 1:-1, 1:-1] = u_ini
PrimVar[2, 1:-1, 1:-1] = v_ini

######First iteration#####
D2Q9.CalEqbDistrFun(c, w, PrimVar, feq)
Grid[:, 1:-1, 1:-1] = feq[:, 1:-1, 1:-1]

iter=0
while(iter<int(SimulTime)):
#while(iter<int(100)):
    iter = iter+1
    print(iter)
    
    D2Q9.Collision(beta, c, w, tau, Grid[:, 1:-1, 1:-1], PrimVar[:, 1:-1, 1:-1], 
          feq[:, 1:-1, 1:-1])
    BCs.PrepareWall(Grid)
    D2Q9.Advection(Grid)
    BCs.WallBcOnGridBounceBack(Grid)
    BCs.MovingWallBcOnGridBounceBack(UtopWall, w, c, PrimVar, Grid)
    

#to plot contours of normal velocity
x = np.linspace(1, Nx, Nx)    
y = np.linspace(1, Ny, Ny)    
X, Y = np.meshgrid(x, y)
VelNorm = np.sqrt(PrimVar[1, 1:-1, 1:-1]*PrimVar[1, 1:-1, 1:-1] + PrimVar[2, 1:-1, 
                 1:-1]*PrimVar[2, 1:-1, 1:-1])
cp=plt.contour(X, Y, np.transpose(VelNorm))
plt.colorbar(cp)
    