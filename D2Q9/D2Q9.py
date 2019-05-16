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

def CalEqbDistrFun(vel, Q, Lambda, feq):
    c = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0], [0.0, 0.0, 
                  1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
    w = np.array([4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 
                  1.0/36.0, 1.0/36.0, 1.0/36.0])
    
    for i in range(Q):
        term1 = (3.0/Lambda**2)*(c[0][i]*vel[0]+c[1][i]*vel[1])
        term2 = (9.0/(2.0*Lambda**4))*(c[0][i]*vel[0]+c[1][i]*vel[1])**2
        term3 = - (3.0/(2.0*Lambda**2))*(vel[0]*vel[0]+vel[1]*vel[1])
        feq[i] = w[i]*rho*(1.0+term1+term2+term3)


def GetDenFromFeq(feq_mom):
    '''Thin function computes density (rho) by taking equlibrium distribution fun-
    ction (feq) as an input'''
    return feq_mom.sum(axis=0)
    
    
    

def GetVelFromFeq(feq_mom, vel):
    '''Thin function computes velocities (both u and v) by taking equlibrium distribution fun-
    ction (feq), particle velocity(c) as an input'''
    c_mom = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0], [0.0, 0.0, 
                  1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
    vel1=np.zeros(np.shape(vel[0]))
    vel2=np.zeros(np.shape(vel[1]))
    for i in range(Q):
        #vel_mom[0][:][:]+=feq_mom[i]*c_mom[0][i]
        #vel_mom[1][:][:]+=feq_mom[i]*c_mom[1][i]
        vel1+=feq_mom[i]*c_mom[0][i]
        vel2+=feq_mom[i]*c_mom[1][i]
        
    vel[0] = vel1
    vel[1] = vel2
    

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

#initialization
rho.fill(rho_ini)
#setting u velocity to 1 and v velocity to zero
vel[0,-1,:] = 1.0
vel[1,:,:] = 0.0

######First iteration#####
CalEqbDistrFun(vel, Q, Lambda, feq)
f = feq

    
iter=0
while(iter<10):
    iter = iter+1
    print(iter)
    
    fnew = f - (1/tau)*(f-feq)
    f = fnew
    WallBcOnGridBounceBack(f)
    MovingWallBcOnGridBounceBack(f,1.0)
    rho=GetDenFromFeq(feq)
    GetVelFromFeq(feq, vel)
    CalEqbDistrFun(vel, Q, Lambda, feq)

    
    
