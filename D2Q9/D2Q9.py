#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:27:01 2019

@author: maruthinh
"""

import numpy as np

def CalEqbDistrFun(PartVel, Weights, PrimVars,  Feq):
    '''to compute distribution function. Takes parti
    cle velocity, weights, and Primitive varibales 
    (rho, u, v) and gives Feq (equlibrium distributi
    on as output. PrimVars and Feq can be 1D array 
    or matrix (nD array)'''
    for i in range(9):
        UdotC = PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]
        Usq = PrimVars[1]*PrimVars[1]+PrimVars[2]*PrimVars[2]
       
        Feq[i] = Weights[i]*PrimVars[0]*(1.0 + 3.0*UdotC - 
           1.5*Usq + 4.5*UdotC**2 + 4.5*UdotC**3 - 4.5*Usq*UdotC )
              
def GetMoments(PartVel, Grids, PrimVars):
    '''Get (rho, u, p) by taking moment
    of distribution fucntion. Takes part
    icle velocity, and Grids (distributi
    on fun) as input and gives PrimVars
    (rho, u, v) as output'''
    vel1=np.zeros(np.shape(PrimVars[1]))
    vel2=np.zeros(np.shape(PrimVars[2]))
    
    PrimVars[0] = Grids.sum(axis=0)
    
    for i in range(9):
        vel1+=Grids[i]*PartVel[0][i]
        vel2+=Grids[i]*PartVel[1][i]
 
    PrimVars[1] = vel1/PrimVars[0]
    PrimVars[2] = vel2/PrimVars[0]
    
def Collision(Beta, PartVel, Weights, Tau, Grids, PrimVars, Feq):
    '''computes collision term of the LBM update 
    formula, i.e., f*=f+(alpha*beta)*(feq-f). Tak
    es beta, Particle Vel, Weights, Grid, PrimVar
    and Feq as inputs. This gives updated Grid as
    output'''
    GetMoments(PartVel, Grids, PrimVars)
    CalEqbDistrFun(PartVel, Weights, PrimVars, Feq)
    
    alpha = 2.0
    for i in range(9):
        Grids[i]+=alpha*Beta*(Feq[i]-Grids[i])    

def Advection(Grids):
    '''Computes Grids after streaming
    '''
    len0, lenx, leny = np.shape(Grids)
    
    for j in range(1, leny-1):
        for i in range(1, lenx-1):
            Grids[3, i, j] = Grids[3, i+1, j]        
            Grids[4, i, j] = Grids[4, i, j+1]        
            Grids[7, i, j] = Grids[7, i+1, j+1]        
            Grids[8, i, j] = Grids[8, i-1, j+1]        
            
    for j in range(leny-2, 0, -1):
        for i in range(lenx-2, 0, -1):        
            Grids[1, i, j] = Grids[1, i-1, j]        
            Grids[2, i, j] = Grids[2, i, j-1]        
            Grids[5, i, j] = Grids[5, i-1, j-1]        
            Grids[6, i, j] = Grids[6, i+1, j-1]       