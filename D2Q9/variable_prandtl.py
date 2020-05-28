#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:18:25 2020

@author: nhmaruthi
"""

import numpy as np
import D2Q9
import sys

def feq_iso(PartVel, Weights, PrimVars,  Feq):
    '''to compute distribution function. Takes parti
    cle velocity, weights, and Primitive varibales 
    (rho, u, v) and gives Feq (equlibrium distributi
    on as output. PrimVars and Feq can be 1D array 
    or matrix (nD array). THis is taken from 
    Quasi-Equilibrium Lattice Boltzmann Method, arXiv, 2007'''
    
    # Usq = PrimVars[1]*PrimVars[1]+PrimVars[2]*PrimVars[2]
    
    if(PrimVars[0].any()<=0):
        sys.exit('negative density')
        
    if(PrimVars[0].any()<=0 ):
        sys.exit('negative density')
        
    for i in range(9):
        # UdotC = PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]
        Csq = PartVel[0][i]**2+PartVel[1][i]**2
        
        pbyrho = PrimVars[3]/PrimVars[0]
        term0 = (1-pbyrho)
        term1 = PrimVars[0]*term0**2.0
        term2 = (pbyrho/(2.0*term0))**Csq
        term3 = (4.0*pbyrho**2+Csq*(1.0-3.0*pbyrho))/(2.0*term0)
        term4 = PartVel[0][i]*PrimVars[0]*PrimVars[1]+PartVel[1][i]*PrimVars[0]*PrimVars[2]
        term5 = (PrimVars[0]*PrimVars[1])**2+(PrimVars[0]*PrimVars[1])**2
        
        if(term0.any()<=0):
            sys.exit('negative term0')
            
        if(term1.any()<=0):
            sys.exit('negative term1')
                        
        if(term2.any()<=0):
            sys.exit('negative term2')
        
        if(term3.any()<=0):
            sys.exit('negative term3')
            
        # if(term4.any()<=0):
            # sys.exit('negative term4')
            
        # Feq[i] = term1*term2*(1.0 + (term4/PrimVars[3]) + (1.0/(2.0*PrimVars[3]**2))*(term4**2 - term5*term3))
        Feq[i] = term1*term2*(1.0 + (term4/PrimVars[3]) + (1.0/(2.0*PrimVars[3]**2))*(term4**2 - term5*term3))
        
        if(Feq[i].any()<=0):
            sys.exit('negative f')
           
            
def collision(dt, Omega1, Tau1, Tau2, PartVel, Weights, Grids, PrimVars, Feq, Fquasi_eq):
    '''computes collision term of the LBM update 
    formula, i.e., f*=f+(alpha*beta)*(feq-f). Tak
    es beta, Particle Vel, Weights, Grid, PrimVar
    and Feq as inputs. This gives updated Grid as
    output'''
    
    # D2Q9.CalEqbDistrFun(PartVel, Weights, PrimVars, Feq)
    feq_iso(PartVel, Weights, PrimVars, Feq)
    
    quasi_feq(dt, Tau2, PartVel, Weights, PrimVars, Grids, Feq, Fquasi_eq)
    for i in range(9):
        Grids[i]+=Omega1*((Tau2/Tau1)*Feq[i]+((Tau2-Tau1)/Tau2)*Fquasi_eq[i]-Grids[i])
        # Grids[i]+=Omega1*(Feq[i]-Grids[i])
        
    
def GetMoments(PartVel, Grids, PrimVars):
    '''Get (rho, u, p) by taking moment
    of distribution fucntion. Takes part
    icle velocity, and Grids (distributi
    on fun) as input and gives PrimVars
    (rho, u, v) as output'''
    vel1=np.zeros(np.shape(PrimVars[1]))
    vel2=np.zeros(np.shape(PrimVars[2]))
    E=np.zeros(np.shape(PrimVars[3]))
    
    PrimVars[0] = Grids.sum(axis=0)
    
    for i in range(9):
        Csq = PartVel[0][i]**2+PartVel[1][i]**2
        vel1+=Grids[i]*PartVel[0][i]
        vel2+=Grids[i]*PartVel[1][i]
        E+=Grids[i]*Csq
        
    PrimVars[1] = vel1/PrimVars[0]
    PrimVars[2] = vel2/PrimVars[0]
    
    jsq = (PrimVars[0]*PrimVars[1])**2+(PrimVars[0]*PrimVars[2])**2
    
    PrimVars[3] = 0.5*(E-(jsq/PrimVars[0]))
    

def quasi_feq(dt, Tau2, PartVel, Weights, PrimVars, Grids, Feq, Fquasi_eq):
    
    Vsq = PrimVars[1]**2+PrimVars[2]**2
    
    for i in range(9):
        Csq = PartVel[0][i]**2+PartVel[1][i]**2
        Amm_11 = PrimVars[0]**2*Feq[i]
        Amm_12 = (PrimVars[0]*PrimVars[0]*PrimVars[1]*Feq[i])
        Amm_13 = (PrimVars[0]*PrimVars[0]*PrimVars[2]*Feq[i])
        Amm_14 = (PrimVars[0]*PrimVars[0]*PrimVars[3]*Feq[i])

        Amm_21 = Amm_12
        Amm_22 = (PrimVars[0]*PrimVars[1])**2*Feq[i]
        Amm_23 = PrimVars[0]*PrimVars[1]*PrimVars[0]*PrimVars[2]*Feq[i]
        Amm_24 = PrimVars[0]*PrimVars[1]*PrimVars[3]*Feq[i]
        
        Amm_31 = Amm_13
        Amm_32 = Amm_23
        Amm_33 = (PrimVars[0]*PrimVars[2])**2*Feq[i]
        Amm_34 = PrimVars[0]*PrimVars[2]*PrimVars[3]*Feq[i]
        
        Amm_41 = Amm_14
        Amm_42 = Amm_24
        Amm_43 = Amm_34
        Amm_44 = PrimVars[3]**2*Feq[i]
        
        q1 = (PartVel[0][i]-PrimVars[1])*(Csq+Vsq-2*(PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]))*Grids[i]
        q2 = (PartVel[1][i]-PrimVars[2])*(Csq+Vsq-2*(PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]))*Grids[i]
    
        q1_eq = (PartVel[0][i]-PrimVars[1])*(Csq+Vsq-2*(PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]))*Feq[i]
        q2_eq = (PartVel[1][i]-PrimVars[2])*(Csq+Vsq-2*(PartVel[0][i]*PrimVars[1]+PartVel[1][i]*PrimVars[2]))*Feq[i]
    
        N1 = q1+(dt/Tau2)*(q1-q1_eq)
        N2 = q2+(dt/Tau2)*(q2-q2_eq)
        
        N1_1 = (1.0-(dt/(2.0*Tau2+dt)))*N1+(dt/(2.0*Tau2+dt))*q1_eq
        N1_2 = (1.0-(dt/(2.0*Tau2+dt)))*N2+(dt/(2.0*Tau2+dt))*q2_eq
        
        Amn_11 = PrimVars[0]*N1_1;             Amn_12 = PrimVars[0]*N1_2        
        Amn_21 = PrimVars[0]*PrimVars[1]*N1_1; Amn_22 = PrimVars[0]*PrimVars[1]*N1_2
        Amn_31 = PrimVars[0]*PrimVars[2]*N1_1; Amn_32 = PrimVars[0]*PrimVars[2]*N1_2
        Amn_41 = PrimVars[3]*N1_1;             Amn_42 = PrimVars[3]*N1_2
        
        Ann_11 = N1_1**2; Ann_12 = N1_1*N1_2
        Ann_21 = Ann_12; Ann_22 = N1_2**2
        
        A1 = np.array([np.sum(Amm_11), np.sum(Amm_12), np.sum(Amm_13), np.sum(Amm_14), np.sum(Amn_11), np.sum(Amn_12)])
        A2 = np.array([np.sum(Amm_21), np.sum(Amm_22), np.sum(Amm_23), np.sum(Amm_24), np.sum(Amn_21), np.sum(Amn_22)])
        A3 = np.array([np.sum(Amm_31), np.sum(Amm_32), np.sum(Amm_33), np.sum(Amm_34), np.sum(Amn_31), np.sum(Amn_32)])
        A4 = np.array([np.sum(Amm_41), np.sum(Amm_42), np.sum(Amm_43), np.sum(Amm_44), np.sum(Amn_41), np.sum(Amn_42)])
        A5 = np.array([np.sum(Amn_11), np.sum(Amm_21), np.sum(Amn_31), np.sum(Amn_41), np.sum(Ann_11), np.sum(Ann_12)])
        A6 = np.array([np.sum(Amn_12), np.sum(Amm_22), np.sum(Amn_32), np.sum(Amn_42), np.sum(Ann_21), np.sum(Ann_22)])
        
        A = np.array([A1, A2, A3, A4, A5, A6])
        print(A)
        A = np.linalg.inv(A)
            
        MN = np.array([0,0,0,0,N1-q1_eq, N2-q2_eq])
        
        Lambda = np.matmul(A,MN)
        
        Fquasi_eq[i]=Feq[i]*(1.0+Lambda[0]*PrimVars[0]+Lambda[1]*PrimVars[0]*PrimVars[1]+Lambda[2]*PrimVars[0]*PrimVars[2]+Lambda[3]*PrimVars[3]+Lambda[4]*q1)
        
        
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    