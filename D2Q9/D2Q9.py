#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:27:01 2019

@author: maruthinh
"""

import numpy as np

def CalEqbDistrFun(PrimVar, feq, c, w):
      for i in range(9):
          UdotC = c[0][i]*PrimVar[1, 1:-1, 1:-1]+c[1][i]*PrimVar[2, 1:-1, 1:-1]
          Usq = PrimVar[1, 1:-1, 1:-1]*PrimVar[1, 1:-1, 1:-1]+PrimVar[2, 1:-1, 
                       1:-1]*PrimVar[2, 1:-1, 1:-1]
       
          feq[i, 1:-1, 1:-1] = w[i]*PrimVar[0, 1:-1, 1:-1]*(1.0 + 3.0*UdotC - 
             1.5*Usq + 4.5*UdotC**2 + 4.5*UdotC**3 - 4.5*Usq*UdotC )

#def GetDenFromFeq(feq_mom):
#    '''Thin function computes density (rho) by taking equlibrium distribution fun-
#    ction (feq) as an input'''
#    return feq_mom.sum(axis=0)
#    
#    
#def GetVelFromFeq(feq_mom, vel):
#    '''Thin function computes velocities (both u and v) by taking equlibrium distribution fun-
#    ction (feq), particle velocity(c) as an input'''
#    c_mom = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0], [0.0, 0.0, 
#                  1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
#    Q = np.size(c_mom[0])
#    vel1=np.zeros(np.shape(vel[0]))
#    vel2=np.zeros(np.shape(vel[1]))
#    for i in range(Q):
#        #vel_mom[0][:][:]+=feq_mom[i]*c_mom[0][i]
#        #vel_mom[1][:][:]+=feq_mom[i]*c_mom[1][i]
#        vel1+=feq_mom[i]*c_mom[0][i]
#        vel2+=feq_mom[i]*c_mom[1][i]
#        
#    vel[0] = vel1
#    vel[1] = vel2
    
def GetMoments(c, Grid, PrimVar):
    
    PrimVar[:, 1:-1, 1:-1] = 0.0
    
    vel1=np.zeros(np.shape(PrimVar[1, 1:-1, 1:-1]))
    vel2=np.zeros(np.shape(PrimVar[2, 1:-1, 1:-1]))
    
    PrimVar[0, 1:-1, 1:-1] = Grid[:, 1:-1, 1:-1].sum(axis=0)
    
    for i in range(9):
        vel1+=Grid[i, 1:-1, 1:-1]*c[0][i]
        vel2+=Grid[i, 1:-1, 1:-1]*c[1][i]
 
    PrimVar[1, 1:-1, 1:-1] = vel1
    PrimVar[2, 1:-1, 1:-1] = vel2
    
    
    