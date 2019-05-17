#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:27:01 2019

@author: maruthinh
"""

import numpy as np

def CalEqbDistrFun(rho, vel, Q, Lambda, feq):
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
    Q = np.size(c_mom[0])
    vel1=np.zeros(np.shape(vel[0]))
    vel2=np.zeros(np.shape(vel[1]))
    for i in range(Q):
        #vel_mom[0][:][:]+=feq_mom[i]*c_mom[0][i]
        #vel_mom[1][:][:]+=feq_mom[i]*c_mom[1][i]
        vel1+=feq_mom[i]*c_mom[0][i]
        vel2+=feq_mom[i]*c_mom[1][i]
        
    vel[0] = vel1
    vel[1] = vel2