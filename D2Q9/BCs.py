#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:47:32 2019

@author: maruthinh
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