#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:47:32 2019

@author: maruthinh
"""
import numpy as np
import D2Q9

def WallBcOnGridBounceBack(Grid):
    
    #bounce back from right wall
    Grid[3, -2, 2:-2] = Grid[1, -1, 2:-2]
    Grid[6, -2, 2:-2] = Grid[8, -1, 2:-2]
    Grid[7, -2, 2:-2] = Grid[5, -1, 2:-2]
    
    
    #left wall 
    Grid[1, 1, 2:-2] = Grid[3, 0, 2:-2]
    Grid[8, 1, 2:-2] = Grid[6, 0, 2:-2]
    Grid[5, 1, 2:-2] = Grid[7, 0, 2:-2]
    
    #bottom wall
    Grid[2, 2:-2, 1] = Grid[4, 2:-2, 0]
    Grid[6, 2:-2, 1] = Grid[8, 2:-2, 0]
    Grid[5, 2:-2, 1] = Grid[7, 2:-2, 0]
    
    
    
def MovingWallBcOnGridBounceBack(PrimVar, Grid, Uwall, c, w, feq):    
    '''implemented from Q. Zou, X.He: Physc Fluids, 1997'''
    
    #first set like a solid wall
    Grid[4, 2:-2, -2] = Grid[2, 2:-2, -1]
    Grid[8, 2:-2, -2] = Grid[6, 2:-2, -1]
    Grid[7, 2:-2, -2] = Grid[5, 2:-2, -1]
    
    #now apply moving wall boundary conditions
    PrimVar[0] = 1.0
    PrimVar[1, :, -2] = Uwall
    PrimVar[2] = 0.0
    
    D2Q9.CalEqbDistrFun(PrimVar, feq, c, w)
    #print(PrimVar,feq)
    denominator = np.zeros(np.shape((feq[0, :, -1])))
    factor = np.zeros(np.shape((Grid[0, :, -1])))
    denominator =1.0/(feq[7, 2:-2, -2]+feq[4, 2:-2, -1]+feq[8, 2:-2, -2])
    
    factor = Grid[2, 2:-2, -1] + Grid[6, 2:-2, -1] + Grid[5, 2:-2, -1]
    
    Grid[4, 2:-2, -2] = factor*denominator*feq[4, 2:-2, -2]
    Grid[8, 2:-2, -2] = factor*denominator*feq[8, 2:-2, -2]
    Grid[7, 2:-2, -2] = factor*denominator*feq[7, 2:-2, -2]
    #print(Grid)
    
def BcOnCorners(Grid, feq):
    #bottom left corner
    Grid[2, 1, 1] = Grid[4, 1, 0]
    Grid[5, 1, 1] = Grid[7, 1, 0]
    Grid[1, 1, 1] = Grid[3, 1, 0]
    Grid[8, 1, 1] = Grid[6, 1, 0]
    Grid[6, 1, 1] = Grid[8, 1, 0]
    
    
    #bottom right corner
    Grid[2, -2, 1] = Grid[2, -2, 0]
    Grid[6, -2, 1] = Grid[8, -2, 0]
    Grid[3, -2, 1] = Grid[1, -2, 0]
    Grid[5, -2, 1] = Grid[7, -2, 0]
    Grid[7, -2, 1] = Grid[5, -2, 0]
    
    #top left corner
    Grid[5, 1, -2] = Grid[7, 1, -1]
    Grid[7, 1, -2] = Grid[5, 1, -1]
    Grid[1, 1, -2] = Grid[3, 1, -1]
    Grid[8, 1, -2] = Grid[6, 1, -1]
    Grid[4, 1, -2] = Grid[2, 1, -1]
    
    rhoTemp = 0
    rhoTemp = sum(Grid[:,1,-2])
    Grid[:, 1, -2] = rhoTemp*feq[:, 1, -2]
    
    #top right corner
    Grid[8, -2, -2] = Grid[6, -2, -1]
    Grid[6, -2, -2] = Grid[8, -2, -1]
    Grid[3, -2, -2] = Grid[1, -2, -1]
    Grid[7, -2, -2] = Grid[5, -2, -1]
    Grid[4, -2, -2] = Grid[2, -2, -1]
    
    rhoTemp = 0
    rhoTemp = sum(Grid[:,-2,-2])
    Grid[:,-2, -2] = rhoTemp*feq[:,-2,-2]
    
    
    
    
    
    
    