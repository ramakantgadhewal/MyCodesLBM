#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:47:32 2019

@author: maruthinh
"""
import numpy as np
import D2Q9

def PrepareWall(Grids):
   '''Values from interior domain are cop
   ied to ghost nodes for applyging bound
   ary conditions. So that we will have 
   some relevant values at the ghose poi
   nts'''
    #copy from right wall to the right ghost
    Grids[:,-1,:] = Grids[:,-2,:]
    
    #copy from left wall to the left ghost
    Grids[:,0,:] = Grids[:,1,:]
    
    #copy from bottom wall to the bottom ghost
    Grids[:,:,-1] = Grids[:,:,-2]
    
    #copy from top wall to the top ghost
    Grids[:,:,0] = Grids[:,:,1]
    
def WallBcOnGridBounceBack(Grid):
    '''To apply wall boundary con
    ditions. Using bounce back method'''
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
    
    
    
def MovingWallBcOnGridBounceBack(Uwall, Weights, PartVel, PrimVars, Grids):    
    '''moving wall boundary conditions. Implemented from Q. Zou, X.He: Phy
    sc Fluids, 1997'''
    
    #first set like a solid wall
    Grids[4, 2:-2, -2] = Grids[2, 2:-2, -1]
    Grids[8, 2:-2, -2] = Grids[6, 2:-2, -1]
    Grids[7, 2:-2, -2] = Grids[5, 2:-2, -1]
    
    #now apply moving wall boundary conditions
    PrimVarTemp = np.array([1, Uwall, 0])
    
    feqTemp = np.zeros(9)
    D2Q9.CalEqbDistrFun(PartVel, Weights, PrimVarTemp, feqTemp)
    #print("priting in the BCs feq", feqTemp)

    denominator =1.0/(feqTemp[7]+feqTemp[4]+feqTemp[8])
    
    factor = Grids[2, 2:-2, -1] + Grids[6, 2:-2, -1] + Grids[5, 2:-2, -1]
    
    Grids[4, 2:-2, -2] = factor*denominator*feqTemp[4]
    Grids[8, 2:-2, -2] = factor*denominator*feqTemp[8]
    Grids[7, 2:-2, -2] = factor*denominator*feqTemp[7]
    #print(Grid)
    
        #bottom left corner
    Grids[2, 1, 1] = Grids[4, 1, 0]
    Grids[5, 1, 1] = Grids[7, 1, 0]
    Grids[1, 1, 1] = Grids[3, 1, 0]
    Grids[8, 1, 1] = Grids[6, 1, 0]
    Grids[6, 1, 1] = Grids[8, 1, 0]
    
    #bottom right corner
    Grids[2, -2, 1] = Grids[4, -2, 0]
    Grids[6, -2, 1] = Grids[8, -2, 0]
    Grids[3, -2, 1] = Grids[1, -2, 0]
    Grids[5, -2, 1] = Grids[7, -2, 0]
    Grids[7, -2, 1] = Grids[5, -2, 0]
    
    #top left corner
    Grids[5, 1, -2] = Grids[7, 1, -1]
    Grids[7, 1, -2] = Grids[5, 1, -1]
    Grids[1, 1, -2] = Grids[3, 1, -1]
    Grids[8, 1, -2] = Grids[6, 1, -1]
    Grids[4, 1, -2] = Grids[2, 1, -1]
    
    rhoTemp = 0
    rhoTemp = sum(Grids[:,1,-2])
    #print("rhoTemp=",rhoTemp)
    Grids[:, 1, -2] = feqTemp*rhoTemp
    
    #top right corner
    Grids[8, -2, -2] = Grids[6, -2, -1]
    Grids[6, -2, -2] = Grids[8, -2, -1]
    Grids[3, -2, -2] = Grids[1, -2, -1]
    Grids[7, -2, -2] = Grids[5, -2, -1]
    Grids[4, -2, -2] = Grids[2, -2, -1]
    
    rhoTemp = 0
    rhoTemp = sum(Grids[:,-2,-2])
    Grids[:,-2, -2] = feqTemp*rhoTemp
    
    
    
    
    
    
    
    
    