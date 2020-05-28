#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 09:05:03 2020

@author: nhmaruthi
"""
import numpy as np


def vel_err(u, v, u0, v0):
    '''to compute error of velocity magnitude
    at the center of the domain'''
    e1=0.0; e2=0.0
    e1=np.sum(np.sqrt((u-u0)**2+(v-v0)**2))
    e2=np.sum(np.sqrt(u**2+v**2))
    print('error is = ', e1/e2)
    return e1/e2    
    
    