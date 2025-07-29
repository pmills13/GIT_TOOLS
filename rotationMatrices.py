# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 06:25:34 2024

@author: mfund
"""

import numpy as np

def rotx(theta_deg):
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    
    rotx = np.array([[1,0,0],[0,c,-s],[0,s,c]])
    
    return rotx

def roty(theta_deg):
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    
    roty = np.array([[c,0,s],[0,1,0],[-s,0,c]])
    
    return roty

def rotz(theta_deg):
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    
    rotz = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    
    return rotz