#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:56:56 2022

@author: e.akkouche
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from CONSTANTS import *

lin_l,lin_up,col_l,col_up = [25,90,80,220]

def plot_band_figures(X,Yh,VZ,band,line):
    
    plt.imshow(X[band,:,:])
    plt.show()
    plt.imshow(VZ[band,:,:])
    plt.show()
    plt.imshow(Yh[band,:,:])
    plt.show()
    plt.plot(X[band,line,:])
    plt.show()
    plt.plot(VZ[band,line,:])
    plt.show()
    plt.plot(Yh[band,line,:])
    plt.show()
    
def get_VZ_J(V,Z,mu):
    
    Z = fits.getdata(SAVE2+'Zoptim_mu_'+str(mu)+'.fits')
    J = np.load(SAVE2+'Jmu_'+str(mu)+'.npy')
    
    VZ = np.dot(V,np.reshape(Z,(Z.shape[0],Z.shape[1]*Z.shape[2])))
    VZ = np.reshape(VZ,(VZ.shape[0],Z.shape[1],Z.shape[2]))
    
    return VZ,J