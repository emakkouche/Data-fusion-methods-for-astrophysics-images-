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
    
def get_VZ_J(V,mu):
        
    Z = fits.getdata(SAVE2+'Zoptim_mu_'+str(mu)+'.fits')
    J = np.load(SAVE2+'Jmu_'+str(mu)+'.npy')
    
    VZ = np.dot(V,np.reshape(Z,(Z.shape[0],Z.shape[1]*Z.shape[2])))
    VZ = np.reshape(VZ,(VZ.shape[0],Z.shape[1],Z.shape[2]))
    
    return VZ,J

def get_VZ_true(file_data):
    
    M = fits.getdata(file_data+'M_1.fits').T
    A = fits.getdata(file_data+'A.fits')
    
    A = A[:,:,START:START+NBCOL_X]
    
    lh, m = M.shape
    n, p, q = A.shape
    
    VZ_true = np.reshape(np.dot(M, np.reshape(A, (n, p*q))), (lh,p, q))
    
    return VZ_true

def get_result(limit,mu,band,line,flag):
    
    Yh = fits.getdata(HS_IM)
    V = fits.getdata(V_acp)
    
    VZ_true = get_VZ_true(DATA) 
    VZ_recov,critJ = get_VZ_J(V,mu)
    
    VZ_true = VZ_true[:,limit[0]:limit[1],limit[2]:limit[3]]
    VZ_recov = VZ_recov[:,limit[0]:limit[1],limit[2]:limit[3]]
    Yh = Yh[:,limit[0]:limit[1],limit[2]:limit[3]]
    
    plot_band_figures(VZ_true, Yh, VZ_recov, band, line)
   
def compute_norm(v_true, v, z_true, z, mean):
    lh, lint = v.shape
    lh, lint_ = v_true.shape
    px, py = z.shape[1:]
    z = np.reshape(z, (lint, px*py))
    z_true = np.reshape(z_true, (lint_, px*py))
    norm = 0
    for l in range(lh):
        band = np.dot(v[l], z) + mean[l]
        band_true = np.dot(v_true[l], z_true)
        norm += np.sum((band-band_true)**2)
    return norm

def compute_norm_true(v_true, z_true):
    lh, lint = v_true.shape
    px, py = z_true.shape[1:]
    z_true = np.reshape(z_true, (lint, px*py))
    norm = 0
    for l in range(lh):
        band_true = np.dot(v_true[l], z_true)
        norm += np.sum(band_true**2)
    return norm
  
#Limite des lignes et colonnes de l'image à extraire
limit = [25,90,80,220]
mu = 9
band = 100
line = 40
flag = True

get_result(limit, mu, band, line, flag)

    norm_true = compute_norm_true(v_true, z_true)
    norm_diff = compute_norm(v_true, v, z_true, z, mean)
    snr = 10*np.log10(norm_true*(norm_diff)**-1)
    print('Global SNR = '+str(snr))



    
    
    