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

def plot_band_figures(X,Yh,VZ,limit,band,line,path_file):
    
    #------- Affichage de l'image full ------------
    plt.imshow(X[band,:,:]);
    plt.savefig(path_file+'VZ_true_full.eps', format='eps')
    plt.show()
    
    plt.imshow(VZ[band,:,:]);
    plt.savefig(path_file+'VZ_recover_full.eps', format='eps')
    plt.show()
    
    plt.imshow(Yh[band,:,:]);
    plt.savefig(path_file+'Yh_full.eps', format='eps')
    plt.show()
    
    #------- Coupe horizontale de l'image ------------
    
    plt.plot(X[band,line,:])
    plt.savefig(path_file+'VZ_true_horiz_section.eps', format='eps')
    plt.show()
    
    plt.plot(VZ[band,line,:])
    plt.savefig(path_file+'VZ_recover_horiz_section.eps', format='eps')
    plt.show()
    
    plt.plot(Yh[band,line,:])
    plt.savefig(path_file+'Yh_horiz_section.eps', format='eps')
    plt.show()
    
    #------- Affichage de l'image zoom ------------  
    X = X[:,limit[0]:limit[1],limit[2]:limit[3]]
    VZ = VZ[:,limit[0]:limit[1],limit[2]:limit[3]]
    Yh = Yh[:,limit[0]:limit[1],limit[2]:limit[3]]
    
    plt.imshow(X[band,:,:])
    plt.savefig(path_file+'VZ_true.eps', format='eps')
    plt.show()
    
    plt.imshow(VZ[band,:,:])
    plt.savefig(path_file+'VZ_recover.eps', format='eps')
    plt.show()
    
    plt.imshow(Yh[band,:,:])
    plt.savefig(path_file+'Yh.eps', format='eps')
    plt.show()
    

# def plot_TV_figures(Band):
def prod_VZ(V,Z):
    
    VZ = np.dot(V,np.reshape(Z,(Z.shape[0],Z.shape[1]*Z.shape[2])))
    VZ = np.reshape(VZ,(VZ.shape[0],Z.shape[1],Z.shape[2]))
    
    return VZ  
    
def get_VZ_J(V,mu):
        
    Z = fits.getdata(SAVE2+'Zoptim_mu_'+str(mu)+'.fits')
    J = np.load(SAVE2+'Jmu_'+str(mu)+'.npy')
    
    VZ = np.dot(V,np.reshape(Z,(Z.shape[0],Z.shape[1]*Z.shape[2])))
    VZ = np.reshape(VZ,(VZ.shape[0],Z.shape[1],Z.shape[2]))
    
    return VZ,Z,J

def get_VZ_true(file_data):
    
    M = fits.getdata(file_data+'M_1.fits').T
    A = fits.getdata(file_data+'A.fits')
    
    A = A[:,:,START:START+NBCOL_X]
    
    lh, m = M.shape
    n, p, q = A.shape
    
    VZ_true = np.reshape(np.dot(M, np.reshape(A, (n, p*q))), (lh,p, q))
    
    return VZ_true,M,A

def get_result(limit,mu,band,line,path_file):
    
    Yh = fits.getdata(HS_IM)
    V = fits.getdata(V_acp)
    mean = fits.getdata(DATA+'mean.fits')
    
    VZ_true,V_true,Z_true = get_VZ_true(DATA) 
    VZ_recov,Z,critJ = get_VZ_J(V,mu)
   
    # VZ_true = VZ_true[:,limit[0]:limit[1],limit[2]:limit[3]]
    # VZ_recov = VZ_recov[:,limit[0]:limit[1],limit[2]:limit[3]]
    # Yh = Yh[:,limit[0]:limit[1],limit[2]:limit[3]]
    
    plot_band_figures(VZ_true, Yh, VZ_recov,limit, band, line,path_file)
    
    return V_true,Z_true,V,Z,mean
   
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
    return norm,band

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
mus = 10**np.linspace(-1, 3, 10)
mu = 9
band = 100
line = 40
position = (51,21)
SNR = []

V_true,Z_true,V,Z,mean = get_result(limit, mu, band, line,SAVE_IMG)

norm_true = compute_norm_true(V_true, Z_true)

for k in range(len(mus)):
    
    Zmu = fits.getdata(SAVE2+'Zoptim_mu_'+str(k)+'.fits')
    norm_diff,band = compute_norm(V_true, V, Z_true, Zmu, mean)
    SNR.append(10*np.log10(norm_true*(norm_diff)**-1))
    
    
plt.plot(mus,SNR);plt.xscale('log')
plt.savefig(SAVE_IMG+'SNR_sobolev.eps', format='eps')    
fname = SAVE2+'SNR_sobolev'
np.save(fname,SNR)
    

#coordoonées des points où le spectre est correcte:
#plt.plot(VZ[:,40,100])
#plt.plot(VZtrue[:,40,100])
    
    