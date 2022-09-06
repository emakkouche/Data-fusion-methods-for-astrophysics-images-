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
    fig,ax =plt.subplots()
    imgplot = plt.imshow(X[band,:,:])
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true_full.eps', format='eps')
    plt.savefig(path_file+'VZ_true_full.pdf', format='pdf')
    plt.savefig(path_file+'VZ_true_full')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(VZ[band,:,:]);
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover_full.eps', format='eps')
    plt.savefig(path_file+'VZ_recover_full.pdf', format='pdf')
    plt.savefig(path_file+'VZ_recover_full')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(Yh[band,:,:]);
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh_full.eps', format='eps')
    plt.savefig(path_file+'Yh_full.pdf', format='pdf')
    plt.savefig(path_file+'Yh_full')
    plt.show()
    
    #------- Affichage spectre de l'image ------------
    fig,ax =plt.subplots()
    imgplot = plt.plot(X[:,40,100])
    ax.set_title('Référence',fontweight='bold')
    ax.set_xlabel('$\lambda$',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_true_spectre.eps', format='eps')
    plt.savefig(path_file+'VZ_true_spectre.pdf', format='pdf')
    plt.savefig(path_file+'VZ_true_spectre')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.plot(VZ[:,40,100])
    ax.set_title('Sobolev',fontweight='bold')
    ax.set_xlabel('$\lambda$',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_recov_spectre.eps', format='eps')
    plt.savefig(path_file+'VZ_recov_spectre.pdf', format='pdf')
    plt.savefig(path_file+'VZ_recov_spectre')
    plt.show()
    
    #------- Coupe horizontale de l'image ------------
    fig,ax = plt.subplots()
    imgplot = plt.plot(X[band,line,:])
    ax.set_title('Référence',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_true_horiz_section.eps', format='eps')
    plt.savefig(path_file+'VZ_true_horiz_section.pdf', format='pdf')
    plt.savefig(path_file+'VZ_true_horiz_section')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(VZ[band,line,:])
    ax.set_title('Sobolev',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_recover_horiz_section.eps', format='eps')
    plt.savefig(path_file+'VZ_recover_horiz_section.pdf', format='pdf')
    plt.savefig(path_file+'VZ_recover_horiz_section')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(Yh[band,line,:])
    ax.set_title('Observée',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'Yh_horiz_section.eps', format='eps')
    plt.savefig(path_file+'Yh_horiz_section.pdf', format='pdf')
    plt.savefig(path_file+'Yh_horiz_section')
    plt.show()
    
    #------------ subplot ------------------
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 2, 1)
    # imgplot = plt.imshow(X[band,:,:]);
    # ax.set_title('Référence')
    # ax = fig.add_subplot(1, 2, 2)
    # imgplot = plt.imshow(VZ[band,:,:]);
    # ax.set_title('Reconstruite')
    
    plt.plot(X[band,line,:],'r',label = 'Référence')
    plt.plot(VZ[band,line,:],'b',label = 'Sobolev')
    plt.ylabel('Intensité',fontweight='bold')
    plt.legend()
    plt.savefig(path_file+'superposition_horiz_coupe.eps', format='eps')
    plt.savefig(path_file+'superposition_horiz_coupe.pdf', format='pdf')
    plt.savefig(path_file+'superposition_horiz_coupe')
    plt.show()
    
    #------- Affichage de l'image zoom ------------  
    X = X[:,limit[0]:limit[1],limit[2]:limit[3]]
    VZ = VZ[:,limit[0]:limit[1],limit[2]:limit[3]]
    Yh = Yh[:,limit[0]:limit[1],limit[2]:limit[3]]
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(X[band,:,:])
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true.eps', format='eps')
    plt.savefig(path_file+'VZ_true.pdf', format='pdf')
    plt.savefig(path_file+'VZ_true')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(VZ[band,:,:])
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover.eps', format='eps')
    plt.savefig(path_file+'VZ_recover.pdf', format='pdf')
    plt.savefig(path_file+'VZ_recover')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(Yh[band,:,:])
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh.eps', format='eps')
    plt.savefig(path_file+'Yh.pdf', format='pdf')
    plt.savefig(path_file+'Yh')
    plt.show()
    

# def plot_TV_figures(Band):
def prod_VZ(V,Z):
    
    VZ = np.dot(V,np.reshape(Z,(Z.shape[0],Z.shape[1]*Z.shape[2])))
    VZ = np.reshape(VZ,(VZ.shape[0],Z.shape[1],Z.shape[2]))
    
    return VZ 

def add_mean_VZ(VZ,mean):
    
    for k in range(VZ.shape[0]):
        
        VZ[k,:,:] += mean[k]
    
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
    
    VZ_recov = add_mean_VZ(VZ_recov, mean)
   
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
mus = 10**np.linspace(3,5,10)
mu = 3
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

fig,ax =plt.subplots()
img_plot = plt.plot(mus,SNR);plt.xscale('log')
ax.set_ylabel('SNR(dB)',fontweight='bold')
ax.set_xlabel('µ',fontweight='bold')
plt.savefig(SAVE_IMG+'SNR_sobolev.eps', format='eps')
plt.savefig(SAVE_IMG+'SNR_sobolev.pdf', format='pdf') 
plt.savefig(SAVE_IMG+'SNR_sobolev')  

   
fname = SAVE2+'SNR_sobolev'
np.save(fname,SNR)
    

#coordoonées des points où le spectre est correcte:
#plt.plot(VZ[:,40,100])
#plt.plot(VZtrue[:,40,100])
    
    