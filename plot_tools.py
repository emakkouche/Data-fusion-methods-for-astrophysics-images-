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

def plot_band_figures(X,Yh,VZ,critJ,limit,band,line,path_file):
    
    #------- Affichage de l'image full ------------
    fig,ax =plt.subplots()
    imgplot = plt.imshow(X[band,:,:])
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true_full.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_full.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_full',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(VZ[band,:,:]);
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover_full.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_full.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_full',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(Yh[band,:,:]);
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh_full.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh_full.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh_full',bbox_inches='tight')
    plt.show()
    
    #------- Affichage spectre de l'image ------------
    fig,ax =plt.subplots()
    imgplot = plt.plot(X[:,40,100])
    ax.set_title('Référence',fontweight='bold')
    ax.set_xlabel('$\lambda$',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_true_spectre.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_spectre.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_spectre',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.plot(VZ[:,40,100])
    ax.set_title('Sobolev',fontweight='bold')
    ax.set_xlabel('$\lambda$',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_recov_spectre.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recov_spectre.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recov_spectre',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.plot(Yh[:,40,100])
    ax.set_title('Observée',fontweight='bold')
    ax.set_xlabel('$\lambda$',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'Yh_spectre.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh_spectre.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh_spectre',bbox_inches='tight')
    plt.show()
    
    
    #------- Coupe horizontale de l'image ------------
    fig,ax = plt.subplots()
    imgplot = plt.plot(X[band,line,:])
    ax.set_title('Référence',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_true_horiz_section.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_horiz_section.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_horiz_section',bbox_inches='tight')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(VZ[band,line,:])
    ax.set_title('Sobolev',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_recover_horiz_section.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_horiz_section.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_horiz_section',bbox_inches='tight')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(Yh[band,line,:])
    ax.set_title('Observée',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'Yh_horiz_section.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh_horiz_section.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh_horiz_section',bbox_inches='tight')
    plt.show()
    
    
    plt.plot(X[band,line,:],'r',label = 'Référence')
    plt.plot(VZ[band,line,:],'b',label = 'Sobolev')
    plt.ylabel('Intensité',fontweight='bold')
    plt.legend()
    plt.savefig(path_file+'superposition_horiz_coupe.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'superposition_horiz_coupe.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'superposition_horiz_coupe',bbox_inches='tight')
    plt.show()
    
    #------- Affichage de l'image zoom ------------  
    X = X[:,limit[0]:limit[1],limit[2]:limit[3]]
    VZ = VZ[:,limit[0]:limit[1],limit[2]:limit[3]]
    Yh = Yh[:,limit[0]:limit[1],limit[2]:limit[3]]
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(X[band,:,:])
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(VZ[band,:,:])
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(Yh[band,:,:])
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh',bbox_inches='tight')
    plt.show()
    
    
    fig,ax = plt.subplots()
    imgplot = plt.plot(critJ)
    ax.set_xlabel('Itération',fontweight='bold')
    ax.set_ylabel('Energie',fontweight='bold')
    plt.savefig(path_file+'CritJ.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'CritJ.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'CritJ',bbox_inches='tight')
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
    
    plot_band_figures(VZ_true, Yh, VZ_recov,critJ,limit, band, line,path_file)
    
    return Yh,V_true,Z_true,V,Z,mean
   
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

def search_max(v, z):
    max_ = 0
    lh, lint = v.shape
    px, py = z.shape[1:]
    z = np.reshape(z, (lint, px*py))
    for l in range(lh):
        band = np.dot(v[l], z)
        max_band = np.max(band)
        if max_band > max_:
            max_ = max_band
    return max_

def compute_error(v_true, v, z_true, z, mean):
    
    lh, lint = v.shape
    px, py = z.shape[1:]
    # Compute sam
    print('Compute SAM ...')
    sam = np.zeros((px, py))
    for i in range(px):
        for j in range(py):
            spec_true = np.dot(v_true, z_true[:, i, j])
            spec = np.dot(v, z[:, i, j]) + mean
            sam[i, j] = compute_sam(spec_true, spec)
    print('Mean SAM = '+str(np.mean(sam)))
    print('Max SAM = '+str(np.max(sam)))
    
def compute_sam(s, s_hat):
    return np.arccos(np.dot(s, s_hat.T)*(np.linalg.norm(s)*np.linalg.norm(s_hat))**-1)

def calc_sigma(Lm, snr=50):
    return np.linalg.norm(Lm)**2*10**(-0.1*snr)*(1/np.prod(Lm.shape))

#Limite des lignes et colonnes de l'image à extraire
limit = [25,90,80,220]
mus = 10**np.linspace(3,5,10)
mu = 7
band = 100
line = 40
position = (51,21)
trunc = 100
SNR = []
PSNR = []

Yh,V_true,Z_true,V,Z,mean = get_result(limit, mu, band, line,SAVE_IMG)

norm_true = compute_norm_true(V_true, Z_true)

max_true = search_max(V_true, Z_true)

SAM = compute_error(V_true, V, Z_true, Z, mean)

S = mycheck_pca.my_choose_subspace(HS_IM)

sigma = calc_sigma(Yh,snr=60)

lh = V.shape[0]
px,py = Z.shape[1:]

for k in range(len(mus)):
    
    Zmu = fits.getdata(SAVE2+'Zoptim_mu_'+str(k)+'.fits')
    norm_diff,band = compute_norm(V_true, V, Z_true, Zmu, mean)
    PSNR.append(10*np.log10(max_true**2*(norm_diff)**-1*lh*px*py))
    SNR.append(10*np.log10(norm_true*(norm_diff)**-1))
    
fig,ax =plt.subplots()
img_plot = plt.plot(mus,PSNR);plt.xscale('log')
ax.set_ylabel('PSNR(dB)',fontweight='bold')
ax.set_xlabel('µ',fontweight='bold')
plt.savefig(SAVE_IMG+'PSNR_sobolev.eps', format='eps',bbox_inches='tight')
plt.savefig(SAVE_IMG+'PSNR_sobolev.pdf', format='pdf',bbox_inches='tight')
plt.savefig(SAVE_IMG+'PSNR_sobolev',bbox_inches='tight')


fig,ax =plt.subplots()
img_plot = plt.plot(mus,SNR);plt.xscale('log')
ax.set_ylabel('SNR(dB)',fontweight='bold')
ax.set_xlabel('µ',fontweight='bold')
plt.savefig(SAVE_IMG+'SNR_sobolev.eps', format='eps',bbox_inches='tight')
plt.savefig(SAVE_IMG+'SNR_sobolev.pdf', format='pdf',bbox_inches='tight')
plt.savefig(SAVE_IMG+'SNR_sobolev',bbox_inches='tight')

       
fname = SAVE2+'PSNR_sobolev'
np.save(fname,PSNR)

fname = SAVE2+'SNR_sobolev'
np.save(fname,SNR)

fig,ax =plt.subplots()
img_plot = plt.plot(np.sqrt(S[:]),'*')
ax.set_xlabel('Indice',fontweight='bold')
ax.set_ylabel('Valeur propre',fontweight='bold')
plt.savefig(SAVE_IMG+'ACP_tronc.eps', format='eps',bbox_inches='tight')
plt.savefig(SAVE_IMG+'ACP_tronc.pdf', format='pdf',bbox_inches='tight')
plt.savefig(SAVE_IMG+'ACP_tronc',bbox_inches='tight')

fname = SAVE2+'Singular_value'
np.save(fname,S)


#coordoonées des points où le spectre est correcte:
#plt.plot(VZ[:,40,100])
#plt.plot(VZtrue[:,40,100])
    
    