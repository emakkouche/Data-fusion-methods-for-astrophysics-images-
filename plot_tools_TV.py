#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 15:18:37 2022

@author: e.akkouche
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from CONSTANTS import *

def plot_band_figures_TV(X,Yh,VZ,limit,line,path_file):
    
    #------- Affichage de l'image full ------------
    fig,ax =plt.subplots()
    imgplot = plt.imshow(X)
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true_full_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_full_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_full_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(VZ);
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover_full_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_full_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_full_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    imgplot = plt.imshow(Yh);
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh_full.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh_full.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh_full',bbox_inches='tight')
    plt.show()
    
    #------- Affichage spectre de l'image ------------
    # fig,ax =plt.subplots()
    # imgplot = plt.plot(X[:,40,100])
    # ax.set_title('Référence',fontweight='bold')
    # ax.set_xlabel('$\lambda$',fontweight='bold')
    # ax.set_ylabel('Intensité',fontweight='bold')
    # plt.savefig(path_file+'VZ_true_spectre_TV.eps', format='eps',bbox_inches='tight')
    # plt.savefig(path_file+'VZ_true_spectre_TV.pdf', format='pdf',bbox_inches='tight')
    # plt.savefig(path_file+'VZ_true_spectre_TV',bbox_inches='tight')
    # plt.show()
    
    # fig,ax =plt.subplots()
    # imgplot = plt.plot(VZ[:,40,100])
    # ax.set_title('Sobolev',fontweight='bold')
    # ax.set_xlabel('$\lambda$',fontweight='bold')
    # ax.set_ylabel('Intensité',fontweight='bold')
    # plt.savefig(path_file+'VZ_recov_spectre_TV.eps', format='eps',bbox_inches='tight')
    # plt.savefig(path_file+'VZ_recov_spectre_TV.pdf', format='pdf',bbox_inches='tight')
    # plt.savefig(path_file+'VZ_recov_spectre_TV',bbox_inches='tight')
    # plt.show()
    
    # fig,ax =plt.subplots()
    # imgplot = plt.plot(Yh[:,40,100])
    # ax.set_title('Observée',fontweight='bold')
    # ax.set_xlabel('$\lambda$',fontweight='bold')
    # ax.set_ylabel('Intensité',fontweight='bold')
    # plt.savefig(path_file+'Yh_spectre.eps', format='eps',bbox_inches='tight')
    # plt.savefig(path_file+'Yh_spectre.pdf', format='pdf',bbox_inches='tight')
    # plt.savefig(path_file+'Yh_spectre',bbox_inches='tight')
    # plt.show()
    
    
    #------- Coupe horizontale de l'image ------------
    fig,ax = plt.subplots()
    imgplot = plt.plot(X[line,:])
    ax.set_title('Référence',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_true_horiz_section_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_horiz_section_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_horiz_section_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(VZ[line,:])
    ax.set_title('Sobolev',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'VZ_recover_horiz_section_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_horiz_section_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_horiz_section_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax = plt.subplots()
    imgplot =plt.plot(Yh[line,:])
    ax.set_title('Observée',fontweight='bold')
    ax.set_ylabel('Intensité',fontweight='bold')
    plt.savefig(path_file+'Yh_horiz_section.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh_horiz_section.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh_horiz_section',bbox_inches='tight')
    plt.show()
    
    
    plt.plot(X[line,:],'r',label = 'Référence')
    plt.plot(VZ[line,:],'b',label = 'Sobolev')
    plt.ylabel('Intensité',fontweight='bold')
    plt.legend()
    plt.savefig(path_file+'superposition_horiz_coupe.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'superposition_horiz_coupe.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'superposition_horiz_coupe',bbox_inches='tight')
    plt.show()
    
    #------- Affichage de l'image zoom ------------  
    X = X[limit[0]:limit[1],limit[2]:limit[3]]
    VZ = VZ[limit[0]:limit[1],limit[2]:limit[3]]
    Yh = Yh[limit[0]:limit[1],limit[2]:limit[3]]
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(X[:,:])
    ax.set_title('Référence',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_true_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_true_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(VZ:,:])
    ax.set_title('Sobolev',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'VZ_recover_TV.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_TV.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'VZ_recover_TV',bbox_inches='tight')
    plt.show()
    
    fig,ax =plt.subplots()
    img_plot = plt.imshow(Yh[:,:])
    ax.set_title('Observée',fontweight='bold')
    plt.axis('off')
    plt.savefig(path_file+'Yh.eps', format='eps',bbox_inches='tight')
    plt.savefig(path_file+'Yh.pdf', format='pdf',bbox_inches='tight')
    plt.savefig(path_file+'Yh',bbox_inches='tight')
    plt.show()
    

plot_band_figures_TV(X, Yh, VZ, limit, line, path_file)
mus = np.linspace(0,3,10)
fname = 'result_sobolev_30_08'+'SNR_TV'
SNR = np.load(fname,err)


fig,ax =plt.subplots()
img_plot = plt.plot(mus,SNR);plt.xscale('log')
ax.set_ylabel('SNR(dB)',fontweight='bold')
ax.set_xlabel('µ',fontweight='bold')
plt.savefig(SAVE_IMG+'SNR_sobolev.eps', format='eps')
plt.savefig(SAVE_IMG+'SNR_sobolev.pdf', format='pdf') 
plt.savefig(SAVE_IMG+'SNR_sobolev') 



