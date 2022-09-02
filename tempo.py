#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:12:22 2022

@author: e.akkouche
"""
from time import time
import tools
import numpy as np
import matplotlib.pyplot as plt

# from Objfun import compute_symmpad_per_band
from astropy.io import fits
from CONSTANTS import *
from simulate_HS_MS import reshape_psf


import warnings
warnings.filterwarnings('ignore')

def snr(x, y):
    """
    snr - signal to noise ratio

       v = snr(x,y);

     v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )

       x is the original clean signal (reference).
       y is the denoised signal.

    Copyright (c) 2014 Gabriel Peyre
    """

    return 20 * np.log10(np.linalg.norm(x) / np.linalg.norm(x - y))

def get_spa_bandpsf_hs(band, sigma=0):
    # Get a spatial Point Spread Function at wavelength number 'band'. It has to be calculated in the Fourier domain and saved in advance with webbpsf.
    g_ = fits.getdata(PSF+'H_fft.fits')[:, band]
    k, m, n = g_.shape
    g = g_[0]+g_[1]*1.j
    return np.reshape(g, (m, n))

"""--------Post GD utilities--------"""
def postprocess(z, lacp):
    
    z = np.reshape(z, (lacp, NR, NC))
    z = np.fft.ifft2(z, norm='ortho')
    z = np.real(tools._centered(z[:, :-2, :-2], (lacp, NR-2*fact_pad-2, NC-2*fact_pad-2)))
    return z


def save_Zoptim(z, filename):
    
    hdu = fits.PrimaryHDU(z)
    hdu.writeto(filename, overwrite=True)
    
def recover_X(V,Z):
    
    Lband,Lin,Col = Z.shape
    Z = np.reshape(Z,(Lband,Lin*Col))
    VZ = np.dot(V,Z)
    
    return np.reshape(VZ,(V.shape[0],Lin,Col))
    
def get_Xtrue(file_data,Band):
    
    M = fits.getdata(file_data+'M_1.fits').T
    A = fits.getdata(file_data+'A.fits')
    
    A = A[:,:,START:START+NBCOL_X]
    
    n, p, q = A.shape
    return np.reshape(np.dot(M[Band], np.reshape(A, (n, p*q))), (p, q))


"""--------------------------------"""
def GradI(I):
    
    Gv = I[:,1:,:] - I[:,:-1,:]
    Gh = I[:,:,1:] - I[:,:,:-1]
    
    return(Gv,Gh)
    

"""--------------------------------"""

def preprocess_D(lin,col):
  
    Dx = np.zeros((lin,col))
    Dy = np.zeros((lin,col))

    Dx[0, 0] = 1
    Dx[0, 1] = -1
    Dy[0, 0] = 1
    Dy[1, 0] = -1

    Dx = np.fft.fft2(Dx)
    Dy = np.fft.fft2(Dy)

    return (Dx, Dy)

"""-------Calculer la FFT de Y par bande----- """ 
    
def process_Yh(Y):
    
    print('\n****** Processing Yh ******\n')

    Lh,Lin,Col = Y.shape
    
    M,N = compute_symmpad_per_band(Y[0],fact_pad).shape
    
    Yfft= np.zeros((Lh,M,N),dtype=complex)
    
    for k in range(Lh):
        
        Ypad = compute_symmpad_per_band(Y[k],fact_pad) 
        Yfft[k] = np.fft.fft2(Ypad,norm='ortho')
    
    print('\n****** Processing Yh Done ******\n')
    
    return Yfft

"""-------Calculer la FFT de Z par bande----- """ 
def process_Z(Z):
    
    print('\n****** Compute fft Z ******\n')
    
    Z = np.fft.fft2(tools.compute_symmpad_3d(Z,fact_pad),axes=(1, 2),norm='ortho')
    
    print('\n****** Compute fft Z Done ******\n')
    
    return Z

"""-------Calculer VZ(f)----- """
""" Il faut que Col(X)>Col(PSF)"""    
def process_VZ(V,Z):
    
    print('\n****** Processing VZ ******\n')
    
    #Lband,Lin,Col = Z.shape
    
    #Z = np.reshape(Z,(Lband,Lin*Col))
    
    VZ = np.dot(V,Z)
       
    print('\n****** Processing VZ Done ******\n')
    
    # return np.reshape(VZ,(V.shape[0],Z.shape[1],Z.shape[2]))
    return VZ

"""-------Reshaper PSF pour toutes les  bandes----- """
def process_PSF(VZ,sigma = 0):
    
    print('\n****** Processing PSF ******\n')
    
    Lband,Lin,Col = VZ.shape
    
    PSF = np.zeros((Lband,Lin,Col),dtype = complex)
    
    #Hreshaped = np.zeros((Lin,Col),dtype = complex)
    
    #H = get_spa_bandpsf_hs(1, sigma)
    
    #L,C = H.shape
    
    #Hamming = np.tile(np.hamming(C),(Lin,1))
    
    for k in range(Lband):
        
        H = get_spa_bandpsf_hs(k, sigma)
        
        Hcrop = reshape_psf(H, VZ[0])
        
        #Hshift = (np.fft.fftshift(H))
    
        #ind = np.unravel_index(np.argmax(Hshift, axis=None), Hshift.shape)
        
        #RLbnd = ind[0] - Lin//2
        #RUbnd = ind[0] + Lin//2
        # CLbnd = ind[1] - C//2
        # CUbnd = ind[1] + C//2
        #Nshift = (Col - C)//2
        
        #Hreshaped[:,Nshift:(Nshift+C)] = Hshift[RLbnd:RUbnd,:] * Hamming

        PSF[k] = Hcrop
        
    
    print('\n****** Processing PSF Done ******\n')
    
    return PSF
# def process_PSF(VZ,sigma = 0):
    
#     print('\n****** Processing PSF ******\n')
    
#     Lband,Lin,Col = VZ.shape
    
#     PSF = np.zeros((Lband,Lin,Col),dtype = complex)
    
#     Hreshaped = np.zeros((Lin,Col),dtype = complex)
    
#     H = get_spa_bandpsf_hs(1, sigma)
    
#     L,C = H.shape
    
#     Hamming = np.tile(np.hamming(C),(Lin,1))
    
#     for k in range(Lband):
        
#         H = get_spa_bandpsf_hs(k, sigma)
        
#         Hshift = (np.fft.fftshift(H))
    
#         ind = np.unravel_index(np.argmax(Hshift, axis=None), Hshift.shape)
        
#         RLbnd = ind[0] - Lin//2
#         RUbnd = ind[0] + Lin//2
#         # CLbnd = ind[1] - C//2
#         # CUbnd = ind[1] + C//2
#         Nshift = (Col - C)//2
        
#         Hreshaped[:,Nshift:(Nshift+C)] = Hshift[RLbnd:RUbnd,:] * Hamming

#         PSF[k] = np.fft.fftshift(Hreshaped)
        
#         Hreshaped *= 0
        
#         print('k = '+ str(k))
    
#     print('\n****** Processing PSF Done ******\n')
    
#     return PSF

def precompute_VtYH(V,H,Y):
    
    print('\n****** Precompute VtYH ******\n')
    t1 = time()
    #compute Vt(Y.*H)
    Lh,Mh,Nh = H.shape
    Ly,My,Ny = Y.shape
    
    res = np.dot(V.T,np.reshape(H,(Lh,Mh*Nh))*np.reshape(Y,(Ly,My*Ny)))
    
    t2 = time()
    print('VtYH Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')
    print('\n****** Precompute VtYH Done ******\n')
    
    return res

def precompute_VtDvH2(V,Lacp,H):
    
    print('\n****** Precompute VtDvH2 ******\n')
    #compute Vt*Diag(Vl)*H**2 for l =1,...,Lacp 
    #Vtranpose : LacpxLh
    #H = PSF**2 : Lhxpm
    #Vtranspose * Diag(Vl) * H**2 : Lacpxpm
    #return matrix : LacpxLacpxpm
    t1 = time()
    Lband,Lin,Col = H.shape
    
    H = np.reshape(H,(Lband,Lin*Col))
    
    H2 = H**2
    
    maxH2 = np.max(np.abs(H2))
    precomputed_term = np.zeros((Lacp,Lacp,Lin*Col))
    
    for k in range(Lacp):
        
        precomputed_term[k] = np.dot(V.T*V[:,k],H2)
    t2 = time()
    print('VtDvH2 Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')
    print('\n****** Precompute VtDvH2 Done ******\n')
        
    return precomputed_term,maxH2


def compute_VtVZH2(repZ,precomp_term):
    #compute Z .* VtDvH2
    #faire la somme des matrices sur l'axe 0 (sur les Lacp bandes)
    
    return np.sum(repZ*precomp_term,0)

"""------ Precompute terms-------"""
def preprocessing(Yh,V,Z,Lacp,mean):
    #calculer les termes en rouge et jaune
    print('\n****** Preprocessing ******\n')
    
    t1 = time()
    
    Yfft = np.fft.fft2(tools.compute_symmpad_3d(Yh,fact_pad),axes=(1, 2),norm='ortho')
    
    dim = Yfft.shape
    
    Zfft = np.fft.fft2(tools.compute_symmpad_3d(Z,fact_pad),axes=(1, 2),norm='ortho')
    
    N = Zfft.shape[1]*Zfft.shape[2]
    mean_ = np.zeros((mean.shape[0], N))
    mean_[:, 0] = mean*np.sqrt(N)
    
    mean_[:, 0] = mean_[:, 0]*tools.get_h_mean()
      
    Yfft = np.reshape(Yfft,(Yfft.shape[0],Yfft.shape[1]*Yfft.shape[2])) - mean_
    
    Yfft = np.reshape(Yfft,dim)
    
    #Y = process_Yh(Yh)
    
    # VZ = process_VZ(V,Z)
    
    # H = process_PSF(VZ)
    
    # Zfft = process_Z(Z)
    
    H = process_PSF(Yfft)
    
    fname = SAVE2+'M_fft_PSF_croped'
    np.save(fname,H)
    
    VtDvH2,maxH2 = precompute_VtDvH2(V,Lacp,H)    #precomputed_term
    
    VtYH = precompute_VtYH(V,H,Yfft)        #term2
    
    D = preprocess_D(NR,NC)
    
    t2 = time()
    
    print('\n****** Preprocessing Done ******\n')
    return Yfft,Zfft,H,VtDvH2,VtYH,maxH2,D,np.round((t2-t1))

"""-------Repliquer les lignes de Z----- """
def replicate_Z(Lacp,Zfft):
    Lband,Lin,Col = Zfft.shape
    
    Zfft = np.reshape(Zfft,(Lband,Lin*Col))
    #faire l'opération pour toutes les bandes a la fois
    #taille de l'input Z est de Lacpxpm
    #taille sortie LacpxLacpxpm
    return np.tile(np.asarray(Zfft),Lacp).reshape(-1,Lacp,Zfft.shape[-1]) #reshape en Nblin = NbrepeatxNbColZ le long de l'axe 0


"""
"-------Calculer la valeur du critère au point Z(f)---------
"        0.5*||Yh(f) - (VZ(f)).*Hpsf||^2
"        J(Z(f)) = sum(i,j) |[Yh(f) - (VZ(f)).*Hpsf](i,j)|^2
"------------------------------------------------------------
"""    
def CritJ(Y,V,Z,H,D,mu,sigma = 0):
    
    print('\n****** Eval CritJ(Z) ******\n')
    
    #dim(Z) = dim(Y) : Lbandx(LinxCol)
    
    # Lband,Lin,Col = Y.shape
    
    # Y = np.reshape(Y,(Lband,Lin*Col))
    
    # VZ = process_VZ(V, Z)
    
    # VZ = np.dot(V,Z)
    
    # VZH = VZ*H
    
    # VZH = np.reshape(VZH,(Y.shape[0],Y.shape[1]))
    
    Ztemp = Z.copy()
    
    Ztemp = np.reshape(Ztemp,(Lacp,NR,NC))
    
    residu = 0.5 * (np.linalg.norm(Y-np.dot(V,Z)*H)**2)
    
    regul = np.sum((D[0]*Ztemp)**2+(D[1]*Ztemp)**2)
    
    #J = 0.5 * (np.linalg.norm(Y-np.dot(V,Z)*H)**2) + mu * np.sum((D[0]*np.reshape(Z,(Lacp,NR,NC)))**2+(D[1]*np.reshape(Z,(Lacp,NR,NC)))**2)
    """mettre Dx Dy en facteur commun"""
    # J = 0.5 * (np.linalg.norm(Y-np.dot(V,Z)*H)**2) + mu * np.linalg.norm((D[0]+D[1])*np.reshape(Z,(Lacp,NR,NC)))

    # print('\n****** Eval CritJ(Z) Done ******\n')
    return (residu + mu * regul),residu,regul


"""-------Evaluer le gradient au point Z(f)---------"""
def GradJ(Zfft,Lacp,term2,precomp_term,D,mu):
    
    #compute Vt[VZ.*H2] - Vt(Y.*H)
    
    #repZ = replicate_Z(Lacp, Zfft)
    # Lband,Lin,Col = Zfft.shape
    # Zfft = np.reshape(Zfft,(Lband,Lin*Col))
    
    # décommenter ça repZ = np.tile(np.asarray(Zfft),Lacp).reshape(-1,Lacp,Zfft.shape[-1])
    
    #term1 = compute_VtVZH2(repZ, precomp_term)
    """Je peux eviter de faire repZ en faisant Z = Lsub x 1 x pm .* lsub x pm"""
    """nr nc prcq je sais que Zfft est en lin+padding col+padding"""
    term1 = np.sum(np.reshape(Zfft,(Lacp,1,NR*NC))*precomp_term,0)
    
    Zfft = np.reshape(Zfft,(Lacp,NR,NC))
    
    return (term1 - term2 + 2 * mu * np.reshape(np.conj(D[0]) * Zfft * D[0] + np.conj(D[1]) * Zfft * D[1],(Lacp,NR*NC)))

    """ Essayer A*A = |A|^2"""
    # return (term1 - term2 + 2 * mu * np.reshape(Zfft * np.abs(D[0])**2 + Zfft * np.abs(D[1])**2,(Lacp,NR*NC)))
    
    # return (term1 - term2 + 2 * mu * np.reshape(Zfft * (np.conj(D[0])*D[0]+np.conj(D[1])*D[1]),(Lacp,NR*NC)))

def GD(Y,V,Z,H,Lacp,term2,precomp_term,D,mu,maxH2):
    print('\n***** GRADIENT DESCENT ALGORITHM *****\n')
    
    #GradJ = Lacpxpm
    #Z(f) = LacpxLinxCol
    
    t1 = time()
    
    time_iter = []
    
    NB_ITER = 1000         #Nbr max d'itérations
    EPS_J = 1e-5        #Seuil de variation du critere
    STEP = 1e-6         #Pas d'incrémentation 
    
    k = 0
    
    Zk = Z.copy()
    
    J,residu,regul = CritJ(Y,V,Zk,H,D,mu,sigma = 0)
    
    J_Zk = [J]
    
    GradJ_Zk = GradJ(Zk,Lacp,term2,precomp_term,D,mu)
    
    Zk_Old = 100 * (Zk)
    
    J,residu,regul = CritJ(Y,V,Zk_Old,H,D,mu,sigma = 0)
    
    J_ZkOld = J
    
    print(str(k)+' -- J(Z) = '+str(J_Zk[-1]))
    
    while (k < NB_ITER) and np.abs((J_ZkOld - J_Zk[-1])/J_ZkOld)> EPS_J: #while k<NB_ITER and (norm(Zold-Z) > EPS_Z or norm(J_Zold-J_Z) > EPS_J)
        
        tstart = time()
        
        k += 1
        
        J_ZkOld = J_Zk[-1]

        Zk =  Zk + STEP * -GradJ_Zk
        
        J,residu,regul = CritJ(Y,V,Zk,H,D,mu,sigma = 0)
        
        J_Zk.append(J)
        
        GradJ_Zk = GradJ(Zk,Lacp,term2,precomp_term,D,mu)
        
        tend = time()
        
        time_iter.append(np.round((tend-tstart)))
        print(str(k)+' -- J(Z) = '+str(J_Zk[-1]))
        print('Iteration Computation time : '+str(np.round((tend-tstart)/60))+'min '+str(np.round((tend-tstart)%60))+'s.')
    

    t2 = time()
    print('GD Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')
    print('\n***** GRADIENT DESCENT Done *****\n')
    
    fname = SAVE2+'residu_mu_'+str((mu))
    np.save(fname,residu)

    fname = SAVE2+'regul_mu_'+str((mu))
    np.save(fname,regul)
        
    return Zk,J_Zk,STEP,time_iter,np.round((t2-t1))
    
    
