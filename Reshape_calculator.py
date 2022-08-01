#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:10:53 2022

@author: e.akkouche
"""

import numpy as np
# import tools
from time import time


from grad import grad
from div import div
from myproduce_HS_MS import add_noise_nocorr

import matplotlib.pyplot as plt
from astropy.io import fits
from tools import compute_symmpad_3d, _centered, compute_symmpad
from skimage.transform import resize
from scipy.signal import convolve2d
from skimage.restoration import denoise_tv_chambolle
# from pandeia.engine.instrument_factory import InstrumentFactory
from CONSTANTS import *

import warnings
warnings.filterwarnings('ignore')


def get_spa_bandpsf_hs(band, sigma=0):
    # Get a spatial Point Spread Function at wavelength number 'band'. It has to be calculated in the Fourier domain and saved in advance with webbpsf.
    g_ = fits.getdata(PSF+'M_fft.fits')[:, band]
    k, m, n = g_.shape
    g = g_[0]+g_[1]*1.j
    return np.reshape(g, (m, n))

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

# (2*fact_pad+2) % d = 0; (fact_pad//d + 1) % 2 = 0
# Number = 147456
# nr = 90
# nc = 1
# d = 3

# Ndim = 500
# Maxpad = 41
# fact_pad = 29

# Divider1 = []
# Divider2 = []

# vect = np.array([0],dtype = 'float')
# vect2 = np.array([0],dtype = 'float')


# for nc in range(Ndim):
    
#     for k in range(Maxpad):
        
#         if( (nc % d == 0) and ((2*k+2) % d ==0) and ((k//d + 1) % 2 == 0)):
            
#             Divider1.append(nc)
#             Divider2.append(k)
            
            
# for i in range(Ncol):
    
#     #if( ((2*i+2) % d == 0) and ((i//d + 1) % 2 ==0 )):
#     #    print(i)
#     #Diff = Ncol - (i+1)
    
#     if( Diff %(2) == 0):
#         Divider.append(i+1)
#     vect = np.append(vect,i+1)
    
    
    
# vect = Eigvalue
    
# vect = vect ** 2
    
# sum_2 = sum(vect)

#vect = vect / sum_2
    
#cumul_var = np.cumsum(vect)




# def tester(a,b):
    
#     L = 3
#     Y = np.zeros((L,a,b))
#     W = np.zeros((L,2*a,2*b))
    
#     for k in range (L):
#         Y[k] = np.random.rand(a,b)
#         W[k] = np.random.rand(2*a,2*b)
        
#     return W,Y


#Q,QQ = tester(3,2)

#L,M,N = QQ.shape

# def compute_symmpad(X, npix):
#     """
#     Compute a symmetric padding of npix and add 2 rows and 2 columns of zero padding (for aliasing)
#     """
#     M = X.shape[0]
#     N = X.shape[1]
#     symm_pad = np.zeros((M+2*npix, N+2*npix), dtype='float64')
#     symm_pad[npix:M+npix, npix:N+npix] = X

#     symm_pad[npix:M+npix, 0:npix] = X[:, npix-1::-1]
#     symm_pad[0:npix, npix:N+npix] = X[npix-1::-1]
#     symm_pad[M+npix:, npix:N+npix] = X[M-1:M-1-npix:-1]
#     symm_pad[npix:M+npix, N+npix:] = X[:, N-1:N-1-npix:-1]

#     symm_pad[0:npix, 0:npix] = symm_pad[2*npix-1:npix-1:-1, 0:npix]
#     symm_pad[M+npix:M+2*npix, 0:npix] = symm_pad[M+npix-1:M-1:-1, 0:npix]
#     symm_pad[0:npix, N+npix:N+2*npix] = symm_pad[2*npix-1:npix-1:-1, N+npix:N+2*npix]
#     symm_pad[M+npix:M+2*npix, N+npix:N+2*npix] = symm_pad[M+npix-1:M-1:-1, N+npix:N+2*npix]

#     symm_pad_ = np.zeros((M+2*npix+2, N+2*npix+2), dtype='float64')
#     symm_pad_[:M+2*npix, :N+2*npix] = symm_pad
#     return symm_pad_  

# def compute_symmpad_3d(X, npix):
#     L, M, N = X.shape
#     Y = np.zeros((L, M+2*npix+2, N+2*npix+2))
#     for i in range(L):
#         Y[i] = compute_symmpad(X[i], npix)
#     return Y
     
# matrice = np.random.randint(1, 10, size=(1,7), dtype=int)


# fact_pad = 3

# symmetric_pad = compute_symmpad_3d(matrice, fact_pad)   

# Lin = 150
# Col = 414
# sigma = 0
# lh = 4974


# Band = 4900

# H = get_spa_bandpsf_hs(Band, sigma=0)

# Hifft = np.fft.ifft2(H)

# plt.imshow(np.log(np.abs(Hifft))) #donne 4 coins verts et image blue (np.abs) et image blanche ave (np.real)
# plt.imshow(np.log(np.abs(np.fft.fftshift(Hifft))))  marche aussi 
# plt.imshow(np.fft.fftshift(np.log(np.real(Hifft)))) marche aussi mais fond blanc




# def preprocess_D(sh):
#     m, n = sh
#     Dx = np.zeros(sh)
#     Dy = np.zeros(sh)

#     Dx[0, 0] = 1
#     Dx[0, 1] = -1
#     Dy[0, 0] = 1
#     Dy[1, 0] = -1

#     Dx = np.fft.fft2(Dx)
#     Dy = np.fft.fft2(Dy)

#     return (Dx, Dy)


# sh = nr,nc

# D = preprocess_D(sh)


# flag = 1
# for i in range(lh):
#         # Compute the 3D scene at a wavelength 
#         # Get spatial PSF at that band
#         k = i *50
#         H = get_spa_bandpsf_hs(k, sigma)
#         M = get_spa_bandpsf_hs_M(k, sigma)
        
#         #Initialize arrays
#         if( flag == 1):
#             Hreshaped = np.zeros((Lin,Col),dtype = complex)
            
#             #lin,col = X.shape
#             #Y = np.zeros((lh,lin,col))
            
#             #Generate Hamming window Nb_coeff = Xcol repeat Xlin 
#             Hamming = np.tile(np.hamming(384),(150,1))
#             flag = 0
#         '''---------------------------------------'''
#         #Get the indexes of the maximal element of |H|
#         Hshift = (np.fft.fftshift(H))             #A Supprimer
#         Mshift = (np.fft.fftshift(M))  
#         ind = np.unravel_index(np.argmax(Hshift, axis=None), Hshift.shape)
        
#         #Set the lower and upper bounds of the area to be extracted 
#         RLbnd = ind[0] - NbpixR//2
#         RUbnd = ind[0] + NbpixR//2
#         CLbnd = ind[1] - NbpixC//2
#         CUbnd = ind[1] + NbpixC//2
        
#         #Reshape PSF
        
#         Hreshaped[:,Nshift:(Nshift+NbpixC)] = Hshift[RLbnd:RUbnd,:] #Hreshape ==> 174x984 Hshift 174x384
#         Hreshaped[:,Nshift:(Nshift+NbpixC)] *= Hamming
#         Hreshaped = np.fft.fftshift(Hreshaped)
        
#         plt.imshow(np.abs(Hshift))
#         plt.show()
#         plt.imshow(np.abs(Mshift))
#         #Convolve with PSF without subsample
#         plt.show()
#         Hreshaped *= 0
        
#         print('i = '+ str(i))  


# V = fits.getdata(V_acp)
# Z = fits.getdata(Z_acp)

# Lband,Lin,Col = Z.shape
# Z = np.reshape(Z,(Lband,Lin*Col))
# t1 = time()
# VZ = np.dot(V,Z)
# t2 = time()
# print('Cg Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')


# Band = V.shape[1]
# Stock = Hreshaped = np.zeros(VZ.shape,dtype = float)
# t1 = time()
# for k in range(Band):
    
#     Stock += np.dot(V[:,k],Z[k,:]) 
    
# t2 = time()

# print('Cg Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')


"""boucle for plus lente"""
# Z1 = np.reshape(Z,(Lacp,Z.shape[1]*Z.shape[2]))
# VZ = np.dot(V,Z1)

# X = np.zeros((VZ.shape[0],VZ.shape[1]),dtype=float)

# X = np.outer(V[:,0],Z1[0,:])

# for k in range(V.shape[1]-2):
    
#     X += np.outer(V[:,k+1],Z1[k+1,:])


# plt.plot(np.arange(1,51),Eigvalue[0:50])
# plt.yscale('log')
# plt.title('Valeurs propres')
# plt.xlabel('Indice de la composante principale')
# plt.savefig(SAVE+'eigenvalues_pca_urgent_50.png')



# plt.plot(np.arange(1,51),cumul_var[:50])
# plt.yscale('log')
# plt.title('Pourcentage inertie expliquée')
# plt.xlabel('Dimension du sous espace')
# # plt.savefig(SAVE+'Cumulative_variance_pca_urgent2.png')
# # plt.show()    

# denois_img = denoise_tv_chambolle(Yh, weight = 1,n_iter_max= 1000,multichannel= True)
    
# eps = 1e-6


# evalGradJ = GradJ(Zfft2+eps,Lacp,T2,T1,D,mu)
# print('Eval GradJ(Z+eps) = ',evalGradJ)

# evalGradJ2 = (CritJ(Yfft,V,Zfft2,H,D,mu,sigma = 0) - CritJ(Yfft,V,Zfft,H,D,mu,sigma = 0)) / eps

# print('Eval GradJ(Z+eps) = ',evalGradJ2)


#####################Version qui fonctionne####################

# Hreshaped = np.zeros((90,354),dtype = complex)
# H = get_spa_bandpsf_hs(4000, sigma=0)
# Hamming = np.tile(np.hamming(H.shape[1]),(Yh.shape[0],1))
# Hshift = (np.fft.fftshift(H))             #A Supprimer
# ind = np.unravel_index(np.argmax(Hshift, axis=None), Hshift.shape)


# RLbnd = ind[0] - 90//2
# RUbnd = ind[0] + 90//2
# CLbnd = ind[1] - 354//2
# CUbnd = ind[1] + 354//2
        
#         #Reshape PSF
        
# # Hreshaped[:,Nshift:(Nshift+NbpixC)] = Hshift[RLbnd:RUbnd,CLbnd:CUbnd] #Hreshape ==> 174x984 Hshift 174x384
# Hreshaped = Hshift[RLbnd:RUbnd,CLbnd:CUbnd]
# Hreshaped[:,Nshift:(Nshift+NbpixC)] *= Hamming
# Hreshaped = np.fft.fftshift(Hreshaped)


# Phi = lambda x,h: np.real(np.fft.ifft2(np.fft.fft2(x) * h))
# repeat3 = lambda x,k: resize( np.repeat( x, k, axis=1), [90, 354, k])
# epsilon = 0.4*1e-2
# Lambda = 0.2
# y = Yh[4000,:,:]
# tau = 1.9 / ( 1 + Lambda * 8 / epsilon)
# Xestim = y
# niter = 2000
# E = np.zeros((niter,1))
# for i in np.arange(0, niter):
#     # Compute the gradient of the smoothed TV functional.
#     Gr = grad(Xestim)
#     d = np.sqrt(epsilon**2 + np.sum(Gr**2, axis=2))
#     G = -div(Gr[:,:,0] / d,Gr[:,:,1] / d )
#     # step
#     e = Phi(Xestim,Hreshaped)-y
#     Xestim = Xestim - tau*( Phi(e,Hreshaped) + Lambda*G)
#     # energy
#     E[i] = 1/2*np.linalg.norm(e.flatten())**2 + Lambda*sum(d.flatten())
# # display energy
# clf;
# plt.plot(E)
# axis('tight')
# xlabel('Iteration #')
# ylabel('Energy')
############################################################

"""---------------------Sans padding bordures--------------------"""
"""REMARQUE!!! dim(H) > dim(Y)"""
# Band = 1000

# Yh = fits.getdata(HS_IM)

# M = fits.getdata(DATA+'M_1.fits').T
# A = fits.getdata(DATA+'A.fits')

# A = A[:,:,Start:Start+NbcolMA]

# n, p, q = A.shape

# Xinit = np.reshape(np.dot(M[Band], np.reshape(A, (n, p*q))), (p, q))

# Nband,Lin,Col = Yh.shape

# H = get_spa_bandpsf_hs(Band, sigma=0)

# Hifft = np.fft.ifft2(H,norm = 'ortho')
# Hifft = np.fft.ifft2(H)

        
# Hifft = np.fft.ifftshift(Hifft)
    
# ind = np.unravel_index(np.argmax(Hifft, axis=None), Hifft.shape)
        
# RLbnd = ind[0] - Lin//2
# RUbnd = ind[0] + Lin//2
# CLbnd = ind[1] - Col//2
# CUbnd = ind[1] + Col//2

# Nshift = (Col - H.shape[1])//2

# Hreshaped[:,Nshift:(Nshift+H.shape[1])] = Hifft[RLbnd:RUbnd,:] 

# Hreshape = np.fft.fftshift(Hifft[RLbnd:RUbnd,CLbnd:CUbnd])

# Hreshape /= np.sum(Hreshape)

# Yblur = np.real(np.fft.ifft2(np.fft.fft2(Hreshape)*np.fft.fft2(Xinit,norm = 'ortho'), norm='ortho'))

# Yblur +=  np.sqrt(0.1)*np.random.randn(Lin,Col) 

# Yh[Band,:,:] = Yblur

# Yblur = Yh[Band,:,:]

# PSF[k] = np.fft.fftshift(Hreshaped)
"""-----------------------------------------------------"""

# Phi = lambda x,h: np.real(np.fft.ifft2(np.fft.fft2(x,norm = 'ortho') * np.fft.fft2(h,norm = 'ortho'),norm = 'ortho'))

# repeat3 = lambda x,k: resize( np.repeat( x, k, axis=1), [90, 354, k])

# dotp = lambda x,y: np.sum(x.flatten()*y.flatten())
# a = np.random.randn(90,354)
# b = np.random.randn(90,354,2)
# dotp(grad(a),b) + dotp(a,div(b[:,:,0],b[:,:,1]))

"""J'ai test de mettre + tau * (Phi+lambda*G) marche pas"""
"""J'ai test division par 1 au lieu de d marche mais pas bien"""
"""J'ai test eps 0.1*1e-2 pas bien """

"""PSF M.ftis"""
"""eps 0.1*e-1 0.2*e-1  0.2*e-2 0.2*e-7  0.2*1e-4 0.2 ne fonctionne pas"""
"""si j'enleve les norm = ortho dans Phi, j'aurai des norm = inf et la boucle s'arrete a qqlq itérations"""


""" Ne fonctionne pas
    eps 0.4*1e-2
    lambda_list(np.log(np.linspace(1,3,20)))
    Xbest  = np.zeros
    Xestim = Yblur * randn
"""

""" Fonctionne mal fonction cout en ligne (peut être il faut augmenter le eps)
    eps 0.4*1e-2
    lambda_list(np.log(np.linspace(1,3,20)))
    Xbest  = np.zeros
    Xestim = Yblur * random
"""

"""eps = 0.2-1e-1
   lambda_list = np.log(np.linspace(1,3,20))
   niter = 4000
   Xbest = np.zeros
   Xestim = Yh[Band,:,:] * np.random.rand(Lin,Col) 
   
   RMQ : SNR qui se stabilise à 7
         Mauvaise reconstruction de l'img'
         Fonction cout diminue peu

"""

"""eps = 0.2-1e-1
   lambda_list = np.log(np.linspace(1,3,20))
   niter = 4000
   Xbest = np.zeros
   Xestim = np.random.rand(Lin,Col) 
   
   RMQ : Energie décroit mais pas suffisament
         SNR décroit et n'est que 1.
         L'image restitué Xbest est comme Yh mais avec du bruit
         
"""
# epsilon = 0.2*1e-1
# # lambda_list = np.linspace(0,2,20)
# lambda_list = np.log(np.linspace(1,3,20))

# tau = 1.9 / ( 1 + max(lambda_list) * 8 / epsilon)

# niter = 4000

# E = np.zeros((niter,1))

# err = np.zeros((len(lambda_list),1))

# Xbest = Yh[Band,:,:]
# Xbest = Yblur
# Xbest = np.zeros(Xinit.shape)

# for k in np.arange(0,len(lambda_list)):
    
    # Xestim = Yh[Band,:,:]
    # Xestim = Yh[Band,:,:] * np.random.rand(Lin,Col) 
    # Xestim = np.zeros(Xinit.shape)
    # Xestim = np.random.rand(Lin,Col) 
    
    # Xestim_prev = Xestim*100
    
    # i = 0
    
    # Lambda = lambda_list[k]
    
    # while i < niter and np.linalg.norm(Xestim - Xestim_prev)/np.linalg.norm(Xestim_prev) > 1e-5 :
    
    #     # Compute the gradient of the smoothed TV functional.
    #     Gr = grad(Xestim)
    #     d = np.sqrt(epsilon**2 + np.sum(Gr**2, axis=2))
    #     # d = np.sqrt(epsilon**2 + np.linalg.norm(Gr))
    #     G = -div(Gr[:,:,0] / d,Gr[:,:,1] / d )
        
    #     # step
    #     e = Phi(Xestim,Hreshape) - Yh[Band,:,:]
    
    #     Xestim_prev = Xestim
    
    #     Xestim = Xestim - tau*( Phi(e,Hreshape) + Lambda*G)
    #     # energy
    #     E[i] = 1/2*np.linalg.norm(e.flatten())**2 + Lambda*np.sum(d.flatten())
    
    #     i = i+1
    
    # err[k] = snr(Xinit,Xestim)
    
    # if err[k] > snr(Xinit,Xbest):
    #     Xbest = Xestim
        
    
# while i < niter and np.linalg.norm(Xestim - Xestim_prev)/np.linalg.norm(Xestim_prev) > 1e-10 :
    
#     # Compute the gradient of the smoothed TV functional.
#     Gr = grad(Xestim)
#     d = np.sqrt(epsilon**2 + np.sum(Gr**2, axis=2))
#     G = -div(Gr[:,:,0] / d,Gr[:,:,1] / d )
#     # step
#     e = Phi(Xestim,Hreshape) - Yh[124,:,:]
    
#     Xestim_prev = Xestim
    
#     Xestim = Xestim - tau*( Phi(e,Hreshape) + Lambda*G)
#     # energy
#     E[i] = 1/2*np.linalg.norm(e.flatten())**2 + Lambda*np.sum(d.flatten())
    
#     i = i+1
# display energy

# plt.plot(E)
# plt.axis('tight')
# plt.xlabel('Iteration #')
# plt.ylabel('Energie')

# plt.imshow(Xinit)
# plt.title('Référence')

# plt.imshow(Yblur)
# plt.title('Observée')

# plt.imshow(Yh[Band,:,:])
# plt.title('Observée')

# plt.imshow(Xestim)
# plt.title('Deconvolution - Variation totale')


# plt.plot(lambda_list,err)
# plt.axis('tight')
# plt.xlabel('\lambda (échelle log)') 
# plt.ylabel('SNR')

# plt.imshow(Xbest)
# plt.title('Variation totale SNR = '+str(round(np.max(err),2))+'dB')    


# plt.plot(Xinit[63,:], 'r',label = 'Référence') # plotting t, a separately 
# plt.plot(Xbest[63,:], 'b',label = 'Variation totale') # plotting t, b separately 
# plt.legend(loc="upper right")


# plt.plot(Yblur[63,:],'g')
# plt.show()
    
# print('Finish')

"""------------Version fonctionnelle V2 12 Juillet------"""
Band = 4000

Yh = fits.getdata(HS_IM)

M = fits.getdata(DATA+'M_1.fits').T
A = fits.getdata(DATA+'A.fits')

A = A[:,:,START:START+NBCOL_X]

n, p, q = A.shape

Xinit = np.reshape(np.dot(M[Band], np.reshape(A, (n, p*q))), (p, q))

Nband,Lin,Col = Yh.shape

H = get_spa_bandpsf_hs(Band, sigma=0)

Hifft = np.fft.ifft2(H)
        
Hifft = np.fft.ifftshift(Hifft)
    
ind = np.unravel_index(np.argmax(Hifft, axis=None), Hifft.shape)
        
RLbnd = ind[0] - Lin//2
RUbnd = ind[0] + Lin//2
CLbnd = ind[1] - Col//2
CUbnd = ind[1] + Col//2

Nshift = (Col - H.shape[1])//2

# Hreshaped[:,Nshift:(Nshift+H.shape[1])] = Hifft[RLbnd:RUbnd,:] 

Hreshape = np.fft.fftshift(Hifft[RLbnd:RUbnd,CLbnd:CUbnd])

Hreshape /= np.sum(Hreshape)

Yblur = np.real(np.fft.ifft2(np.fft.fft2(Hreshape)*np.fft.fft2(Xinit,norm = 'ortho'), norm='ortho'))

Yblur +=  np.sqrt(0.1)*np.random.randn(Lin,Col) 

Yh[Band,:,:] = Yblur

# Yblur = Yh[Band,:,:]

# PSF[k] = np.fft.fftshift(Hreshaped)

"""----------"""
Phi = lambda x,h: np.real(np.fft.ifft2(np.fft.fft2(x,norm = 'ortho') * np.fft.fft2(h,norm = 'ortho'),norm = 'ortho'))

repeat3 = lambda x,k: resize( np.repeat( x, k, axis=1), [90, 354, k])

"""-------"""
epsilon = 0.2*1e-1
# lambda_list = np.linspace(0,2,20)
lambda_list = np.log(np.linspace(1,3,20))

tau = 1.9 / ( 1 + max(lambda_list) * 8 / epsilon)

niter = 4000

E = np.zeros((niter,1))

err = np.zeros((len(lambda_list),1))

# Xbest = Yh[Band,:,:]
# Xbest = Yblur
Xbest = np.zeros(Xinit.shape)

for k in np.arange(0,len(lambda_list)):
    
    # Xestim = Yh[Band,:,:]
    # Xestim = Yh[Band,:,:] * np.random.rand(Lin,Col) 
    # Xestim = np.zeros(Xinit.shape)
    Xestim = np.random.rand(Lin,Col) 
    
    Xestim_prev = Xestim*100
    
    i = 0
    
    Lambda = lambda_list[k]
    
    while i < niter and np.linalg.norm(Xestim - Xestim_prev)/np.linalg.norm(Xestim_prev) > 1e-5 :
    
        # Compute the gradient of the smoothed TV functional.
        Gr = grad(Xestim)
        d = np.sqrt(epsilon**2 + np.sum(Gr**2, axis=2))
        # d = np.sqrt(epsilon**2 + np.linalg.norm(Gr))
        G = -div(Gr[:,:,0] / d,Gr[:,:,1] / d )
        
        # step
        e = Phi(Xestim,Hreshape) - Yh[Band,:,:]
    
        Xestim_prev = Xestim
    
        Xestim = Xestim - tau*( Phi(e,Hreshape) + Lambda*G)
        # energy
        E[i] = 1/2*np.linalg.norm(e.flatten())**2 + Lambda*np.sum(d.flatten())
    
        i = i+1
    
    err[k] = snr(Xinit,Xestim)
    
    if err[k] > snr(Xinit,Xbest):
        Xbest = Xestim
"""-------------Version fonctionnelle---------"""
# Band = 3541

# Yh = fits.getdata(HS_IM)

# M = fits.getdata(DATA+'M_1.fits').T
# A = fits.getdata(DATA+'A.fits')

# A = A[:,:,Start:Start+NbcolMA]

# n, p, q = A.shape

# Xinit = np.reshape(np.dot(M[Band], np.reshape(A, (n, p*q))), (p, q))

# Nband,Lin,Col = Yh.shape

# H = get_spa_bandpsf_hs(Band, sigma=0)

# Hifft = np.fft.ifft2(H)
        
# Hifft = np.fft.ifftshift(Hifft)
    
# ind = np.unravel_index(np.argmax(Hifft, axis=None), Hifft.shape)
        
# RLbnd = ind[0] - Lin//2
# RUbnd = ind[0] + Lin//2
# CLbnd = ind[1] - Col//2
# CUbnd = ind[1] + Col//2

# Nshift = (Col - H.shape[1])//2

# # Hreshaped[:,Nshift:(Nshift+H.shape[1])] = Hifft[RLbnd:RUbnd,:] 

# Hreshape = np.fft.fftshift(Hifft[RLbnd:RUbnd,CLbnd:CUbnd])

# # PSF[k] = np.fft.fftshift(Hreshaped)
# """-----------------------------------------------------"""

# Phi = lambda x,h: np.real(np.fft.ifft2(np.fft.fft2(x,norm = 'ortho') * np.fft.fft2(h,norm = 'ortho')))

# repeat3 = lambda x,k: resize( np.repeat( x, k, axis=1), [90, 354, k])

# # dotp = lambda x,y: np.sum(x.flatten()*y.flatten())
# # a = np.random.randn(90,354)
# # b = np.random.randn(90,354,2)
# # dotp(grad(a),b) + dotp(a,div(b[:,:,0],b[:,:,1]))


# epsilon = 0.4*1e-2

# tau = 1.9 / ( 1 + Lambda * 8 / epsilon)

# niter = 2000
# E = np.zeros((niter,1))

# lambda_list = np.linspace(0,1,20)
# err = np.zeros((len(lambda_list),1))

# Xbest = Yh[Band,:,:]

# for k in np.arange(0,len(lambda_list)):
    
#     Xestim = Yh[Band,:,:]

#     Xestim_prev = Xestim*100
    
#     i = 0
    
#     Lambda = lambda_list[k]
    
#     while i < niter and np.linalg.norm(Xestim - Xestim_prev)/np.linalg.norm(Xestim_prev) > 1e-10 :
    
#         # Compute the gradient of the smoothed TV functional.
#         Gr = grad(Xestim)
#         d = np.sqrt(epsilon**2 + np.sum(Gr**2, axis=2))
#         G = -div(Gr[:,:,0] / d,Gr[:,:,1] / d )
#         # step
#         e = Phi(Xestim,Hreshape) - Yh[Band,:,:]
    
#         Xestim_prev = Xestim
    
#         Xestim = Xestim - tau*( Phi(e,Hreshape) + Lambda*G)
#         # energy
#         E[i] = 1/2*np.linalg.norm(e.flatten())**2 + Lambda*np.sum(d.flatten())
    
#         i = i+1
    
#     err[k] = snr(Xinit,Xestim)
    
#     if err[k] > snr(Xinit,Xbest):
#         Xbest = Xestim


# def estimated_autocorrelation(x):
#     """
#     http://stackoverflow.com/q/14297012/190597
#     http://en.wikipedia.org/wiki/Autocorrelation#Estimation
#     """
#     n = len(x)
#     variance = x.var()
#     x = x-x.mean()
#     r = np.correlate(x, x, mode = 'full')[-n:]
#     assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in range(n)]))
#     result = r/(variance*(np.arange(n, 0, -1)))
#     return result

