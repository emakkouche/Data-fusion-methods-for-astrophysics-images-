#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 16:48:01 2022

@author: e.akkouche
"""
import matplotlib.pyplot as plt 
from scipy import linalg
from sklearn.utils.extmath import svd_flip
import numpy as np
from numpy import matlib


import warnings

from astropy.io import fits
from tools import compute_symmpad_3d, _centered, compute_symmpad
from skimage.transform import resize
from scipy.signal import convolve2d
# from pandeia.engine.instrument_factory import InstrumentFactory
from CONSTANTS import *

warnings.filterwarnings('ignore')
    
def my_choose_subspace(args):
    """
    Plots and saves PCA eigenvalues of the covariance matrix to choose the spectral subspace dimension.
    """
    #Ym, Yh = get_hsms(MS_IM, HS_IM) 
    Yh = fits.getdata(HS_IM)
    
    #Yh = Yh*(d**2*FLUXCONV_NC)
    
    # PCA
    L_h, S_hx, S_hy = Yh.shape
    X = np.reshape(Yh.copy(), (L_h, S_hx*S_hy)).T
    L, M = X.shape
    X_mean = np.mean(X, axis=0)
    X -= X_mean
    U, S, V = linalg.svd(X, full_matrices=False)

    return S

def my_compute_pca(X, nb_comp):
    L, M = X.shape
    X_mean = np.mean(X, axis=0)
#    X_mean=np.zeros(M)
    X -= X_mean
    U, S, V = linalg.svd(X, full_matrices=False)
    U, V = svd_flip(U, V)
    S = S[:nb_comp]
#    Z=U[:,:nb_comp]*S
#    V=V[:nb_comp]
    Z = U[:, :nb_comp]*(S**(1/2))
    V = np.dot(np.diag(S**(1/2)), V[:nb_comp])
    return V.T, Z.T, X_mean


def my_pca_nirspec(Yns, nb_comp):
    
    print('\n********** Compute PCA **********\n')
    L_h, S_hx, S_hy = Yns.shape
    X = np.reshape(Yns.copy(), (L_h, S_hx*S_hy))
    V, Z, X_mean = my_compute_pca(X.T, nb_comp)
    Z = np.reshape(Z, (nb_comp, S_hx, S_hy))
    
    print('\n********** PCA Done !!! **********\n')
    return V, Z, X_mean


def check_pca(V, Z, X_mean, Yns):
    print('\n********** PCA Checking **********\n')
    L, M, N = Z.shape
    Z = np.reshape(Z, (L, M*N))
    Yns = np.reshape(Yns, (Yns.shape[0], M*N))
    X = np.dot(V, Z)+np.matlib.repmat(X_mean, M*N, 1).T
    error = np.linalg.norm(Yns-X)/(L*M*N)
    
    print('\n********** PCA Checking Done !!! **********\n')
    return np.reshape(X, (Yns.shape[0], M, N)), error

#####COMMENTAIRES####
'''
j'ai modifié la fonction choose_subspace car je ne prend pas en compte les dégradation spectrales'
commenter gethsms car 
    je n'ai besoin que de l'image Yhs
    pas de resize de Yhs en fonction de Yms
    pas besoin de sig2hs et sig2ms

'''