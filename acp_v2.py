#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 13:54:27 2018

@author: cguillot3

Ref 1 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
“Simulated JWST datasets for multispectral and hyperspectral image fusion”
The Astronomical Journal, vol. 160, no. 1, p. 28, Jun. 2020.

Ref 2 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
"Hyperspectral and Multispectral Image Fusion Under Spectrally Varying Spatial Blurs – Application to High Dimensional Infrared Astronomical Imaging"
IEEE Transactions on Computatonal Imaging, vol.6, Sept. 2020.

This code implents PCA performed on the HS image for spectral dimension reduction.
"""

import numpy as np
from scipy import linalg
from sklearn.utils.extmath import svd_flip


def compute_pca(X, nb_comp):
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


def pca_nirspec(Yns, nb_comp):
    L_h, S_hx, S_hy = Yns.shape
    X = np.reshape(Yns.copy(), (L_h, S_hx*S_hy))
    V, Z, X_mean = compute_pca(X.T, nb_comp)
    Z = np.reshape(Z, (nb_comp, S_hx, S_hy))
    return V, Z, X_mean


# def check_pca(V, Z, X_mean, Yns):
#     L, M, N = Z.shape
#     Z = np.reshape(Z, (L, M*N))
#     Yns = np.reshape(Yns, (Yns.shape[0], M*N))
#     X = np.dot(V, Z)+np.matlib.repmat(X_mean, M*N, 1).T
#     error = np.linalg.norm(Yns-X)/(L*M*N)
#     return np.reshape(X, (Yns.shape[0], M, N)), error
