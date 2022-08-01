# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:51:36 2019

@author: cguillot3

Ref 1 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
“Simulated JWST datasets for multispectral and hyperspectral image fusion”
The Astronomical Journal, vol. 160, no. 1, p. 28, Jun. 2020.

Ref 2 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
"Hyperspectral and Multispectral Image Fusion Under Spectrally Varying Spatial Blurs – Application to High Dimensional Infrared Astronomical Imaging"
IEEE Transactions on Computatonal Imaging, vol.6, Sept. 2020.

Modify paths and constants in the 'CONSTANTS.py' file, if necessary.
Run this code to test fusion.
"""

import fusion
import errors
import sys
from CONSTANTS import *
import numpy as np
from astropy.io import fits

from sparse_preprocess import get_hsms, get_throughput
from sparse_preprocess import get_hsms
from scipy import linalg

import warnings
warnings.filterwarnings('ignore')


def main(args):
    return fusion.fusion_reginf(lsub, MS_IM, HS_IM)


def choose_subspace(args):
    """
    Plots and saves PCA eigenvalues of the covariance matrix to choose the spectral subspace dimension.
    """
    #### Preprocessing of HS and MS images, operators and regularization
    Ym, Yh, tabwave, sig2 = get_hsms(MS_IM, HS_IM)
    
    #Lh = get_throughput(tabwave, 'f170lp', 'g235h')
    Lh = fits.getdata(LH)
    
    # Perform PCA on the HS image
    # Reshapes
    l, m, n = Yns.shape
    Yh_ = np.reshape(np.dot(np.diag(Lh**-1), np.reshape(Yh, (l, m*n))), (l, m, n))
    # PCA
    L_h, S_hx, S_hy = Yh_.shape
    X = np.reshape(Yh_.copy(), (L_h, S_hx*S_hy)).T
    L, M = X.shape
    X_mean = np.mean(X, axis=0)
    X -= X_mean
    U, S, V = linalg.svd(X, full_matrices=False)
    plt.semilogy(S)
    plt.title('Eigenvalues -- PCA')
    plt.savefig(SAVE+'eigenvalues_pca.png')
    plt.show()
    return S


if __name__ == "__main__":
    #image = main(sys.argv)
    #eigen = choose_subspace(args)
    print("Argument list ")