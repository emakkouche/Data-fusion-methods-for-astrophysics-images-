#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 13:21:34 2019

@author: claireguilloteau

Ref 1 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
“Simulated JWST datasets for multispectral and hyperspectral image fusion”
The Astronomical Journal, vol. 160, no. 1, p. 28, Jun. 2020.

Ref 2 : C. Guilloteau, T. Oberlin, O. Berné, É. Habart, and N. Dobigeon
"Hyperspectral and Multispectral Image Fusion Under Spectrally Varying Spatial Blurs – Application to High Dimensional Infrared Astronomical Imaging"
IEEE Transactions on Computatonal Imaging, vol.6, Sept. 2020.

This code implents the fusion process of infrared astronomical hyperspectral and multispectral images as described in [1,2].
"""
import numpy as np
import numpy.matlib
from astropy.io import fits
from time import time
from tools import _centered
from mysparse import get_linearsyst_reginf
import errors
from CONSTANTS import *

import warnings
warnings.filterwarnings('ignore')



""" ******************* Regularization term gradient computation ************************ """

def compute_Areg(lacp, mu1, D, Wd, z):
    #### Smooth spatial regularization with weights.
    # Compute AregZ = W^2 . (D^T D Z)
    z = np.reshape(z, (lacp, Wd[0].shape[1], Wd[0].shape[2]))
    gx = np.fft.fft2(np.real(np.fft.ifft2(z*D[0]))*Wd[0]**2, norm='ortho')*np.conj(D[0])
    gy = np.fft.fft2(np.real(np.fft.ifft2(z*D[1]))*Wd[1]**2, norm='ortho')*np.conj(D[1])
    return 2*mu1*np.reshape(gx+gy, np.prod(gx.shape))


""" ******************* CONJUGATE GRADIENT ROUTINE ************************ """


def cg(lacp, mu1, A, D, Wd, b, c, z, save_it=False):
    #### Conjugate gradient iterations #
    print('--- CONJUGATE GRADIENT ALGORITHM ---')
    t1 = time()
    nb_it = 0
    z0 = z.copy()
    # Initialize

    ########## Control procedure ##########
    # print('NANs in Az, b and A_reg z : '+str(np.sum(np.isnan(A.dot(z))))+' ; '+str(np.sum(np.isnan(b)))+' ; '+str(np.sum(np.isnan(compute_Areg(lacp, mu1, D, Wd, z)))))
    #######################################

    r = A.dot(z)+b+compute_Areg(lacp, mu1, D, Wd, z)
    p = -r.copy()
    # Objective function
    obj = [0.5*np.dot(np.conj(z).T, A.dot(z))+np.dot(np.conj(b.T), z)+c+0.5*np.dot(np.conj(z).T, compute_Areg(lacp, mu1, D, Wd, z))]
    print(str(nb_it)+' -- Objective function value : '+str(obj[nb_it]))
    
    ########## Control procedure ##########
    # if np.isnan(obj):
    #     print('Objective value is NAN :')
    #     print('NANs in 0,5* zt A z : '+str(np.sum(np.isnan(0.5*np.dot(np.conj(z).T, A.dot(z))))))
    #     print('NANs in bt z : '+str(np.sum(np.isnan(np.dot(np.conj(b.T), z)))))
    #     print('NANs in 0,5* zt A_reg z : '+str(np.sum(np.isnan(0.5*np.dot(np.conj(z).T, compute_Areg(lacp, mu1, D, Wd, z))))))
    ######################################

    # Stopping criterion
    stop = obj[nb_it]
    # Iterations
    while stop > 1e-8 and nb_it < 2000:
        areg = compute_Areg(lacp, mu1, D, Wd, p)
        alpha = np.dot(np.conj(r).T, r)*(np.dot(np.conj(p).T, A.dot(p))+np.dot(np.conj(p).T, areg))**-1
        z = z+alpha*p
        r_old = r.copy()
        r = r+alpha*(A.dot(p)+areg)
        beta = np.dot(np.conj(r).T, r)/np.dot(np.conj(r_old).T, r_old)
        p = -r+beta*p
        nb_it += 1
        obj.append(0.5*np.dot(np.conj(z).T, A.dot(z))+np.dot(np.conj(b.T), z)+c+0.5*np.dot(np.conj(z).T, compute_Areg(lacp, mu1, D, Wd, z)))
        print(str(nb_it)+' -- Objective function value : '+str(obj[nb_it]))
        stop = (obj[-2]-obj[-1])/obj[-2]
        if save_it:
            hdu = fits.PrimaryHDU(postprocess(z, lacp))
            hdu.writeto(SAVE+'z_'+str(nb_it+1)+'.fits')

    t2 = time()
    print('Cg Computation time : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')
    return z, obj


############################## POST-PROCESSING ###################################


def postprocess(z, lacp):
    # Reshape and ifft on the solution of the conjugate gradient.
    z = np.reshape(z, (lacp, nr, nc))
    z = np.fft.ifft2(z, norm='ortho')
    z = np.real(_centered(z[:, :-2, :-2], (lacp, nr-2*fact_pad, nc-2*fact_pad)))
    return z


def save(z, filename):
    # Save solution
    hdu = fits.PrimaryHDU(z)
    hdu.writeto(filename+'z_opt.fits', overwrite=True)


def calc_mu_auto(Am, Ah, bm, bh, cm, ch, zopt, filename):
    # Select the regularization parameter
    sm = 0.5*np.dot(np.conj(zopt).T, Am.dot(zopt)) + np.dot(np.conj(bm).T, zopt) + cm
    # print('0.5*np.dot(np.conj(zopt).T, Am.dot(zopt)) = '+str(0.5*np.dot(np.conj(zopt).T, Am.dot(zopt))))
    # print('np.dot(np.conj(bm).T, zopt) = '+str(np.dot(np.conj(bm).T, zopt)))
    # print('cm = '+str(cm))
    # print('sm = '+str(sm))
    sh = (0.5*np.dot(np.conj(zopt).T, Ah.dot(zopt)) + np.dot(np.conj(bh).T, zopt) + ch)/9
    # print('0.5*np.dot(np.conj(zopt).T, Ah.dot(zopt)) = '+str((0.5*np.dot(np.conj(zopt).T, Ah.dot(zopt)))))
    # print('np.dot(np.conj(bh).T, zopt) = '+str(np.dot(np.conj(bh).T, zopt)))
    # print('ch = '+str(ch))
    # print('sh = '+str(sh))
    hdu = fits.PrimaryHDU([np.real(sm), np.real(sh)])
    hdu.writeto(filename+'_smsh_.fits', overwrite=True)


############################# GETTERS #########################################


def get_vz_true(filename=DATA):
    #### Get the ground truth
    #
    m = fits.getdata(filename+'M_1_conv.fits')
    a = fits.getdata(filename+'A.fits')
    # A = fits.getdata('/usr/local/home/cguillot3/Fusion_dataset/dataset1/A.fits')
    # A = fits.getdata('/usr/local/home/cguillot3/Fusion_dataset/Hst/A_hst.fits')[:, 100:400, 100:400]
    # A = fits.getdata(filename+'A_hst.fits')[:, 100:400, 100:400]
    # for i in range(4):
        # Amin = np.min(A[i])
        # Amax = np.max(A[i]-Amin)
        # A[i] = (A[i]-Amin)/Amax
    return m, a


def get_v_mean(filename=SAVE2):
    #### Get the spectra matrix learned from the HS image, and the mean spectrum.
    return fits.getdata(filename+'V_1_mjy.fits'), fits.getdata(filename+'mean_1_mjy.fits')


############################### MAIN FUSION PROCEDURE ############################


def fusion_reginf(lacp=4, MS_IM=MS_IM, HS_IM=HS_IM):
    # Set regularization parameters
    mus1 = 10**np.linspace(0, 2, 10)
    # Time
    t1 = time()
    # Preprocessing
    Am, Ah, bm, bh, cm, ch, z, D, Wd = get_linearsyst_reginf(lacp, MS_IM, HS_IM)
    # Get the ground truth and spectra matrix
    v, mean = get_v_mean()
    v_true, z_true = get_vz_true()
    # For each reg. parameter, solve the problem

    ####### Control procedure ####### 
    # fname = SAVE2+'_control_0'
    # calc_mu_auto(Am, Ah, bm, bh, cm, ch, z.copy(), fname)
    #################################

    for mu1 in mus1:
        print('----------------- MU1 '+str(mu1)+' -----------------')
        # Solve
        zf, obj = cg(lacp, mu1, Am+Ah, D, Wd, bm+bh, cm+ch, z.copy(), save_it=False)
        # Post-process and save
        zf_ = postprocess(zf, lacp)
        fname = SAVE2+'_full_'+str(mu1)
        save(zf_, fname)
        # Compute and save errors
        # errors.compute_errors(v_true, v, z_true, zf_, mean, fname)
        calc_mu_auto(Am, Ah, bm, bh, cm, ch, zf, fname)
    t2 = time()
    print('******************************************')
    print('******* TOTAL COMPUTATION TIME : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')
    return zf
