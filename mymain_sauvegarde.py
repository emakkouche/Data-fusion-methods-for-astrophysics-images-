#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:26:35 2022

@author: e.akkouche
"""

# import myproduce_HS_MS
# import mycheck_pca
# import precompute

import matplotlib.pyplot as plt

import myproduce_HS_MS
import mycheck_pca
from tempo import *
from astropy.io import fits
from tools import *

import matplotlib.pyplot as plt

from CONSTANTS import *

import warnings
warnings.filterwarnings('ignore')

"""--------------------Generate image X and Yh--------------------"""
# Yh = myproduce_HS_MS.main(DATA,DATA)
#--------------------Check generated images --------------------
#mycheck_pca.mycheck_HS_MS_img(HS_IM,Band)
#--------------------PCA --------------------
# Eigvalue,cumul_var = mycheck_pca.my_choose_subspace(HS_IM)

# plt.semilogy(Eigvalue)
# # plt.yscale('log')
# plt.title('Eigenvalues -- ACP')
# plt.savefig(SAVE+'eigenvalues_pca.png')
# plt.show()

# plt.plot(cumul_var[:Lacp])
# plt.title('Cumulative Variance Explained -- PCA')
# plt.xlabel('PC index')
# plt.savefig(SAVE+'Cumulative_variance_pca.png')
# plt.show()
"""--------------------Total Variance PCA--------------------"""
# Yh = fits.getdata(HS_IM)
# Lacp = 10

# V, Z, mean = mycheck_pca.my_pca_nirspec(Yh,Lacp)

# # print('\n********** Saving PCA **********\n')

# hdu = fits.PrimaryHDU(Z)
# hdu.writeto(DATA+'Z.fits', overwrite=True)

# hdu = fits.PrimaryHDU(V)
# hdu.writeto(DATA+'V.fits', overwrite=True)

# hdu = fits.PrimaryHDU(mean)
# hdu.writeto(DATA+'mean.fits', overwrite=True)

# print('\n********** Saving Done !! **********\n')

# V = np.load(DATA+'V.npy')
# Z = np.load(DATA+'Z.npy')
"""--------------------Check PCA--------------------"""
# Yh = fits.getdata(HS_IM)
# V = fits.getdata(V_acp)
# Z = fits.getdata(Z_acp)

# Xback,error = mycheck_pca.check_pca(V, Z, mean, Yh)

#hdu = fits.PrimaryHDU(np.real(np.fft.ifft2(Difference, axes=(1, 2), norm='ortho')))
#hdu.writeto(DATA+'Difference.fits', overwrite=True)

# Difference = Objfun.compute_diff(Yh, V, Z)
# print('CritJ = ',tempo.test(Yh,V,Z));

"""--------------------Denoising--------------------"""
Yh = fits.getdata(HS_IM)
V = fits.getdata(V_acp)
Z = fits.getdata(Z_acp)

B,L,C = Z.shape
Z = np.random.rand(B,L,C)
"""--------------------Preprocessing--------------------"""
Yfft,Zfft,H,T1,T2,maxH2,D = preprocessing(Yh,V,Z,Lacp)

Lband,Lin,Col = Yfft.shape
Yfft = np.reshape(Yfft,(Lband,Lin*Col))

H = np.reshape(H,(Lband,Lin*Col))

Lband,Lin,Col = Zfft.shape
Zfft = np.reshape(Zfft,(Lband,Lin*Col))
"""--------------------Gradient Descent--------------------"""
"""faire boucle ici"""
mu = 5
Zoptim,J_Zoptim,step = GD(Yfft,V,Zfft,H,Lacp,T2,T1,D,mu,maxH2)


"""--------------------Postprocess & saving--------------------"""
Zifft = postprocess(Zoptim,Lacp)

save_Zoptim(Zifft, SAVE)

# Zoptim = fits.getdata(SAVE_PATH)

Band = 4900
Xestim = recover_X(V, Zifft)

Xtrue = get_Xtrue(DATA,Band)

"""--------------------Affichage--------------------"""
plt.imshow(Xestim[Band,:,:])
plt.show()
plt.imshow(Yh[Band,:,:])
plt.show()
plt.imshow(Xtrue)
plt.show()

# np.save(SAVE+'J_Zoptim',J_Zoptim)
# plt.semilogy(J)
plt.plot(J_Zoptim)
plt.yscale('log')
plt.title('Evolution de la fonction coût '+'(pas = ' +str(step)+')')
plt.xlabel('Itération')
plt.savefig(SAVE+'Evolution_fonction_coût_'+'pas_' +str(step))
plt.show()

# JZe_7 = np.load(SAVE+'J_Zoptim_1e7.npy')
# plt.plot(JZe_7,"-b", label="pas = 1e-7")
# plt.plot(J_Zoptim, "-r", label="pas = 1e-6")
# plt.legend(loc="upper right")
# plt.title('Evolution de la fonction coût au cours des itérations')
# plt.xlabel('Itération')
# plt.savefig(SAVE+'Evolution_fonction_coût_2_pas')