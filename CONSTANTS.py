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

This code is used to store constants and paths for data fusion.
"""

"""
Paths for data and saving files
"""
#DATA = 'Data/'
DATA = '/media/e.akkouche/DONNEES/Fusion_data/Data/'
#SAVE = 'Expe/'
SAVE = '/media/e.akkouche/DONNEES/Fusion_data/Save_result/'
SAVE2 = '/media/e.akkouche/DONNEES/Fusion_data/result_sobolev_02_09_Mfft/'
SAVE_IMG = '/media/e.akkouche/DONNEES/Fusion_data/result_img/'
PSF = DATA+'PSF/'
SAVE_PATH = SAVE+'z_optim.fits'


"""
File name and path for hyperspectral and multispectral images
"""
V_acp = DATA+'V.fits'
Z_acp = DATA+'Z.fits'

MS_IM = DATA+'Ym_highsnr_hst.fits'
HS_IM = DATA+'Yh_highsnr_hst.fits'
# MS_IM = SAVE+'Ym_lowsnr_10.fits'
# HS_IM = SAVE+'Yh_lowsnr_10.fits'
# MS_IM = SAVE+'Ym_highsnr_psfsig-1.fits'
# HS_IM = SAVE+'Yh_highsnr_psfsig-1.fits'

"""
Dimensions de reshape des PSF Dim(MA) = 90x900 fact_pad = 41
IL faut que nc > NbpixC pour Nshift>0 donc je peux ajouter des zeros a gauche et droite pour lisser ensuite avec Hamming
"""
START = 500
NBCOL_X = 354
NbpixR =  150# = nr 
NbpixC = 384  #Nb col de PSF
Nshift = (414 - NbpixC)//2 #414 Nb col de MA apres padding

"""
Subspace dimension
"""
lsub = 4
Lacp = 10

"""
Downsampling factor (integer)
"""
d = 3

"""
High spatial resolution image size with padding and padding factor
Example :
Original MS band size : 300x300
fact_pad : 41
nr = 300 + 2*fact_pad + 2 
nc = 300 + 2*fact_pad + 2
Conditions (to chose fact_pad) : nr % d = 0; nc % d = 0; (2*fact_pad+2) % d = 0; (fact_pad//d + 1) % 2 = 0
"""

"""
Origiginal HS size : 90x354
"""
fact_pad = 29
NR =150 #90+2*29+2  =150 nrow of MA after padding
NC =414 #354+2*29+2  = 414 ncol of MA after padding

"""
Conversion constants : mJy/arcsec**2 to mJy/pixel (depends on the observation instruments)
"""
FLUXCONV_NC = 0.031**2
FLUXCONV_NS = 0.01

"""
*******************************************************************************************************
CONSTANTS FOR NOISE GENERATION, only used to generate/simulate HS and MS data
*******************************************************************************************************
"""

"""
Path for H2RG correlation matrix for noise generation (only used to generate/simulate HS and MS data)
"""
CORR = DATA+'h2rg_corr.fits'

"""
Constants for noise generation
"""
EXPTIME_NC = 128.841 # Exposure time in seconds for the MS image
SIG2READ_NC = 2*(16.2)**2/2 # Noise constant in e-/s for the MS image
# SIG2READ_NC = 29000
NFRAMES_NC = 2 # Number of frames during the MS image acquisition
EXPTIME_NS = 257.682 # Exposure time in seconds for the HS image
SIG2READ_NS = 1.6*(6/88)**2 # Noise constant in e-/s for the HS image
# SIG2READ_NS = 500
NFRAMES_NS = 5 # Number of frames during the HS image acquisition
GAIN_NC = 2 # Gain in e-/ph for MS instrument
GAIN_NS = 1 # Gain in e-/ph for HS instrument