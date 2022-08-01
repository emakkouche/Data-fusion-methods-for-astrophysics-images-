
import warnings
import numpy as np
import numpy.matlib
from astropy.io import fits
from skimage.transform import resize
from time import time
import acp_v2 as pca
from CONSTANTS import *
import errors
# import RGB

warnings.filterwarnings('ignore')


# Code in python for HS ans MS image fusion performed with the Brovey method.
# Ref :

def get_vz_true(filename=DATA):
    m = fits.getdata(filename+'M_1_conv.fits')
    # a = fits.getdata(filename+'A.fits')
    # a = fits.getdata('/usr/local/home/cguillot3/Fusion_dataset/dataset1/A.fits')
    A = fits.getdata('/usr/local/home/cguillot3/Fusion_dataset/Hst/A_hst.fits')[:, 100:400, 100:400]
    # A = fits.getdata(filename+'A_hst.fits')[:, 100:400, 100:400]
    for i in range(4):
        Amin = np.min(A[i])
        Amax = np.max(A[i]-Amin)
        A[i] = (A[i]-Amin)/Amax
    return m, A


def get_v_mean(sigma, filename=SAVE2):
    #### Get the spectra matrix learned from the HS image, and the mean spectrum.
    return fits.getdata(filename+'V_1_mjy.fits'), fits.getdata(filename+'mean_1_mjy.fits')


########################## INITIALIZE IMAGES ############################

def get_hsms(MS_IM, HS_IM, lacp):
    # Get HS and MS images from files MS_IM and HS_IM
    Ync = fits.getdata(MS_IM)
    Yns = fits.getdata(HS_IM)
    l, m, n = Ync.shape
    l_, m_, n_ = Yns.shape
    # Get HS spectral operator Lh and correct HS image
    Lh = fits.getdata(DATA+'pce_conv_1.fits')
    Yns = np.reshape(np.dot(np.diag(Lh**-1), np.reshape(Yns, (l_, m_*n_))), (l_, m_, n_))
    Yns = resize(Yns, (l_, m//3, n//3), order=3, mode='symmetric')*(0.93**2)
    # Denoise HS image with PCA and interpolate with bicubic splines
    Yns = denoising(Yns, lacp)
    Yns = resize(Yns, (l_, m, n), order=3, mode='symmetric')
    # Get MS spectral degradation operator Lm
    Lm = fits.getdata('/usr/local/home/cguillot3/Fourier_Res/Weighted_Reg/Lm.fits')
    return Ync, Yns, Lm


########################## DENOISING #####################################

def denoising(Yns, n_comp):
    # PCA performing, only n_comp are kept
    print(' PCA on the HS image : ')
    t1 = time()
    l, m, n = Yns.shape
    V, Z, mean = pca.pca_nirspec(Yns, n_comp)
    Z = np.reshape(Z, (n_comp, m*n))
    mean = np.matlib.repmat(mean, m*n, 1).T
    # Back into image domain
    Yns_ = np.reshape(np.dot(V, Z)+mean, (l, m, n))
    t2 = time()
    print(str(t2-t1)+'s.')
    print('------------------------------------------------------------------')
    return Yns_


########################## BROVEY ROUTINE ##################################

def brovey(MS_IM, HS_IM, lacp):
    # Get images and operators
    Ym, Yh, Lm = get_hsms(MS_IM, HS_IM, lacp)
    # Reshape 3d to 2d procedures
    l, m, n = Yh.shape
    Yh = np.reshape(Yh, (l, m*n))
    Ym = np.reshape(Ym, (Ym.shape[0], m*n))
    lm = Lm.shape[0]
    # Spatial details matrix I
    t1 = time()
    I = lm**(-1)*np.sum(Ym*(np.dot(Lm, Yh)+1e-10)**(-1), axis=0)
    # Save I
    # hdu = fits.PrimaryHDU(np.reshape(I, (m, n)))
    # hdu.writeto('Brovey/brovey_I0.fits')
    # Compute solution by injecting details I into interpolated HS image
    X_hat = Yh*I
    t2 = time()
    print('Brovey : '+str(t2-t1)+'s.')
    # Save solution
    # RGB.create_rgb_(X_hat, 'brovey0')
    hdu = fits.PrimaryHDU(np.reshape(X_hat, (l, m, n)))
    hdu.writeto('Expe_hst/brovey_solution.fits', overwrite=True)

    v, mean = get_v_mean(-1)
    v_true, z_true = get_vz_true()

    errors.compute_errors(v_true, v, z_true, X_hat, mean, 'Expe_hst/brovey_errors')
