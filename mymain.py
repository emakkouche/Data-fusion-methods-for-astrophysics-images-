#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:26:35 2022

@author: e.akkouche
"""

# import myproduce_HS_MS
import mycheck_pca
# import precompute

import simulate_HS_MS

import matplotlib.pyplot as plt

import myproduce_HS_MS
import mycheck_pca


from tempo import *
from astropy.io import fits
from tools import *

from CONSTANTS import *

import warnings
warnings.filterwarnings('ignore')

"""--------------------Generate image X and Yh--------------------"""
#Yh = myproduce_HS_MS.main(DATA,DATA)
#Yh = simulate_HS_MS.main(DATA, DATA,sigma2 = 0.1)
"""--------------------Choose subspace dimension------------------"""
#S = mycheck_pca.my_choose_subspace(HS_IM)

#Affichage complet
# plt.semilogy(S);
# plt.ylabel('Valeurs propres')
# plt.xlabel('Indice de la composante principale')
# plt.show()
    
#Affichage tronqué
# plt.semilogy(S[:nbeig]);
# plt.xticks(np.arange(len(S[:nbeig])), np.arange(1, len(S[:nbeig])+1))
# plt.ylabel('Valeurs propres')
# plt.xlabel('Indice de la composante principale')

#plt.savefig(SAVE+'eigenvalues_pca.png')

#Afficher variance cumulee
# plt.plot(cumul_var)
# plt.title('Cumulative Explained Variance -- PCA')
# plt.xlabel('Dimensions du sous-espace')
"""--------------------Total Variance PCA----------------------"""
# Yh = fits.getdata(HS_IM)
Lacp = 10

#V, Z, mean = mycheck_pca.my_pca_nirspec(Yh,Lacp)

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
Yh = fits.getdata(HS_IM)
V = fits.getdata(V_acp)
Z = fits.getdata(Z_acp)

#Xback,error = mycheck_pca.check_pca(V, Z, mean, Yh)

#hdu = fits.PrimaryHDU(np.real(np.fft.ifft2(Difference, axes=(1, 2), norm='ortho')))
#hdu.writeto(DATA+'Difference.fits', overwrite=True)

# Difference = Objfun.compute_diff(Yh, V, Z)
# print('CritJ = ',tempo.test(Yh,V,Z));


"""-----Init et Stockage----"""
Jmu = []
Regmu = []
J2mu = []
#Time
t1 = time()

"""--------------------Denoising--------------------"""
#Yh = fits.getdata(HS_IM)
#V = fits.getdata(V_acp)
#Z = fits.getdata(Z_acp)

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
mus = np.log(np.linspace(1,3,10))

err = np.zeros((len(mus),1))

Band = 1000

m = 0

Xbest = Yh[Band,:,:]

for mu in mus:
    
        print('----------------- MU '+str(mu)+' -----------------')

        Zoptim,J_Zoptim,step = GD(Yfft,V,Zfft,H,Lacp,T2,T1,D,mu,maxH2)
        
        # Jmu.append(np.linalg.norm(Yfft-np.dot(V,Zoptim)*H)**2)
        # J2mu.append(J_Zoptim)
        # Regmu.append(np.linalg.norm((D[0]+D[1])*np.reshape(Zoptim,(Lacp,NR,NC)))**2)

        Zifft = postprocess(Zoptim,Lacp)
        
        Xestim = recover_X(V, Zifft)
        
        Xtrue = get_Xtrue(DATA,Band)
        
        err[m] = snr(Xtrue,Xestim[Band,:,:])
        
        if err[m] > snr(Xtrue,Xbest):
            Xbest = Xestim[Band,:,:]
     
        m = m+1
        
        fname = SAVE2+'Zoptim_mu_'+str(m)+'.fits'
        save_Zoptim(Zifft, fname)
        
        fname = SAVE2+'J_mu_'+str(m)
        np.save(fname,J_Zoptim)
 
      
fname = SAVE2+'SNR_sobolev'
np.save(fname,err)


t2 = time()


print('******************************************')
print('******* TOTAL COMPUTATION TIME : '+str(np.round((t2-t1)/60))+'min '+str(np.round((t2-t1)%60))+'s.')


plt.plot(mus,err)
plt.axis('tight')
plt.xlabel('\lambda (échelle log)') 
plt.ylabel('SNR')
plt.show()

plt.imshow(Xbest)
plt.title('Variation totale SNR = '+str(round(np.max(err),2))+'dB') 


"""----------------------------"""
# for mu in mus:
    
#     Zmu = fits.getdata(SAVE2+'_full_'+str(mu)+'z_opti.fits')
    
#     Zmu = np.fft.fft2(tools.compute_symmpad_3d(Zmu,fact_pad),axes=(1, 2),norm='ortho')
    
#     Zmu = np.reshape(Zmu,(Zmu.shape[0],Zmu.shape[1]*Zmu.shape[2]))
    
#     Jmu.append(np.linalg.norm(Yfft-np.dot(V,Zmu)*H)**2)
    
#     Regmu.append(np.linalg.norm((D[0]+D[1])*np.reshape(Zmu,(Lacp,NR,NC))))
    
"""--------------------Postprocess & saving--------------------"""
# save_Zoptim(Zifft, SAVE)

# Zoptim = fits.getdata(SAVE_PATH)

# Band = 4900

# Xestim = recover_X(V, Zifft)

# Xtrue = get_Xtrue(DATA,Band)

"""--------------------Affichage--------------------"""
# plt.imshow(Xestim[Band,:,:])
# plt.show()
# plt.imshow(Yh[Band,:,:])
# plt.show()
# plt.imshow(Xtrue)
# plt.show()

# np.save(SAVE+'J_Zoptim',J_Zoptim)
# plt.semilogy(J)


# plt.plot(J_Zoptim)
# plt.yscale('log')
# plt.title('Evolution de la fonction coût '+'(pas = ' +str(step)+')')
# plt.xlabel('Itération')
# plt.savefig(SAVE+'Evolution_fonction_coût_'+'pas_' +str(step))
# plt.show()

# JZe_7 = np.load(SAVE+'J_Zoptim_1e7.npy')
# plt.plot(JZe_7,"-b", label="pas = 1e-7")
# plt.plot(J_Zoptim, "-r", label="pas = 1e-6")
# plt.legend(loc="upper right")
# plt.title('Evolution de la fonction coût au cours des itérations')
# plt.xlabel('Itération')
# plt.savefig(SAVE+'Evolution_fonction_coût_2_pas')

# plt.plot(Regmu,Jmu)
# plt.title('Courbe en L')
# plt.xlabel('Regul Sobolev(mu)')
# plt.ylabel('Attache aux données(mu)')

# plt.imshow(Xestim[Band,:,:])
# plt.title('Sobolev SNR = '+str(round(snr(Xtrue,Xestim[Band,:,]),2))+'dB')  


#Zmu = fits.getdata(SAVE2+'_full_4.641588833612778z_opti.fits')
#Xmu = recover_X(V, Zmu)

#plt.imshow(Xmu[Band,:,:])
#plt.plot(Xmu[Band,63,:])

#plt.plot(Xtrue[63,:], 'r',label = 'Référence') # plotting t, a separately 
#plt.plot(Xmu[Band,63,:], 'b',label = 'Sobolev') # plotting t, b separately 
#plt.legend(loc="upper right")


#plt.plot()