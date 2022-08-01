#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:19:55 2022

@author: e.akkouche
"""

def produce_HS_nir_bis(M, A, tabwave, fname, dname, snr=50, sigma=0, d=0.31):
    """
    Produce the HS image (each band at a time)
    """
    
    flag = 1
    
    Lh = fits.getdata(DATA+'Lh.fits')

    # Get spectral PSF
    L = get_spec_psf()
    
    # Shapes
    lh, m = M.shape
    
    A = A[:,:,Start:Start+NbcolMA]
    
    n, p, q = A.shape
    p_ = int(p//(1/d)+1)
    q_ = int(q//(1/d)+1)
    # Initialize HS image
    #Y = np.zeros((lh, p_, q_))
    Y = np.zeros((lh,p,q))
    
    
    for i in range(lh):
        # Compute the 3D scene at a wavelength 
        X = compute_symmpad(np.reshape(np.dot(M[i], np.reshape(A, (n, p*q))), (p, q)), fact_pad)
        # Get spatial PSF at that band
        H = get_spa_bandpsf_hs(i, sigma)
        
        #Initialize arrays
        if( flag == 1):
            Hreshaped = np.zeros(X.shape,dtype = complex)
            
            #lin,col = X.shape
            #Y = np.zeros((lh,lin,col))
            
            #Generate Hamming window Nb_coeff = Xcol repeat Xlin 
            Hamming = np.tile(np.hamming(H.shape[1]),(X.shape[0],1))
            flag = 0
        '''---------------------------------------'''
        #Get the indexes of the maximal element of |H|
        Hshift = (np.fft.fftshift(H))             #A Supprimer
        ind = np.unravel_index(np.argmax(Hshift, axis=None), Hshift.shape)
        
        #Set the lower and upper bounds of the area to be extracted 
        RLbnd = ind[0] - NbpixR//2
        RUbnd = ind[0] + NbpixR//2
        CLbnd = ind[1] - NbpixC//2
        CUbnd = ind[1] + NbpixC//2
        
        #Reshape PSF
        
        Hreshaped[:,Nshift:(Nshift+NbpixC)] = Hshift[RLbnd:RUbnd,:] #Hreshape ==> 174x984 Hshift 174x384
        Hreshaped[:,Nshift:(Nshift+NbpixC)] *= Hamming
        Hreshaped = np.fft.fftshift(Hreshaped)
          
        #Convolve with PSF without subsample
        Y[i] = np.real(_centered(np.fft.ifft2(Hreshaped*np.fft.fft2(X, norm='ortho'), norm='ortho')[:-2, :-2], (p, q)))
    
        Hreshaped *= 0
        
        print('i = '+ str(i))  
        '''---------------------------------------'''
        
    print('\n********** Simulation Done  !!**********\n')

    return Y