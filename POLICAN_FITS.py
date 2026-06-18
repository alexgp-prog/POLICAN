# -*- coding: utf-8 -*-
"""
Created on Sat Sep 13 10:57:22 2025

@author: AleX
"""
#### ARCHIVO DE GENERACIÓN DE FITS
import PolarStokes as stokes
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
# import matplotlib.colors as colors
import numpy as np
from astropy.io import fits



# Nombre del objeto y archivo, se introduce la dirección de los archivos.
obj = 'M1_H'
root = 'HM1H'
date = '2024_12_01_02'
# cap = 'D:/Alex/INAOE/DOCTORADO/Observaciones/Reducciones/M1_Anterior/2017_02_08_09_J/'

# date = '2024_11_04_05/2024_11_08_09'
cap = '' # 'D:/Alex/INAOE/DOCTORADO/Observaciones/Reducciones/Observacion_Nov_2022/2022_11_22_23/M1_H/M1_OMAR/' 

### PARAMETROS A CAMBIAR
sigmask = False
savefits = True
save = False

# Nombre del archivo
file1 = cap+root+'00_1b.fits'
file2 = cap+root+'22_1b.fits'
file3 = cap+root+'45_1b.fits'
file4 = cap+root+'67_1b.fits'

# EMPIEZA PROGRAMA

SRA = 1
bz = 1
bzs = 1
sigsnr = 0.01
sigstokes = 0.2
ZP = 20.6
doubleI = False

###
Nx=0
Ny=0
N = []
plt.close('all')
plt.style.use('default')#'default')

#Lectura de los archivos.
im1, head2 = stokes.openfits(file1)
im2, head1 = stokes.openfits(file2)
im3, head3 = stokes.openfits(file3)
im4, head4 = stokes.openfits(file4)

n1 = head1['NCOMBINE']
n2 = head2['NCOMBINE']
n3 = head3['NCOMBINE']
n4 = head4['NCOMBINE']

exptime = head1['EXPTIME']
print('=========================================================')
print('n1: ',n1)
print('n2: ',n2)
print('n3: ',n3)
print('n4: ',n4)
print('=========================================================')

# ANGULO SOLIDO DEL PIXEL COMPUESTO
PS = 0.32 #arcosegundos
PS = (1*0.32/3600)*(np.pi/180) # Radianes
dSolid = PS*PS #(1-np.cos(PS))
print('Angulo solido de 1 PIXEL COMPUESTO: ', dSolid, 'str')
print('=========================================================')

Sz = np.size(im1,axis = 0)
Nx = int(np.size(im1,axis = 1)/2)
Ny = int(np.size(im1,axis = 0)/2)



im00,im22,im45,im67 = stokes.centro(Sz, Nx, Ny, im1, im2, im3, im4, SRA)
I_imagen = (im00 + im22 + im45 + im67)/2

im00,im22,im45,im67,im00sd,im22sd,im45sd,im67sd,Fzx,Fzy = stokes.binning(bz,im00,
                                                                          im22,im45,
                                                                          im67,n1,n2,
                                                                          n3,n4)
I1 = (im00 + im22 + im45 + im67)/2
n = n1+n2+n3+n4

im00,im22,im45,im67,im00sd,im22sd,im45sd,im67sd,snr00,snr22,snr45,snr67 = stokes.umbral(bz,im00, im22, 
                                                            im45, im67, 
                                                            im00sd, im22sd, 
                                                            im45sd, im67sd, 
                                                            sigsnr,sigmask,
                                                            exptime,ZP,
                                                            n,n,n,n)

I,Isd,Qf,Qfsd,Uf,Ufsd,P,Psd,Pbin,A,Asd,Abin = stokes.stokes(im00, im22, 
                                                            im45, im67, 
                                                            im00sd, im22sd, 
                                                            im45sd, im67sd,
                                                            Fzx,Fzy,bzs,
                                                            sigstokes,
                                                            sigmask, doubleI)

stokes.figure(I,'Stokes I',imstd=I,z=1,save=save, comap='inferno')
stokes.figure(Isd,'sigma_I',z=0,imstd=Isd,comap='nipy_spectral',save=save)
stokes.figure(Qf,'Stokes Q',z=0,comap='nipy_spectral',save=save)
stokes.figure(Qfsd,'sigma_Q',z=0,imstd=Qfsd,comap='nipy_spectral',save=save)
stokes.figure(Uf,'Stokes U',z=0,comap='nipy_spectral',save=save)
stokes.figure(Ufsd,'sigma_U',z=0,imstd=Qfsd,comap='nipy_spectral',save=save)
stokes.figure(P,'Polarization fraction',z=0,comap='jet',save=save)
stokes.figure(Psd,'sigma_P',z=0,comap='jet',save=save)

if (savefits):
    
    # head1['RA'] = '15:26:30.3'#'15:01:17.2'#'15:26:30.3'
    # head1['DEC'] = '+54:37:06'#'+63:57:12.9' #'+54:37:06'
    head1['SIGUNIT']='MJy/sr'
    # hdu = fits.PrimaryHDU(data = I_imagen,header = head1)
    # hdu.writeto(cap + root + '_I.fits',overwrite=True)
    # hdu.writeto(cap + root + '_I.fits',overwrite=True)
    
    head1['CRPIX1'] = head1['CRPIX1']/(bz)
    head1['CRPIX2'] = head1['CRPIX2']/(bz)
    head1['CD1_1'] =  head1['CD1_1']*bz
    head1['CD2_2'] =  head1['CD2_2']*bz
    
    hdu = fits.PrimaryHDU(data = [],header = head1)
    
    hduI = fits.ImageHDU(data = I,header = head1)
    hduIsd = fits.ImageHDU(data = Isd,header = head1)
    
    hduQ = fits.ImageHDU(data = Qf * I,header = head1)
    hduQsd = fits.ImageHDU(data = Qfsd,header = head1)
    # hdu.writeto(cap + root + '_Q.fits',overwrite=True)
    
    hduU = fits.ImageHDU(data = Uf * I,header = head1)
    hduUsd = fits.ImageHDU(data = Ufsd,header = head1)
    # hdu.writeto(cap + root + '_U.fits',overwrite=True)
    
    hduP = fits.ImageHDU(data = P,header = head1)
    hduPsd = fits.ImageHDU(data = Psd,header = head1)
    # hdu.writeto(cap + root + '_P.fits',overwrite=True)
    
    hduA = fits.ImageHDU(data = A,header = head1)
    hduAsd = fits.ImageHDU(data = Asd,header = head1)
    # hdu.writeto(cap + root + '_A.fits',overwrite=True)
    
    hdulIQU = fits.HDUList([hdu, hduI,hduIsd,hduQ,hduQsd,hduU,hduUsd,hduP,hduPsd,
                            hduA,hduAsd])
    
    hdulIQU[0].name = 'PRIMARY'
    hdulIQU[1].name = 'STOKES I'
    hdulIQU[2].name = 'SIGMA I'
    
    hdulIQU[3].name = 'STOKES Q'
    hdulIQU[4].name = 'SIGMA Q'
    
    hdulIQU[5].name = 'STOKES U'
    hdulIQU[6].name = 'SIGMA U'
    
    hdulIQU[7].name = 'POL_DEGREE'
    hdulIQU[8].name = 'SIGMA POL DEGREE'
    
    hdulIQU[9].name = 'POL_ANGLE'
    hdulIQU[10].name = 'SIGMA POL_ANGLE'
    
    print('SAVE IN: ', cap + root + '_POLICAN_TEST.fits')
    
    hdulIQU.writeto(cap + root + '_POLICAN_TEST.fits',overwrite=False,output_verify='warn')


plt.show()