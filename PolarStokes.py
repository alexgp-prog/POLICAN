import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import matplotlib.colors as colors
from photutils.centroids import centroid_com#, centroid_sources
from photutils.aperture import CircularAperture,EllipticalAperture, CircularAnnulus
from photutils.aperture import aperture_photometry, ApertureStats
import os
import pandas as pd
# from scipy import ndimage

def openfits(archivo):
    hdul = fits.open(archivo)
    hdul.info()
    imagen = hdul[0].data
    head = hdul[0].header
    hdul.close()
    return imagen, head

def figure(imagen,titulo,imstd=0,z=1,comap='inferno',save=False): 
    fig, ax = plt.subplots()
    ax.set_title(titulo)
    if(z == 1):
        im1 = ax.imshow(imagen,cmap=comap,
                    origin = 'lower',
                    vmax=2.0, vmin=0)
        print('DEBUG')
    if(z == 2): 
        im1 = ax.imshow(imagen,cmap=comap,
                    origin = 'lower',
                    norm=colors.LogNorm(vmax=np.max(imagen)))
    else: 
        im1 = ax.imshow(imagen,cmap=comap,
                    origin = 'lower')
    ax.set_xlabel('pixeles')
    ax.set_ylabel('pixeles')
    fig.colorbar(im1,ax=ax)
    ax.grid(False)
    if (save == True):
        plt.savefig(titulo+'.png', dpi=150)
   
def figure_hist(imagen,titulo,imstd=0,z=1,comap='inferno'):
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15,8))
    ax[0].set_title(titulo)
    if(z == 1):
        im1 = ax[0].imshow(imagen,cmap=comap,
                    origin = 'lower',
                    vmax=3*imstd.std(), 
                    vmin=-imstd.std())
    if(z == 2): 
        im1 = ax[0].imshow(imagen,cmap=comap,
                    origin = 'lower',
                    norm=colors.LogNorm(vmax=np.max(imstd)))
    else: 
        im1 = ax[0].imshow(imagen,cmap=comap,
                    origin = 'lower')
    ax[0].set_xlabel('pixeles')
    ax[0].set_ylabel('pixeles')
    fig.colorbar(im1,ax=ax[0],location='right', anchor=(0.3, 0.5), shrink=0.5)
    ax[0].grid(False)

    ax[1].hist(imagen.flatten(),log = True, bins = 100, color = 'blue')
    ax[1].axvline(np.nanmean(imagen), color = 'red', linestyle = '--', 
            linewidth = 2, 
            label = '$ mean $ = ' + str(round(np.nanmean(imagen),3)))
    ax[1].axvline(np.nanmedian(imagen), color = 'green', linestyle = '--', 
            linewidth = 2, 
            label = '$ median $ = ' + str(round(np.nanmedian(imagen),3)))
    ax[1].axvline(3*np.nanstd(imagen), color = 'yellow', linestyle = '--', 
            linewidth = 2, 
            label = '$ 3 std $ = ' + str(round(3*np.nanstd(imagen),3)))
    
    ax[1].axvline(-3*np.nanstd(imagen), color = 'yellow', linestyle = '--', 
            linewidth = 2)
    
    ax[1].legend()
    
def contour(array, arrayct, titulo, level = [0, 10, 20, 30, 40, 50, 60], comap = 'nipy_spectral'):
    x = np.linspace(0,np.size(array, axis = 1)-1, np.size(array, axis = 1))
    y = np.linspace(0,np.size(array, axis = 0)-1, np.size(array, axis = 0))
    X,Y = np.meshgrid(x,y)

    fig, ax = plt.subplots(figsize = (10,10))
    ax.set_title(titulo)
    contours = ax.contour(X,Y,arrayct, level, colors = 'black')
                          # extent = [RA_init, RA_fin, DEC_init, DEC_fin])
    ax.clabel(contours, inline = True, fontsize = 10)
    im = ax.imshow(array, origin = 'lower', cmap = comap, 
                     alpha = 0.8)
                    # extent = [RA_init, RA_fin, DEC_init, DEC_fin])
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('MJy/sr', rotation = 270,fontsize = 10)
    ax.set_xlabel('pixeles',fontsize = 10)
    ax.set_ylabel('pixeles',fontsize = 10)

def centro(Sz,Nx, Ny, im1, im2, im3, im4, SRA):
    print('Coordenadas seleccionadas: ', Nx, ' ', Ny)
    cxi = 0
    cyi = 0
    cxs = np.size(im1,axis=1)
    cys = np.size(im1,axis=0)
    
    print(im1.shape)
    
    cmaxx = False
    cmaxy = False
    
    if ((Nx - int(Sz/2))>cxi): 
        cxi = Nx - int(Sz/2)
        cmaxx = True
    else: 
        cxs = cxs + np.abs((Nx - int(Sz/2)))
        cmaxx = False
        
    if ((Nx + int(Sz/2))<cxs):
        cxs = Nx + int(Sz/2)
        cmaxx = True
    else: 
        cxi = cxi - np.abs((Nx + int(Sz/2) - cxs))
        cmaxx = False
    
    if ((Ny - int(Sz/2))>cyi): 
        cyi = Ny - int(Sz/2)
        cmaxy = True        
    else: 
        cys = cys + np.abs((Ny - int(Sz/2))) 
        cmaxy = False
    
    if ((Ny + int(Sz/2))<cys):
        cys = Ny + int(Sz/2)
        cmaxy = True
    else: 
        cyi = cyi - np.abs((Ny + int(Sz/2) - cys))
        cmaxy = False
        
    print('cxi,cxs: ',cxi,cxs)
    print('cyi,cys: ',cyi,cys)
    print('cmaxx: ', cmaxx)
    print('cmaxy: ', cmaxx)
    
    im00 = im1[cyi:cys,cxi:cxs]
    im22 = im2[cyi:cys,cxi:cxs]
    im45 = im3[cyi:cys,cxi:cxs]
    im67 = im4[cyi:cys,cxi:cxs]
    
    if((cmaxx and cmaxy)):# and (cyi != 0 and cxi !=0 )):
        co00i = 0
        co00j = 0
        for j in range(np.size(im00,axis = 1)): 
            for i in range (np.size(im00,axis = 0)):
                if (im00[i,j] == np.max(im00)): 
                    co00i = i
                    co00j = j
        
        Nx00 = Nx -int(Sz/2) + co00j
        Ny00 = Ny -int(Sz/2) + co00i
        
        cxi = Nx00 - int(Sz/2)
        cxs = Nx00 + int(Sz/2)
        
        cyi = Ny00 - int(Sz/2)
        cys = Ny00 + int(Sz/2)
        
        im00 = im1[cyi:cys,cxi:cxs]
        im22 = im2[cyi:cys,cxi:cxs]
        im45 = im3[cyi:cys,cxi:cxs]
        im67 = im4[cyi:cys,cxi:cxs]   
        
        print('cxi,cxs: ',cxi,cxs)
        print('cyi,cys: ',cyi,cys)
        print("\n COORDENADAS CORREGIDAS:\n ", Nx00, Ny00)   
    print('=========================================================')
    print('Dimensiones: ', im00.shape)
    print('Punto máximo: ', im00.max())
    print('=========================================================')
    return im00, im22, im45, im67

def binning(bz, im00, im22, im45, im67,n00,n22,n45,n67, debias = True):
    
    Fzx = int(np.size(im00,axis=1)/bz)
    Fzy = int(np.size(im00,axis=0)/bz)
    
    print('x: ',Fzx,', y: ',Fzy)
    print('=========================================================')
    
    im00[np.isnan(im00)] = 0
    im22[np.isnan(im22)] = 0
    im45[np.isnan(im45)] = 0
    im67[np.isnan(im67)] = 0
    
    im00_minor = im00[im00<0].mean()
    im22_minor = im22[im22<0].mean()
    im45_minor = im45[im45<0].mean()
    im67_minor = im67[im67<0].mean()
    
    if(debias):
        all_minor = np.array([im00_minor, im22_minor, im45_minor, im67_minor])
        print('Minor values: ', all_minor)
        print('Minor value: ', all_minor.min())
        
        im00+=np.abs(all_minor.min())
        im22+=np.abs(all_minor.min())
        im45+=np.abs(all_minor.min())
        im67+=np.abs(all_minor.min())

    
    im00b = np.zeros([Fzy,Fzx])
    im22b = np.zeros([Fzy,Fzx])
    im45b = np.zeros([Fzy,Fzx])
    im67b = np.zeros([Fzy,Fzx])
    
    im00sd = np.zeros([Fzy,Fzx])
    im22sd = np.zeros([Fzy,Fzx])
    im45sd = np.zeros([Fzy,Fzx])
    im67sd = np.zeros([Fzy,Fzx])
    
    
    for i in range(0,Fzy):
        for j in range(0,Fzx): 
            im00b[i,j] = im00[i*bz:(i+1)*bz,j*bz:(j+1)*bz].mean()
            im22b[i,j] = im22[i*bz:(i+1)*bz,j*bz:(j+1)*bz].mean()
            im45b[i,j] = im45[i*bz:(i+1)*bz,j*bz:(j+1)*bz].mean()
            im67b[i,j] = im67[i*bz:(i+1)*bz,j*bz:(j+1)*bz].mean()
            
            
            im00sd[i,j] = im00[i*bz:(i+1)*bz,j*bz:(j+1)*bz].std()
            im22sd[i,j] = im22[i*bz:(i+1)*bz,j*bz:(j+1)*bz].std()
            im45sd[i,j] = im45[i*bz:(i+1)*bz,j*bz:(j+1)*bz].std()
            im67sd[i,j] = im67[i*bz:(i+1)*bz,j*bz:(j+1)*bz].std()
            
    
    print('Promedio 00', im00b.mean())
    print('Promedio 22', im22b.mean())
    print('Promedio 45', im45b.mean())
    print('Promedio 67', im67b.mean())
    print('=========================================================')
    
    print('Promedio 00sd', im00sd.mean())
    print('Promedio 22sd', im22sd.mean())
    print('Promedio 45sd', im45sd.mean())
    print('Promedio 67sd', im67sd.mean())
    print('=========================================================')
    
    im00sd = np.sqrt(np.abs(im00b)/(5.8*n00)) + im00sd.mean()#[np.abs(im00b)<10].mean()
    im22sd = np.sqrt(np.abs(im22b)/(5.8*n22)) + im22sd.mean()#[np.abs(im22b)<10].mean()
    im45sd = np.sqrt(np.abs(im45b)/(5.8*n45)) + im45sd.mean()#[np.abs(im45b)<10].mean()
    im67sd = np.sqrt(np.abs(im67b)/(5.8*n67)) + im67sd.mean()#[np.abs(im67b)<10].mean()
    
    print('Promedio 00sd final', im00sd.mean())
    print('Promedio 22sd final', im22sd.mean())
    print('Promedio 45sd final', im45sd.mean())
    print('Promedio 67sd final', im67sd.mean())
    print('=========================================================')
    print('Dimensiones: ', im00.shape)
    print('=========================================================')
    
    
    return im00b,im22b,im45b,im67b, im00sd,im22sd,im45sd,im67sd,Fzx,Fzy

def std(im00,im22,im45,im67,sig): 
    sky00 = im00[np.abs(im00)<10]
    sky22 = im22[np.abs(im22)<10]
    sky45 = im45[np.abs(im45)<10]
    sky67 = im67[np.abs(im67)<10]
    
    # skymean00 = sky00.mean()
    # skymean22 = sky22.mean()
    # skymean45 = sky45.mean()
    # skymean67 = sky67.mean()
    
    skystd00 = sky00.std()
    skystd22 = sky22.std()
    skystd45 = sky45.std()
    skystd67 = sky67.std()
    
    if (skystd00 >= skystd22 and skystd00 >= skystd45
        and skystd00 >= skystd67): 
        umbral_std = skystd00
    elif (skystd22 >= skystd00 and skystd22 >= skystd45
        and skystd22 >= skystd67): 
        umbral_std = skystd22
    elif (skystd45 >= skystd22 and skystd45 >= skystd00
        and skystd45 >= skystd67): 
        umbral_std = skystd45
    elif (skystd67 >= skystd22 and skystd67 >= skystd45
        and skystd67 >= skystd00): 
        umbral_std = skystd67
        
    umbral_sd = sig#*umbral_std
    
    print('umbral: ', umbral_sd)
    
    im00[np.abs(im00)<umbral_sd] = umbral_std
    im22[np.abs(im00)<umbral_sd] = umbral_std
    im45[np.abs(im00)<umbral_sd] = umbral_std
    im67[np.abs(im00)<umbral_sd] = umbral_std
    
    im00[np.abs(im22)<umbral_sd] = umbral_std
    im22[np.abs(im22)<umbral_sd] = umbral_std
    im45[np.abs(im22)<umbral_sd] = umbral_std
    im67[np.abs(im22)<umbral_sd] = umbral_std
    
    im00[np.abs(im45)<umbral_sd] = umbral_std
    im22[np.abs(im45)<umbral_sd] = umbral_std
    im45[np.abs(im45)<umbral_sd] = umbral_std
    im67[np.abs(im45)<umbral_sd] = umbral_std
    
    im00[np.abs(im67)<umbral_sd] = umbral_std
    im22[np.abs(im67)<umbral_sd] = umbral_std
    im45[np.abs(im67)<umbral_sd] = umbral_std
    im67[np.abs(im67)<umbral_sd] = umbral_std
    
    
    
    return im00,im22,im45,im67

def umbral(bz,im00,im22,im45,im67,im00sd,im22sd,im45sd,im67sd,sig,sigmask
           ,exptime,ZP,
           n1,n2,n3,n4): 
    
        
    umbral_std = 0
    
    snr00 = im00/im00sd
    snr22 = im22/im22sd
    snr45 = im45/im45sd
    snr67 = im67/im67sd
    
    snr00[snr00<0] = 0
    snr22[snr22<0] = 0
    snr45[snr45<0] = 0
    snr67[snr67<0] = 0
    
    print('UMBRAL: ',umbral_std)
    print('=========================================================')
       
    print('SEÑALES A RUIDO PROMEDIO: ')
    print(snr00.mean())
    print(snr22.mean())
    print(snr45.mean())
    print(snr67.mean())
    print('=========================================================')
    
    ind00 = snr00 < sig#*snr00.mean()
    ind22 = snr22 < sig#*snr22.mean()
    ind45 = snr45 < sig#*snr45.mean()
    ind67 = snr67 < sig#*snr67.mean()
    
    zerosQ = ind00 & ind45
    zerosU = ind22 & ind67
    zeros = zerosQ | zerosU
    
    # figure(zerosQ, 'Ubicación zeros Q', z = 0,comap='jet')
    # figure(zerosU, 'Ubicación zeros U', z = 0,comap='jet')
    # figure(zeros, 'Ubicación zeros', z = 0,comap='jet')
    
    if(sigmask == True):
    
        im00[im00<0] = umbral_std
        im22[im22<0] = umbral_std
        im45[im45<0] = umbral_std
        im67[im67<0] = umbral_std
        
        
        im00[zeros] = umbral_std
        im22[zeros] = umbral_std
        im45[zeros] = umbral_std
        im67[zeros] = umbral_std
    
    # ANGULO SOLIDO DEL PIXEL COMPUESTO
    Sz = np.size(im00, axis=0)
    PS = 0.32 #arcosegundos
    PS = (bz*0.32/3600)*(np.pi/180) # Radianes
    dSolid = PS*PS #(1-np.cos(PS))
    print('Angulo solido: ', dSolid, 'str')
    print('=========================================================')
    
    # TRANSFORMANDO DE ADUs a Flujo
    
    # Minst00 = -2.5*np.log10(im00/(exptime*n1))
    # Minst22 = -2.5*np.log10(im22/(exptime*n2))
    # Minst45 = -2.5*np.log10(im45/(exptime*n3))
    # Minst67 = -2.5*np.log10(im67/(exptime*n4))
    
    # Minst00sd = -2.5*np.log10(im00sd/(exptime*n1))
    # Minst22sd = -2.5*np.log10(im22sd/(exptime*n2))
    # Minst45sd = -2.5*np.log10(im45sd/(exptime*n3))
    # Minst67sd = -2.5*np.log10(im67sd/(exptime*n4))
    
    # Mag00 = Minst00 + ZP
    # Mag22 = Minst22 + ZP
    # Mag45 = Minst45 + ZP
    # Mag67 = Minst67 + ZP
    
    # Mag00sd = Minst00sd + ZP
    # Mag22sd = Minst22sd + ZP
    # Mag45sd = Minst45sd + ZP
    # Mag67sd = Minst67sd + ZP
    
    # CONVERSION A FLUJO erg cm⁻² s⁻¹ A⁻¹
    # CUENTAS POR SEGUNDO
    # im00 = (10**((Mag00-ZP)/-2.5))
    # im22 = (10**((Mag22-ZP)/-2.5))
    # im45 = (10**((Mag45-ZP)/-2.5))
    # im67 = (10**((Mag67-ZP)/-2.5))
    
    # im00sd = (10**((Mag00sd-ZP)/-2.5))
    # im22sd = (10**((Mag22sd-ZP)/-2.5))
    # im45sd = (10**((Mag45sd-ZP)/-2.5))
    # im67sd = (10**((Mag67sd-ZP)/-2.5))
    
    im00 = (im00/(exptime*n1))
    im22 = (im22/(exptime*n2))
    im45 = (im45/(exptime*n3))
    im67 = (im67/(exptime*n4))
    
    im00sd = (im00sd/(exptime*n1))
    im22sd = (im22sd/(exptime*n2))
    im45sd = (im45sd/(exptime*n3))
    im67sd = (im67sd/(exptime*n4))
    
    print('=========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'cuentas por segundo')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'cuentas por segundo')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'cuentas por segundo')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'cuentas por segundo')
    
    # ERGIOS POR SEGUNDO
    im00 = im00*((2.998e10)*(6.626e-27))/(1.65e-4)
    im22 = im22*((2.998e10)*(6.626e-27))/(1.65e-4)
    im45 = im45*((2.998e10)*(6.626e-27))/(1.65e-4)
    im67 = im67*((2.998e10)*(6.626e-27))/(1.65e-4)
    
    im00sd = im00sd*((2.998e10)*(6.626e-27))/(1.65e-4)
    im22sd = im22sd*((2.998e10)*(6.626e-27))/(1.65e-4)
    im45sd = im45sd*((2.998e10)*(6.626e-27))/(1.65e-4)
    im67sd = im67sd*((2.998e10)*(6.626e-27))/(1.65e-4)
    
    print('=========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'ergios s⁻¹')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'ergios s⁻¹')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'ergios s⁻¹')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'ergios s⁻¹')
    
    
    # ERGIOS POR SEGUNDO POR CM2
   
    im00 = im00/((Sz*18.5e-4)**2)
    im22 = im22/((Sz*18.5e-4)**2)
    im45 = im45/((Sz*18.5e-4)**2)
    im67 = im67/((Sz*18.5e-4)**2)
    
    im00sd = im00sd/((Sz*18.5e-4)**2)
    im22sd = im22sd/((Sz*18.5e-4)**2)
    im45sd = im45sd/((Sz*18.5e-4)**2)
    im67sd = im67sd/((Sz*18.5e-4)**2)
    print('=========================================================')
    print('Area: ', (Sz*18.5e-4)**2, 'cm⁻²')
    print('=========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'ergios s⁻¹ cm⁻²')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'ergios s⁻¹ cm⁻²')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'ergios s⁻¹ cm⁻²')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'ergios s⁻¹ cm⁻²')
    
    # ERGIOS POR SEGUNDO POR CM2 POR ANGSTROM
    
    L = 16500
    
    im00 = im00 / L
    im22 = im22 / L
    im45 = im45 / L
    im67 = im67 / L
    
    im00sd = im00sd / L
    im22sd = im22sd / L
    im45sd = im45sd / L
    im67sd = im67sd / L
    
    print('=========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'ergios s⁻¹ cm⁻² A⁻¹')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'ergios s⁻¹ cm⁻² A⁻¹')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'ergios s⁻¹ cm⁻² A⁻¹')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'ergios s⁻¹ cm⁻² A⁻¹')
    
    # CONVERSION A JANSKIS POR UNIDAD DE LONGITUD DE ONDA erg cm⁻² s⁻¹ A⁻¹
    # 1 Jy = 1.1011937557392103e-11 erg cm⁻² s⁻¹ A⁻¹
    # im00 = im00/1.1011937557392103e-11 #/dSolid/1000
    # im22 = im22/1.1011937557392103e-11 #/dSolid/1000
    # im45 = im45/1.1011937557392103e-11 #/dSolid/1000
    # im67 = im67/1.1011937557392103e-11 #/dSolid/1000
    
    # im00sd = im00sd/1.1011937557392103e-11 #/dSolid/1000
    # im22sd = im22sd/1.1011937557392103e-11 #/dSolid/1000
    # im45sd = im45sd/1.1011937557392103e-11 #/dSolid/1000
    # im67sd = im67sd/1.1011937557392103e-11 #/dSolid/1000
    
    # CONVERSION A FLUJO POR UNIDAD DE FRECUENCIA 1e-22 erg cm⁻² s⁻¹ Hz⁻¹
    im00 = (((1.65e-4)**2)/(2.998e10))*im00/10e-8
    im22 = (((1.65e-4)**2)/(2.998e10))*im22/10e-8
    im45 = (((1.65e-4)**2)/(2.998e10))*im45/10e-8
    im67 = (((1.65e-4)**2)/(2.998e10))*im67/10e-8
    
    im00sd = (((1.65e-4)**2)/(2.998e10))*im00sd/10e-8
    im22sd = (((1.65e-4)**2)/(2.998e10))*im22sd/10e-8
    im45sd = (((1.65e-4)**2)/(2.998e10))*im45sd/10e-8
    im67sd = (((1.65e-4)**2)/(2.998e10))*im67sd/10e-8
    
    print('==========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'ergios s⁻¹ cm⁻² Hz⁻¹')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'ergios s⁻¹ cm⁻² Hz⁻¹')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'ergios s⁻¹ cm⁻² Hz⁻¹')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'ergios s⁻¹ cm⁻² Hz⁻¹')
    
    im00 = im00/(1e-23)
    im22 = im22/(1e-23)
    im45 = im45/(1e-23)
    im67 = im67/(1e-23)
    
    im00sd = im00sd/(1e-23)
    im22sd = im22sd/(1e-23)
    im45sd = im45sd/(1e-23)
    im67sd = im67sd/(1e-23)
    
    print('=========================================================')
    print('PROMEDIO DE FLUJO 00: ', im00.mean(),'Jy')
    print('PROMEDIO DE FLUJO 22: ', im22.mean(),'Jy')
    print('PROMEDIO DE FLUJO 45: ', im45.mean(),'Jy')
    print('PROMEDIO DE FLUJO 67: ', im67.mean(),'Jy')
    
    im00 = im00/dSolid/1e6
    im22 = im22/dSolid/1e6
    im45 = im45/dSolid/1e6
    im67 = im67/dSolid/1e6
    
    im00sd = im00sd/dSolid/1e6
    im22sd = im22sd/dSolid/1e6
    im45sd = im45sd/dSolid/1e6
    im67sd = im67sd/dSolid/1e6
    print('=========================================================')
    print('PROMEDIO DE BRILLO SUPERFICIAL 00: ', im00.mean(),'MJy/sr')
    print('PROMEDIO DE BRILLO SUPERFICIAL 22: ', im22.mean(),'MJy/sr')
    print('PROMEDIO DE BRILLO SUPERFICIAL 45: ', im45.mean(),'MJy/sr')
    print('PROMEDIO DE BRILLO SUPERFICIAL 67: ', im67.mean(),'MJy/sr')
    print('=========================================================')
    print('UNIDADES EN MJy/sr')
           
    
    return im00, im22, im45, im67, im00sd,im22sd,im45sd,im67sd, snr00,snr22,snr45,snr67

def stokes(im00,im22,im45,im67,im00sd,im22sd,im45sd,im67sd,Fzx,Fzy,bzs,sig,sigmask, 
           doubleI = True):
    
    I=(im00+im22+im45+im67)/2
    Isd = np.sqrt(im00sd**2 + im22sd**2 + im45sd**2 + im67sd**2)/2
    
    # #Procesamiento Parametros U y Q
    if(doubleI): 
        Q=(im00-im45)/(2*I)
        U=(im22-im67)/(2*I)
    else: 
        Q=(im00-im45)/(I)
        U=(im22-im67)/(I)
    
    Qsd = np.sqrt(((im00sd**2+im45sd**2)/I**2)+((Q/I)*Isd)**2)
    Usd = np.sqrt(((im22sd**2+im67sd**2)/I**2)+((U/I)*Isd)**2)
    #Parametros Q y U corregidos pol. instrumental
    #Angulo de offset y eficiencia
    theta = 139*(np.pi)/180
    etha = 0.963
    Qeq = (Q*np.cos(2*theta)-U*np.sin(2*theta))/etha
    Ueq = (U*np.cos(2*theta)+Q*np.sin(2*theta))/etha
    Qeqsd = np.sqrt((Qsd*(np.cos(2*theta))/etha)**2 + (Usd*(np.sin(2*theta))/etha)**2 )
    Ueqsd = np.sqrt((Usd*(np.cos(2*theta))/etha)**2 + (Qsd*(np.sin(2*theta))/etha)**2 )
    
    # Parametros instrumentales
    Qinst = -0.0050
    Uinst = 0.0012
    # PARAMETROS Q y U finales
    Qf = Qeq - Qinst
    Uf = Ueq - Uinst
    
    Qfsd = np.sqrt(Qeqsd**2 + 0.002**2)
    Ufsd = np.sqrt(Ueqsd**2 + 0.002**2)
    
    # LIMPIEZA INICIAL
    # EXCLUYENDO VALORES MAYORES A 1 (NO SENTIDO FÍSICO)
    Qf[np.abs(Qf)>1]=np.nan
    Uf[np.abs(Uf)>1]=np.nan
    # EXCLUYENDO VALORES MAYORES A 6 SIGMA (REMOCIÓN DEL RUIDO EN IMAGEN)
    # Qf[np.abs(Qf)>6*np.nanstd(Qf)]=np.nan
    # Uf[np.abs(Uf)>6*np.nanstd(Uf)]=np.nan
    # EXCLUYENDO VALORES MAYORES A 1 EN MAPAS DE PORPAGACION DE ERRORES
    Qfsd[(Qfsd)>1]=1
    Ufsd[(Ufsd)>1]=1
    
    # ENMASCARAMIENTO POR PROPAGACIÓN DE ERRORES
    if(sigmask):
        # ENMASCARANDO IMAGEN POR VARIACIÓN DE LA PROPAGACIÓN DE ERRORES
        # EXCLUYENDO VALORES CUYA PROPAGACIÓN DE ERRORES SEA MAYOR A UN CIERTO 
        # SIGMA DEFINIDO
        Qf[(Qfsd)>sig]=np.nan
        Uf[(Ufsd)>sig]=np.nan
        # ENMASCARANDO IMAGEN DE SIGMA Q Y U POR VARIACIÓN DE 
        # LA PROPAGACIÓN DE ERRORES
        # EXCLUYENDO VALORES CUYA PROPAGACIÓN DE ERRORES SEA MAYOR A UN CIERTO 
        # SIGMA DEFINIDO
        Qfsd[(Qfsd)>sig]=sig
        Ufsd[(Ufsd)>sig]=sig
    
    
    P=100*(np.sqrt((Qf**2)+(Uf**2)))
    A=0.5*((np.arctan(np.abs(Uf)/np.abs(Qf))))*(180/np.pi)
    # Peq = 100*(np.sqrt((Qeq**2)+(Ueq**2)))
    Psd = ((100**2)/P)*np.sqrt((Qf**2)*(Qfsd**2)+(Uf**2)*(Ufsd**2))
    Asd = np.sqrt((-(Uf*Qfsd)/((Qf**2)+(Uf**2)))**2 + ((Qf*Ufsd)/((Qf**2)+(Uf**2)))**2)
    
    for i in range (Fzy): 
        for j in range (Fzx):
            if(P[i,j]>Psd[i,j]):    
                P[i,j]=np.sqrt(P[i,j]**2 - Psd[i,j]**2)
    P[P>100] = np.nan
    P[P == 0] = np.nan
    # P[np.isnan(P) == True] = 0
            
    # Binning
    Pbin = np.zeros([Fzy,Fzx])
    Abin = np.zeros([Fzy,Fzx])

    for i in range(0,Fzy-bzs,bzs):
        for j in range (0,Fzx-bzs,bzs):
            Pbin[i+bzs][j+bzs] = P[i+bzs][j+bzs]
            Abin[i+bzs][j+bzs] = A[i+bzs][j+bzs] 
            
    # for i in range (Fzy): 
    #     for j in range (Fzx):
    #         if (Qf[i][j] >= 0.0) and (Uf[i][j] >= 0.0): 
    #             A[i][j] = 90 - A[i][j]
    #             Abin[i][j] = 90 - Abin[i][j]
    #         if (Qf[i][j] <= 0.0) and (Uf[i][j] >= 0.0):
    #             A[i][j] = A[i][j]
    #             Abin[i][j] = Abin[i][j]
    #         if (Qf[i][j] <= 0.0) and (Uf[i][j] <= 0.0): 
    #             A[i][j] = 180 - A[i][j]
    #             Abin[i][j] = 180 - Abin[i][j]
    #         if (Qf[i][j] >= 0.0) and (Uf[i][j] <= 0.0):
    #             A[i][j] = 90 + A[i][j]
    #             Abin[i][j] = 90 + Abin[i][j]
                
    print('ESQUEMA PROPIO')
    for i in range (Fzy): 
        for j in range (Fzx):
            if (Qf[i,j] >= 0.0) and (Uf[i,j] >= 0.0): 
                A[i][j] = 90 + A[i][j]
                Abin[i][j] = 90 + Abin[i][j]
            if (Qf[i,j] <= 0.0) and (Uf[i,j] >= 0.0):
                A[i][j] = 180 - A[i][j]
                Abin[i][j] = 180 - Abin[i][j]
            if (Qf[i,j] <= 0.0) and (Uf[i,j] <= 0.0): 
                A[i][j] = A[i][j]
                Abin[i][j] = Abin[i][j]
            if (Qf[i,j] >= 0.0) and (Uf[i,j] <= 0.0):
                A[i][j] = 90 - A[i][j]
                Abin[i][j] = 90 - Abin[i][j]
                
    # for i in range (Fzy): 
    #     for j in range (Fzx):
    #         if(A[i,j]>=90):
    #             A[i,j]=A[i,j]-90
    #             Abin[i,j]=Abin[i,j]-90
    #         elif(A[i,j]<90):
    #             Abin[i,j] = Abin[i,j]+90
                
                
    return I,Isd,Qf,Qfsd,Uf,Ufsd,P,Psd,Pbin,A,Asd,Abin

def graf_vector(obj,Sz, Sz2, bz, mask, If, P,Pbin, A, Abin,RA=0,DEC=0,save=False):
    
    Fz = int(Sz/bz)
    Fz2 = int(Sz2/bz)
    PSoagh = 0.32
    arcsec2 = Sz2*PSoagh/2
    Lz2 = np.linspace(-arcsec2, arcsec2, Fz2)
    X3,Y3 = np.meshgrid(Lz2, Lz2)
    
    Iff = np.zeros([Fz2,Fz2])
    Pf = np.zeros([Fz2,Fz2])
    Af = np.zeros([Fz2,Fz2])
    Pv = np.zeros([Fz2,Fz2])
    Av = np.zeros([Fz2,Fz2])
    #Valores de los vectores
    xv=np.zeros([Fz2,Fz2])
    yv=np.zeros([Fz2,Fz2])
    # for j in range(Fz2): 
    #     for i in range (Fz2):
    #         l = i + int(Fz/2)-int(Fz2/2)
    #         k = j + int(Fz/2)-int(Fz2/2)

    #         if (mask == True):
    #             Pv[i,j]=Pbin[l,k]
    #             Av[i,j]=Abin[l,k]
    #         else: 
    #             Pv[i,j]=P[l,k]
    #             Av[i,j]=A[l,k]
            
    #         # if (Pv[i][j]!= 0): 
    #         #     Pv[i][j] = 6
    #         # if (Pv[i][j]<5 and Pv[i][j]>0): 
    #             # Pv[i][j] = 10
    #         # elif (Pv[i][j]<10 and Pv[i][j]>=5):
    #         #     Pv[i][j] = 10
                
    #         xv[i][j]=Pv[i][j]*np.cos(Av[i][j]*np.pi/180)
    #         yv[i][j]=Pv[i][j]*np.sin(Av[i][j]*np.pi/180)
    #         Pf[i][j]=P[l][k]
    #         Af[i][j]=A[l][k]
    #         Iff[i][j]= If[l][k]
    
    
    print('=========================================================')
    print("\n Polarización Máxima:\n", np.nanmax(P))
    Pmax = int(input('Introduzca escala: '))
    print('=========================================================')
    
    if (Pmax<=10 and Pmax >0): 
        Pc = int(10*(Pmax/10))
    elif(Pmax>10):
        Pc = (10*int(Pmax/10))
    else:
        Pc = 30
    print("\n Polarización escala:\n", Pc)
    print('=========================================================')
    sc = 10*(Pc)
    # if (sc < 300): sc = sc + 100
    if (sc < 200): sc = sc + 20
    
    for j in range(Fz2): 
        for i in range (Fz2):
            l = i + int(Fz/2)-int(Fz2/2)
            k = j + int(Fz/2)-int(Fz2/2)

            if (mask == True):
                Pv[i,j]=Pbin[l,k]
                Av[i,j]=Abin[l,k]
            else: 
                Pv[i,j]=P[l,k]
                Av[i,j]=A[l,k]
            
            # if (Pv[i][j]!= 0): 
            #     Pv[i][j] = 6
            # if (Pv[i][j]<5 and Pv[i][j]>0): 
                # Pv[i][j] = 10
            # elif (Pv[i][j]<10 and Pv[i][j]>=5):
            #     Pv[i][j] = 10
            if (Pmax == 0 and Pv[i][j] > 0):
                xv[i][j]=10*np.cos(Av[i][j]*np.pi/180)
                yv[i][j]=10*np.sin(Av[i][j]*np.pi/180)
            else: 
                xv[i][j]=Pv[i][j]*np.cos(Av[i][j]*np.pi/180)
                yv[i][j]=Pv[i][j]*np.sin(Av[i][j]*np.pi/180)
            Pf[i][j]=P[l][k]
            Af[i][j]=A[l][k]
            Iff[i][j]= If[l][k]

    # plt.style.use('dark_background')#vmin=np.min(Iff),
    print('ESCALA VECTORES: ', sc)
    Iff_imagen = np.copy(Iff)
    
    fig, ax0 = plt.subplots(figsize=(10,10))
    ax0.set_title(obj + ' (I) + POLARIZACIÓN')
    im = ax0.imshow(Iff_imagen, cmap='Greys_r',origin = 'lower',
                    vmax=3*Iff.std(), vmin=0, 
                    extent = [-arcsec2, arcsec2, -arcsec2, arcsec2])#,norm=colors.LogNorm(vmax=np.max(Iff)))
                    # vmax=2*Iff.std(), vmin=10)#-Iff.std())
    # im = ax0.pcolor(X3,Y3,Iff, vmax = np.max(Iff/4), cmap='magma')
    ax0.grid(False)
    cb = fig.colorbar(im, ax = ax0)
    cb.set_label('I (MJy/sr)', rotation = 270)
    ax0.quiver(X3,Y3,xv,yv, scale=sc, color='red', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
    ax0.set_xlabel('arcsec')
    ax0.set_ylabel('arcsec')
    if (Pmax != 0):
        ax0.text(arcsec2 - (arcsec2/4),-arcsec2 + (arcsec2/6), '%i'% Pc +'%',color = 'red', size = 15)
        ax0.axhline(y = -arcsec2 + (arcsec2/10), xmin = 0.86, xmax = 0.94, color = 'red', linestyle = '-')
    if (save == True):
        plt.savefig('Stokes_I.png', dpi=150)
        
    #Grafica de vectores, grado de polarización
    fig, ax1 = plt.subplots(figsize=(10,10))
    ax1.set_title(obj + ' POLARIZACIÓN')
    # im1 = ax1.pcolor(X3,Y3,Pf, cmap='jet', norm=colors.LogNorm(vmax=np.max(Pf)))
    im1 = ax1.imshow(Pf, cmap='jet', origin='lower',
                     extent = [-arcsec2, arcsec2, -arcsec2, arcsec2])#,vmax = Pmax,vmin = 0)
    ax1.grid(False)
    cb = plt.colorbar(im1, ax = ax1)
    cb.set_label('%P', rotation = 270)
    ax1.quiver(X3,Y3,xv,yv, scale=sc, color='black', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
    ax1.set_xlabel('arcsec')
    ax1.set_ylabel('arcsec')
    if (Pmax != 0):
        ax1.text(arcsec2 -  (arcsec2/4), -arcsec2 + (arcsec2/6), '%i'% Pc +'%',color = 'black', size = 15)
        ax1.axhline(y = -arcsec2 + (arcsec2/10), xmin = 0.86, xmax = 0.94, color = 'black', linestyle = '-')
    if (save == True):
        plt.savefig('Polarizacion.png', dpi=150)
    
    #Grafica de vectores, grado de polarización
    fig, ax1 = plt.subplots(figsize=(10,10))
    ax1.set_title(obj + ' POLARIZACIÓN')
    # im1 = ax1.pcolor(X3,Y3,Pf, cmap='jet', norm=colors.LogNorm(vmax=np.max(Pf)))
    im1 = ax1.imshow(Af, cmap='jet',origin='lower',vmax = 180,vmin = 0,
                     extent = [-arcsec2, arcsec2, -arcsec2, arcsec2])
    ax1.grid(False)
    cb = plt.colorbar(im1, ax = ax1)
    cb.set_label('A°', rotation = 270)
    # ax1.quiver(X3,Y3,xv,yv, scale=sc, color='black', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
    ax1.set_xlabel('arcsec')
    ax1.set_ylabel('arcsec')
    # ax1.text(arcsec2 -  (arcsec2/4), -arcsec2 + (arcsec2/6), '%i'% Pc +'%',color = 'black', size = 15)
    # ax1.axhline(y = -arcsec2 + (arcsec2/10), xmin = 0.86, xmax = 0.94, color = 'black', linestyle = '-')
    if (save == True):
        plt.savefig('Angulo.png', dpi=150)
    # fig, ax1 = plt.subplots()#figsize=(10,10))
    # ax1.set_title(obj + ' POLARIZACIÓN')
    # # im1 = ax1.pcolor(X3,Y3,Pf, cmap='jet', norm=colors.LogNorm(vmax=np.max(Pf)))
    # im1 = ax1.contourf(X3,Y3,Pf, cmap='jet')
    # ax1.grid(False)
    # cb = plt.colorbar(im1, ax = ax1)
    # cb.set_label('%P', rotation = 270)
    # # ax1.quiver(X3,Y3,xv,yv, scale=sc, color='black', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
    # ax1.set_xlabel('arcsec')
    # ax1.set_ylabel('arcsec')
    # # ax1.text(arcsec2 -  (arcsec2/4), -arcsec2 + (arcsec2/6), '%i'% Pc +'%',color = 'black', size = 15)
    # # ax1.axhline(y = -arcsec2 + (arcsec2/10), xmin = 0.86, xmax = 0.94, color = 'black', linestyle = '-')
    
    # fig, ax1 = plt.subplots()#figsize=(10,10))
    # ax1.set_title(obj + ' POLARIZACIÓN')
    # # im1 = ax1.pcolor(X3,Y3,Pf, cmap='jet', norm=colors.LogNorm(vmax=np.max(Pf)))
    # im1 = ax1.contour(X3,Y3,Pf,2, cmap='jet')
    # ax1.grid(False)
    # cb = plt.colorbar(im1, ax = ax1)
    # cb.set_label('%P', rotation = 270)
    # # ax1.quiver(X3,Y3,xv,yv, scale=sc, color='black', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
    # ax1.set_xlabel('arcsec')
    # ax1.set_ylabel('arcsec')
    
    return X3,Y3,xv,yv, Pf, Af, sc

def photometry(Aperture,bz,r,a,b,t,im,exptime,ZP,n): 
    
    x1,y1=centroid_com(im)
    Sz = np.size(im, axis = 0)
    # DEFINIENDO LAS APERTURAS 
    if (Aperture == 1):
        aperture1 = CircularAperture((x1,y1), r)

    elif (Aperture == 2):
        aperture1 = EllipticalAperture((x1,y1), a, b, theta=t)
    
    # ANGULO SOLIDO DE LA APERTURA
    PS = 0.32 #arcosegundos
    print('=========================================================')
    print('Radio de la apertura: ', r, 'pixeles')
    print('Radio de la apertura: ', PS, 'arcsec')
    print('Diametro de la apertura: ', 2*PS, 'arcsec')
    print('=========================================================')
    dS = (bz*0.32/3600)*(np.pi/180) # Radianes (bz*0.32/3600)*(np.pi/180)
    dSolid = PS*PS #2*np.pi*(1-np.cos(PS))
    print('=========================================================')
    print('Angulo solido de la apertura: ', dSolid, 'str')
    print('=========================================================')
    # ESTADISTICAS DE LA APERTURAS
    aperstats = ApertureStats(im, aperture1) 
    ap_mean = aperstats.mean
    ap_sig = aperstats.std
    ap_sum = aperstats.sum
    phot_table = aperture_photometry(im, aperture1)
    print(phot_table)
    print('=========================================================')
    print('SUMA APERTURA',ap_sum/(exptime*n))#*dSolid)
    print('=========================================================')
    # MAGNITUD INSTRUMENTAL
    
    Minst = -2.5*np.log10(ap_sum/(exptime*n))
    Mag = Minst + ZP
    F = 10**((Mag-ZP)/-2.5)
    F = ap_sum/(exptime*n)
    print('FLUJO: ', F, 'cuentas por segundo')
    print('=========================================================')
    
    F = 10**((Mag-ZP)/-2.5)*((2.998e10)*(6.626e-27))/(1.65e-4)
    print('FLUJO: ', F, 'ergios s⁻¹')
    print('=========================================================')
    
    A = (Sz*18.5e-4)**2#np.pi*((18.5e-4)*r*bz)**2
    print('Area: ', A, 'cm²')
    print('=========================================================')
    F = F/A
    print('FLUJO: ', F, 'ergios s⁻¹ cm⁻²')
    print('=========================================================')
    
    F = F/(16500)
    print('FLUJO: ', F, 'ergios s⁻¹ cm⁻² A⁻¹')
    print('=========================================================')
    
    F = (((1.65e-4)**2)/(2.998e10))*F/(10e-8)
    print('MAGNITUD INSTRUMENTAL: ', Minst)
    print('MAGNITUD CORREGIDA: ', Mag)
    print('FLUJO: ', F, "erg cm⁻² s⁻¹ Hz⁻¹")
    print('=========================================================')
    
    # CONVERSION A JANSKYS
    F = F/(1e-23)
    print('FLUJO: ', F, "Jy")
    print('=========================================================')
    
    plt.figure()
    plt.title("Fotometría",fontsize = 10)
    plt.imshow(im, origin = 'lower',cmap = 'inferno',vmax=3*im.std(), vmin=10)
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 10)
    aperture1.plot(color='white', lw=2, label='Photometry aperture')
    plt.xlabel('pixeles',fontsize = 10)
    plt.ylabel('pixeles',fontsize = 10)
    # plt.legend(labelcolor = 'white')
    plt.text(x = 1, y = 1,s = 'F = ' + str(round(F,3)) + ' Jy'
             ,color = 'white', size = 15)
    
    return F

def validation_photometry(Aperture,bz,r,a,b,t,im,exptime,ZP,n): 
     
    x1,y1=centroid_com(im)
     
    # DEFINIENDO LAS APERTURAS 
    if (Aperture == 1):
        aperture1 = CircularAperture((x1,y1), r)
    elif (Aperture == 2):
        aperture1 = EllipticalAperture((x1,y1), a, b, theta=t)
    
    # ANGULO SOLIDO DE LA APERTURA
    Sz = np.size(im, axis = 0)
    PS = bz*0.32 #arcosegundos
    print('=========================================================')
    print('Radio de la apertura: ', r, 'pixeles')
    print('Radio de la apertura: ', PS/2, 'arcsec')
    print('Diametro de la apertura: ', PS, 'arcsec')
    print('=========================================================')
    PS = (PS/3600)*(np.pi/180) # Radianes
    dSolid = PS*PS #2*np.pi*(1-np.cos(PS))
    print('=========================================================')
    print('Angulo solido de la apertura: ', dSolid, 'sr')
    print('=========================================================')
    # ESTADISTICAS DE LA APERTURAS
    aperstats = ApertureStats(im, aperture1) 
    ap_mean = aperstats.mean
    ap_sig = aperstats.std
    ap_sum = aperstats.sum
    phot_table = aperture_photometry(im, aperture1)
    print(phot_table)
    print('=========================================================')
    print('SUMA APERTURA',ap_sum,'MJy/sr')#*dSolid)
    print('=========================================================')
    F = ap_sum * dSolid * 1e6
    print('FLUJO: ', F, "Jy")
    print('=========================================================')
    
    plt.figure()
    plt.title("Fotometría",fontsize = 10)
    plt.imshow(im, origin = 'lower',cmap = 'inferno',vmax=3*im.std(), vmin=0)
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('MJy/sr', rotation = 270,fontsize = 10)
    aperture1.plot(color='white', lw=2, label='Photometry aperture')
    plt.xlabel('pixeles',fontsize = 10)
    plt.ylabel('pixeles',fontsize = 10)
    # plt.legend(labelcolor = 'white')
    plt.text(x = 1, y = 1,s = 'F = ' + str(round(F,3)) + ' Jy'
             ,color = 'white', size = 15)
    
    return F

def photometry_stokes(Aperture,r,a,b,t,rin, rout, im00, im22, im45, im67,archivo_csv,objeto,filtro): 
    
    x1, y1 = centroid_com(im00)
    x2, y2 = centroid_com(im22)
    x3, y3 = centroid_com(im45)
    x4, y4 = centroid_com(im67)
    # positions = [(x1,y1), (x2,y2),(x3,y3), (x4,y4)]
    
    # DEFINIENDO LAS APERTURAS 
    if (Aperture == 1):
        aperture1 = CircularAperture((x1,y1), r)
        aperture2 = CircularAperture((x2,y2), r)
        aperture3 = CircularAperture((x3,y3), r)
        aperture4 = CircularAperture((x4,y4), r)
    elif (Aperture == 2):
        aperture1 = EllipticalAperture((x1,y1), a, b, theta=t)
        aperture2 = EllipticalAperture((x2,y2), a, b, theta=t)
        aperture3 = EllipticalAperture((x3,y3), a, b, theta=t)
        aperture4 = EllipticalAperture((x4,y4), a, b, theta=t)
    
    # ESTADISTICAS DE LAS APERTURAS
    aperstats00 = ApertureStats(im00, aperture1) 
    aperstats22 = ApertureStats(im22, aperture2) 
    aperstats45 = ApertureStats(im45, aperture3) 
    aperstats67 = ApertureStats(im67, aperture4) 
    
    # ap_mean00 = aperstats00.mean
    # ap_mean22 = aperstats22.mean
    # ap_mean45 = aperstats45.mean
    # ap_mean67 = aperstats67.mean
    
    ap_sig00 = aperstats00.std
    ap_sig22 = aperstats22.std
    ap_sig45 = aperstats45.std
    ap_sig67 = aperstats67.std
    
    # DEFINIENDO EL ANNULUS
    annulus_aperture1 = CircularAnnulus((x1,y1), r_in=rin, r_out=rout)
    annulus_aperture2 = CircularAnnulus((x2,y2), r_in=rin, r_out=rout)
    annulus_aperture3 = CircularAnnulus((x3,y3), r_in=rin, r_out=rout)
    annulus_aperture4 = CircularAnnulus((x4,y4), r_in=rin, r_out=rout)
    
    # ESTADISTICAS DEL ANNULUS
    aperstats1 = ApertureStats(im00, annulus_aperture1) 
    aperstats2 = ApertureStats(im22, annulus_aperture2) 
    aperstats3 = ApertureStats(im45, annulus_aperture3) 
    aperstats4 = ApertureStats(im67, annulus_aperture4) 
    
    bkg_mean1 = aperstats1.mean
    bkg_mean2 = aperstats2.mean
    bkg_mean3 = aperstats3.mean
    bkg_mean4 = aperstats4.mean
    
    bkg_sig1 = aperstats1.std
    bkg_sig2 = aperstats2.std
    bkg_sig3 = aperstats3.std
    bkg_sig4 = aperstats4.std
    
    # AREA DEL APERTURA
    aperture_area1 = aperture1.area_overlap(im00)
    aperture_area2 = aperture2.area_overlap(im22)
    aperture_area3 = aperture3.area_overlap(im45)
    aperture_area4 = aperture4.area_overlap(im67)
    
    total_bkg1 = bkg_mean1 * aperture_area1
    total_bkg2 = bkg_mean2 * aperture_area2
    total_bkg3 = bkg_mean3 * aperture_area3
    total_bkg4 = bkg_mean4 * aperture_area4
    
    total_bkg_sig1 = bkg_sig1 * aperture_area1
    total_bkg_sig2 = bkg_sig2 * aperture_area2
    total_bkg_sig3 = bkg_sig3 * aperture_area3
    total_bkg_sig4 = bkg_sig4 * aperture_area4
    

    phot_table1 = aperture_photometry(im00, aperture1)
    phot_table2 = aperture_photometry(im22, aperture2)
    phot_table3 = aperture_photometry(im45, aperture3)
    phot_table4 = aperture_photometry(im67, aperture4)
    
    # print(phot_table1['aperture_sum'])
    # print(phot_table2['aperture_sum'])
    # print(phot_table3['aperture_sum'])
    # print(phot_table4['aperture_sum'])
    
    # print(total_bkg1)
    # print(total_bkg2)
    # print(total_bkg3)
    # print(total_bkg4)
    
    a00 = (phot_table1['aperture_sum'])# - total_bkg1)
    a22 = (phot_table2['aperture_sum'])# - total_bkg2)
    a45 = (phot_table3['aperture_sum'])# - total_bkg3)
    a67 = (phot_table4['aperture_sum'])# - total_bkg4)
    
    # a00 = (ap_mean00 - total_bkg1)
    # a22 = (ap_mean22 - total_bkg2)
    # a45 = (ap_mean45 - total_bkg3)
    # a67 = (ap_mean67 - total_bkg4)
    
    a00sd = ap_sig00 - total_bkg_sig1
    a22sd = ap_sig22 - total_bkg_sig2
    a45sd = ap_sig45 - total_bkg_sig3
    a67sd = ap_sig67 - total_bkg_sig4
    
    # MAGNITUD
    
    print(a00)
    print(a22)
    print(a45)
    print(a67)
    
    I=(a00+a22+a45+a67)/2
    Isd = np.sqrt(a00sd**2 + a22sd**2 + a45sd**2 + a67sd**2)/2
    
    # #Procesamiento Parametros U y Q
    Q=(a00-a45)/(I)
    U=(a22-a67)/(I)
    
    # Parametros instrumentales
    Qinst = -0.0050
    Uinst = 0.0012
    
    Qsd = np.sqrt(((a00sd**2+a45sd**2)/I**2)+((Q/I)*Isd)**2)
    Usd = np.sqrt(((a22sd**2+a67sd**2)/I**2)+((U/I)*Isd)**2)
    
    #Parametros Q y U corregidos pol. instrumental
    #Angulo de offset y eficiencia
    theta = 139*(np.pi)/180
    etha = 0.963
    Qeq = (Q*np.cos(2*theta)-U*np.sin(2*theta))/etha
    Ueq = (U*np.cos(2*theta)-Q*np.sin(2*theta))/etha
    Qeqsd = np.sqrt((Qsd*(np.cos(2*theta))/etha)**2 + (Usd*(np.sin(2*theta))/etha)**2 )
    Ueqsd = np.sqrt((Usd*(np.cos(2*theta))/etha)**2 + (Qsd*(np.sin(2*theta))/etha)**2 )
    Qf = Qeq - Qinst
    Uf = Ueq - Uinst
    Qsd = np.sqrt(Qeqsd**2 + 0.002**2)
    Usd = np.sqrt(Ueqsd**2 + 0.002**2)
       
    
    P = 100*np.sqrt(Qf**2 + Uf**2)
    if (P < 0): P = 0
    Ain = (0.5*(np.arctan(abs(Uf)/abs(Qf))))*(180/np.pi)
    
    Peq=100*(np.sqrt((Qeq**2)+(Ueq**2)))
    if (Peq != 0):
        Psd = (100/Peq)*np.sqrt((Qf**2)*(Qsd**2)+(Uf**2)*(Usd**2))
    else:
        Psd = 0
        
    if (Psd!=0):
        if ((P/Psd)>=1.0): 
            P= np.sqrt(P**2 - Psd**2)
        else:
            P = 0
        
    print("Q: ",Q)
    print("U: ",U)
    
    if (Qf >= 0.0) and (Uf >= 0.0): 
        A = 90 - Ain
    elif (Qf <= 0.0) and (Uf >= 0.0):
        A = Ain
    elif (Qf <= 0.0) and (Uf <= 0.0): 
        A = 180 - Ain
    elif (Qf >= 0.0) and (Uf <= 0.0):
        A = 90 + Ain
    
    if (A >= 90):
        Adata = A - 90
    elif(A < 90):
        Adata = A + 90
    
        
    plt.figure()
    plt.title("00",fontsize = 20)
    plt.imshow(im00, cmap = 'inferno', vmax=im00.flatten().std(), vmin=-im00.flatten().std())
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
    aperture1.plot(color='white', lw=2, label='Photometry aperture')
    annulus_aperture1.plot(color='red', lw=2,label='Background annulus')
    plt.xlabel('pixeles',fontsize = 20)
    plt.ylabel('pixeles',fontsize = 20)
    plt.legend(labelcolor = 'white')
    
    plt.show()

    plt.figure()
    plt.title("22",fontsize = 20)
    plt.imshow(im22, cmap = 'inferno', vmax=im22.flatten().std(), vmin=-im22.flatten().std())
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
    aperture2.plot(color='white', lw=2, label='Photometry aperture')
    annulus_aperture2.plot(color='red', lw=2,label='Background annulus')
    plt.xlabel('pixeles',fontsize = 20)
    plt.ylabel('pixeles',fontsize = 20)
    plt.legend(labelcolor = 'white')
    
    plt.show()

    plt.figure()
    plt.title("45",fontsize = 20)
    plt.imshow(im45, cmap = 'inferno', vmax=im45.flatten().std(), vmin=-im45.flatten().std())
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
    aperture3.plot(color='white', lw=2, label='Photometry aperture')
    annulus_aperture3.plot(color='red', lw=2,label='Background annulus')
    plt.xlabel('pixeles',fontsize = 20)
    plt.ylabel('pixeles',fontsize = 20)
    plt.legend(labelcolor = 'white')
    
    plt.show()

    plt.figure()
    plt.title("67",fontsize = 20)
    plt.imshow(im67, cmap = 'inferno', vmax=im67.flatten().std(), vmin=-im67.flatten().std())
    plt.grid(False)
    cb = plt.colorbar()
    cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
    aperture4.plot(color='white', lw=2, label='Photometry aperture')
    annulus_aperture4.plot(color='red', lw=2,label='Background annulus')
    plt.xlabel('pixeles',fontsize = 20)
    plt.ylabel('pixeles',fontsize = 20)
    plt.legend(labelcolor = 'white')
    
    
    print("Parametros de Stokes")
    print("I: ",I)
    print("Q: ",Qf)
    print("U: ",Uf)
    print("P: ",P)
    print("Adata: ",Adata)
    
    datos = {'Objeto': objeto, 'filtro':filtro,
             '00':a00,'22':a22,'45':a45,'67':a67,
             'I':I,'Q':Q,'U':U,
             'P':P,'A':Adata}
    
    data = pd.DataFrame(datos, index = [0])
    
    
    for root, dirs, files in os.walk('.'):
        if archivo_csv in files:
            find = True
            break
            # return os.path.join(root, archivo_cvs)
        else:
            find = False
            # return 0
    if(find):
        poldata = pd.read_csv(archivo_csv)
        poldatan = pd.concat([poldata,data],ignore_index=True)
        poldatan.head()
        poldatan.to_csv(archivo_csv, index = False)
        print('ACTUALIZANDO')
    else:
        data.to_csv(archivo_csv)
        print('GENERANDO NUEVO')



#CODIGO RECICLADO

# if (im00[im00<10].std() >= im22[im22<10].std() 
#     and im00[im00<10].std() >= im45[im45<10].std()
#     and im00[im00<10].std() >= im67[im67<10].std()):
#     umbral_std = im00[im00<10].std()
    
# elif (im22[im22<10].std() >= im00[im00<10].std()
#     and im22[im22<10].std() >= im45[im45<10].std()
#     and im22[im22<10].std() >= im67[im67<10].std()):
#     umbral_std = im22[im22<10].std()
    
# elif (im45[im45<10].std() >= im22[im22<10].std() 
#     and im45[im45<10].std() >= im00[im00<10].std()
#     and im45[im45<10].std() >= im67[im67<10].std()):
#     umbral_std = im45[im45<10].std()
    
# elif (im67[im67<10].std() >= im22[im22<10].std() 
#     and im67[im67<10].std() >= im45[im45<10].std()
#     and im67[im67<10].std() >= im00[im00<10].std()):
#     umbral_std = im67[im67<10].std()

 # im00[im00/im00sd < sig*snr00.mean()] = umbral_std
 # im22[im22/im22sd < sig*snr22.mean()] = umbral_std
 # im45[im45/im45sd < sig*snr45.mean()] = umbral_std
 # im67[im67/im67sd < sig*snr67.mean()] = umbral_std
 
 # if (im00[im00<10].std() <= im22[im22<10].std() 
 #     and im00[im00<10].std() <= im45[im45<10].std()
 #     and im00[im00<10].std() <= im67[im67<10].std()):
 #     umbral_std = im00[im00<10].std()
     
 # elif (im22[im22<10].std() <= im00[im00<10].std()
 #     and im22[im22<10].std() <= im45[im45<10].std()
 #     and im22[im22<10].std() <= im67[im67<10].std()):
 #     umbral_std = im22[im22<10].std()
     
 # elif (im45[im45<10].std() <= im22[im22<10].std() 
 #     and im45[im45<10].std() <= im00[im00<10].std()
 #     and im45[im45<10].std() <= im67[im67<10].std()):
 #     umbral_std = im45[im45<10].std()
     
 # elif (im67[im67<10].std() <= im22[im22<10].std() 
 #     and im67[im67<10].std() <= im45[im45<10].std()
 #     and im67[im67<10].std() <= im00[im00<10].std()):
 #     umbral_std = im67[im67<10].std()
 
 # Estadistica de Parámetros de Stokes
 # Qq1, Qq3 = np.percentile(Qf.flatten(),[25,75])  
 # Uq1, Uq3 = np.percentile(Uf.flatten(),[25,75])
 # Qiqd = Qq3 - Qq1
 # Uiqd = Uq3 - Uq1
 
 # print('Qq1, Qq3: ', Qq1,Qq3)
 # print('Uq1, Uq3: ', Uq1,Uq3)
 # print('Qiqd: ', Qiqd)
 # print('Uiqd: ', Uiqd)
 
 # Qli = Qq1 - 2.5*Qiqd
 # Qls = Qq3 + 2.5*Qiqd
 
 # Uli = Uq1 - 2.5*Uiqd
 # Uls = Uq3 + 2.5*Uiqd
 
 # print('Qli, Qls: ', Qli,Qls)
 # print('Uli, Qls: ', Uli,Qls)
 
 # Qf[Qf<Qli] = np.nan
 # Uf[Uf<Uli] = np.nan
 
 # Qf[Qf>Qls] = np.nan
 # Uf[Uf>Uls] = np.nan
 
