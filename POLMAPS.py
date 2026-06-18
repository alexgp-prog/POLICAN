import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS as WCS_astropy
from astropy.coordinates import SkyCoord
from astropy.wcs.wcsapi import SlicedLowLevelWCS
import numpy as np
#from mpdaf.obj import Image, WCS

name = 'M82/HM82H_POLICAN_NEW.fits'
obj = 'M82 H ($1.65 \mu m$)'
ra_obj = [9,55,52.430]
dec_obj = [69,40,46.93]

wsz = 120
snr_thres = 0.04
sep_vector = 3

def coords_hours_to_deg(ra_in, dec_in):

    print('RA_IN: ', ra_in)
    print('DEC_IN: ', dec_in)
    ra_hours = ra_in[0]
    ra_min = ra_in[1]
    ra_sec = ra_in[2]

    dec_deg = dec_in[0]
    dec_min = dec_in[1]
    dec_sec = dec_in[2]

    ra_f = (ra_hours + ra_min*(60/3600) + ra_sec/3600)*(360/24)
    dec_f = np.abs(dec_deg) + dec_min*(60/3600) + dec_sec/3600

    if dec_deg<1:
        dec_f *=-1

    print('RA_OUT: ', ra_f)
    print('DEC_OUT: ', dec_f)
    return ra_f, dec_f



ra, dec = coords_hours_to_deg(ra_obj, dec_obj)

hdul = fits.open(name)
hdul.info()

# Cargando los datos a cada arreglo
head_polican = hdul[0].header
headI = hdul[1].header
dataI = hdul[1].data
sigI = hdul[2].data
dataQ = hdul[3].data
sigQ = hdul[4].data
dataU = hdul[5].data
sigU = hdul[6].data
P = hdul[7].data
sigP = hdul[8].data
A = hdul[9].data
sigA = hdul[10].data

wcs = WCS_astropy(headI)
print(wcs)

coords = SkyCoord(ra = ra, dec = dec, unit = 'deg')
print(coords)
### Conversion de coordenadas a pixeles para cada uno de los casos. 
wr = wcs.world_to_pixel(coords)
print('===============================')
print(wr)
print('===============================')

I_visual = np.copy(dataI[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])
Q_visual = np.copy(dataQ[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])
U_visual = np.copy(dataU[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])

sigQ_visual = np.copy(sigQ[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])
sigU_visual = np.copy(sigU[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])
P_visual = np.copy(P[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])
A_visual = np.copy(A[int(wr[1]) - int(wsz/2): int(wr[1]) + int(wsz/2), int(wr[0]) - int(wsz/2): int(wr[0]) + int(wsz/2)])

P_visual[sigQ_visual > snr_thres] = 0
P_visual[sigU_visual > snr_thres] = 0
P_visual[np.isnan(P_visual)] = 0
A_visual[np.isnan(A_visual)] = 0

## Corrected coordinates

slices = [slice(int(wr[1]) - int(wsz/2), int(wr[1]) + int(wsz/2)), 
          slice(int(wr[0]) - int(wsz/2), int(wr[0]) + int(wsz/2))]
subwcs = SlicedLowLevelWCS(wcs, slices=slices)

### PARAMETROS DE STOKES
fig = plt.figure()
fig.suptitle(obj + 'Parámetros de Stokes')
ax0 = plt.subplot(131, projection=subwcs, slices=('x', 'y'))
ax0.set_title('Stokes I')
im = ax0.imshow(I_visual, origin = 'lower')
cb = fig.colorbar(im, ax = ax0, location = 'bottom', cmap = 'jet')
cb.set_label('MJy/sr')
ax0.set_xlabel('RA (J2000)')
ax0.set_ylabel('DEC (J2000)')

ax1 = plt.subplot(132, projection=subwcs, slices=('x', 'y'))
ax1.set_title('Stokes Q')
im = ax1.imshow(Q_visual * I_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax1, location = 'bottom')
cb.set_label('MJy/sr')
ax1.set_xlabel('RA (J2000)')
ax1.set_ylabel('DEC (J2000)')

ax2 = plt.subplot(133, projection=subwcs, slices=('x', 'y'))
ax2.set_title('Stokes U')
im = ax2.imshow(U_visual * I_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax2, location = 'bottom')
cb.set_label('MJy/sr')
ax2.set_xlabel('RA (J2000)')
ax2.set_ylabel('DEC (J2000)')

### SIGMA PARAMETROS DE STOKES
fig = plt.figure()
ax1 = plt.subplot(121, projection=subwcs, slices=('x', 'y'))
im = ax1.imshow(sigQ_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax1, location = 'bottom')
cb.set_label('%')
ax1.set_xlabel('RA (J2000)')
ax1.set_ylabel('DEC (J2000)')

ax2 = plt.subplot(122, projection=subwcs, slices=('x', 'y'))
im = ax2.imshow(sigU_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax2, location = 'bottom')
cb.set_label('%')
ax2.set_xlabel('RA (J2000)')
ax2.set_ylabel('DEC (J2000)')

### GRADO Y ÁNGULO DE POLARIZACIÓN
fig = plt.figure()
ax1 = plt.subplot(121, projection=subwcs, slices=('x', 'y'))
ax1.set_title('P')
im = ax1.imshow(P_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax1, location = 'bottom')
cb.set_label('P%')
ax1.set_xlabel('RA (J2000)')
ax1.set_ylabel('DEC (J2000)')

ax2 = plt.subplot(122, projection=subwcs, slices=('x', 'y'))
ax2.set_title('A')
im = ax2.imshow(A_visual, origin = 'lower', cmap = 'jet')
cb = fig.colorbar(im, ax = ax2, location = 'bottom')
cb.set_label('A (°)')
ax2.set_xlabel('RA (J2000)')
ax2.set_ylabel('DEC (J2000)')

x_text = 0
y_text = 0
xv = np.zeros([wsz,wsz])
yv = np.zeros([wsz,wsz])

scale = int(np.nanmax(P_visual)/10)*10
if(np.nanmax(P_visual)<10): 
    scale = int(np.nanmax(P_visual))
scale_asig = False
for i in range(0, np.size(P_visual,axis = 0), sep_vector): 
    for j in range(0, np.size(P_visual, axis = 1), sep_vector): 
        xv[i,j] = P_visual[i,j]*np.cos(A_visual[i,j]*np.pi/180)
        yv[i,j] = P_visual[i,j]*np.sin(A_visual[i,j]*np.pi/180)

        if((i >= int((wsz*0.15)) and j >= wsz - int((wsz*0.15))) and scale_asig==False): 
            xv[i,j] = scale*np.cos(0)
            yv[i,j] = scale*np.sin(0)
            print('SCALE VECTOR SETTED AT: ', i, j)
            x_text = j
            y_text = i
            scale_asig = True
            
fig = plt.figure()


ax1 = plt.subplot(121, projection=subwcs, slices=('x', 'y'))
ax1.set_title(obj +' OAGH-POLICAN.')
im1 = ax1.imshow(I_visual, origin='lower', cmap='Greys_r', vmax = 60)

#ax1.contour(P_visual, levels = 5, colors = 'black', alpha = 0.5)
ax1.quiver(xv,yv, color='red', scale = 5*np.nanmax(P_visual),
          width = 0.002, headwidth = 0, 
          pivot='middle', minlength = 0.01)
ax1.text(x = x_text - x_text*0.07, y = y_text - y_text*0.45, 
         s = str(scale) + ' %', size = 15, color = 'red')

#ax1.text(x = 0, y =0, 
#         s = 'PRELIMINARY', size = 20)
ax1.set_xlabel('RA (J2000)')
ax1.set_ylabel('DEC (J2000)')
cb = fig.colorbar(im1, ax=ax1, location ='bottom')
cb.set_label('MJy/sr')

ax2 = plt.subplot(122, projection=subwcs, slices=('x', 'y'))
ax2.set_title(obj +' OAGH-POLICAN.')
im2 = ax2.imshow(P_visual, origin='lower', cmap='jet_r')

ax2.contour(P_visual, levels = 5, colors = 'black', alpha = 0.5)
ax2.quiver(xv,yv, color='black', scale = 6*np.nanmax(P_visual),
          width = 0.002, headwidth = 0, 
          pivot='middle', minlength = 0.01)
ax2.text(x = x_text - x_text*0.07, y = y_text - y_text*0.45, 
         s = str(scale) + ' %', size = 15)

#ax1.text(x = 0, y =0, 
#         s = 'PRELIMINARY', size = 20)
ax2.set_xlabel('RA (J2000)')
ax2.set_ylabel('DEC (J2000)')
cb = fig.colorbar(im2, ax=ax2, location = 'bottom')
cb.set_label('P (%)')

plt.show()
