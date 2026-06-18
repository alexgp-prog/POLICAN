# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 14:48:26 2024

@author: AleX
"""
### PStokes MINI v1.0
### PROGRAMA BÁSICO DE ANÁLISIS DE LOS PARÁMETROS DE STOKES

import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import matplotlib.colors as colors
import numpy as np
from astropy.io import fits

plt.close('all')
## PARÁMETROS A MODIFICAR POR PARTE DEL USUARIO
cap = '//wsl.localhost/Ubuntu/home/alejandro/alex/Reducciones/Observaciones_2024/2024_12_01_02/M82_H/' # Carpeta donde se encuentra el archivo, si se está en ella dejar vacío. '/home/alejandro/alex/Reducciones/Observaciones_2024/2024_12_01_02/M82_H/'
root_name ='HM82H01' # raiz del nombre de los objetos 
paso_vector = 8     # Espacio entre vectores, para mejor visualización
wSz = 300           # Tamaño de la ventana deseado en pixeles
umbral = 600        # Umbral de intensidad para tomar valores relevantes, en ADUs
mask = np.nan       # Valor de enmascaramiento

x_offset = 0        # Offset en x
y_offset = -100     # Offset en y


# INICIA PROGRAMA--------------------------------------------
## SE RECOMIENDA NO MODIFICAR MÁS ALLÁ DE ESTE PUNTO
# Nombre del archivo
file1 = cap + root_name+'00_1b.fits'
file2 = cap + root_name+'22_1b.fits'
file3 = cap + root_name+'45_1b.fits'
file4 = cap + root_name+'67_1b.fits'

def openfits(filename): 

    hdul = fits.open(filename)
    hdul.info()
    data = hdul[0].data
    
    return data

## LECTURA DE ARCHIVOS
im00 = openfits(file1)
im22 = openfits(file2)
im45 = openfits(file3)
im67 = openfits(file4)
## CÁLCULO DE LOS PARÁMETROS DE STOKES

I=(im00+im22+im45+im67)/2
dataI = np.copy(I)
Qi=(im00-im45)/(I)
Ui=(im22-im67)/(I)

## MATRÍZ DE ROTACIÓN

#Angulo de offset y eficiencia
theta = 139*(np.pi)/180
etha = 0.963
Qeq = (Qi*np.cos(2*theta)-Ui*np.sin(2*theta))/etha
Ueq = (Ui*np.cos(2*theta)+Qi*np.sin(2*theta))/etha
# Parametros instrumentales
Qinst = -0.0050
Uinst = 0.0012
# PARAMETROS Q y U finales
Q = Qeq - Qinst
U = Ueq - Uinst

### PARÁMETROS DE STOKES FINALES Y CÁLCULO DE GRADO Y ÁNGULO DE POLARIZACIÓN


Q[I<umbral] = mask
U[I<umbral] = mask
I[I<umbral] = mask

Q[np.abs(Q)>1] = mask
U[np.abs(Q)>1] = mask

Q[np.abs(U)>1] = mask
U[np.abs(U)>1] = mask

print(I.shape)
print(Q.shape)
print(U.shape)

Szx = np.size(I, axis = 1)
Szy = np.size(I, axis = 0)

list_stokes = [dataI,Q,U]
list_stokes_title = ['Stokes I', 'Stokes Q', 'Stokes U']
list_maxval = [3000, Q.max()/2,U.max()/2]
list_minval = [0, Q.min()/2,U.min()/2]
list_labels = ['ADUs', 'Q/I', 'U/I']

fig, ax = plt.subplots(nrows = 1, ncols = 3, sharey = True)

for i in range(0,3):

    ax[i].set_title(list_stokes_title[i])
    im = ax[i].imshow(list_stokes[i], origin = 'lower', cmap = 'jet', 
                      vmax = list_maxval[i], vmin = list_minval[i])
    cb = fig.colorbar(im, ax =  ax[i], location = 'bottom')

    cb.set_label(list_labels[i])
    ax[i].set_xlabel('Pixels')
    if(i<1):
        ax[i].set_ylabel('Pixels')
    
    ax[i].set_xlim([int(Szx/2) - int(wSz/2)-x_offset, int(Szx/2) + int(wSz/2)-x_offset])
    ax[i].set_ylim([int(Szy/2) - int(wSz/2)-y_offset, int(Szy/2) + int(wSz/2)-y_offset])
        
        
P = 100*np.sqrt(Q**2+U**2)
A = 0.5*np.arctan(np.abs(U)/np.abs(Q)) * (180/np.pi)


for j in range (np.size(A, axis = 0)):
    for i in range (np.size(A, axis = 1)):
        if (Q[i,j] >= 0.0) and (U[i,j] >= 0.0): 
            A[i][j] = 90 + A[i][j]
        if (Q[i,j] <= 0.0) and (U[i,j] >= 0.0):
            A[i][j] = 180 - A[i][j]
        if (Q[i,j] <= 0.0) and (U[i,j] <= 0.0): 
            A[i][j] = A[i][j]
        if (Q[i,j] >= 0.0) and (U[i,j] <= 0.0):
            A[i][j] = 90 - A[i][j]

A[P>100] = mask
P[P>100] = mask

A[I<umbral] = mask
P[I<umbral] = mask


fig, ax = plt.subplots(nrows = 1, ncols = 2)
pol_list = [P,A]
pol_list_title = ['Polarization Fraction', 'Polarization Angle']
pol_label = ['P(%)', 'A(°)']
pol_vmax = [np.nanmax(P)/2, np.nanmax(A)]
for i in range(2):
    ax[i].set_title(pol_list_title[i])
    im = ax[i].imshow(pol_list[i], origin = 'lower', cmap = 'jet', vmax = pol_vmax[i])
    cb = fig.colorbar(im, ax = ax[i], location = 'bottom')
    cb.set_label(pol_label[i])
    ax[i].set_xlabel('Pixels')
    ax[i].set_ylabel('Pixels')
    
    ax[i].set_xlim([int(Szx/2) - int(wSz/2)-x_offset, int(Szx/2) + int(wSz/2)-x_offset])
    ax[i].set_ylim([int(Szy/2) - int(wSz/2)-y_offset, int(Szy/2) + int(wSz/2)-y_offset])
    
x = np.linspace(0, np.size(P, axis = 1)-1, np.size(P, axis = 1))
y = np.linspace(0, np.size(P, axis = 0)-1, np.size(P, axis = 0))
X,Y = np.meshgrid(x,y)
    
xv = np.zeros([np.size(P, axis = 1),np.size(P, axis = 1)])
yv = np.zeros([np.size(P, axis = 0),np.size(P, axis = 0)])
    
for i in range(0,np.size(A, axis = 0), paso_vector):
    for j in range(0,np.size(A, axis = 1), paso_vector):
        
        xv[i,j] = P[i,j]*np.cos(A[i,j]*np.pi/180)
        yv[i,j] = P[i,j]*np.sin(A[i,j]*np.pi/180)
        
sc = round(np.nanmax(P)/10 + 3)*10

print('POLARIZACIÓN MÁXIMA: ', np.nanmax(P))
# print('ESCALA DE POLARIZACIÓN: ', sc)
print('ESCALA DE POLARIZACIÓN: ', ((sc/10) - 3)*10)

fig, ax = plt.subplots() 
ax.set_title('Polarization Fraction w/Vectors')
im = ax.imshow(P, origin = 'lower', cmap = 'jet', vmax = np.nanmax(P)/2)
cb = fig.colorbar(im, ax = ax)
ax.quiver(X,Y,xv,yv,
              color='black', scale = sc,
          width = 0.004, headwidth = 0, 
          pivot='middle', minlength = 0.001)

cb.set_label('P (%)')
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')

ax.set_xlim([int(Szx/2) - int(wSz/2)-x_offset, int(Szx/2) + int(wSz/2)-x_offset])
ax.set_ylim([int(Szy/2) - int(wSz/2)-y_offset, int(Szy/2) + int(wSz/2)-y_offset])

ax.text(int(Szx/2) + int(wSz/2)-x_offset - 50, int(Szy/2) - int(wSz/2)-y_offset + 10, 
        str(((sc/10) - 3)*10) +'%',color = 'black')
ax.axhline(y = int(Szy/2) - int(wSz/2)-y_offset + 25, xmin = 0.86, xmax = 0.94, 
           color = 'black', linestyle = '-')

fig, ax = plt.subplots() 
ax.set_title('Stokes I w/Vectors')
im = ax.imshow(dataI, origin = 'lower', cmap = 'inferno', vmax = 3000)
cb = fig.colorbar(im, ax = ax)
ax.quiver(X,Y,xv,yv,
              color='white', scale = sc,
          width = 0.004, headwidth = 0, 
          pivot='middle', minlength = 0.001)

cb.set_label('ADUs')
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')

ax.set_xlim([int(Szx/2) - int(wSz/2)-x_offset, int(Szx/2) + int(wSz/2)-x_offset])
ax.set_ylim([int(Szy/2) - int(wSz/2)-y_offset, int(Szy/2) + int(wSz/2)-y_offset])

ax.text(int(Szx/2) + int(wSz/2)-x_offset - 50, int(Szy/2) - int(wSz/2)-y_offset + 10, 
        str(((sc/10) - 3)*10) +'%',color = 'white')
ax.axhline(y = int(Szy/2) - int(wSz/2)-y_offset + 25, xmin = 0.86, xmax = 0.94, 
           color = 'white', linestyle = '-')
        
    
            
plt.show()

# Qp = Q>=0.0
# Qm = Q<=0.0
# Up = U>=0.0
# Um = U<=0.0

# A[Qp & Up] = 90 + A
# A[Qm & Up] = 180 - A
# A[Qm & Um] = A
# A[Qp & Um] = 90 - A

# cap = '/home/alejandro/alex/Reducciones/Observaciones_2024/2024_12_01_02/M82_H/'
    
