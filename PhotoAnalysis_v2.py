import matplotlib.pyplot as plt
import numpy as np
from photutils.centroids import centroid_sources, centroid_com
from photutils.aperture import CircularAperture,EllipticalAperture, CircularAnnulus
from photutils.aperture import aperture_photometry, ApertureStats
from astropy.io import fits

def graph(im1):
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))
    im = ax[0].imshow(im1, origin = 'lower', vmax = im1.flatten().std(), vmin = -im1.flatten().std())
    fig.colorbar(im, ax = ax[0])

    ax[1].hist(im1.flatten(),log = True, bins = 100, color = 'blue')
    ax[1].axvline(-im1.flatten().std(), color = 'red', linestyle = '--', 
            linewidth = 2, 
            label = '$ -\sigma $ = ' + str(round(im1.flatten().std(),3)))
    ax[1].axvline(im1.flatten().std(), color = 'green', linestyle = '--', 
            linewidth = 2, 
            label = '$ \sigma $ = ' + str(round(im1.flatten().std(),3)))
    ax[1].legend()

def openfits(archivo):
    hdul = fits.open(archivo)
    hdul.info()
    imagen = hdul[0].data
    head = hdul[0].header
    hdul.close()
    return imagen, head

#%%//wsl.localhost/Ubuntu/home/alejandro/alex/Reducciones/2023_12_12_13/HD38_H/HD38CH0122_1b.fits
# Nombre del objeto y archivo 
plt.close('all')

#obj = 'NGC2068'
# D:/Alex/INAOE/MAESTRIA/Reducciones/NGC2068/Red2_15/HN2068H0100_1b.fits
#root = 'HN2068H01'
#cap = 'D:/Alex/INAOE/MAESTRIA/Reducciones/NGC2068/Red2_15/' 

obj = 'NGC 2068 H'
root = 'HN2068H01'
cap = '//wsl.localhost/Ubuntu/home/alejandro/alex/Reducciones/Observacion_Nov_2025/2025_11_24_25/NGC2068/'

file1 = cap+root+'00_1b.fits'
file2 = cap+root+'22_1b.fits'
file3 = cap+root+'45_1b.fits'
file4 = cap+root+'67_1b.fits'

im1,head1 = openfits(file1)
im2,head2 = openfits(file2)
im3,head3 = openfits(file3)
im4,head4 = openfits(file4)

all_minor = np.array([im1[im1<0].mean(), 
                      im2[im1<0].mean(), 
                      im3[im1<0].mean(), 
                      im4[im1<0].mean()])

print(head1)
print(all_minor)
print(all_minor.min())

if(all_minor.min()<0):
   im1+=np.abs(all_minor.min())
   im2+=np.abs(all_minor.min())
   im3+=np.abs(all_minor.min())
   im4+=np.abs(all_minor.min())
    
exptime = head1['EXPTIME']
UT = head1['DATE-OBS']

print(type(exptime))

n1 = head1['NCOMBINE']
n2 = head2['NCOMBINE']
n3 = head3['NCOMBINE']
n4 = head4['NCOMBINE']

sizey = np.size(im1, axis = 0)
sizex = np.size(im1, axis = 1)

mean = im1.flatten().mean()
standard = im1.flatten().std()
im11 = im1[0:int(sizey/2), 0:int(sizex/2)]
im12 = im1[int(sizey/2):sizey, 0:int(sizex/2)]
im13 = im1[int(sizey/2):sizey, int(sizex/2):sizex]
im14 = im1[0:int(sizey/2), int(sizex/2):sizex]

# graph(im11)
# graph(im12)
# graph(im13)
# graph(im14)

standard1 = im11.flatten().std()
standard2 = im12.flatten().std()
standard3 = im13.flatten().std()
standard4 = im14.flatten().std()

#%%
# Inicia procesamiento
plt.close('all')
im1stars = np.zeros([sizey,sizex])
sec = 40
iterations = 10
sig = 1.5

threshold1 = sig*standard1
threshold2 = sig*standard2
threshold3 = sig*standard3
threshold4 = sig*standard4

if (threshold1 >= threshold2 and threshold1 >= threshold3 and threshold1 >= threshold4):
    threshold = threshold1
elif (threshold2 >= threshold1 and threshold2 >= threshold3 and threshold2 >= threshold4):
    threshold = threshold2
elif (threshold3 >= threshold1 and threshold3 >= threshold2 and threshold3 >= threshold4):
    threshold = threshold3
elif (threshold4 >= threshold1 and threshold4 >= threshold2 and threshold4 >= threshold3):
    threshold = threshold4
    
x = []
y = []
xy = []
c = 0

for i in range(0,np.size(im1, axis = 0)):
    for j in range(0,np.size(im1, axis = 1)): 
        if(im1[i,j] >= threshold and (i >= sec and i <= sizey-sec) 
           and (j >= sec and j <= sizex-sec)):
            im1stars[i,j] = im1[i,j]
            if (c == 1):
                if (i > y[-1] + 10 or i < y[-1] - 10):
                    x.append(j)
                    y.append(i)
            else:
                x.append(j)
                y.append(i)
                c = 1
           
fig, ax = plt.subplots()
    
ax.imshow(im1stars,origin = 'lower', 
          interpolation='nearest', cmap = 'viridis', 
          vmin = -standard, 
          vmax = standard)
# ax.scatter(x, y, marker='+', s=80, color='red')
for i in range(iterations):
    x,y = centroid_sources(im1stars, xpos = x, ypos = y, 
                          box_size=21,
                          centroid_func=centroid_com)

ax.scatter(x, y, marker='+', s=80, color='green')
#%%
plt.close('all')

for i in range(len(x)): 
    xy.append((x[i],y[i]))

rad = 10

aperture1 = CircularAperture(xy, r = rad)
aperture2 = CircularAperture(xy, r = rad)
aperture3 = CircularAperture(xy, r = rad)
aperture4 = CircularAperture(xy, r = rad)

# ESTADISTICAS DE LAS APERTURAS
aperstats00 = ApertureStats(im1, aperture1) 
aperstats22 = ApertureStats(im2, aperture2) 
aperstats45 = ApertureStats(im3, aperture3) 
aperstats67 = ApertureStats(im4, aperture4) 

ap_mean00 = aperstats00.mean
ap_mean22 = aperstats22.mean
ap_mean45 = aperstats45.mean
ap_mean67 = aperstats67.mean

ap_sig00 = aperstats00.std
ap_sig22 = aperstats22.std
ap_sig45 = aperstats45.std
ap_sig67 = aperstats67.std

phot_table1 = aperture_photometry(im1, aperture1)
phot_table2 = aperture_photometry(im2, aperture2)
phot_table3 = aperture_photometry(im3, aperture3)
phot_table4 = aperture_photometry(im4, aperture4)   

a00 = (phot_table1['aperture_sum'])
a22 = (phot_table2['aperture_sum'])
a45 = (phot_table3['aperture_sum'])
a67 = (phot_table4['aperture_sum'])

a00sd = np.sqrt(np.abs(a00)/(5.8*n1)) + ap_sig00
a22sd = np.sqrt(np.abs(a22)/(5.8*n2)) + ap_sig22
a45sd = np.sqrt(np.abs(a45)/(5.8*n3)) + ap_sig45
a67sd = np.sqrt(np.abs(a67)/(5.8*n4)) + ap_sig67

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
Ueq = (U*np.cos(2*theta)+Q*np.sin(2*theta))/etha
Qeqsd = np.sqrt((Qsd*(np.cos(2*theta))/etha)**2 + (Usd*(np.sin(2*theta))/etha)**2 )
Ueqsd = np.sqrt((Usd*(np.cos(2*theta))/etha)**2 + (Qsd*(np.sin(2*theta))/etha)**2 )
Qf = Qeq - Qinst
Uf = Ueq - Uinst
Qsd = np.sqrt(Qeqsd**2 + 0.002**2)
Usd = np.sqrt(Ueqsd**2 + 0.002**2)
   

P = 100*np.sqrt(Qf**2 + Uf**2)
Psd = (100**2/P)*np.sqrt((Qf**2)*(Qsd**2)+(Uf**2)*(Usd**2))
for i in range(len(P)):
    if (P[i] > Psd[i]): 
        P[i] = np.sqrt(P[i]**2 - Psd[i]**2)

print(P)
print(Psd)

Ain = (0.5*(np.arctan(np.abs(Uf)/np.abs(Qf))))*(180/np.pi)
Aor = np.zeros([len(Ain)])

# for i in range(len(Ain)):
#     if (Qf[i] > 0.0) and (Uf[i] > 0.0): 
#         Ain[i] = 90 - Ain[i]
#     elif (Qf[i] <= 0.0) and (Uf[i] > 0.0):
#         Ain[i] = Ain[i]
#     elif (Qf[i] <= 0.0) and (Uf[i] <= 0.0): 
#         Ain[i] = 180 - Ain[i]
#     elif (Qf[i] > 0.0) and (Uf[i] <= 0.0):
#         Ain[i] = 90 + Ain[i]
    
#     if (Ain[i] >= 90):
#         Aor[i] = Ain[i] - 90
#     elif(Ain[i] < 90):
#         Aor[i] = Ain[i] + 90
        

print('ESQUEMA PROPIO')
for i in range (len(Ain)):
    if (Qf[i] >= 0.0) and (Uf[i] >= 0.0): 
            Ain[i] = 90 + Ain[i]
    if (Qf[i] <= 0.0) and (Uf[i] >= 0.0):
            Ain[i] = 180 - Ain[i]
    if (Qf[i] <= 0.0) and (Uf[i] <= 0.0): 
            Ain[i] = Ain[i]
    if (Qf[i] >= 0.0) and (Uf[i] <= 0.0):
            Ain[i] = 90 - Ain[i]
    
    # Aor = Ain
    if (Ain[i] >= 90):
        Aor[i] = Ain[i] - 90
    elif(Ain[i] < 90):
        Aor[i] = Ain[i] + 90
        
Asd = 28.65*(Psd/P)#np.sqrt(((-Uf*Qsd)/((Qf**2)+(Uf**2)))**2 + ((Qf*Usd)/((Qf**2)+(Uf**2)))**2)

Peq=100*(np.sqrt((Qeq**2)+(Ueq**2)))

print(Qf)
print(Uf)
print(P)
print(Ain)
print(Aor)


# print("Ain: ",Ain)
#%%
plt.figure()
plt.title(obj + ' (I) + POLARIZACIÓN (t:'+str(exptime)+' s, r: '+ str(rad)+  'pix ap). UT: ' + UT,fontsize = 10)
plt.imshow(im1, origin = 'lower', cmap = 'jet', vmax = 6*im1.flatten().std(), vmin = 0)
plt.grid(False)
cb = plt.colorbar()
cb.set_label('I (ADUs)', rotation = 270,fontsize = 10)
aperture1.plot(color='white', lw=2, label='Photometry aperture')
# annulus_aperture1.plot(color='red', lw=2,label='Background annulus')
plt.xlabel('pixeles',fontsize = 10)
plt.ylabel('pixeles',fontsize = 10)
# plt.legend(labelcolor = 'white')
s = 0
for i,j in zip(x,y):
    plt.text(i+25,j,'P: '+ str(round(P[s],2)) + '$\pm$' + str(round(Psd[s],2)),
             color='white', fontsize = 10)
    plt.text(i+25,j-20,'A: '+str(round(Aor[s],2))+ '$\pm$' + str(round(Asd[s],2)),
             color='white',fontsize = 10)
    s=s+1
    

#%%]
# GRAFICACION

PSoagh = 0.32
arcsec1 = sizex*PSoagh/2
arcsec2 = sizey*PSoagh/2
Lz1 = np.linspace(-arcsec1, arcsec1,sizex)
Lz2 = np.linspace(-arcsec2, arcsec2,sizey)
X,Y = np.meshgrid(Lz1, Lz2)

Pm = np.zeros([sizey,sizex])
Am = np.zeros([sizey,sizex])
xv=np.zeros([sizey,sizex])
yv=np.zeros([sizey,sizex])

for i in range(sizey): 
    for j in range(sizex): 
        for ob in range(len(P)):
            if (i == round(y[ob]) and j == round(x[ob])):
                Pm[i,j] = P[ob]
                Am[i,j] = Ain[ob]
                break
            
        xv[i][j]=Pm[i][j]*np.cos(Am[i][j]*np.pi/180)
        yv[i][j]=Pm[i][j]*np.sin(Am[i][j]*np.pi/180)

Pmax = np.mean(P)
print("\n Polarización Máxima:\n", Pmax)
if (Pmax<10): 
    Pc = int(10*(Pmax/10))
else:
    Pc = (10*int(Pmax/10))
print("\n Polarización escala:\n", Pc)
sc = 10*(Pc)
#%%

fig, ax0 = plt.subplots()
ax0.set_title(obj + ' (I) + POLARIZACIÓN',fontsize = 10)
im = ax0.imshow(im1,vmax = im1.std(), vmin = -im1.std(), cmap='Greys_r'
                ,extent = [-arcsec1, arcsec1, -arcsec2, arcsec2], origin = 'lower')
# im = ax0.pcolor(X3,Y3,Iff, vmax = np.max(Iff/4), cmap='magma')
ax0.grid(False)
cb = fig.colorbar(im, ax = ax0)
cb.set_label('I(ADUs)', rotation = 270,fontsize = 10)
ax0.quiver(X,Y,xv,yv, scale = sc ,color='red', width = 0.003, headwidth = 0, pivot='middle', minlength = 0.01)
ax0.set_xlabel('arcsec',fontsize = 10)
ax0.set_ylabel('arcsec',fontsize = 10)
ax0.text(arcsec2 - (arcsec2/4),-arcsec2 + (arcsec2/6), '%i'% Pc +'%',color = 'red', size = 15)
ax0.axhline(y = -arcsec2 + (arcsec2/10), xmin = 0.86, xmax = 0.94, color = 'red', linestyle = '-')

#%%
### GENERACIÓN DE DATAFRAME

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


#I = (im1 + im2 + im3 + im4)/2
#Q = (im1 - im3)
#U = (im2 - im4)

#hduI = fits.PrimaryHDU(data = I,header = head1)
#hduI.writeto(cap + root + '_POLICAN.fits',overwrite=True,output_verify='warn')

plt.show()

# def graph(im1):
#     fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (15,10))
#     im = ax[0].imshow(im1, origin = 'lower', vmax = im1.flatten().std(), vmin = -im1.flatten().std())
#     fig.colorbar(im, ax = ax[0])

#     ax[1].hist(im1.flatten(),log = True, bins = 100, color = 'blue')
#     ax[1].axvline(-im1.flatten().std(), color = 'red', linestyle = '--', 
#             linewidth = 2, 
#             label = '$ -\sigma $ = ' + str(round(im1.flatten().std(),3)))
#     ax[1].axvline(im1.flatten().std(), color = 'green', linestyle = '--', 
#             linewidth = 2, 
#             label = '$ \sigma $ = ' + str(round(im1.flatten().std(),3)))
#     ax[1].legend()

#if(im11[i,j] >= threshold1 and (i >= sec and i <= sizey-sec) 
   # and (j >= sec and j <= sizex-sec)):

#GRAFICAS

# plt.figure()
# plt.title("00",fontsize = 20)
# plt.imshow(im1, origin = 'lower', cmap = 'inferno', vmax = im1.flatten().std(), vmin = -im1.flatten().std())
# plt.grid(False)
# cb = plt.colorbar()
# cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
# aperture1.plot(color='white', lw=2, label='Photometry aperture')
# # annulus_aperture1.plot(color='red', lw=2,label='Background annulus')
# plt.xlabel('pixeles',fontsize = 20)
# plt.ylabel('pixeles',fontsize = 20)
# # plt.legend(labelcolor = 'white')

# plt.show()

# plt.figure()
# plt.title("22",fontsize = 20)
# plt.imshow(im1, origin = 'lower', cmap = 'inferno', vmax = im2.flatten().std(), vmin = -im2.flatten().std())
# plt.grid(False)
# cb = plt.colorbar()
# cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
# aperture2.plot(color='white', lw=2, label='Photometry aperture')
# # annulus_aperture2.plot(color='red', lw=2,label='Background annulus')
# plt.xlabel('pixeles',fontsize = 20)
# plt.ylabel('pixeles',fontsize = 20)
# # plt.legend(labelcolor = 'white')

# plt.show()

# plt.figure()
# plt.title("45",fontsize = 20)
# plt.imshow(im3, origin = 'lower', cmap = 'inferno', vmax = im3.flatten().std(), vmin = -im3.flatten().std())
# plt.grid(False)
# cb = plt.colorbar()
# cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
# aperture3.plot(color='white', lw=2, label='Photometry aperture')
# # annulus_aperture3.plot(color='red', lw=2,label='Background annulus')
# plt.xlabel('pixeles',fontsize = 20)
# plt.ylabel('pixeles',fontsize = 20)
# # plt.legend(labelcolor = 'white')

# plt.show()

# plt.figure()
# plt.title("67",fontsize = 20)
# plt.imshow(im4, origin = 'lower', cmap = 'inferno', vmax = im4.flatten().std(), vmin = -im4.flatten().std())
# plt.grid(False)
# cb = plt.colorbar()
# cb.set_label('I (ADU\'s)', rotation = 270,fontsize = 20)
# aperture4.plot(color='white', lw=2, label='Photometry aperture')
# # annulus_aperture4.plot(color='red', lw=2,label='Background annulus')
# plt.xlabel('pixeles',fontsize = 20)
# plt.ylabel('pixeles',fontsize = 20)
# # plt.legend(labelcolor = 'white')


# plt.show()

# for i in range(0,int(sizey/2)):
#     for j in range(0,int(sizex/2)): 
#         if(im11[i,j] >= threshold):
#             im1stars[i,j] = im11[i,j]
#             if (c == 1):
#                 if (i > y[-1] + 10 or i < y[-1] - 10):
#                     x.append(j)
#                     y.append(i)
#             else:
#                 x.append(j)
#                 y.append(i)
#                 c = 1
#         if(im12[i,j] >= threshold):
#             im1stars[i+int(sizey/2),j] = im12[i,j]
#             if (c == 1):
#                 if (i + int(sizey/2) > y[-1] + 10 or i + int(sizey/2) < y[-1] - 10):
#                     x.append(j)
#                     y.append(i + int(sizey/2))
#             else:
#                 x.append(j)
#                 y.append(i + int(sizey/2))
#                 c = 1
#         if(im13[i,j] >= threshold):
#             im1stars[i+int(sizey/2),j+int(sizex/2)] = im13[i,j]
#             if (c == 1):
#                 if (i + int(sizey/2)> y[-1] + 10 or i + int(sizey/2) < y[-1] - 10):
#                     x.append(j + int(sizex/2))
#                     y.append(i + int(sizey/2))
#             else:
#                 x.append(j + int(sizex/2))
#                 y.append(i + int(sizey/2))
#                 c = 1
#         if(im14[i,j] >= threshold):
#             im1stars[i,j+int(sizex/2)] = im14[i,j]
#             if (c == 1):
#                 if (i > y[-1] + 10 or i < y[-1] - 10):
#                     x.append(j + int(sizex/2))
#                     y.append(i)
#             else:
#                 x.append(j + int(sizex/2))
#                 y.append(i)
#                 c = 1
