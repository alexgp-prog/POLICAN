# CAMBIO DE NOMBRE DE LOS ARCHIVOS OBSERVACIONALES
# V1. 
# SE IMPORTA LA LIBRERIA NECESARIA
import os 
from astropy.io import fits

def openfits(archivo):
    hdul = fits.open(archivo)
    head = hdul[0].header
    hdul.close()
    return head

# Base del nombre anterior, ejemplo, si el archivo se llama 
# "N2068H1010001_I.fits" solamente incluir el "N2068H101".
nombre_anterior = "N2071H04"
cap = '//wsl.localhost/Ubuntu/home/alejandro/alex/Reducciones/Observacion_Dic_2023/2023_12_15_16/'
# De igual forma solamente se pone la base del nombre que se quiere
# Utilizar, ya que el programa agregará el ángulo correspondiente a la
# imagen.
nombre_nuevo = "HN2071H03"

i00 = 16
i22 = 16
i45 = 16
i67 = 16

for i in range (2,62):
    file = cap+ nombre_anterior+str(i).zfill(4)+'_I.fits'
    head = openfits(file)
    print(head['POL_ANG'])
    
    if (head['POL_ANG'] == 0):
        filenn = cap+ nombre_nuevo+str(0).zfill(2)+str(i00).zfill(4)+'_I.fits'
        i00+=1
    if (head['POL_ANG'] == 22.5):
        filenn = cap+ nombre_nuevo+str(22).zfill(2)+str(i22).zfill(4)+'_I.fits'
        i22+=1
    if (head['POL_ANG'] == 45):
        filenn = cap+ nombre_nuevo+str(45).zfill(2)+str(i45).zfill(4)+'_I.fits'
        i45+=1
    if (head['POL_ANG'] == 67.5):
        filenn = cap+ nombre_nuevo+str(67).zfill(2)+str(i67).zfill(4)+'_I.fits'
        i67+=1
    print(file) 
    print(filenn)   
    # Renombramiento de archivos.
    os.rename(file,filenn)

# FIN DEL PROGRAMA.
