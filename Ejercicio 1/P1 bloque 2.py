#https://www.webelements.com/polonium/crystal_structure.html
# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import *
import os

#####Defninimos Matrices y vectores que utilizaremos
Matxyz = []; xmat = []; ymat = []; zmat = []

###############################################################################
#                               FUNCIONES
###############################################################################

'''
Definimos las funciones para crear la celda unidad en donde tendremos nx=ny=nz=2
'''

def CSUnitCell(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    nx = 2; ny = 2; nz = 2
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z])
                xmat.append(x)
                ymat.append(y)
                zmat.append(z)
    return Matxyz


def BCCUnitCell(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    nx = 2; ny = 2; nz = 2
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;y=0.5*Lattb;z=0.5*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    return Matxyz


def FCCUnitCell(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    nx = 2; ny = 2; nz = 2
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;y=0.5*Lattb;z=0.0*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;y=0.0*Lattb;z=0.5*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.0*Latta;y=0.5*Lattb;z=0.5*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;y=0.5*Lattb;z=1.0*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;y=1.0*Lattb;z=0.5*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=1.0*Latta;y=0.5*Lattb;z=0.5*Lattc
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    return Matxyz


def HCPUnitCell(Latta, Lattb):
    Latta=Latta #(dos r)
    Lattb=Lattb #(< dos r)
    Lattc=Latta #(< dos r)
    
    x=0.0;         y=0.0;         z=0.5*Lattb 
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.5*Latta;  y=0.0;         z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.25*Latta;  y=0.5*Lattc;   z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.25*Latta;   y=0.5*Lattc;   z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;   y=0.0;         z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.25*Latta;   y=-0.5*Lattc;  z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.25*Latta;  y=-0.5*Lattc;  z=0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    
    x=0.0;         y=-0.25*Lattb;         z=0.0 
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.25*Latta;  y=0.25*Lattb;         z=0.0 
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.25*Latta;   y=0.25*Lattb;         z=0.0 
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    
    x=0.0;         y=0.0;         z=-0.5*Lattb 
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.5*Latta;  y=0.0;         z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.25*Latta;  y=0.5*Lattc;   z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.25*Latta;   y=0.5*Lattc;   z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.5*Latta;   y=0.0;         z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=0.25*Latta;   y=-0.5*Lattc;  z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    x=-0.25*Latta;  y=-0.5*Lattc;  z=-0.5*Lattb
    Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    
    return Matxyz

'''
Definimos las funciones para crear el cristal, para estas definiremos los parametros de red
y colocaremos los minimos atomos para poder generar la estructura deseada sin huecos o atomos 
repetidos al desplazar los ejes.
'''

def CScrystal(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z])
                xmat.append(x)
                ymat.append(y)
                zmat.append(z)
    return Matxyz


def BCCcrystal(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(0.5+i)*Latta;y=(0.5+j)*Lattb;z=(0.5+k)*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    return Matxyz


def FCCcrystal(Latta):
    Latta=Latta
    Lattb=Latta
    Lattc=Latta
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                x=i*Latta;y=j*Lattb;z=k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(0.5+i)*Latta;y=(0.5+j)*Lattb;z=(0.0+k)*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(0.5+i)*Latta;y=(0.0+j)*Lattb;z=(0.5+k)*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(0.0+i)*Latta;y=(0.5+j)*Lattb;z=(0.5+k)*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
    return Matxyz

'''
En el caso del HCP es un caso especial en el que necesitamos introducir dos terminos el mu y el lam, 
con el fin de que no se repitan los atomos al desplazarlo. Siendo el lam el avance hacia la derecha/izq. 
y el mu de arriba/abajo. Simulando asi un panal de abejas.
Los atomos escogidos son los minimos para generar la estructura HCP sin repetir posiciones o generar huecos
'''
def HCPcrystal(Latta,Lattb):
    Latta=Latta #(dos r)
    Lattb=Lattb #(< dos r)
    Lattc=Lattb #(< dos r)
    mu = 0.5*Lattb
    lam = 0.75*Latta
    for i in np.arange(nx):
        for j in np.arange(ny):
            for k in np.arange(nz):
                
                x=(0.0)*Latta + i*lam + j*lam;    y=(-0.25)*Lattb - i*mu + j*mu;  z=(0.0)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(-0.25)*Latta + i*lam + j*lam;  y=(0.25)*Lattb - i*mu + j*mu;   z=(0.0)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(0.25)*Latta + i*lam + j*lam;   y=(0.25)*Lattb - i*mu + j*mu;   z=(0.0)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                
                x=(0.0)*Latta + i*lam + j*lam;    y=(0.0)*Lattb - i*mu + j*mu;     z=(-0.5)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(-0.5)*Latta + i*lam + j*lam;   y=(0.0)*Lattb - i*mu + j*mu;    z=(-0.5)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                x=(-0.25)*Latta + i*lam + j*lam;  y=(0.5)*Lattb - i*mu + j*mu;    z=(-0.5)*Lattc + k*Lattc
                Matxyz.append([x,y,z]); xmat.append(x); ymat.append(y); zmat.append(z)
                
    return Matxyz

###############################################################################
#                               INPUTS
###############################################################################

'''
Nuestros inputs seran: 
    si queremos celda unidad o cristal
    que tipo de estructura
    el elemento del que queremos 
    y en caso de cristal, cuantas celdas queremos en cada posicion
'''

print("Selecciona si se desea la celda unidad o el cristal")
print("*****************")
print("(a). celda")
print("(b). cristal")
print("*****************")
display_type = input('Opción elegida: ')
print('Selecciona la celdad unidad')
print("*****************")
print("(1). SC")
print("(2). FCC")
print("(3). BCC")
print("(4). HCP")
print("*****************")
m = input('Opción elegida: ')

if m=="1":
    coloret = 'blue'
    print("CS")
    print("Selecciona un elemento")
    print("*****************")
    print("Po")
    print("*****************")
    Element = input('Escribe el elemento:')
    if Element=="Po":
       Latta=0.3359
       print('Parámetro celda unidad=',Latta,'nm')
       print("*****************")
       if display_type == 'a':
           CSUnitCell(Latta)   
       elif display_type == 'b':
           nx = int(input('Numero de celdas en la componente x: '))
           ny = int(input('Numero de celdas en la componente y: '))
           nz = int(input('Numero de celdas en la componente z: '))
           CScrystal(Latta)
                 
    else:
       print('El elemento no estaba en la tabla. Debe incluirlo manualmentes introducir su parámetro de red')
      
  
                  
elif m== "2":
    coloret = 'orange'
    print("FCC")
    print("Selecciona un elemento")
    print("*****************")
    print("Xe")
    print("*****************")
    Element = input('Escribe el elemento:')
    if Element=="Xe":
       Latta = 0.624
       print('Parámetro celda unidad=',Latta,'nm')
       print("*****************")
       if display_type == 'a':
           FCCUnitCell(Latta) 
       elif display_type == 'b':
            nx = int(input('Numero de celdas en la componente x: '))
            ny = int(input('Numero de celdas en la componente y: '))
            nz = int(input('Numero de celdas en la componente z: '))
            FCCcrystal(Latta)
    else:
       print('El elemento no estaba en la tabla. Debe incluirlo manualmentes introducir su parámetro de red')

elif m== "3":
    coloret = 'green'
    print("BCC")
    print("Selecciona un elemento")
    print("*****************")
    print("V")
    print("*****************")
    Element = input('Escribe el elemento:')
    if Element=="V":
       Latta=0.3359
       print('Parámetro celda unidad=',Latta,'nm')
       print("*****************")
       if display_type == 'a':
           BCCUnitCell(Latta) 
       elif display_type == 'b':
           nx = int(input('Numero de celdas en la componente x: '))
           ny = int(input('Numero de celdas en la componente y: '))
           nz = int(input('Numero de celdas en la componente z: '))
           BCCcrystal(Latta)
    else:
       print('El elemento no estaba en la tabla. Debe incluirlo manualmentes introducir su parámetro de red')

elif m== "4":
    coloret = 'gray'
    print("HCP")
    print("Selecciona un elemento")
    print("*****************")
    print("Co")
    print("*****************")
    Element = input('Escribe el elemento:')
    if Element=="Co":
       Latta=0.2507
       Lattc = 0.4069
       print('Parámetro celda unidad en a=',Latta,'nm')
       print('Parámetro celda unidad en c=',Lattc,'nm')
       print("*****************")
       #Lattb = np.sqrt(3)/2*Latta
    
       if display_type == 'a':
           HCPUnitCell(Latta, Lattc) 
       elif display_type == 'b':
           nx = int(input('Numero de celdas en la componente x: '))
           ny = int(input('Numero de celdas en la componente y: '))
           nz = int(input('Numero de celdas en la componente z: '))
           HCPcrystal(Latta, Lattc)
    else:
       print('El elemento no estaba en la tabla. Debe incluirlo manualmentes introducir su parámetro de red')
      
else:
    print("ERROR FATAL")

###############################################################################
#                               OUTPUTS
###############################################################################


out=str(len(Matxyz)) + '\n' + "Electrodo de " + str(m) + "\n"
for i in range(0,len(Matxyz)):   
    xx = xmat[i]
    yy = ymat[i]
    zz = zmat[i]
    
    out+= str(Element) +"         "+ str(xx)  +"          " + str(yy)  +"          "  + str(zz) + '\n'
        
script_directory = os.path.dirname(__file__)
fn = os.path.join(script_directory, display_type + Element + m + '.xyz')
with open(fn, "w") as f:
    f.write(out)
f.close()

print("Número total de posiciones generadas: "+str(len(Matxyz)))

###############################################################################
#                               GRAFICA
###############################################################################
 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xmat, ymat, zmat, s=50, zdir='z',c=coloret, marker='o')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.set_aspect('equal')
#plt.savefig('Cristal hcp.png')
plt.show()