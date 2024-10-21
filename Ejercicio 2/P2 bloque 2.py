#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 16:52:02 2023

@author: gabriel
"""

import numpy as np
import matplotlib.pyplot as plt



#Definimos las constantes que vamos a necesitar

GM = 6.6738e-11 * 1.9891e30   #podemos definirlas juntos 
mt = 5.9e24
h = 3600                #nuestro paso de tiempo es de 1h
t0 = 0.
tf = 5*365*24*60*60


#Definimos una funcion que nos devuelva la ecuacion de mov. a resolver y f1 para RK2 y RK4(facilidad)

def f(x,y):
    r = np.sqrt(x**2 + y**2)
    return (-GM*x/r**3 , -GM*y/r**3)

def f1(r):
    return -GM*r/np.dot(r,r)**1.5


#Definimos las listas que vamos a necesitar

t_an = np.linspace(0,5,5*365*24 + 1)
t_s = np.linspace(t0,tf+h,h)
N = len(t_an)
x , y , r = np.empty(N) , np.empty(N) , np.empty(N)
vx , vy , v = np.empty(2*N) , np.empty(2*N) , np.empty(N) #Es 2N porque los contamos semienteros
Ep , Ec  = np.empty(N) , np.empty(N)


#Condiciones iniciales, llenanmos las listas anteriores con los primeros valores

x[0], y[0]  = 1.4719e11  , 0.   #Consideramos que en el punto mas cercano (perihelio)
vx[0], vy[0]  = 0. , 3.0287e4   


r[0] = np.sqrt(x[0]**2+y[0]**2)
v[0] = np.sqrt(vx[0]**2+vy[0]**2)

Ep[0] = -GM*mt/r[0]
Ec[0] = 0.5*mt*v[0]**2

#La velocidad la contamos por semienteros, entonces calculamos el primer semientero
#Tenemos en cuenta que necesitaremos h/2

vx[1] = vx[0] + h/2 * f(x[0],y[0])[0] 
vy[1] = vy[0] + h/2 * f(x[0],y[0])[1]

#Ahora haremos el bucle para llenar todas las listas anteriores

for i in range(1,N): #empezamos en 1 porque el 0 es la condicion inicial
    #Llenamos x,y, y radio
    x[i] = x[i-1] + h * vx[2*i-1]
    y[i] = y[i-1] + h * vy[2*i-1]  
    r[i] = np.sqrt(x[i]**2 + y[i]**2)
    
    kx = h*f(x[i],y[i])[0]
    ky = h*f(x[i],y[i])[1]

    vx[2*i] = vx[2*i-1] + kx*0.5
    vy[2*i] = vy[2*i-1] + ky*0.5
    
    
    vx[2*i+1] = vx[2*i-1] + kx
    vy[2*i+1] = vy[2*i-1] + ky
    v[i] = np.sqrt(vx[2*i]**2+vy[2*i]**2)
    
    
    Ep[i] = -GM*mt/r[i]
    Ec[i] = 0.5 * mt * v[i]**2

# =============================================================================
#                           GRAFICAS VERLET
# =============================================================================

''' APARTADO a'''

plt.figure(1)
plt.title('Distancia al sol $(r)$ en función del tiempo')
plt.xlabel('$t (años)$')
plt.ylabel('$r (m)$')
plt.plot(t_an,r)
#plt.savefig('Distancia Verlet.png')
plt.show()


'''APARTADO b'''

plt.figure(2)
plt.title('Trayectoria x = f(y) ')
plt.ylabel('$x (m)$')
plt.xlabel('$y (m)$')
plt.plot(y,x,c='r')
#plt.savefig('Trayectoria Verlet.png')
plt.show()


'''APARTADO c'''

plt.figure(3)
plt.title('Energías')
plt.plot(t_an,Ep,c='b')
plt.plot(t_an,Ec,c='y')
plt.plot(t_an,Ec+Ep,c='g')
#plt.savefig('Energias Verlet.png')
plt.show()


# =============================================================================
#                               RUNGE-KUTA 2
# =============================================================================

def RK2():
    #definimos una matriz para rellenar con x e y
    r=np.zeros((N,2))      
    v=np.zeros((N,2))
    #condiciones iniciales
    r[0]=x[0],y[0]
    v[0]=vx[0],vy[0]
    #aplicamos el algoritmo de RK2
    for i in range(N-1):
        k1=h*v[i]
        l1=h*f1(r[i])
        k2=h*(v[i]+l1/2)
        l2=h*f1(r[i]+k1/2)
        
        r[i+1]=r[i]+k2
        v[i+1]=v[i]+l2
        
    R=np.sqrt(r[:,0]**2+r[:,1]**2)
        
    Ec=0.5*mt*(v[:,0]**2+v[:,1]**2)
    Ep=-GM*mt/R
    E=Ec+Ep
    
    return r,v,R,Ec,Ep,E

 # =============================================================================
 #                               RUNGE-KUTA 4
 # =============================================================================
   
                           
def RK4():
    #definimos una matriz para rellenar con x e y
    r=np.zeros((N,2))
    v=np.zeros((N,2))
    #condiciones iniciales
    r[0]=x[0],y[0]
    v[0]=vx[0],vy[0]
    #aplicamos el algoritmo de RK4
    for i in range(N-1):
        k1=h*v[i]
        l1=h*f1(r[i])
        k2=h*(v[i]+l1/2)
        l2=h*f1(r[i]+k1/2)
        k3=h*(v[i]+l2/2)
        l3=h*f1(r[i]+k2/2)
        k4=h*(v[i]+l3)
        l4=h*f1(r[i]+k3)
        
        r[i+1]=r[i]+1/6*(k1+2*k2+2*k3+k4)
        v[i+1]=v[i]+1/6*(l1+2*l2+2*l3+l4)
        
    R=np.sqrt(r[:,0]**2+r[:,1]**2)
        
    Ec=0.5*mt*(v[:,0]**2+v[:,1]**2)
    Ep=-GM*mt/R
    E=Ec+Ep
    
    return r,v,R,Ec,Ep,E

# =============================================================================
#                           GRAFICAS RK2
# =============================================================================

r_RK2,v_RK2,R_RK2,Ec_RK2,Ep_RK2,E_RK2=RK2()

''' APARTADO a'''

plt.figure(1)
plt.title('Distancia al sol $(r)$ en función del tiempo')
plt.xlabel('$t (años)$')
plt.ylabel('$r (m)$')
plt.plot(t_an,R_RK2)
#plt.savefig('Distancia RK2.png')
plt.show()


'''APARTADO b'''

plt.figure(2)
plt.title('Trayectoria x = f(y) ')
plt.ylabel('$x (m)$')
plt.xlabel('$y (m)$')
plt.plot(r_RK2[::2,0],r_RK2[::2,1])
#plt.savefig('Trayectoria RK2.png')
plt.show()


'''APARTADO c'''

plt.figure(3)
plt.title('Energías')
plt.plot(t_an,Ep_RK2,c='b')
plt.plot(t_an,Ec_RK2,c='y')
plt.plot(t_an,E_RK2,c='g')
#plt.savefig('Energias RK2.png')
plt.show()


# =============================================================================
#                           GRAFICAS RK4
# =============================================================================

r_RK4,v_RK4,R_RK4,Ec_RK4,Ep_RK4,E_RK4=RK4()

''' APARTADO a'''

plt.figure(1)
plt.title('Distancia al sol $(r)$ en función del tiempo')
plt.xlabel('$t (años)$')
plt.ylabel('$r (m)$')
plt.plot(t_an,R_RK4)
#plt.savefig('Distancia RK4.png')
plt.show()


'''APARTADO b'''

plt.figure(2)
plt.title('Trayectoria x = f(y) ')
plt.ylabel('$x (m)$')
plt.xlabel('$y (m)$')
plt.plot(r_RK4[::2,0],r_RK4[::2,1])
#plt.savefig('Trayectoria RK4.png')
plt.show()


'''APARTADO c'''

plt.figure(3)
plt.title('Energías')
plt.plot(t_an,Ep_RK4,c='b')
plt.plot(t_an,Ec_RK4,c='y')
plt.plot(t_an,E_RK4,c='g')
#plt.savefig('Energias RK4.png')
plt.show()

########COMPARACION ENERGIAS###########

plt.figure()
plt.title('Energía')
plt.plot(t_an,Ec+Ep,label='$E_{Verlet}$')
plt.plot(t_an,E_RK2,label='$E_{RK2}$')
plt.plot(t_an,E_RK4,label='$E_{RK4}$')
plt.xlabel('t (años)')
plt.ylabel('Energía (J)')
plt.legend(loc='best')
#plt.savefig('Comparacion energias.png')
plt.show()


