import numpy as np
import matplotlib.pyplot as plt
from random import random,randrange,choice

# =============================================================================
#                                INPUTS
# =============================================================================

N = int(input('Número de partículas '))
T = float(input('Temperatura '))
pas = int(input('Pasos del Monte Carlo '))

'''
Queremos crear un programa el cual eliga de una particula de forma aleatoria y 
la direccion en x,y,z y realice el cambio de +1 o -1, el cual aceptamos o rechazamos
considerando el algoritmo de metropolis
'''

n = np.ones([N,3],int) #Creamos nx,ny,nz para los atomos deseados
Et = 3/2*N*np.pi**2    #Energia para N particulas no interactuantes 
E = []                 #Creamos lista de energias a la que añadimos los valores aceptados

for paso in range(pas):
    
    i = randrange(N)
    j = randrange(3)
    k = choice([0,1])   #Elige entra la probabilidad que le damos

    if k == 0.: 
        dn = 1
        dE = ((np.pi**2)/(2))*(2*n[i,j]+1)
    if k == 1.: 
        dn = -1
        dE = ((np.pi**2)/(2))*(-2*n[i,j]+1)
    

    if n[i,j]>1 or dn==1:
        if random() < np.exp(-dE/T):
            n[i,j] += dn
            Et+=dE
        E.append(Et)


plt.figure()
t = np.linspace(0,pas,len(E))
plt.plot(t,E)
plt.xlabel('Pasos del Monte Carlo')
plt.ylabel('Energía')
#plt.savefig('N = 10')
plt.show()

# =============================================================================
#                   EXTRA(EXPOSICION NO PROGRAMA FINAL)
# =============================================================================

'''
Como extra utilizaré el método montecarlo para conseguir que el Cu (cuya estructura es la fcc)
se organice en una fcc, en donde supuestamente tendrá su mínimo de energía.
Es decir, creo un programa el cual nos pregunte el número de particulas que queremos soltar,
el lado del cubo de la simulación, la temperatura a la que queremos trabajar y el numero
de pasos de Montecarlo deseados.

Utilizaremos el potencial de L-J usado en las practicas anteriores y necesitaremos definir el parametro
de red del Cu.

De esta esta manera, cabe esperar que, dado que la estructura de menor energia para el Cu es la fcc,
obtengamos efectivamente una fcc.
'''

###############################################################################
#                               FUNCIONES
###############################################################################
def LennardJones(r):
    sigma = 2.3151; eps = 0.167; n = 12; m = 6
    return 4*eps*( (sigma/r)**n - (sigma/r)**m )


def CalcR(r):
    n = len(r)
    R = np.empty((n, n))
    for i in range(n):
        R[i] = np.sqrt( (r[:,0]-r[:,0][i])**2 + (r[:,1]-r[:,1][i])**2 + (r[:,2]-r[:,2][i])**2)
    return R

def CalcV(r, factor_de_corte=3):
    sigma=2.3151
    R = CalcR(r)
    R_filtered = np.logical_and(R < factor_de_corte * sigma, R > 0)
    V_red = 0.5 * np.sum(LennardJones(R[R_filtered]))
    return V_red

def random_pick_vect(N, dim):
    v = np.random.random_sample(dim)-0.5
    v = v/np.sqrt(np.sum(v**2))
    return  np.random.randint(0,N), v

def random_move(r, step, dim, L):
    while True:
        i, vect = random_pick_vect(len(r), dim)
        if (r[i]+vect*step < [L,L,L]).any() or (r[i]+vect*step > [0,0,0]).any():
             r[i] = r[i]+vect*step
             return r

def Metropolis(r, T, E_i):
    global kb, dim, step
    rj = random_move(r.copy(), step, dim,L)
    E_j = CalcV(rj)
    if E_j <= E_i:
        Pr = 1
    else:
        Pr = np.exp(-(E_j-E_i)/(kb*T))
    if np.random.choice((True, False), p = [Pr, 1-Pr]):
        return rj, E_j
    return r, E_i


###############################################################################
#                               PARÁMETROS
###############################################################################
kb = 8.6181024e-5
dim = 3
a = 3.603
step = a/100
L = 2*a
sigma=2.3151
eps=0.167
n=12
m=6

###############################################################################
#                               INPUTS
###############################################################################
L = a*int(float(input('Lado del cubo de simulación:   ')))
n = int(float(input('Defina el número de partículas:   ')))
T = float(input('Defina la temperatura del sistema (Kelvin):    '))
n_p = int(float(input('Defina el número de pasos de montecarlo necesarios:   ')))
r = a*np.random.random_sample((n,dim))

###############################################################################
#                               GRÁFICAS
###############################################################################

Energy = np.empty(n_p+1)
Energy[0] = CalcV(r)


fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Posición inicial.')
ax.scatter(r[:,0], r[:,1], r[:,2], c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

marker = int(n_p/20)
print('Progreso de la simulación: |', end = '')
for i in range(n_p):
    r, Energy[i+1] = Metropolis(r, T, Energy[i])
    if i%marker == 0:
        print('<3', end = '')
print('|  Completada')
Energy /= n
print(f'Energía final:  {Energy[-1]}')

fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
ax.set_title('Posición final')
ax.scatter(r[:,0], r[:,1], r[:,2], c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

plt.figure()
plt.title(f'Metrópolis con {n_p} pasos, {n} partículas y T = {T}')
plt.plot(range(n_p+1), Energy, color = 'deepskyblue', linewidth = 0.7)
plt.xlabel('Nº de pasos')
plt.ylabel('Energía')
plt.show()













