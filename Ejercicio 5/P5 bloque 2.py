import numpy as np
import matplotlib.pyplot as plt

"""
Queremos calcular la fuerza con un potencial de pares
"""

###############################################################################
#                               FUNCIONES
###############################################################################

#implementamos el programa que genera la red fcc (ejercicio 1 Bloque II)
def cristal_fcc(nx, ny, nz):
    positions = []
    count = 0
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                positions.append((x, y, z))
                positions.append((x + 0.5, y + 0.5, z))
                positions.append((x, y + 0.5, z + 0.5))
                positions.append((x + 0.5, y, z + 0.5))
                count += 4
    return positions, count

# graficamos el cristal
def plot_3d_positions(positions, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs, ys, zs = zip(*positions)
    ax.scatter(xs, ys, zs, marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title(title)
    plt.show()
   
#Definimos el potencial de L-J

def V(r):
    return 4*eps*( (sigma/r)**n - (sigma/r)**m )

#Escribimos la derivada del potencial de L-J

def derV(r):   
    return 4*eps*( (m*sigma**m)/r**(m+1) - (n*sigma**n)/r**(n+1) ) 

'''Queremos una funcion que calcule la energía potencial entre todas las partículas en un cristal 
en función de sus posiciones espaciales y condiciones de contorno, con ciertas consideraciones de distancia y potencial, y 
almacena estas energías en una lista llamada E'''

def energia(cc, positions , V ,derV, x,y, z, alat): #-> radio de corte : 3 sigmas
    Ei=0
    Ep=0
    E=[]    #para rellenar con las soluciones
    F = np.zeros([len(positions),3])
    F_lista=[]

    if cc == 'F': #SUPERFICIES LIBRES: No imponemos ninguna condición especial a las partículas del borde de la caja de simulación
    #La energía de cada partícula dependerá de su distancia a la superficie
        for i in range(len(positions)): #itera sobre cada partícula i en el cristal 
            Ei = 0
            for j in range(len(positions)):    #itera sobre cada partícula j del cristal para calcular la energía potencial 
                #calculamos las distancias entre partículas
                dx = x[i]-x[j]
                dy = y[i]-y[j]
                dz = z[i]-z[j]
                r = np.sqrt(dx**2+dy**2+dz**2)   #módulo vector posición -> distancia entre partículas
                #tenemos en cuenta el radio de corte
                if r==0:
                    V_ij= 0
                elif r <= 3*sigma:  #condiciones de corto
                    V_ij=V(r)
                    Ep += 0.5*V_ij #añado el potencial a la energia
                    #añado la fuerza ejercida por cada átomo
                    F[i][0] += -0.5* derV(r)*dx/r
                    F[i][1] += -0.5* derV(r)*dy/r
                    F[i][2] += -0.5* derV(r)*dz/r
                    Ft = np.round(np.sqrt(F[i][0]**2+F[i][1]**2+F[i][2]**2),decimals=2)

                else:
                    V_ij = 0
                    
                Ei +=0.5*V_ij 
                
            F_lista.append(Ft)
            E.append(Ei)    #añadimos la energía individual de cada partícula a la lista E
        return np.sum(E), E , F_lista, F
        
#CONDICIONES PERIODICAS: cuando una partícula está cerca de los límites de la caja de simulación consideramos las partículas vecinas que podrían estar en la caja repetida  
    if cc == 'P':
        
        for i in range(len(positions)):
            
            Ei = 0
            for j in range(len(positions)):
                #diferencia entre coordenadas de todas las partículas respecto de la partícula i-ésima
                dx = x[i]-x[j]
                dy = y[i]-y[j]
                dz = z[i]-z[j]
                
                #Si una partícula está más allá de la mirad de la caja en alguna de las tres dimensiones
                #se ajusta para que esté en la caja principal sumando la longitud total de la caja
                #ajuste para condiciones periódicas en x
                if dx >0:
                    if dx > nx*alat/2:
                        dx-=nx*alat
                if dx<0:
                    if dx < -nx*alat/2:
                        dx+=nx*alat 
                
                #ajuste para condiciones periódicas en y
                if dy>0:
                    if dy > ny*alat/2:
                        dy-=ny*alat
                if dy<0:
                    if dy < -ny*alat/2:
                        dy+=ny*alat

                #ajuste para condiciones periódicas en z    
                if dz>0:
                    if dz > nz*alat/2:
                        dz-=nz*alat
                if dz<0:
                    if dz < -nz*alat/2:
                        dz+=nz*alat


                r = np.sqrt(dx**2+dy**2+dz**2)  #módulo distancia
                
                if r==0:
                    V_ij= 0
                elif r <= 3*sigma:
                    V_ij=V(r)
                    Ep += 0.5*V_ij #añado el potencial a la energia
                    #añado la fuerza ejercida por cada átomo
                    F[i][0] += -0.5* derV(r)*dx/r
                    F[i][1] += -0.5* derV(r)*dy/r
                    F[i][2] += -0.5* derV(r)*dz/r
                    
                    Ft = np.round(np.sqrt(F[i][0]**2+F[i][1]**2+F[i][2]**2),decimals=2)

                else:
                    V_ij = 0
                    
                Ei += 0.5*V_ij
            
            F_lista.append(Ft)
            E.append(np.round(Ei,decimals=13))
            
        return np.sum(E), E, F_lista, F

        
    else:
        print('Condiciones de contorno incorrectas.')


###############################################################################
#                               PARÁMETROS
###############################################################################


#Definimos los parametros que vamos a necesitar para simular el cobre
sigma = 2.3151 #Distancia (finita) en la que el potencial entre partículas es cero en armstrong
n = 12.   #n y m para potencial L-J
m = 6.
eps = 0.167  #Profundidad del potencial en eV
alat = 3.603 #lattice parameter (armstrong)

###############################################################################
#                               INPUTS
###############################################################################

nx = int(input("Número de celdas unidad en la dirección x: "))
ny = int(input("Número de celdas unidad en la dirección y: "))
nz = int(input("Número de celdas unidad en la dirección z: "))
cc = input("Condiciones de contorno (F o P): ")

###############################################################################
#                        CÁLCULO ENERGÍA Y FUERZA
###############################################################################

#llamamos a la función y a la gráfica
positions, count=cristal_fcc(nx,ny,nz)
plot_3d_positions(positions, 'Red Cristalina (fcc).')

#listas de posiciones de las partículas
x = []
y = []
z = []

for i in range(len(positions)):
    x.append(alat*positions[i][0])
    y.append(alat*positions[i][1])
    z.append(alat*positions[i][2])


energia, E , F_lista, F = energia(cc, positions, V,derV, x,y, z, alat)
N = [nx,ny,nz]
nAt = int(4*N[0]*N[1]*N[2])

###############################################################################
#                               OUTPUTS
###############################################################################

print('La energia potencial es', np.round(np.sum(E),4))
print('Energia por atomo = ', np.round(np.sum(E)/nAt),4)

###############################################################################
#                               GRAFICAS
###############################################################################

#Grafica de las energias

fig = plt.figure()  #crea una figura en blanco
ax = fig.add_subplot(111,projection='3d')   #gráfico tridimensional

col = ax.scatter(x, y, z,c=E,cmap='magma',lw=5) #grafico de dispersión 3D. el color de cada punto se determina
                                                #según los valores de E
fig.colorbar(col,label='Energía potencial')     #añado una barra de color y una etiqueta
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()


#Grafica de las fuerzas

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

col = ax.scatter(x, y, z,c=F_lista,cmap='magma',lw=4)
fig.colorbar(col,label='Fuerza')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.savefig('Fuerza P D2.png')
plt.show()


#Grafica de los vectores de las fuerzas

#Creamos las listas de los vectores
U = []
W = []
H = []
#Añadimos valores
for i in range(len(F)):
    U.append(F[i,0])  
    W.append(F[i,1])
    H.append(F[i,2])

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
#col = ax.scatter(x, y, z,c=F_lista,cmap='inferno',lw=3)
ax.quiver(x, y, z,U,W,H,color='purple', normalize = True)
#fig.colorbar(col,label='Fuerza')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.savefig(' Vectores Fuerza P D2.png')
plt.show()




