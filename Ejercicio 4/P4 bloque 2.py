import numpy as np
import matplotlib.pyplot as plt

"""
Queremos generar una red fcc y calcular la energia total del sistema
considerando un potencial de interaccion de Lennard-Jones
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

def V_LJ(r, eps, sigma):
    return 4*eps*( (sigma/r)**n - (sigma/r)**m )


'''Queremos una funcion que calcule la energía potencial entre todas las partículas en un cristal 
en función de sus posiciones espaciales y condiciones de contorno, con ciertas consideraciones de distancia y potencial, y 
almacena estas energías en una lista llamada E'''

def energia(cc, positions , V , x,y, z, alat): #-> radio de corte : 3 sigmas
    Ei=0
    Ep=0
    E=[]    #para rellenar con las soluciones
    if cc == 'F': #SUPERFICIES LIBRES: No imponemos ninguna condición especial a las partículas del borde de la caja de simulación
    #La energía de cada partícula dependerá de su distancia a la superficie
        for i in range(len(positions)): #itera sobre cada partícula i en el cristal 
            Ei = 0
            for j in range(len(positions)):    #itera sobre cada partícula j del cristal para calcular la energía potencial 
                #calculamos las distancias entre partículas
                dx = x[j]-x[i]
                dy = y[j]-y[i]
                dz = z[j]-z[i]
                r = np.sqrt(dx**2+dy**2+dz**2)   #módulo vector posición -> distancia entre partículas
                #tenemos en cuenta el radio de corte
                if r==0:
                    V_ij= 0
                elif r <= 3.5*sigma:  #condiciones de corto
                    V_ij=V(r,eps,sigma)
                    Ep += 0.5*V_ij #añado el potencial a la energia
                else:
                    V_ij = 0
                    
                Ei +=0.5*V_ij  
                
            E.append(Ei)    #añadimos la energía individual de cada partícula a la lista E
        return np.sum(E), E
        
#CONDICIONES PERIODICAS: cuando una partícula está cerca de los límites de la caja de simulación consideramos las partículas vecinas que podrían estar en la caja repetida  
    if cc == 'P':
        
        for i in range(len(positions)):
            
            Ei = 0
            for j in range(len(positions)):
                #diferencia entre coordenadas de todas las partículas respecto de la partícula i-ésima
                dx = x[j]-x[i]
                dy = y[j]-y[i]
                dz = z[j]-z[i]
                
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
                elif r <= 3.5*sigma:
                    V_ij=V(r,eps,sigma)
                    Ep += 0.5*V_ij #añado el potencial a la energia
                else:
                    V_ij = 0
                    
                Ei += 0.5*V_ij
                
            E.append(Ei)
            
        return np.sum(E), E

        
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
#                           CÁLCULO ENERGÍA
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

#utilizando el de Lennard-Jonnes
energia, E = energia(cc, positions, V_LJ, x,y, z, alat)
#utilizando el potencial de 
N = [nx,ny,nz]
nAt = int(4*N[0]*N[1]*N[2])

###############################################################################
#                               OUTPUTS
###############################################################################

print('La energia potencial es',np.round(np.sum(E),decimals=4))
print('Energia por atomo = ',np.round(np.sum(E)/nAt, decimals=4))

###############################################################################
#                               GRAFICAS
###############################################################################

    #Energia por partícula
fig = plt.figure()  #crea una figura en blanco
ax = fig.add_subplot(111,projection='3d')   #gráfico tridimensional

col = ax.scatter(x, y, z,c=E,cmap='magma',lw=5) #grafico de dispersión 3D. el color de cada punto se determina
                                                #según los valores de E
fig.colorbar(col,label='Energía potencial')     #añado una barra de color y una etiqueta
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.title('Energía')
#plt.savefig('Energía NPerio D2.png')
plt.show()

#Potencial de Lennard-Jones
lennard_jones= plt.figure()
r=np.linspace(0.9, 3, 1000)
# Parámetros del potencial
epsilon = 1.0  # Profundidad del potencial
sigma1 = 1.0    # Distancia finita a la cual el potencial es nulo
VLJ=V_LJ(r,epsilon,sigma1)
plt.plot(r,VLJ)   
plt.xlabel('Distancia (r)')
plt.ylabel('Potencial')
plt.title('Potencial de Lennard-Jones')
#plt.savefig('Potencial L-J.png')
plt.show()

