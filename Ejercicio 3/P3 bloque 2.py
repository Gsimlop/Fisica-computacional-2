import numpy as np
import matplotlib.pyplot as plt
import numpy.random as npr

'''
El objetivo de esta práctica es escribir un programa Monte Carlo para integrar funciones.
'''
###############################################################################
#                               EXTRA Funcion
###############################################################################

def f(x):
    return np.exp(-x**2)

###############################################################################
#                               Funcion pedida
###############################################################################

#Definimos la función que queremos integral
'''
def f(x):
    return np.sin(1/(x*(2-x)))**2
'''


#Ahoradefinimos los limites de integracion 

h,l = 1,0 # maximo y minimo del rectangulo
a,b = -2,2 # intervalo de integracion

N = 10000

#Definimos los puntos aleatorios en el dominio

random_x = np.random.uniform(a, b, N)
random_y = np.random.uniform(l, h, N)

#El valor de la integral vendra dado por el area del rectangulo que creamos mult. por la densidad de pts. que caen dentro

x = np.linspace(l,h,N)
a_rec = (h-l) * (b-a)

#Creamos las condiciones para puntos aceptados y no aceptados

cond = random_y <= f(random_x)
condno = random_y > f(random_x)

#Creamos las listas para los valores aceptados y no aceptados

x_acep, y_acep = random_x[cond], random_y[cond]
x_nacep, y_nacep = random_x[condno], random_y[condno]

sumay = np.sum(cond) #Creamos un array de valores booleanos que entran dentro de la funcion

valor_anal = 0.9953 * 2       #Valor analitico calc.
valor = a_rec * (sumay/N) 
err = np.sqrt((valor*(a_rec-valor))/N)

print('La integral vale ',valor,'+-',err)

#Queremos representar el valor de la integral (con su error) en función del numero de pasos de Monte Carlo
#Calculamos el valor de la integral en funcion de N

n = 10 # valor del intervalo de pasos 

integrales = []
err = []
err_max = []
err_min = []

while n < N:
    
    bajo = 0
    for i in range(0,n):
        if random_y[i]<f(random_x[i]):
            bajo+=1
            
    n_p = bajo/n # el porcentaje de puntos que caen bajo la funcion
    area_n = (h-l)*(b-a) # area del rectangulo
    valor_n = n_p * area_n #el valor del area de la funcion dentro del rectangulo [a,b]x[l,h]
    err_n = np.sqrt((valor_n*(area_n-valor_n))/N)
    
    
    integrales.append(valor_n)
    err.append(err_n)
    err_max.append(valor_n + err_n)
    err_min.append(valor_n - err_n)
    n+=10


#Graficamos apartado 1

t = np.linspace(a,b,N)

plt.figure(1)
plt.title('$f(x)$')
plt.plot(t,f(t))
#plt.savefig('Gráfica f2.png')
plt.show()

#Graficamos apartado 2

t2 = np.linspace(0,N,len(integrales))

plt.figure(2)
plt.title('Convergencia de I al aumentar pasos en el MC ')
plt.plot(t2,integrales,c='black',label = 'I')
plt.plot(t2,err_max,c='red',label = 'I $\pm$ error')
plt.plot(t2,err_min,c='red')
plt.legend()
#plt.savefig('Evolucion valor integral f2.png')
plt.show()

#Graficamos apartado 3

plt.figure()
plt.plot(t2,err)
plt.xlabel('N')
plt.ylabel('error')
plt.ylim(bottom=0.008)
#plt.savefig('Evolucion error f2.png')
plt.show()

#Graficamos apartado 4

plt.figure(3)
plt.title('Puntos generados aleatoriamente(Monte Carlo)')
plt.plot(x_acep,y_acep,'.',c='g',label='puntos aceptados')
plt.plot(x_nacep,y_nacep,'.',c='r',label='puntos no aceptados')
plt.plot(t,f(t),c='black',lw=2,label='función')
plt.legend()
#plt.savefig('f2 y montecarlo.png')
plt.show()


'''
Los calculos se haran menos precisos dependiendo de lo grande que sea la altura
de nuestro rectangulo (si mantenemos el numero de puntos aleatorio)
Esto es porque hay un area mayor y la probabilidad de que cada punto aleatorio
caiga en la region menor que la funcion es menor.

Si aumentamos el numero de pasos del Monte Carlo, el valor de la integral sera
mucho mas preciso (su error disminuira) ya que tenemos un conjunto de datos 
mayor y nuestro metodo se basa en probabilidad
'''

###############################################################################
#                               EXTRA ESFERA
###############################################################################
D = 10
prandom = np.random.uniform(-1,1,(N,D))  #lo centramos en -1,1

radios = np.zeros(N)

for i in range(N):
    radios[i] = np.linalg.norm(prandom[i])**2
    
cond_r = radios <=1

Nesf = np.sum(cond_r)

Vesf = Nesf/N*2**10

print('Volumen de la hiperesfera para ', N, 'pasos =', Vesf)


