from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D
import math
import random

#unidades SI
_m= 1.
_kg = 1.
_s = 1.
_mm = 1e-3*_m
_gr = 1e-3*_kg


vfx = 5.0*_m/_s #m/s 
vfy = 0.1 *_m/_s#m/s

n = 4   #Numero de particulas

x01 = zeros((n,2))   #matriz posicion de las particulas
v01 = zeros((n,2))   #matriz velocidad de las particulas
for i in range(n):
	x01[i][:2] = array([0.,double(random.randint(0,10)*_mm)])                         #valores iniciales de la posicion
	v01[i][:2] = array([double(random.randint(0,10)),double(random.randint(0,10))])   # y la velocidad de la particula

print("Condiciones iniciales:")
print("Posiociones =",x01)
print("Velocidades =", v01)


#x01 = array([[0.,1.*_mm],[0.,2*_mm]], dtype=double)   # matriz posicion inicial [x,y] particula [1, 2, 3, ..., n]
#v01 = array([[1.,1.],[1.,2.]], dtype=double)   # matriz velocidad inicial [vx,vy] particula [1, 2, 3, ..., n]
#x02 = array([0.,1.*_mm], dtype=double)
#v02 = array([1.,2.], dtype=double)

#xi = x0   #zeros (2, dtype= double) #posicion actual
#vi = v0   #zeros (2, dtype= double)	#velocidad actual
#xim1 = zeros (2, dtype= double)	#posicion siguiente
#vim1 = zeros (2, dtype= double)	#velocidad siguiente

g = 9.81*_m/_s**2
d = 1 * _mm
rho_particula = 2650.*_kg/(_m**3)
rho_agua = 1000.*_kg/(_m**3)
Cd = 0.47  

A = pi*(d/2)**2   #Area transversal de particula
V = (4./3.)*pi*(d/2)**3   #Volumen de particula
m = rho_particula*V   #masa de la particula, grano de arena



dt = 0.001*_s #paso de tiempo
tmax = 2 *_s#tiempo maximo de simulacion
ti = 0.*_s#tiempo actual

W = array([0, -m*g])   #Vector peso de la particula
fB = array([0, rho_agua*V*g])   #Vetcor Empuje

t = arange(0,tmax,dt)
Nt = len(t)
#Nt = int32(2*tmax/ dt) #numero de espacios de tiempo, tiene que ser entero porque hare un for

norm = lambda v: sqrt(dot(v,v))

k_penal = lambda v: 1000*0.5*Cd*rho_agua*A*norm(v)/(1*_mm)

vfx = lambda xy: 1./0.4*log(30*(xy+1./30))
choque = lambda x1,x2: sqrt((x1[0]-x2[0])**2+(x1[1]-x2[1])**2)   #Distancia entre 2 particulas
def angulo(x1,x2):   #Angulo de choque entre 2 particulas
	if abs(x1[0]-x2[0]) > 0:
		return arctan(abs(x1[1]-x2[1])/abs(x1[0]-x2[0])) 
	else:
		return 0.


def particula(z,t):

	zp = zeros(len(z))
	for p1 in range((len(z))/4):
		xi = z[p1*4:2+p1*4]
		vi = z[2+p1*4:4+p1*4]
		vf = array([vfx(xi[1]),vfy])
		vrel = vf-vi
		fD = (0.5*Cd*rho_agua*norm(vrel)*A)*vrel
		Fi = W + fD + fB

		if xi[1]<0:
			Fi[1] += -k_penal(vi)*xi[1]

		for p2 in range((len(z))/4):   #Revisamos si hay choque entre la particula p1 y alguna otra particula
			if p1 != p2:
				x2 = z[p2*4:2+p2*4]    #Posicion y velocidad de la 2da particula
				v2 = z[2+p2*4:4+p2*4]
				xdif = choque(xi,x2)   #Distancia entre los centros de las particulas
				if xdif < 1 * _mm:     #Revisamos si hay choque revisando si la diferencia de posiciones es menor a 1 diametro de particula
					vdif = vi-v2
					Fct = abs(-k_penal(vdif)*xdif)   #La misma idea de choque contra el suelo
					theta = angulo(xi,x2)            #Angulo con que chocan las particulas
					if xi[0] > x2[0]:                #Revisamos posicion de la particula con respecto a la que choco para ver si
						Fcx = Fct*cos(theta)         # acelera o desacelera en el eje x e y
					else:
						Fcx = -Fct*cos(theta)
					if xi[1] >= x2[1]:
						Fcy = Fct*cos(theta)
					else:
						Fcy = -Fct*cos(theta)
					Fc = [Fcx, Fcy]
					Fi += Fc
					print ("choque")

		zp[p1*4:2+p1*4] = vi      #Guardamos la derivada de la posicion con respecto al tiempo
		zp[2+p1*4:4+p1*4] = Fi/m  #Guardamos la derivada de la velocidad con respecto al tiempo

	return zp   #Retornamos el vector derivada


from scipy.integrate import odeint

z0 = zeros(4*n)      #Vector de condiciones iniciales para la integracion, guardamos de la forma
for p in range(n):   # z0 = [x1,y1,vx1,vy1,x2,y2,vx2,vy2,x3,y3,vx3,vy3,....,xn,yn,vxn,vyn]
    z0[p*4:2+p*4] = x01[p]
    z0[2+p*4:4+p*4] = v01[p]
    
print z0
#z0[:2] = x01
#z0[2:4] = v01
#z0[4:6] = x02
#z0[6:8] = v02


z = odeint(particula, z0, t)


#x1 = z[:,:2]
#v1 = z[:,2:4]
#x2 = z[:,4:6]
#v2 = z[:,6:8]

figure()
for p in range(n):
    x1 = z[:,p*4:2+p*4]
    plot(x1[:,0],x1[:,1],label="p"+str(p+1))
    
#plot(x1[:,0],x1[:,1],label="x1")
#plot(x2[:,0],x2[:,1],label="x2")
ylim([0,50*_mm])
plt.title("Posicion particulas plano XY")
plt.legend()


figure()
for p in range(n):
    subplot(2,n,p+1)
    x1 = z[:,p*4:2+p*4]
    plot(t,x1[:,0],label="x"+str(p+1))
    plot(t,x1[:,1],label="y"+str(p+1))
    plt.title("Particula "+str(p+1)+'\nPosicion')
    plt.legend()

for p in range(n):
    subplot(2,n,p+n+1)
    v1 = z[:,2+p*4:4+p*4]
    plot(t,v1[:,0],label="vx"+str(p+1))
    plot(t,v1[:,1],label="vy"+str(p+1))
    plt.title("Velocidad")
    plt.legend()

#fig = plt.figure()
#ax = Axes3D(fig)
#for p in range(n):
#	x1 = z[:,p*4:2+p*4]
#	ax.plot(t, x1[:,0],x1[:,1]) 

show()


