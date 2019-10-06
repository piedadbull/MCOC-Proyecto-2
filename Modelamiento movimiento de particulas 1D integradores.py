from matplotlib.pylab import *

#unidades SI

_m= 1.
_kg = 1.
_s = 1.
_mm = 1e-3*_m
_gr = 1e-3*_kg


vfx = 5.0*_m/_s #m/s 
vfy = 0.1 *_m/_s#m/s


x01 = array([[0.,1.*_mm],[0.,1.*_mm]], dtype=double)   # matriz posicion inicial [x,y] particula [1, 2, 3, ..., n]
v01 = array([[1.,1.],[1.,2.]], dtype=double)   # matriz velocidad inicial [vx,vy] particula [1, 2, 3, ..., n]
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

A = pi*(d/2)**2
V = (4./3.)*pi*(d/2)**3
m = rho_particula*V #masa de la particula, grano de arena



dt = 0.001*_s #paso de tiempo
tmax = 2 *_s#tiempo maximo de simulacion
ti = 0.*_s#tiempo actual

W = array([0, -m*g])
fB = array([0, rho_agua*V*g])

t = arange(0,tmax,dt)
Nt = len(t)
#Nt = int32(2*tmax/ dt) #numero de espacios de tiempo, tiene que ser entero porque hare un for

norm = lambda v: sqrt(dot(v,v))

k_penal = lambda v: 1000*0.5*Cd*rho_agua*A*norm(v)/(1*_mm)


def particula(z,t):
	zp = zeros(len(z))
	for p in range((len(z))/4):
		xi = z[p*4:2+p*4]
		vi = z[2+p*4:4+p*4]
		vf = array([vfx,vfy])
		vrel = vf-vi
		fD = (0.5*Cd*rho_agua*norm(vrel)*A)*vrel
		Fi = W + fD + fB

		if xi[1]<0:
			Fi[1] += -k_penal(vi)*xi[1]

		zp[p*4:2+p*4] = vi
		zp[2+p*4:4+p*4] = Fi/m

	return zp


from scipy.integrate import odeint
n = len(x01)   #Numero de particulas
z0 = zeros(4*n)
for p in range(n):
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
    plot(x1[:,0],x1[:,1],label="x"+str(p+1))
    
#plot(x1[:,0],x1[:,1],label="x1")
#plot(x2[:,0],x2[:,1],label="x2")
ylim([0,20*_mm])
plt.legend()


figure()
for p in range(n):
    subplot(2,n,p+1)
    x1 = z[:,p*4:2+p*4]
    plot(t,x1[:,0],label="x"+str(p+1))
    plot(t,x1[:,1],label="y"+str(p+1))
    plt.legend()

for p in range(n):
    subplot(2,n,p+n+1)
    v1 = z[:,2+p*4:4+p*4]
    plot(t,v1[:,0],label="vx"+str(p+1))
    plot(t,v1[:,1],label="vy"+str(p+1))
    plt.legend()

show()



