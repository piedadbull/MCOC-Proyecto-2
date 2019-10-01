# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:52:48 2019

"""

from matplotlib.pylab import *

#unidades base
_m = 1.
_Kg = 1.
_s = 1.
_mm = 1e-3*_m
_gr = 1e-3*_Kg

vfx =5.0 #m/s
vfy =0.0 #m/s

x0 = array([0., 1.], dtype=double)
v0 = array([1., 1.], dtype=double)

xi= x0  #zeros(2,dtype=double)     #posicion actual
vi = v0  #zeros(2,dtype=double)    #velocidad actual
xim1 = zeros(2, dtype=double)         #posicion siguiente
vim1 = zeros(2, dtype=double)         #velocidad siguiente

#masa de la particula
g = 9.81*_m/_s**2
d = 1*_mm
rho = 2700.*_Kg/(_m**3) #suponiendo parecido a una roca
Cd = 0.47 #grag coef; para una particula 
m = rho*(4./3./8.)*pi*(d**3)

#inicializar Euler en x0 

dt = 2e-6*_s         #paso del tiempo
tmax = 1*_s #3*dt  #0.1*_s      #tiempo maximo de simulacion
ti = 0.*_s          #tiempo actual

W = array([0, -m*g])
vf = array([vfx, vfy])

Nt = int32(2*tmax/dt)
x_store = zeros((2,Nt))
v_store = zeros((2,Nt))
t_store = zeros(Nt)

#metodo de euler
i = 0
while ti < tmax:
    
    if i % 100 == 0: #que imprima cada 100 pasos
        print 'ti= ',ti , '|xi|= ', sqrt(dot(xi,xi))
    #print 'xi= ',xi
    #print 'vi= ',vi
    #evaluar v. relativa
    vrel = vf-vi
    norm_vrel = sqrt(dot(vrel,vrel))    #norma de vrel
    
    #evalucar feurzas sobre la particula
    fD=0.5*Cd*norm_vrel*vrel
    Fi = W + fD
   
    print 'Fi= ',Fi
   
    #evaliucar aceleracion
    ai=Fi/m
    print 'ai= ',ai

    #integrar
    xim1 = xi + vi*dt + ai*((dt**2)/2)
    vim1 = vi + ai*dt
    
    #avanzar al siguiente paso
    x_store[:,i] = xi
    v_store[:,i] = vi
    t_store[:] = ti
    ti+=dt
    
    i+=1
    xi = xim1
    vi = vim1

#guardar ultimo paso  
x_store[:,i] = xi   
v_store[:,i] = vi
t_store[i] = ti

print x_store
   
figure()
#plot(x_store0,:i)
plot(x_store[0, :i], x_store[1,:i])
show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    