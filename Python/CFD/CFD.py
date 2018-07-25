# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 09:50:13 2016

@author: Alonyan
"""

#%% importing libraries, 
import numpy                       #here we load numpy
from matplotlib import pyplot      #here we load matplotlib
import time, sys   
#matplotlib inline 


#%% Linear convection
nx = 41  # X steps
dx = 2/(nx-1) #X diff
nt = 25    #t steps
dt = .025  # dt
c = 1      # wavespeed

u = numpy.ones(nx)      #numpy function ones()
u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s Initial conditions
#print(u)

pyplot.plot(numpy.linspace(0,2,nx), u); #draw initial conditions

un = numpy.ones(nx)

for n in range(nt):
    un=u.copy()
    for i in range(1,nx): ## you can try commenting this line and...
        u[i] = un[i]-c*dt/dx*(un[i]-un[i-1]) #difference equation for convection    
pyplot.plot(numpy.linspace(0,2,nx),u);
#%% Nonlinear convection
nx = 41  # X steps
dx = 2/(nx-1) #X diff
nt = 50    #t steps
dt = .01  # dt
c = 1      # wavespeed

u = numpy.ones(nx)      #numpy function ones()
u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s Initial conditions
#print(u)

pyplot.plot(numpy.linspace(0,2,nx), u); #draw initial conditions

un = numpy.ones(nx)

for n in range(nt):
    un=u.copy()
    for i in range(1,nx): ## you can try commenting this line and...
        u[i] = un[i]-un[i]*dt/dx*(un[i]-un[i-1]) #difference equation for convection    
pyplot.plot(numpy.linspace(0,2,nx),u);


#%% Diffusion
D = 1      # diffusion
nx = 41  # X steps
dx = 2/(nx-1) #X diff
nt = 50    #t steps
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = abs(sigma*dx**2/D) #dt is defined using sigma ... more later!


u = numpy.ones(nx)      #numpy function ones()
u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s Initial conditions
#print(u)

pyplot.plot(numpy.linspace(0,2,nx), u); #draw initial conditions

un = numpy.ones(nx)

for n in range(nt):
    un=u.copy()
    for i in range(1,nx-1): ## you can try commenting this line and...
        u[i] = un[i]+D*dt/dx**2*(un[i+1]-2*un[i]+un[i-1]) #difference equation for convection    
pyplot.plot(numpy.linspace(0,2,nx),u);

#%% Burgers equation

nu = 10      # diffusion
nx = 41  # X steps
dx = 2/(nx-1) #X diff
nt = 200    #t steps
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = abs(sigma*dx**2/nu) #dt is defined using sigma ... more later!


u = numpy.ones(nx)      #numpy function ones()
u[.5/dx : 1/dx+1]=2  #setting u = 2 between 0.5 and 1 as per our I.C.s Initial conditions
#print(u)

pyplot.plot(numpy.linspace(0,2,nx), u); #draw initial conditions

un = numpy.ones(nx)

for n in range(nt):
    un=u.copy()
    for i in range(1,nx-1): ## you can try commenting this line and...
        u[i] = un[i]+nu*dt/dx**2*(un[i+1]-2*un[i]+un[i-1])-un[i]*dt/dx*(un[i]-un[i-1]) #difference equation for convection    
pyplot.plot(numpy.linspace(0,2,nx),u);


#%% using sympy
import sympy
from sympy import init_printing
init_printing(use_latex=True) #latex rendering
x, nu, t = sympy.symbols('x nu t') #symbolic representation of x nu t
phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + sympy.exp(-(x-4*t-2*numpy.pi)**2/(4*nu*(t+1)))
phiprime = phi.diff(x)
phiprime

u = -2*nu*(phiprime/phi)+4 #initial condition

from sympy.utilities.lambdify import lambdify
ufunc = lambdify((t, x, nu), u) #lambdify takes a symbolic expression and turns it into a callable


###variable declarations
nx = 101
nt = 100
dx = 2*numpy.pi/(nx-1)
nu = 0.05
dt = dx*nu

x = numpy.linspace(0, 2*numpy.pi, nx)
u = numpy.empty(nx)
un = numpy.empty(nx)
t = 0

u = numpy.asarray([ufunc(t, x0, nu) for x0 in x]) # this is all analytical

pyplot.figure(figsize=(6,4), dpi=100)
pyplot.plot(x,u, marker='o', lw=2)
pyplot.xlim([0,2*numpy.pi])
pyplot.ylim([0,10]);



#%% Burgers equation
u = numpy.asarray([ufunc(t, x0, nu) for x0 in x])
for n in range(nt):
    un = u.copy()
    for i in range(-1, nx-1):
        u[i] = un[i] - un[i] * dt/dx *(un[i] - un[i-1]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])
#        u[0] = un[0] - un[0] * dt/dx * (un[0] - un[-1]) + nu*dt/dx**2*\
#                (un[1]-2*un[0]+un[-1])
#        u[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
#                (un[0]-2*un[-1]+un[-2])

        
u_analytical = numpy.asarray([ufunc(nt*dt, xi, nu) for xi in x])

pyplot.figure(figsize=(6,4), dpi=100)
pyplot.plot(x,u, marker='o', lw=2, label='Computational')
pyplot.plot(x, u_analytical, label='Analytical')
pyplot.xlim([0,2*numpy.pi])
pyplot.ylim([0,10])
pyplot.legend();

#%% 2D Linear Convection
from mpl_toolkits.mplot3d import Axes3D

nx = 81
ny = 81
nt = 100
c = 1
dx = 2/(nx-1)
dy = 2/(ny-1)
sigma = .2
dt = sigma*dx

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx)) ##create a 1xn vector of 1's
un = numpy.ones((ny,nx)) ##

###Assign initial conditions

u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

###Plot Initial Condition
fig = pyplot.figure(figsize=(6,4), dpi=100)          ##the figsize parameter can be used to produce different sized images
ax = fig.gca(projection='3d')                      
X, Y = numpy.meshgrid(x,y)                            
surf = ax.plot_surface(X,Y,u[:])

#%% advance in time

u = numpy.ones((ny,nx))
u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2

for n in range(nt+1): ##loop across number of time steps
    un = u.copy()
    row, col = u.shape
    u[1:,1:]=un[1:,1:]-c*dt/dx*(un[1:,1:]-un[1:,:-1])-c*dt/dy*(un[1:,1:]-un[:-1,1:])
    u[0,:] = 1 #bc
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1

fig = pyplot.figure(figsize=(6,4), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X,Y,u[:])



#%% 2D Non-Linear convection
from mpl_toolkits.mplot3d import Axes3D

nx = 81
ny = 81
nt = 100
c = 1
dx = 2/(nx-1)
dy = 2/(ny-1)
sigma = .2
dt = sigma*dx

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx)) ##create a 1xn vector of 1's
v = numpy.ones((ny,nx))
un = numpy.ones((ny,nx)) ##
vn = numpy.ones((ny,nx))

###Assign initial conditions

u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v[.5/dy:1/dy+1,.5/dx:1/dx+1]=2


#% advance in time

u = numpy.ones((ny,nx))
u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2

for n in range(nt+1): ##loop across number of time steps
    un = u.copy()
    vn = v.copy()
    row, col = u.shape
    u[1:,1:]=un[1:,1:]-vn[1:,1:]*dt/dx*(un[1:,1:]-un[1:,:-1])-un[1:,1:]*dt/dy*(un[1:,1:]-un[:-1,1:])
    v[1:,1:]=vn[1:,1:]-un[1:,1:]*dt/dx*(vn[1:,1:]-vn[1:,:-1])-vn[1:,1:]*dt/dy*(vn[1:,1:]-vn[:-1,1:])

    u[0,:] = 1 #bc
    u[-1,:] = 1
    u[:,0] = 1
    u[:,-1] = 1
    
    v[0,:] = 1 #bc
    v[-1,:] = 1
    v[:,0] = 1
    v[:,-1] = 1

fig = pyplot.figure(figsize=(6,4), dpi=100)
ax = fig.gca(projection='3d')
from matplotlib import cm
surf2 = ax.plot_surface(X,Y,u[:],cmap = cm.viridis)



#%% 2D diffusion
import numpy
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D ##library for 3d projection plots
from matplotlib import cm ##cm = "colormap" for changing the 3d plot color palette
%matplotlib inline

###variable declarations
nx = 31
ny = 31

nu=.05
dx = 2/(nx-1)
dy = 2/(ny-1)
sigma = .25
dt = sigma*dx*dy/nu

x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)

u = numpy.ones((ny,nx)) ##create a 1xn vector of 1's
un = numpy.ones((ny,nx)) ##

###Assign initial conditions

u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

fig = pyplot.figure()
ax = fig.gca(projection='3d')
X,Y = numpy.meshgrid(x,y)
surf = ax.plot_surface(X,Y,u, rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=False)
ax.set_xlim(0,2)
ax.set_ylim(0,2)
ax.set_zlim(1,2.5);

#%% Define 2D diffusion function
def diffuse(nt, u):
    
    for n in range(nt+1): 
        un = u.copy()
        u[1:-1,1:-1]=un[1:-1,1:-1]+nu*dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
                                    nu*dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])
        u[0,:]=1
        u[-1,:]=1
        u[:,0]=1
        u[:,-1]=1

    
    fig = pyplot.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X,Y,u[:], rstride=1, cstride=1, cmap=cm.viridis,
        linewidth=0, antialiased=True)
    ax.set_zlim(1,2.5)
#%%
    f = numpy.ones((ny,nx)) ##create a 1xn vector of 1's
    f[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

    diffuse(10, f)

#%% Define 2D burgers function, for fixed edge boundaries
def burgers(nt, u, v):
    
    for n in range(nt+1): 
        un = u.copy()
        vn = v.copy()
        u[1:-1,1:-1]=un[1:-1,1:-1]+nu*dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
                                    nu*dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])-\
                                    vn[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])-un[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[:-2,1:-1])
                                    
        v[1:-1,1:-1]=vn[1:-1,1:-1]+nu*dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+\
                                    nu*dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])-\
                                    un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])
        
        
        u[0,:]=1
        u[-1,:]=1
        u[:,0]=1
        u[:,-1]=1
        v[0,:]=1
        v[-1,:]=1
        v[:,0]=1
        v[:,-1]=1

    

#%%    
def plot2D(x, y, p):
    fig = pyplot.figure(figsize=(6,4), dpi=100)
    ax = fig.gca(projection='3d')
    X,Y = numpy.meshgrid(x,y)
    surf = ax.plot_surface(X,Y,p[:], rstride=1, cstride=1, cmap=cm.viridis,
            linewidth=0, antialiased=False)
    ax.set_xlim(0,2)
    ax.set_ylim(0,2)
    ax.view_init(30,225)
#%%
###variable declarations
nx = 41
ny = 41
c = 10
dx = 2/(nx-1)
dy = 2/(ny-1)
sigma = .0009
nu = 0.01
dt = sigma*dx*dy/nu


x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
X,Y = numpy.meshgrid(x,y)

u = numpy.ones((ny,nx)) ##create a 1xn vector of 1's
v = numpy.ones((ny,nx))
comb = numpy.ones((ny,nx))
u[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v[.5/dy:1/dy+1,.5/dx:1/dx+1]=2 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2


plot2D(x,y,u)
#%%
burgers(500, u, v)
#%%
plot2D(x,y,u)



#%% 2D Laplace equation
def Laplace2D(p,y,dx,dy,Tol):
    Cost = 1;
    pn = numpy.empty_like(p)
    
    while Cost>Tol:
        pn = p.copy()
        p[1:-1, 1:-1] = (dy**2*(pn[1:-1,2:]+pn[1:-1,0:-2])+\
                            dx**2*(pn[2:,1:-1]+pn[0:-2,1:-1]))/(2*(dx**2+dy**2)) 
                
        p[:,0] = 0       ##p = 0 @ x = 0
        p[:,-1] = y      ##p = y @ x = 2
        p[0,:] = p[1,:]  ##dp/dy = 0 @ y = 0
        p[-1,:] = p[-2,:]    ##dp/dy = 0 @ y = 1
        Cost = (numpy.sum(numpy.abs(p[:])-numpy.abs(pn[:])))\
                        /numpy.sum(numpy.abs(pn[:]))
         
    return p
        
#%%
    ##variable declarations
nx = 21
ny = 21
c = 1
dx = 2/(nx-1)
dy = 2/(ny-1)


##initial conditions
p = numpy.zeros((ny,nx)) ##create a XxY vector of 0's


##plotting aids
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,1,ny)



#%%
p = Laplace2D(p, y, dx, dy, 1e-4)
plot2D(x, y, p)






#%% 2D Poisson equation
def Poisson2D(p,b,dx,dy,Tol):
    Cost = 1;
    pn = numpy.empty_like(p)
    
    while Cost>Tol:
        pn = p.copy()
        p[1:-1, 1:-1] = (dy**2*(pn[1:-1,2:]+pn[1:-1,0:-2])+\
                dx**2*(pn[2:,1:-1]+pn[0:-2,1:-1]))/(2*(dx**2+dy**2))-\
                    (b[1:-1,1:-1]*dx**2*dy**2)/(2*(dx**2+dy**2)) 
                
        p[:,0] = 0       ##p = 0 @ x = 0
        p[:,-1] = 0      ##p = y @ x = 2
        p[0,:] = 0  ##dp/dy = 0 @ y = 0
        p[-1,:] = 0    ##dp/dy = 0 @ y = 1
        Cost = (numpy.sum(numpy.abs(p[:])-numpy.abs(pn[:])))\
                        /numpy.sum(numpy.abs(pn[:]))
         
    return p
    
#%%
    # Parameters
nx = 50
ny = 50
nt  = 100
xmin = 0
xmax = 2
ymin = 0
ymax = 2

dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)

# Initialization
p  = numpy.zeros((ny,nx))
b  = numpy.zeros((ny,nx))
x  = numpy.linspace(xmin,xmax,nx)
y  = numpy.linspace(xmin,xmax,ny)

# Source
b[3*ny/4,nx/4]  = -100
b[ny/4,3*nx/4] = 100
#%%
plot2D(x, y, -b)
#%%
p = Poisson2D(p,b,dx,dy,1e-4)
plot2D(x, y, p)



#%% Cavity flow

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy
%matplotlib inline

nx = 41
ny = 41
nt = 500
nit=50
c = 1
dx = 2/(nx-1)
dy = 2/(ny-1)
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
X,Y = numpy.meshgrid(x,y)

rho = 1
nu = .1
dt = .001

u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx)) 
b = numpy.zeros((ny, nx))

def buildUpB(b, rho, dt, u, v, dx, dy):
    
    b[1:-1,1:-1]=rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))-\
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2-\
                      2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))-\
                      ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)

    return b
    
def presPoisson(p, dx, dy, b):
    pn = numpy.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/\
                        (2*(dx**2+dy**2)) -\
                        dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

        p[:,-1] =p[:,-2] ##dp/dy = 0 at x = 2
        p[0,:] = p[1,:]  ##dp/dy = 0 at y = 0
        p[:,0]=p[:,1]    ##dp/dx = 0 at x = 0
        p[-1,:]=0        ##p = 0 at y = 2
        
    return p
    
def cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))
    
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = buildUpB(b, rho, dt, u, v, dx, dy)
        p = presPoisson(p, dx, dy, b)
        
        u[1:-1,1:-1] = un[1:-1,1:-1]-\
                        un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
                        vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
                        dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])+\
                        nu*(dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
                        dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))

        v[1:-1,1:-1] = vn[1:-1,1:-1]-\
                        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
                        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
                        dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+\
                        nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+\
                        (dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])))

        u[0,:] = 0
        u[:,0] = 0
        u[:,-1] = 0
        u[-1,:] = 3    #set velocity on cavity lid equal to 1
        v[0,:] = 0
        v[-1,:]= 0
        v[:,0] = 0
        v[:,-1] = 0
        
        
    return u, v, p
#%%    
u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx))
b = numpy.zeros((ny, nx))
nt = 700
u, v, p = cavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
fig = pyplot.figure(figsize=(6,4), dpi=100)
pyplot.contourf(X,Y,p,alpha=0.5)    ###plnttong the pressure field as a contour
pyplot.colorbar()
pyplot.contour(X,Y,p)               ###plotting the pressure field outlines
pyplot.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
#pyplot.quiver(X,Y,u,v) ##plotting velocity
pyplot.xlabel('X')
pyplot.ylabel('Y');





#%%
def buildUpB(rho, dt, dx, dy, u, v):
    b = numpy.zeros_like(u)
    b[1:-1,1:-1]=rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))-\
                      ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2-\
                      2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))-\
                      ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)
    
    ####Periodic BC Pressure @ x = 2
    b[1:-1,-1]=rho*(1/dt*((u[1:-1,0]-u[1:-1,-2])/(2*dx)+(v[2:,-1]-v[0:-2,-1])/(2*dy))-\
                    ((u[1:-1,0]-u[1:-1,-2])/(2*dx))**2-\
                    2*((u[2:,-1]-u[0:-2,-1])/(2*dy)*(v[1:-1,0]-v[1:-1,-2])/(2*dx))-\
                    ((v[2:,-1]-v[0:-2,-1])/(2*dy))**2)

    ####Periodic BC Pressure @ x = 0
    b[1:-1,0]=rho*(1/dt*((u[1:-1,1]-u[1:-1,-1])/(2*dx)+(v[2:,0]-v[0:-2,0])/(2*dy))-\
                    ((u[1:-1,1]-u[1:-1,-1])/(2*dx))**2-\
                    2*((u[2:,0]-u[0:-2,0])/(2*dy)*(v[1:-1,1]-v[1:-1,-1])/(2*dx))-\
                    ((v[2:,0]-v[0:-2,0])/(2*dy))**2)
    
    return b

def presPoissPeriodic(p, dx, dy):
    pn = numpy.empty_like(p)
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1,1:-1] = ((pn[1:-1,2:]+pn[1:-1,0:-2])*dy**2+(pn[2:,1:-1]+pn[0:-2,1:-1])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,1:-1]

        ####Periodic BC Pressure @ x = 2
        p[1:-1,-1] = ((pn[1:-1,0]+pn[1:-1,-2])*dy**2+(pn[2:,-1]+pn[0:-2,-1])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,-1]

        ####Periodic BC Pressure @ x = 0
        p[1:-1,0] = ((pn[1:-1,1]+pn[1:-1,-1])*dy**2+(pn[2:,0]+pn[0:-2,0])*dx**2)/\
            (2*(dx**2+dy**2)) -\
            dx**2*dy**2/(2*(dx**2+dy**2))*b[1:-1,0]
        
        ####Wall boundary conditions, pressure
        p[-1,:] =p[-2,:]     ##dp/dy = 0 at y = 2
        p[0,:] = p[1,:]      ##dp/dy = 0 at y = 0
    
    return p

##variable declarations
nx = 41
ny = 41
nt = 10
nit=50 
c = 1
dx = 2/(nx-1)
dy = 2/(ny-1)
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
X,Y = numpy.meshgrid(x,y)


##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions
u = numpy.zeros((ny,nx)) ##create a XxY vector of 0's
un = numpy.zeros((ny,nx)) ##create a XxY vector of 0's

v = numpy.zeros((ny,nx)) ##create a XxY vector of 0's
vn = numpy.zeros((ny,nx)) ##create a XxY vector of 0's

p = numpy.ones((ny,nx)) ##create a XxY vector of 0's
pn = numpy.ones((ny,nx)) ##create a XxY vector of 0's

b = numpy.zeros((ny,nx))

#%%
udiff = 1
stepcount = 0

while udiff > .001:
    un = u.copy()
    vn = v.copy()

    b = buildUpB(rho, dt, dx, dy, u, v)
    p = presPoissPeriodic(p, dx, dy)

    u[1:-1,1:-1] = un[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(un[1:-1,1:-1]-un[1:-1,0:-2])-\
        vn[1:-1,1:-1]*dt/dy*(un[1:-1,1:-1]-un[0:-2,1:-1])-\
        dt/(2*rho*dx)*(p[1:-1,2:]-p[1:-1,0:-2])+\
        nu*(dt/dx**2*(un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])+\
        dt/dy**2*(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1]))+F*dt

    v[1:-1,1:-1] = vn[1:-1,1:-1]-\
        un[1:-1,1:-1]*dt/dx*(vn[1:-1,1:-1]-vn[1:-1,0:-2])-\
        vn[1:-1,1:-1]*dt/dy*(vn[1:-1,1:-1]-vn[0:-2,1:-1])-\
        dt/(2*rho*dy)*(p[2:,1:-1]-p[0:-2,1:-1])+\
        nu*(dt/dx**2*(vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])+\
        (dt/dy**2*(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])))

    ####Periodic BC u @ x = 2     
    u[1:-1,-1] = un[1:-1,-1]-\
        un[1:-1,-1]*dt/dx*(un[1:-1,-1]-un[1:-1,-2])-\
        vn[1:-1,-1]*dt/dy*(un[1:-1,-1]-un[0:-2,-1])-\
        dt/(2*rho*dx)*(p[1:-1,0]-p[1:-1,-2])+\
        nu*(dt/dx**2*(un[1:-1,0]-2*un[1:-1,-1]+un[1:-1,-2])+\
        dt/dy**2*(un[2:,-1]-2*un[1:-1,-1]+un[0:-2,-1]))+F*dt

    ####Periodic BC u @ x = 0
    u[1:-1,0] = un[1:-1,0]-\
        un[1:-1,0]*dt/dx*(un[1:-1,0]-un[1:-1,-1])-\
        vn[1:-1,0]*dt/dy*(un[1:-1,0]-un[0:-2,0])-\
        dt/(2*rho*dx)*(p[1:-1,1]-p[1:-1,-1])+\
        nu*(dt/dx**2*(un[1:-1,1]-2*un[1:-1,0]+un[1:-1,-1])+\
        dt/dy**2*(un[2:,0]-2*un[1:-1,0]+un[0:-2,0]))+F*dt

    ####Periodic BC v @ x = 2
    v[1:-1,-1] = vn[1:-1,-1]-\
    un[1:-1,-1]*dt/dx*(vn[1:-1,-1]-vn[1:-1,-2])-\
        vn[1:-1,-1]*dt/dy*(vn[1:-1,-1]-vn[0:-2,-1])-\
        dt/(2*rho*dy)*(p[2:,-1]-p[0:-2,-1])+\
        nu*(dt/dx**2*(vn[1:-1,0]-2*vn[1:-1,-1]+vn[1:-1,-2])+\
        (dt/dy**2*(vn[2:,-1]-2*vn[1:-1,-1]+vn[0:-2,-1])))

    ####Periodic BC v @ x = 0
    v[1:-1,0] = vn[1:-1,0]-\
        un[1:-1,0]*dt/dx*(vn[1:-1,0]-vn[1:-1,-1])-\
        vn[1:-1,0]*dt/dy*(vn[1:-1,0]-vn[0:-2,0])-\
        dt/(2*rho*dy)*(p[2:,0]-p[0:-2,0])+\
        nu*(dt/dx**2*(vn[1:-1,1]-2*vn[1:-1,0]+vn[1:-1,-1])+\
        (dt/dy**2*(vn[2:,0]-2*vn[1:-1,0]+vn[0:-2,0])))     


    ####Wall BC: u,v = 0 @ y = 0,2
    u[0,:] = 0
    u[-1,:] = 0
    v[0,:] = 0
    v[-1,:]=0
    
    udiff = (numpy.sum(u)-numpy.sum(un))/numpy.sum(u)
    stepcount += 1
    
    
    
fig = pyplot.figure(figsize = (6,4), dpi=100)
pyplot.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3]);

