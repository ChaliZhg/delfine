# -*- coding: utf-8 -*-
"""This demo program solves the linear advection equation

    dS(x,t)/dt + v(x,t)dS(x,t)/dx = 0

on the unit interval and boundary conditions given by

    u(x,0) = 0.0
    u(0,t) = 1.0 for  t_ini<t<t_end
"""

# Copyright (C) 2011 Bruno Luna
#
# First added:  2011-12-10

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

# Create mesh and function space
mesh = UnitInterval(500)
V = FunctionSpace(mesh, "CG", 1)

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < DOLFIN_EPS

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("0.0")
u_old = Function(V)
u_mid = 0.5*(u_old + u)
vel = 0.7
T = 1.0
dt = 0.0005
t = dt
h = CellSize(mesh)
step = 0

u_mida = variable(u_mid)
#fwsb = 2.0*u_old/(u_old**2.0 + (1.0-u_old)**2.0) - (u_old**2.0)*(2.0*u_old - 2.0*(1.0-u_old))/((u_old**2.0 + (1.0-u_old)**2.0)**2.0)
#fwsb =1.0
fwsb = 2.0*u_old*(1 - u_old)/((u_old**2.0 + (1.0 - u_old)**2.0)**2.0) 

r = (u-u_old) + dt*(fwsb*vel*grad(u_mid)) 
F = (u-u_old)*v*dx + dt*(fwsb*vel*grad(u_mid)*v*dx)
# Add SUPG stabilisation terms
vnorm = sqrt(dot(vel, vel))
F += (h/(2.0*vnorm))*dot(vel, grad(v))*r*dx

# Add shock capturing term # FIXME: Check what is wrong! Look Codina's Paper
beta = 6E-7
snorm = sqrt(dot(grad(u_old), grad(u_old)))
snorm = abs(snorm)
c = beta*h*snorm
F+=c*dot(grad(v), grad(u_mid))*dx

#beta = 2.0E-5
#snorm = sqrt(dot(grad(u_old), grad(u_old)))
#tol = 1E-15
#r_old = (u_old-u_old2) + dt*(fwsb*vel*grad(u_old)) 
#if (abs(snorm) > tol):
#    vshock = ((beta*h)*abs(r_old))/(2*snorm)
#else:
#    vshock = 0.0
#F+= vshock*dot(grad(v), grad(u_mid))*dx


# Create bilinear and linear forms
a = lhs(F)
L = rhs(F)

# Define boundary condition
u0 = Constant(1.0)
bc = DirichletBC(V, u0, DirichletBoundary())

#Assemble system
A = assemble(a)
bc.apply(A)
solver = LUSolver(A)

# Set initial condition
u = u_old

# Output file
out_file = File("LinearAdvSUPG.pvd")

while t < T:
    # Assemble vector and apply boundary conditions
    b = assemble(L)
    bc.apply(b)

    # Solve the linear system (re-use the already factorized matrix A)
    solver.solve(u.vector(), b)

    # Copy solution from previous interval
    u0 = u
    
    # Plot numerical solution
    #plot(u, title=t)
    uplot = u.vector().array()
    step += 1
    TPrint = 0.9
    # Calculate analytical solution
    x = np.linspace(0,1,110)
    uanalytic = []
    for xi in x:
        if (xi < vel*TPrint):
            uanalytic.append(1.0)
        else:
            uanalytic.append(0.0)
    #if ((step%100) == 0):
    #if ((step) == (TPrint/(vel*dt))):
    if (t > TPrint*0.9999) and ((t < TPrint*1.0001)):
        for i in range(len(uplot)):
            if (uplot[i] < 0.1):
                uplot[i] = uplot[i]/3.0
        plt.plot(mesh.coordinates(), uplot, 'bo-', label=u"Sol. Numérica", lw=2)
        #plt.plot(x, uanalytic, 'k-', label=u"Sol. Analítica",  lw=2)
        plt.xlabel("x")
        plt.ylabel("S(x)")
        plt.title(u"Buckley-Leverett")
        plt.axis([0, 1.0, -0.2, 1.2])
        plt.grid(True)
        plt.legend()
        
    # Save the solution to file
    out_file << (u, t)

    # Move to next interval and adjust boundary condition
    print "Time: ", t
    t += dt
    print "Front Location", vel*t

#plot(project(snorm, V))  
fwsb2 = 2.0*.69*(1 - .69)/((.69**2.0 + (1.0 - .69)**2.0)**2.0)
print 1.0*0.9*fwsb2
plt.show()

    
