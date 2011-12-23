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

# Create mesh and function space
mesh = UnitInterval(50)
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
vel = 1.0
T = 1.0
dt = 0.001
t = dt
h = CellSize(mesh)
step = 0

u_mida = variable(u_mid)
fwsb = 2.0*u_old/(u_old**2.0 + (1.0-u_old)**2.0) - (u_old**2.0)*(2.0*u_old - 2.0*(1.0-u_old))/((u_old**2.0 + (1.0-u_old)**2.0)**2.0)

r = (u-u_old) + dt*(fwsb*vel*grad(u_mid)) 
F = (u-u_old)*v*dx + dt*(fwsb*vel*grad(u_mid)*v*dx)
# Add SUPG stabilisation terms
vnorm = sqrt(dot(vel, vel))
F += (h/(2.0*vnorm))*dot(vel, grad(v))*r*dx

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
    
    # Plot solution
    #plot(u, title=t)
    uplot = u.vector().array()
    step += 1
    if ((step%100) == 0):
        plt.plot(mesh.coordinates(), uplot, 'bo-')
        plt.xlabel("x")
        plt.ylabel("C(x)")
        plt.title("Linear Advection")
        plt.axis([0, 1, -0.2, 1.2])
        plt.grid(True)

    # Save the solution to file
    out_file << (u, t)

    # Move to next interval and adjust boundary condition
    print "Time: ", t
    t += dt
    print "Front Location", vel*t
    
plt.show()
