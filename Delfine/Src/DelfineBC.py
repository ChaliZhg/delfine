#############################################################
# File: DelfineBC.py
# Function: Defines Boundary Conditions for the Elliptic Assemble
# Author: Bruno Luna
# Date: 29/05/11
# Modifications Date:
# 29/05/11 - Created this module to handle the boundary condition
#
#############################################################
from dolfin import *
from math import *

def u0_boundary(x, on_boundary):
    """Defines of a point x lies on the actual boundary of the domain"""  
    return on_boundary
            
class BCCond(Expression):
    """Defines values for the boundary conditions"""
    def eval (self, values, x):
        pass
    ############### Begin - Heterogeneous Anisotropic Ex.##############
#    def eval(self, values, x):
#        alpha = 1000
#        if (x[0] > 0):
#            values[0] = exp(x[0])*sin(x[1])
#        else:
#            values[0] = (2*sin(x[1]) + cos(x[1]))*alpha*x[0] + sin(x[1])
    ############### End - Heterogeneous Anisotropic Ex.###############
    ############### Begin - Homogeneous Anisotropic Ex.###############
#    def eval (self, values, x):
#        values[0] = exp(x[0]*x[1])
    ############### End - Homogeneous Anisotropic Ex.################
    ############### Begin - Homogeneous Isotropic Ex.#################
    # To use this example, cut the code snippet below and paste it in the AssembleElliptic.py
    # file in the define boundary conditions part (it is also necessary to comment the 2 lines
    # which deal with the general case for the BCs). You have to use the default eval function
    # for BC (the one above which just the 'pass' statement, to avoid erros)
#    # Define general boundary (x=0 or x=1 or y=0 or y=1)
#    def any_boundary(x):
#        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
#              x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
#    
#    # Define boundary (0<x<1 and y=0) - Homo Iso example
#    def dirichlet_boundary1(x):
#        return any_boundary(x) and abs(x[1]) < DOLFIN_EPS
#        
#    # Define boundary (0<x<1 and y=1) - Homo Iso example
#    def dirichlet_boundary2(x):
#        return any_boundary(x) and abs(x[1] - 1) < DOLFIN_EPS
#        
#    # Define boundary condition -  Homo Iso example
#    u0= Expression('cos(pi*x[0])*cos(pi*x[1])') # Analytical solution
#    u01= Expression('cos(pi*x[0])')
#    u02= Expression('- cos(pi*x[0])')
#    
#    bc1 = DirichletBC(V, u01, dirichlet_boundary1)
#    bc2 = DirichletBC(V, u02, dirichlet_boundary2)
#    bc = [bc1,  bc2]
    ############### End - Homogeneous Isotropic Ex.#################
