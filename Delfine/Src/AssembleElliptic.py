#############################################################
# File: AssembleElliptic.py
# Function: Assemble variational problem with Dolfin,
# which is a module of Fenics (http://www.fenicsproject.org/)
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 04/04/11 - Mesh read function transfered to its own class
#
#
#
#
#############################################################
from dolfin import *

class AssembleElliptic:
    """Assembles the elliptic problem variational form with Dolfin."""  
    
    def __init__(self):
        pass
     
    def assemble_withDolfin(self, delfineVar, parameter,  meshData):
        """Assembles the elliptic problem variational form with Dolfin."""  
        
        # Defining linear algebra package to be used (Petsc, uBlas, etc.)
        # Attention: this parameter(s) is different from the parameter
        # variable used by delfine to store the data file informations
        parameters.linear_algebra_backend = "uBLAS"
        
        # Create mesh with dolfin functions or read mesh from file
        mesh = meshData.allData
        
        # Create function space
        V = FunctionSpace(mesh, "CG", 1)
        
        # Define general boundary (x=0 or x=1 or y=0 or y=1)
        def any_boundary(x):
            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
                  x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
    
        # Define boundary (0<x<1 and y=0)
        def dirichlet_boundary1(x):
            return any_boundary(x) and abs(x[1]) < DOLFIN_EPS
            
        # Define boundary (0<x<1 and y=1)
        def dirichlet_boundary2(x):
            return any_boundary(x) and abs(x[1] - 1) < DOLFIN_EPS
            
        # Define boundary condition
        u0= Expression('cos(pi*x[0])*cos(pi*x[1])') # Analytical solution
        u01= Expression('cos(pi*x[0])')
        u02= Expression('- cos(pi*x[0])')
        
        bc1 = DirichletBC(V, u01, dirichlet_boundary1)
        bc2 = DirichletBC(V, u02, dirichlet_boundary2)
        bc = [bc1,  bc2]
        
        # Define permeability matrix
        C = as_matrix(((1.0, 0.0), (0.0, 1.0)))
        
        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)
        f = Expression("2*pow(pi,2)*cos(pi*x[0])*cos(pi*x[1])")
        g = Expression("0.0")
        a = inner(grad(v), C*grad(u))*dx
        L = v*f*dx - g*v*ds
        
        # Assemble matrices and vectors
        A, rhs = assemble_system(a, L, bc)
        
        # Define return parameters
        delfineVar.A = A
        delfineVar.rhs = rhs
        delfineVar.V = V
        delfineVar.u0 = u0
############################################################
