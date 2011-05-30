#############################################################
# File: AssembleElliptic.py
# Function: Assemble variational problem with Dolfin,
# which is a module of Fenics (http://www.fenicsproject.org/)
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 04/04/11 - Mesh read function transfered to its own class
# 26/04/11 - Read subdomains data through mesh_function
# 28/05/11 - Permeability tensor classes moved to DelfineData module
#
#
#############################################################
from dolfin import *
from DelfineData import RockType,  PermeabilityTensor3D, PermeabilityTensor2D
from DelfineBC import BCCond, u0_boundary
from DelfineSource import Source

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
        
        # Mesh data (see DelfineMesh for details)
        mesh = meshData.allData
        dim = parameter.geom.mesh.dim
        order = parameter.geom.mesh.order
        
        # Read subdomains data for permeability
        if (dim == 2):
            K = PermeabilityTensor2D(mesh, parameter)
        elif (dim == 3):
            K = PermeabilityTensor3D(mesh, parameter)
        
        # Create function space
        V = FunctionSpace(mesh, "CG", order)
        
        # Define boundary conditions (except wells)
        u0= BCCond()
        bc = DirichletBC(V, u0, u0_boundary)
        
        # Define source terms (wells, etc.)
        f = Source()
        g = Expression("0.0")
        
        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)
        
        a = inner(grad(v), K*grad(u))*dx
        L = v*f*dx - g*v*ds
        
        # Assemble matrices and vectors
        A, rhs = assemble_system(a, L, bc)

        # Define return parameters
        delfineVar.A = A
        delfineVar.rhs = rhs
        delfineVar.V = V
        delfineVar.u0 = u0
############################################################
