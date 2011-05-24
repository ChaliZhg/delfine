#############################################################
# File: AssembleElliptic.py
# Function: Assemble variational problem with Dolfin,
# which is a module of Fenics (http://www.fenicsproject.org/)
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 04/04/11 - Mesh read function transfered to its own class
# 26/04/11 - Read subdomains data through mesh_function
#
#
#
#############################################################
from dolfin import *
from DelfineData import RockType

from math import *
# FIXME: Merge both PermeabilityTensor classes in a more
# "elegant" one(deal with value_shape), or take them to their own file
class PermeabilityTensor3D(Expression):
    """Defines the 3D permeability tensor for each cell"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K =[]
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            self.DomainID.append(int(self.rocks[i].id))
            self.K.append(self.rocks[i].permeability.K)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        values[0] = self.K[j][0] # Kxx
                        values[1] = self.K[j][1] # Kxy
                        values[2] = self.K[j][2] # Kxz            | Kxx  Kxy  Kxz |       | values[0] values[1] values[2] |
                        values[3] = values[1]   # Kyx     K = | Kyx  Kyy  Kyz  | => | values[3] values[4] values[5] |
                        values[4] = self.K[j][3] # Kyy            | Kzx  Kzy  Kzz  |       | values[6] values[7] values[8] |
                        values[5] = self.K[j][4] # Kyz
                        values[6] = values[2]   # Kzx
                        values[7] = values[5]   # Kzy
                        values[8] = self.K[j][5] # Kzz
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                values[0] = self.K[0][0] # Kxx
                values[1] = self.K[0][1] # Kxy
                values[2] = self.K[0][2] # Kxz            | Kxx  Kxy  Kxz |       | values[0] values[1] values[2] |
                values[3] = values[1]    # Kyx     K = | Kyx  Kyy  Kyz  | => | values[3] values[4] values[5] |
                values[4] = self.K[0][3] # Kyy            | Kzx  Kzy  Kzz  |       | values[6] values[7] values[8] |
                values[5] = self.K[0][4] # Kyz
                values[6] = values[2]    # Kzx
                values[7] = values[5]    # Kzy
                values[8] = self.K[0][5] # Kzz
        else:
            pass
    def value_shape(self):
        return (3, 3)
        
# FIXME: Merge both PermeabilityTensor classes in a more
# "elegant" one(deal with value_shape), or take them to their own file
class PermeabilityTensor2D(Expression):
    """Defines the 2D permeability tensor for each cell"""
    def __init__(self, mesh, parameter):
        self.mesh = mesh
        self.rocks = parameter.phys.rock.rocks
        self.DomainID = []
        self.K =[]
        # Getting rocks IDs and respective K from parameter structure
        for i in range(len(self.rocks)):
            self.DomainID.append(int(self.rocks[i].id))
            self.K.append(self.rocks[i].permeability.K)
    def eval_cell(self, values, x, ufc_cell):
        # Get material indicator(which corresponds to rock ID) from mesh
        mf = self.mesh.data().mesh_function("material_indicators")
        i = ufc_cell.index
        j = 0
        if (type(i) is int):
            if (mf != None):
                for id in self.DomainID:
                    if (mf[i] == id):
                        values[0] = self.K[j][0] # Kxx
                        values[1] = self.K[j][1] # Kxy     K = | Kxx  Kxy | => | values[0] values[1] |
                        values[2] = values[1]   # Kyx            | Kyx  Kyy  |       | values[2] values[3] |      
                        values[3] = self.K[j][2] # Kyy
                    j += 1
            else: 
                # For dolfin-generated meshes or meshes without "material_indicator"
                # This option consider just homogeneous cases
                values[0] = self.K[0][0] # Kxx
                values[1] = self.K[0][1] # Kxy
                values[2] = values[1]    # Kyx
                values[3] = self.K[0][2] # Kyy
        else:
            pass
    def value_shape(self):
        return (2, 2)

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

        # Read subdomains data for permeability with meshfunction
        dim = parameter.geom.mesh.dim
        if (dim == 2):
            K = PermeabilityTensor2D(mesh, parameter)
        elif (dim == 3):
            K = PermeabilityTensor3D(mesh, parameter)
        
        # Create function space
        # FIXME: Define element order as input variable
        V = FunctionSpace(mesh, "CG", 1)

        ############### Begin - Homogeneous Isotropic Ex.###############
#         Define general boundary (x=0 or x=1 or y=0 or y=1)
#        def any_boundary(x):
#            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
#                  x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
#        
#        # Define boundary (0<x<1 and y=0) - Homo Iso example
#        def dirichlet_boundary1(x):
#            return any_boundary(x) and abs(x[1]) < DOLFIN_EPS
#            
#        # Define boundary (0<x<1 and y=1) - Homo Iso example
#        def dirichlet_boundary2(x):
#            return any_boundary(x) and abs(x[1] - 1) < DOLFIN_EPS
#            
#        # Define boundary condition -  Homo Iso example
#        u0= Expression('cos(pi*x[0])*cos(pi*x[1])') # Analytical solution
#        u01= Expression('cos(pi*x[0])')
#        u02= Expression('- cos(pi*x[0])')
#        
#        bc1 = DirichletBC(V, u01, dirichlet_boundary1)
#        bc2 = DirichletBC(V, u02, dirichlet_boundary2)
#        bc = [bc1,  bc2]
#
#        f = Expression("2*pow(pi,2)*cos(pi*x[0])*cos(pi*x[1])") # Homo Iso example
        ############### End - Homogeneous Isotropic Ex.#################
        ############### Begin - Homogeneous Anisotropic Ex.###############
#        #Define general boundary (x=0 or x=1 or y=0 or y=1)
#        def any_boundary(x):
#            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or \
#                x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
#
#        #Define boundary condition -  Homo Aniso example
#        u0= Expression('exp(x[0]*x[1])') # Analytical solution and bc
#        bc = DirichletBC(V, u0, any_boundary)      
#       
#        f = Expression("-2*(1 + pow(x[0],2) + x[0]*x[1] + pow(x[1],2))*exp(x[0]*x[1])") # -  Homo Aniso example
        ############### End - Homogeneous Anisotropic Ex.#################
        ############### Begin - Heterogeneous Anisotropic Ex.#################
        # Define general boundary (x=-1 or x=1 or y=-1 or y=1)

        
        def any_boundary(x):
            return ((x[0] + 1) < DOLFIN_EPS) or (x[0] > 1.0 - DOLFIN_EPS) or \
                ((x[1] + 1) < DOLFIN_EPS) or (x[1] > 1.0 - DOLFIN_EPS)
        
        class BCCond(Expression):
            def eval(self, values, x):
                alpha = 1000
                if (x[0] > 0):
                    values[0] = exp(x[0])*sin(x[1])
                else:
                    values[0] = (2*sin(x[1]) + cos(x[1]))*alpha*x[0] + sin(x[1])
                    
        
        u0= BCCond()
        
        bc = DirichletBC(V, u0, any_boundary)
        
        class Source(Expression):
            def eval(self, values, x):
                alpha = 1000
                if (x[0] > 0):
                    values[0] = -2*alpha*exp(x[0])*cos(x[1])
                else:
                    values[0] = (2*sin(x[1]) + cos(x[1]))*alpha*x[0] + sin(x[1])

        f = Source()
        ############### End - Heterogeneous Anisotropic Ex.#################
        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)
        # f = Expression("...")
        g = Expression("0.0")
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
