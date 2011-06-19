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
from DelfineData import ViscScalar, RelatPermScalar
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

        # Read data for fluid viscosity
        muW = ViscScalar( "water", parameter)
        muO = ViscScalar( "oil", parameter)

        # Read data for relative permeabilities
        Sw = 1.0 # FIXME: This comes from sat. solution
        krW = RelatPermScalar("water", parameter, Sw)
        krO = RelatPermScalar("oil", parameter, Sw)

        # Define partials and total mobilites
        mobW = krW/muW
        mobO = krO/muO
        mob = mobW + mobO

        # Create function space
        V = FunctionSpace(mesh, "CG", order)

        # Define Dirichlet boundary conditions
        u0= BCCond()
        bc = DirichletBC(V, Expression("0.0"), u0_boundary,  "pointwise") 
#        from DelfineBC import u1_boundary
#        u1 = BCCond()
#        bc = DirichletBC(V, Expression("1.0"), u1_boundary,  "pointwise") 

        # Define distributed source terms (not wells)
        f = Source()

        # Define Neumman boundary conditions
        g = Expression("0.0")

        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)

        a = inner(grad(v), K*mob*grad(u))*dx
        L = v*f*dx - g*v*ds

        # Assemble matrices and vectors
        A, rhs = assemble_system(a, L)
        bc.apply(A, rhs)

        # Apply well source/sink conditions
        # TEST ---------------------------------------------------
#        vertexs = vertices(mesh)
#        wells = parameter.geom.bc.wells
#        print "Get Wells:"
#        mf = mesh.data().mesh_function("well_indicators")
#        if (mf != None):
#            for vertex in vertexs:
#                for well in wells:
#                    if (mf[vertex.index()] == int(well.id)):
#                        print "Vertex:",  vertex.index()
#                        print "Vertex Coord:",  vertex.midpoint().x(), vertex.midpoint().y(), vertex.midpoint().z()
#                        print "Well ID:",  well.id
#                        print "Well Value:",  float(well.value)
#                        if (int(well.id) > 300):
#                            p = Point(vertex.midpoint().x(), vertex.midpoint().y())
#                            print vertex.midpoint().x(),  type(vertex.midpoint().y())
#                            q = PointSource(V, p, float(well.value)) 
#                            print float(well.value)
#                            q.apply(rhs)

        
        p = Point(0, 0)
        q = PointSource(V, p,1.0) 
        q.apply(rhs)
        
        # END TEST ------------------------------------------

        # Define return parameters
        delfineVar.A = A
        delfineVar.rhs = rhs
        delfineVar.V = V
        delfineVar.u0 = u0
############################################################



    # Define variational problem
    V = FunctionSpace(mesh, "Lagrange", order)
    v = TestFunction(V)
    p = TrialFunction(V)

    a = inner(grad(v), K*mob*grad(p))*dx
    L = v*f*dx - g*v*ds
    A, rhs = assemble_system(a, L) # A*p=rhs
    
    
    # Project velocity vector on vector function space
    Vg = VectorFunctionSpace(mesh, "Lagrange", order)
    vT = project(-K*mob*grad(p), Vg) - 
        project(-K*mobW*dPc*grad(Sw), Vg)


