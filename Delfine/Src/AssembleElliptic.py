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
# 19/06/11 - Added the possibility to use mixed element method
#
#############################################################
from dolfin import *
from DelfineData import RockType,  PermeabilityTensor3D, PermeabilityTensor2D
from DelfineData import InversePermeabilityTensor2D, InversePermeabilityTensor3D
from DelfineData import ViscScalar, RelatPermScalar
from DelfineBC import BCCond, u0_boundary
from DelfineSource import Source

class AssembleElliptic:
    """Assembles the elliptic problem variational form with Dolfin."""  

    def __init__(self, parameter,  meshData):
        
        # Defining linear algebra package to be used (Petsc, uBlas, etc.)
        # Attention: this parameter(s) is different from the parameter
        # variable used by delfine to store the data file informations
        parameters.linear_algebra_backend = "uBLAS" #FIXME: Verificar se precisa de self
        
        # Mesh data (see DelfineMesh for details)
        self.mesh = meshData.allData
        self.dim = parameter.geom.mesh.dim
        self.order = parameter.geom.mesh.order

        # Read subdomains data for permeability
        if (self.dim == 2):
            self.K = PermeabilityTensor2D(self.mesh, parameter)
            self.invK = InversePermeabilityTensor2D(self.mesh, parameter)
        elif (self.dim == 3):
            self.K = PermeabilityTensor3D(self.mesh, parameter)
            self.invK = InversePermeabilityTensor3D(self.mesh, parameter)

        # Read data for fluid viscosity
        #muW = ViscScalar( "water", parameter)
        #muO = ViscScalar( "oil", parameter)

        # Read data for relative permeabilities
        #Sw = 1.0 # FIXME: This comes from sat. solution
        #Sw = Expression("1.0")
        #krW = RelatPermScalar("water", parameter, Sw)
        #krO = RelatPermScalar("oil", parameter, Sw)

        # Define partials and total mobilites
        #mobW = krW/muW
        #mobO = krO/muO
        #self.mob = mobW + mobO
############################################################
    def Galerkin(self, delfineVar, parameter):
        """Assembles the elliptic problem in a Galerkin variational form with Dolfin.""" 

        # Get data from base initialization
        mesh = self.mesh
        order = self.order
        K = self.K
        mob = self.mob
        
        # Create function space
        V = FunctionSpace(mesh, "CG", order)
        
        # Define Dirichlet boundary conditions
#        u0= BCCond()
#        bc = DirichletBC(V, Expression("0.0"), u0_boundary,  "pointwise") 
#        from DelfineBC import u1_boundary
#        u1 = BCCond()
#        bc = DirichletBC(V, Expression("1.0"), u1_boundary,  "pointwise") 
        u0 = Constant("0.0") # FIXME: Define BC correctly
        
        # Define distributed source terms (not wells)
        f = Source()

        # Define Neumman boundary conditions (vel.n=0 on borders)
        g = Expression("0.0")

        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)

        a = inner(grad(v), K*mob*grad(u))*dx
        L = v*f*dx - g*v*ds

        # Assemble matrices and vectors
        A, rhs = assemble_system(a, L)
        #bc.apply(A, rhs)

        # Apply well source/sink conditions # FIXME: Pressure/Saturation Well still not implemented
        vertexs = vertices(mesh)
        wells = parameter.geom.bc.wells
        mf = mesh.data().mesh_function("well_indicators")
        if (mf != None):
            for vertex in vertexs:
                for well in wells:
                    if (mf[vertex.index()] == int(well.id)):
                        if (int(well.id) > 300): # Source/Sink Well
                            p = Point(vertex.midpoint().x(), vertex.midpoint().y())
                            q = PointSource(V, p, float(well.value)) 
                            q.apply(rhs)
        

        # Define return parameters
        delfineVar.A = A
        delfineVar.rhs = rhs
        delfineVar.V = V
        delfineVar.u0 = u0
        delfineVar.K = K
        delfineVar.mob = mob
############################################################
    def MixedFEM(self, delfineVar, parameter):
        """Assembles the elliptic problem in a Mixed FEM variational form with Dolfin.""" 
        
        # Get data from base initialization
        mesh = self.mesh
        invK = self.invK
        #mob = self.mob
        
        # FIXME: Decide if these parameters are initialized here
        Sw = delfineVar.satW
        muW = ViscScalar( "water", parameter) 
        muO = ViscScalar( "oil", parameter)
        krW = delfineVar.krW
        krO = delfineVar.krO
        krW.Sw = Sw
        krO.Sw = Sw
        
        # Define partials and total mobilites
        mobW = krW/muW
        mobO = krO/muO
        mob = mobW + mobO
        
        # Define function spaces and mixed (product) space
        BDM = FunctionSpace(mesh, "BDM", 1)
        DG = FunctionSpace(mesh, "DG", 0)
        W = BDM * DG
        
        # Define trial and test functions
        (sigma, u) = TrialFunctions(W)
        (tau, v) = TestFunctions(W)
        
        # Define source function
        f = Expression("0.0")
        
        # Define variational form
        a = (dot((invK/mob)*sigma, tau) - div(tau)*u - div(sigma)*v)*dx
        L = - f*v*dx
        
        # Define function G such that vel.n=0 on borders
        class BoundarySource(Expression):
            def __init__(self, mesh):
                self.mesh = mesh
            def eval_cell(self, values, x, ufc_cell):
                cell = Cell(self.mesh, ufc_cell.index)
                n = cell.normal(ufc_cell.local_facet)
                g = 0.0
                values[0] = g*n[0]
                values[1] = g*n[1]
            def value_shape(self):
                return (2,)
        
        G = BoundarySource(mesh)
        
        # Define essential boundary
        def boundary(x, on_boundary):
            return on_boundary
        
        bc = DirichletBC(W.sub(0), G, boundary)
        
        # Assemble system
        A, rhs = assemble_system(a, L)
        bc.apply(A, rhs)
        
        # Apply well source/sink conditions # FIXME: Pressure/Saturation Well still not implemented
        vertexs = vertices(mesh)
        wells = parameter.geom.bc.wells
        mf = mesh.data().mesh_function("well_indicators")
        if (mf != None):
            for vertex in vertexs:
                for well in wells:
                    if (mf[vertex.index()] == int(well.id)):
                        if (int(well.id) > 300): # Source/Sink Well
                            p = Point(vertex.midpoint().x(), vertex.midpoint().y())
                            q = PointSource(W.sub(1), p, -float(well.value)) # Verify minus signal
                            q.apply(rhs)
        
        
        # Compute solution # FIXME: Rename the next 'u' and 'v' variable to avoid confusion with
        # trial and test functions above
        u = Function(W)
        solve(A, u.vector(), rhs)
        
        # Velocity
        v = u.split()[0]
        
        # Force condition \ int p dx = 0
        p = u.split(deepcopy=True)[1]
        int_p = p*dx
        average_p = assemble (int_p, mesh=mesh)
        p_array = p.vector().array() - average_p
        p.vector()[:] = p_array
        
        # Define return parameters
        delfineVar.pressTot = p
        delfineVar.velocity = v
############################################################
