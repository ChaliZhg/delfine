#############################################################
# File: DelfineSat.py
# Function: Calculates the saturation with the results for velocity
# Author: Bruno Luna
# Date: 26/10/11
# Modifications Date:
# 01/11/11 - Added saturation solver using SUPG
#############################################################
from dolfin import *

def boundary_value(n):
    if n < 10:
        return float(n)/10.0
    else:
        return 1.0

def u1_boundary(x, on_boundary):
    """Defines of a point x lies on the actual boundary of the domain"""  
    #return on_boundary
    return ((x[0] == 0.0) and (x[1] == 0.0)) #FIXME: Generalize boundaries
    
class DelfineSat:
    """Assembles the hyperbolic problem variational form with Dolfin."""   

    def __init__(self, parameter,  meshData):
        # Defining linear algebra package to be used (Petsc, uBlas, etc.)
        # Attention: this parameter(s) is different from the parameter
        # variable used by delfine to store the data file informations
        parameters.linear_algebra_backend = "uBLAS" #FIXME: Verificar se precisa de self

        # Mesh data (see DelfineMesh for details)
        self.mesh = meshData.allData
        self.dim = parameter.geom.mesh.dim
        self.order = parameter.geom.mesh.order 

    def SUPG(self, delfineVar, parameter):
        """Assembles the hyperbolic problem in a SUPG variational form with Dolfin."""  
        
        # Get data from base initialization
        mesh = self.mesh
        order = self.order
        K = delfineVar.K
        mob = delfineVar.mob
        
        # Get velocity from elliptic solver
        velElliptic = delfineVar.velocity
        
        # Get elements size
        h = CellSize(mesh)
    
        # Create FunctionSpaces
        Q = FunctionSpace(mesh, "CG", 1)
        V = VectorFunctionSpace(mesh, "CG", 2)
        
        # Create velocity Function from file
        velocity = Function(V);
        velocity = project(velElliptic, V)
        plot(velocity,  title='Test Velocity')
        #interactive()
        
        # Initialise source function and previous solution function
        f  = Constant(0.0)
        u0 = Function(Q)
        
        # Parameters
        T = 5.0
        dt = 0.01
        t = dt
        c = 0.00005
        
        # Test and trial functions
        u, v = TrialFunction(Q), TestFunction(Q)
        
        # Mid-point solution
        u_mid = 0.5*(u0 + u)
        
        # Residual
        r = u-u0 + dt*(dot(velocity, grad(u_mid)) - c*div(grad(u_mid)) - f)
        
        # Galerkin variational problem
        F = v*(u-u0)*dx + dt*(v*dot(velocity, grad(u_mid))*dx + c*dot(grad(v), grad(u_mid))*dx)
        
        # Add SUPG stabilisation terms
        vnorm = sqrt(dot(velocity, velocity))
        F += (h/2.0*vnorm)*dot(velocity, grad(v))*r*dx
        
        # Create bilinear and linear forms
        a = lhs(F)
        L = rhs(F)
        
        # Set up boundary condition
        g = Constant(boundary_value(0) )
        bc = DirichletBC(Q, g, u1_boundary,  "pointwise")
        
        # Assemble matrix
        A = assemble(a)
        bc.apply(A)
        
        # Create linear solver and factorize matrix
        solver = LUSolver(A)
        solver.parameters["reuse_factorization"] = True
        
        # Output file
        out_file = File("Results/supg_saturation.pvd")
        
        # Set intial condition
        u = u0
        
        # Time-stepping
        while t < T:
        
            # Assemble vector and apply boundary conditions
            b = assemble(L)
            bc.apply(b)
        
            # Solve the linear system (re-use the already factorized matrix A)
            solver.solve(u.vector(), b)
        
            # Copy solution from previous interval
            u0 = u
        
            # Plot solution
            plot(u)
        
            # Save the solution to file
            out_file << (u, t)
        
            # Move to next interval and adjust boundary condition
            t += dt
            g.assign( boundary_value( int(t/dt) ) )
