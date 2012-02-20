#############################################################
# File: DelfineSat.py
# Function: Calculates the saturation with the results for velocity
# Author: Bruno Luna
# Date: 26/10/11
# Modifications Date:
# 01/11/11 - Added saturation solver using SUPG
#############################################################
from dolfin import *
from DelfineData import WaterFracFlowDiff, RelatPermScalar

def boundary_value(n):
    if n < 5:
        return float(n)/5.0
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

    def SUPG(self, delfineVar, parameter,  elliptic):
        """Assembles the hyperbolic problem in a SUPG variational form with Dolfin."""  
        
        # Initial value for saturation FIXME: This is a parameter or inital condition
        Sw = Expression("0.0")
        delfineVar.satW = Sw
        fws = WaterFracFlowDiff(parameter, Sw)
        
        krW = RelatPermScalar("water", parameter, Sw)
        krO = RelatPermScalar("oil", parameter,  Sw)
        delfineVar.krW = krW
        delfineVar.krO = krO
        elliptic.MixedFEM(delfineVar, parameter)
        
        # Get data from base initialization
        mesh = self.mesh
        order = self.order
        K = delfineVar.K
        phi = 1.0 # FIXME: Define porosity from parameters
        
        # Get elements size
        h = CellSize(mesh)
    
        # Create FunctionSpaces
        Q = FunctionSpace(mesh, "CG", 1)
        V = VectorFunctionSpace(mesh, "CG", 2)
        
        # Get velocity from elliptic solver
        velElliptic = delfineVar.velocity
        velocity = Function(V);
        velocity = project(velElliptic, V)
        plot(velocity,  title='Initial Velocity')
    
        # Initialise source function and previous solution function
        f  = Constant(0.0)
        u0 = Function(Q)
        
        # Parameters #FIXME: Should be in input file
        T = 5.0
        dt = 0.005
        t = dt
        step = 0
        c = 0.0
        
        # Test and trial functions
        u, v = TrialFunction(Q), TestFunction(Q)
        
        # Mid-point solution
        u_mid = 0.5*(u0 + u)
        
        u_mida = variable(u_mid)
        #fwsb = 2.0*u0/(u0**2.0 + (1.0-u0)**2.0) - (u0**2.0)*(2.0*u0 - 2.0*(1.0-u0))/((u0**2.0 + (1.0-u0)**2.0)**2.0)
        #fwsb =1.0
        fwsb = 2.0*u0*(1 - u0)/((u0**2.0 + (1.0 - u0)**2.0)**2.0) 
        
        # Residual
        r = phi*(u-u0) + dt*((fwsb*dot(velocity, grad(u_mid))) - c*div(grad(u_mid)) - f)
        
        # Galerkin variational problem
        F = v*phi*(u-u0)*dx + dt*(v*fwsb*dot(velocity, grad(u_mid))*dx + c*dot(grad(v), grad(u_mid))*dx)
        
        # Add SUPG stabilisation terms
        vnorm = sqrt(dot(velocity, velocity))
        F += (h/(2.0*vnorm))*dot(velocity, grad(v))*r*dx
        
        # Add shock capturing term # FIXME: Check what is wrong! Look Codina's Paper
#        beta = 2.0
#        snorm = sqrt(dot(grad(u0), grad(u0)))
#        tol = 1E-15
#        if (abs(snorm) > tol):
#            print "Here!!!!!!!!!!!!!!!!!!!"
#            vshock = (beta*h)/(2.0*snorm)#*abs(r))#/(2*snorm)
#            plot(project(vshock, Q),  title='Shock Capturing Term')
#        else:
#            vshock = 0.0
#        F+= vshock*dot(grad(v), grad(u_mid))*dx
        beta = 9E-6
        snorm = sqrt(dot(grad(u0), grad(u0)))
        snorm = abs(snorm)
        c = beta*h*snorm
        F+=c*dot(grad(v), grad(u_mid))*dx
        
        # Create bilinear and linear forms
        a = lhs(F)
        L = rhs(F)
        
        # Set up boundary condition
        g = Constant(boundary_value(0) )
        bc = DirichletBC(Q, g, u1_boundary,  "pointwise")
        
        # Assemble matrix
        A = assemble(a)
        bc.apply(A)
        
        # Create linear solver
        solver = LUSolver(A)
        
        # Output files
        out_file_sat = File("Results/supg_saturation.pvd")
        out_file_vel = File ( "Results/mfem_velocity.pvd" )
        out_file_pressTot = File ( "Results/mfem_pressure.pvd" )
        
        # Set intial condition
        u = u0
        
        # Time-stepping
        while t < T:
            
            print "Time: ",  t,  " of ",  T
            # Assemble vector and apply boundary conditions
            A = assemble(a)
            bc.apply(A)
            solver = LUSolver(A)
            b = assemble(L)
            bc.apply(b)
        
            # Solve the linear system
            solver.solve(u.vector(), b)
            
            #bacalho (Sw<=1.0) and (Sw>0.0)
            u_array = u.vector().array()
            for i in range(len(u_array)):
                if(u_array[i] > 1.0):
                    u_array[i] = 1.0
                elif(u_array[i] < 0.0):
                    u_array[i] = 0.0
            u.vector()[:] = u_array
            
            # Copy solution from previous interval
            u0 = u
        
            # Plot solution
            plot(u,  title='Water Saturation')
            
            # Save the solution to file
            out_file_sat << (u, t)
            delfineVar.satW = u
            
            # Move to next interval and adjust boundary condition
            t += dt
            g.assign( boundary_value( int(t/dt) ) )
            
            # FIXME: Bacalho de Sequential Implicit
            elliptic.MixedFEM(delfineVar, parameter)
            fws.Sw = u
            
            velElliptic = delfineVar.velocity
            velocity = project(velElliptic, V)
            
            out_file_pressTot << delfineVar.pressTot
            out_file_vel << velocity
        
            step += 1
            #if ((step%10) == 0):
                #plot(velocity,  title='Updated Velocity')
