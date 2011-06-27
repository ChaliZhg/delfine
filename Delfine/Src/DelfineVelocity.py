#############################################################
# File: DelfineVelocity.py
# Function: Calculates the velocity with the result for the pressure
# Author: Bruno Luna
# Date: 16/06/11
# Modifications Date:
#
#############################################################
from dolfin import *

class DelfineVelocity:
    """Calculates the velocity from the pressure results"""  

    def __init__(self):
        pass 

    def velCalc(self, delfineVar, parameter,  meshData):
        """Calculates the velocity from the pressure results"""  
    
        mesh = meshData.allData
        dim = parameter.geom.mesh.dim
        order = parameter.geom.mesh.order
        
        # Compute solution with Dolfin internal solver
        # FIXME: Use solution from PyAMG
        V = delfineVar.V
        A = delfineVar.A
        rhs = delfineVar.rhs
        K = delfineVar.K
        mob = delfineVar.mob
        
        w = Function(V)
        solve(A, w.vector(), rhs, 'gmres')
        
        int_p = w*dx
        average_p = assemble (int_p, mesh=mesh)
        p_array = w.vector().array() - average_p
        w.vector()[:] = p_array
        
        plot(w, title='Pressure - Dolfin')
        
        mesh = meshData.allData
        dim = parameter.geom.mesh.dim
        order = parameter.geom.mesh.order
        
        # Compute and plot flux function
        Vg = VectorFunctionSpace(mesh, "CG", order)
        velocity = project(-K*mob*grad(w), Vg)
        plot(velocity,  title='Velocity')
    
        flux_x, flux_y = velocity.split(deepcopy=True) # extract components
        plot(flux_x, title='x-component of flux (-p*grad(u))')
        plot(flux_y, title='y-component of flux (-p*grad(u))')
