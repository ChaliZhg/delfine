#############################################################
# File: DelfinePlot.py
# Function: Print graphics on screen and write files with the results
# Author: Bruno Luna
# Date: 03/02/11
# Modifications Date:
# 
#
#
#
#
#############################################################
from dolfin import *

class DelfinePlot:
    """Plots the results as contour plots and the residuals history as x-y plot"""
    
    def __init__(self):
        pass
    
    def plotResults(self, delfineVar ,  parameter):
        """Plots the results as contour plots and the residuals history as x-y plot
            Input:
            elliptic = AssembleElliptic() [AssembleElliptic.py]
            solver = SolveEqSystem() [SolveEqSystem.py]
        """
        
        # Getting data from elliptic eq. assembler
        A = delfineVar.A
        rhs = delfineVar.rhs
        V = delfineVar.V
        u0 = delfineVar.u0
    
        # Getting data from eq. system solver
        x = delfineVar.x
        residuals = delfineVar.residuals
        solverType = parameter.num.pressSolv.type
        preCondType = parameter.num.pressSolv.preConditioning.type
        
        # Function to convert NumPy array to Dolfin vector
        def as_vector(x):
            v = Vector(len(x))
            v.set_local(x)
            return v
        
        vx = as_vector(x)
        fx = Function(V, vx)
        
        # Save solution in VTK format
        fileVtk = File("Results/poisson.pvd")
        fileVtk << fx
        
        # Save solution in raw text format
        resList = list(residuals)

        # File name definition
        # PS: Name definitions just valid while
        if (preCondType == "none"):
            name = solverType.swapcase() 
        elif (solverType == "none"):
            name = preCondType.swapcase() 
        else:
            name = preCondType.swapcase() + '+' + solverType.swapcase() 
        
        fileTxt = open('Results/residuals_' + name + '.txt',"w")
        
        for res in resList:
            fileTxt.write(str(res) + "\n")
        
        fileTxt.close()
        
        # Calculate max error norm
        u_numerical = fx.vector().array()
        u0_interp = interpolate(u0, V)
        u_exact = u0_interp.vector().array()
        
        diff = abs(u_numerical - u_exact)
        print 'Max Error: %1.2E' % diff.max()
        
        # Plot solution
        plot(fx, interactive=True)
        
        import pylab
        pylab.figure(2)
        pylab.semilogy(residuals)
        pylab.show()
#############################################################
