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
import numpy
import pylab

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
        
        # Plot numerical solution
        plot(fx, interactive=True,  title="Numerical Solution")
        
        # Create array with residuals of the solution
        resList = list(residuals)

        # Name definition for the residuals file
        if (preCondType == "none"):
            name = solverType.swapcase() 
        elif (solverType == "none"):
            name = preCondType.swapcase() 
        else:
            name = preCondType.swapcase() + '+' + solverType.swapcase() 
        
        fileTxt = open('Results/residuals_' + name + '.txt',"w")
        
        # Save residuals of solution in raw text format
        for res in resList:
            fileTxt.write(str(res) + "\n")
        
        fileTxt.close()
        
        # Calculate error norms
        u_numerical = fx.vector().array()
        u0_interp = interpolate(u0, V)
        u_exact = u0_interp.vector().array()
        
        diff = u_numerical - u_exact
        diffMax = abs(diff)
        print 'Error for Max Norm: %1.2E' % diffMax.max()
        E = errornorm(u0, fx, norm_type='l2', degree=3)
        print 'Error norm L2: %1.2E' % E
        
        # Plot analytic solution
        #plot(u0_interp, interactive=True, title="Analytical Solution")
        fileVtk = File("Results/poissonanalytic.pvd")
        fileVtk << u0_interp
        
        # Plot difference between numerical and analytical solution
        e = as_vector(diff)
        fe = Function(V, e)
        #plot(fe, interactive=True,  title="Error")
        fileVtk = File("Results/poissonerror.pvd")
        fileVtk << fe
        
        pylab.figure(2)
        pylab.semilogy(residuals)
        pylab.show()
#############################################################
