#############################################################
# File: Delfine.py
# Function: Class containing the main methods and attributes of the program 
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 04/04/11 - Added function to read mesh data outside the main loop
#
#
#
#
############################################################
import sys
from dolfin import *
from AssembleElliptic import AssembleElliptic
from DelfineData import DelfineData
from DelfineInput import DelfineInput
from DelfineMesh import DelfineMesh
from DelfinePlot import DelfinePlot
from SolveEqSystem import SolveEqSystem

class Delfine:
    """Delfine class"""
    
    def __init__(self, argv=None):
        if (len(argv) == 2):
            self.argv = argv
        else:
            print "Error(1):"
            print "Correct program usage: ./main.py 'caseName' "
            print " "
            sys.exit(1)

    def initialize(self):
        """Initialize the main objects"""
        
        # Object used to store input data
        self.parameter = DelfineData("input")
        
        # Read input data
        self.input = DelfineInput()
        self.input.readData(self.parameter,  self.argv)
        
        # Object used for data transfer among objects
        self.transferData = DelfineData("transfer")
        
        # Defines mesh for dolfin use (dolfin-generated or .xml)
        self.mesh = DelfineMesh(self.parameter)
        
        # Elliptic equations assembler
        self.elliptic = AssembleElliptic()
        
        # Equation system solver
        self.solver = SolveEqSystem()
        
        # Results and residual plotter
        self.outPlot = DelfinePlot()
        
    def drive(self):
        """Manages solution procedure"""
        
        # Assemble elliptic (pressure) equation
        self.elliptic.assemble_withDolfin(self.transferData, self.parameter,  self.mesh)
        
        # Solve equations system provenient from the variational form of the problem
        self.solver.solve_withPyAMG(self.transferData, self.parameter)
        
        # Plot solution, residual and compare numerical to analytical solution
        self.outPlot.plotResults(self.transferData,  self.parameter)
        

#############################################################
