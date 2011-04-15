#############################################################
# File: SolveEqSystem.py
# Function: Solves equation system of the type Ax=b using mainly multigrid methods provenient
# from the pyAMG package (http://code.google.com/p/pyamg/)
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 
#
#
#
#
#############################################################
from scipy.sparse import csr_matrix
from pyamg import smoothed_aggregation_solver
from pyamg.krylov._cg import cg
from pyamg.krylov._gmres import gmres
import sys

class SolveEqSystem:
    """Solves the equation system with pyAMG.""" 
    
    def __init__(self):
        pass
    
    def solve_withPyAMG(self, delfineVar, parameter):
        """Solves the equation system with pyAMG.
            Input:
            elliptic = AssembleElliptic() [AssembleElliptic.py]
        """
        
        # Reads  solver data
        # Alternatives: CG (with or without PreCond) | GMRES (w/ or w/o PC)
        #                     | LU (still pending))
        solverType = parameter.num.pressSolv.type
        tolerance = parameter.num.pressSolv.tolerance
        maxStep = parameter.num.pressSolv.maxNumSteps
        
        # Reads type of preconditioner
        # Alternatives: AMG | ILU | None
        preCondType = parameter.num.pressSolv.preConditioning.type
        
        
        # Getting data from elliptic eq. assembler
        A = delfineVar.A
        rhs = delfineVar.rhs
        
        # Get sparse matrix data
        (row,col,data) = A.data()
        n = A.size(0)
        
        # It was commented the line 126 and 129 of the file compressed.py
        # of the scipy package located in /usr/local/lib/python2.6/dist-packages
        # /scipy/sparse in order to avoid some annoying warning messages.
        # If some strange behaviour takes place in this solver routine, this 
        # information may be important for debugging
        Asp = csr_matrix( (data,col.view(),row.view()), shape=(n,n))
        
        # Get right-hand side vector(rhs) data
        b = rhs.data()
        residuals = []
        
        # Solves equation system
        if (preCondType == "amg"):
            # Using AMG Solver as preconditioner with 'solverType' (cg,gmres) as 
            # accelerator. If ilu is defined as solver, it will use pure AMG without
            # accelaration.
            nCoarse = parameter.num.pressSolv.preConditioning.numCoarseLevel
            ml = smoothed_aggregation_solver(Asp, max_levels=nCoarse, max_coarse=1)
            if (solverType != "lu"):
                # Use CG or GMRES acceleration
                x = ml.solve(b,tol=tolerance, maxiter=maxStep, cycle='V', accel=solverType,residuals=residuals)
            else:
                # No accelaration (stand-alone AMG)
                x = ml.solve(b,tol=tolerance, maxiter=maxStep, cycle='V', residuals=residuals)
            
            print ml
        
        elif (preCondType == "none") :
            # Iterate without preconditioner
            if (solverType == "cg"):
                # Using conventional Conjugate Gradients Solver
                (x, flag) = cg(Asp,b, maxiter=maxStep, tol=tolerance,  residuals=residuals)
            elif (solverType == "gmres"):
                # Using conventional Generalized Minimum Residual Method Solver
                (x, flag) = gmres(Asp,b, maxiter=maxStep, tol=tolerance,  residuals=residuals)
            elif (solverType == "lu"):
                # Using a direct LU solver
                # (still pending, to be done with dolfin solver schema, not pyamg)
                print "Error(7):"
                print "Direct solver still not available, use cg or gmres instead"
                print " "
                sys.exit(1)
                
        elif (preCondType == "ilu"):
                # Using a ILU preconditioner
                # (still pending, to be done with dolfin solver schema, not pyamg)
                print "Error(8):"
                print "ILU Preconditioner still not available, use amg or none instead"
                print " "
                sys.exit(1)
        
        # Print residuals history
        residuals = residuals/residuals[0]
  
        # Define return parameters
        delfineVar.x = x
        delfineVar.residuals = residuals
############################################################
