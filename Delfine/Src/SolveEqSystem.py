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

# For tests only###############
from pyamg.krylov._cg import cg
#########################

class SolveEqSystem:
    """Solves the equation system with pyAMG.""" 
    
    def __init__(self):
        pass
    
    def solve_withPyAMG(self, delfineVar):
        """Solves the equation system with pyAMG.
            Input:
            elliptic = AssembleElliptic() [AssembleElliptic.py]
        """
        
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
        
        # Solves equation system
        # OBS: The max_coarse parameters tells the maximum number of Degrees of
        # freedom in the coarsest mesh, and not the number of coarse levels as would
        # be expected (at least by me!).
        
        # Using AMG Solver
        ml = smoothed_aggregation_solver(Asp,max_coarse=1)
        residuals = []
        x = ml.solve(b,tol=1e-10,accel="cg",residuals=residuals)
        
        # Using conventional Conjugate Gradients Solver (use parameter.num to define it
        # automatically on initialization)
        #(x, flag) = cg(Asp,b, maxiter=200, tol=1e-10,  residuals=residuals)
   
        # Print residuals history
        residuals = residuals/residuals[0]
        print ml
        ##print residuals
        
        # Define return parameters
        delfineVar.x = x
        delfineVar.residuals = residuals
############################################################
