#############################################################
# File: DelfineSource.py
# Function: Defines Source Terms for the Elliptic Assemble 
# Author: Bruno Luna
# Date: 29/05/11
# Modifications Date:
# 29/05/11 - Created this module to handle the source terms
#
#############################################################
from dolfin import *
from math import *

class Source(Expression):
    """Assembles the source terms (Wells) for elliptic problem variational form with Dolfin."""
    def eval (self, values, x):
        pass
    ############### Begin - Heterogeneous Anisotropic Ex.##############  
#    def eval(self, values, x):
#        alpha = 1000
#        if (x[0] > 0):
#            values[0] = -2*alpha*exp(x[0])*cos(x[1])
#        else:
#            values[0] = (2*sin(x[1]) + cos(x[1]))*alpha*x[0] + sin(x[1])
    ############### End - Heterogeneous Anisotropic Ex.###############
    ############### Begin - Homogeneous Anisotropic Ex.###############
#    def eval(self, values, x):
#        values[0]  = -2*(1 + pow(x[0],2) + x[0]*x[1] + pow(x[1],2))*exp(x[0]*x[1]) 
    ############### End - Homogeneous Anisotropic Ex.################
    ############### Begin - Homogeneous Isotropic Ex.#################
#    def eval(self, values, x):
#        values[0]  = 2*pow(pi,2)*cos(pi*x[0])*cos(pi*x[1])
    ############### End - Homogeneous Isotropic Ex.#################
