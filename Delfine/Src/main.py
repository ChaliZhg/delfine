#!/usr/bin/env python
#############################################################
# File: main.py
# Function: DELFINE - Program to solve the porous media flow equations using automated 
# scientific programming (Fenics) with the Finite Element Method (Dolfin) and Multigrid (pyAMG)
# Author: Bruno Luna
# Date: 07/02/11
# Modifications Date:
# 
#
#
#
#
############################################################
import sys
from dolfin import *
from Delfine import Delfine

delfine = Delfine(sys.argv)

# Initialize the solver
delfine.initialize()

# Simulate the porous media problem
delfine.drive()

#############################################################
