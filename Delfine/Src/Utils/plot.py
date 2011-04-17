#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################
# File: plot.py
# Function: Auxiliary script to create nice plot figures of residuals for thesis
# Author: Bruno Luna
# Date: 16/04/11
# Modifications Date:
# 
#
#
#############################################################
import sys
import getopt
import os
import numpy as np
import matplotlib.pyplot as plt

def main(argv):
    "Main function"
    
    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "h", ["help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Get options
    for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()

    # Check that we at least one filename
    if (len(args) == 0):
        usage()
        sys.exit(2)
    
    # Define auxiliary variables
    filelist= args
    format = ['bo-', 'rs-', 'k^-', 'gd-']
    i = 0 # format index
    
    for name in filelist:
        file = open(name, "r")
        t = [ ]
        
        if (i > (len(format)-1)):
            i = 0 # Cycle in format list
        
        for line in file:
            t.append(float(line))
        
        # Blue circles
        plt.plot(t, format[i], label=name.split('.')[0].split('_')[1])
        plt.xlabel(u"Iteração")
        plt.ylabel(u"Resíduo Normalizado")
        plt.title(u"Evolução do Resíduo")
        #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
        plt.semilogy()
        #plt.axis([40, 160, 0, 0.03])
        plt.grid(True)
        plt.legend()
        # Format list advance
        i = i+1
        
    plt.show()
    
def usage():
    "Display usage"
    print """\

Usage: plot.py [OPTIONS]  ... input1.x input2.x input3.x

Options:
 
    -h: Display this text and exit

"""
if __name__ == "__main__":
    main(sys.argv[1:])
