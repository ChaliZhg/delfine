#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################
# File: plot_fw.py
# Function: Auxiliary script to create nice plot figures of fractional flux functions
# Author: Bruno Luna
# Date: 21/12/11
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

def main():
    "Main function"
    
    # Define auxiliary variables
    format = ['b-', 'r-', 'k^-', 'gd-', 'mx-']

    # Define 1-D domain
    s = np.linspace(0,1,110)
    # Define fractional flux functions
    fw = s**2/(s**2 + (1-s)**2)
    dfwds = 2*s*(1-s)/(s**2 + (1-s)**2)**2

    # Plot Fractional Flux
    plt.plot(s, fw, format[0], label='$f_{w}(S_{w})$', lw=2)
    plt.plot(s, dfwds, format[1], label='${d f_{w}}/{d S_{w}}$', lw=2)
    plt.xlabel(u"$S_{w}$",  fontsize = 20)
    plt.ylabel(u"$f_{w}$; ${d f_{w}}/{d S_{w}}$",  fontsize = 20)
    plt.title(u"Fluxo Fracional e Derivada",  fontsize = 20)
    plt.axis([0, 1.0, 0, 2.05])
    plt.grid(True)
    plt.legend(loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='xx-large')    # the legend text fontsize
    plt.show()
    
if __name__ == "__main__":
    main()
