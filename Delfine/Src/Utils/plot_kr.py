#!/usr/bin/env python
# -*- coding: utf-8 -*-
#############################################################
# File: plot_kr.py
# Function: Auxiliary script to create nice plot figures of relative permeability functions
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
    krow = 0.65
    nw = 2
    no = 2
    krw = 0.65*s**nw
    kro = (1-s)**no

    # Plot Fractional Flux
    plt.plot(s, krw, format[0], label=r'$k_{rw}$', lw=2)
    plt.plot(s, kro, format[1], label=r'$k_{ro}$', lw=2)
    plt.xlabel(r"$S_{wn}$",  fontsize = 20)
    plt.ylabel(r"$k_{r}$",  fontsize = 20)
    plt.title(u"Permeabilidades Relativas",  fontsize=20)
    plt.axis([0, 1.0, 0, 1.0])
    plt.grid(True)
    plt.legend(loc=0)
    plt.annotate('$k_{rw}^{o}$', xy=(1.0, 0.65), xytext=(0.8, 0.65),
    arrowprops=dict(facecolor='black', shrink=0.05), fontsize=20
                      )
    ########################################
#    # set some legend properties.  All the code below is optional.  The
#    # defaults are usually sensible but if you need more control, this
#    # shows you how
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
#    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
#    frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend
#    
#    # see text.Text, lines.Line2D, and patches.Rectangle for more info on
#    # the settable properties of lines, text, and rectangles
#    frame.set_facecolor('0.80')      # set the frame face color to light gray
    plt.setp(ltext, fontsize='xx-large')    # the legend text fontsize
#    plt.setp(llines, linewidth=1.5)      # the legend linewidth
    #########################################
    plt.show()
    
if __name__ == "__main__":
    main()
