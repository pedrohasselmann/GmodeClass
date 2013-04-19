#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: central_method.py
'''Neighboring Point Central Method'''

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

from gmode_module import stats, hyp_test, Invert, free
from barycenter import barycenter_hist, barycenter_density
#from plot_module import plot_group, plot_distribution
from file_module import l_to_s
from itertools import imap
from numpy import array, sqrt, amax, amin, median, diagonal, diagflat, eye, ndenumerate
from numpy import any as aany
from numpy import sum as asum
from numpy import round as arround
from plot_module import plot_clump

#from Gmode_v3 import make_dir
#import os
#pathjoin = os.path.join
  
############ Neighboring Central Method: Recognition of Classes and Classification ###########
 
def classifying(q1, vlim, minlim, grid, design, data, devt, report):

    N = data.shape[0]    # Sample Number
    M = data.shape[1]

#_______________________________Whitenning the sample__________________________

    #elems = data
    data = data/devt

#________________________________Barycenter______________________________________

    #from subroutines_gmode import barycenter_fast
    #seed = barycenter_fast(50, design, data, devt)

    #seed = barycenter_hist(grid, design, data)
    seed = barycenter_density(data, grid, amax(data, axis=0), amin(data, axis=0))

    if len(seed) > 2:
       report.append("Barycenter: "+l_to_s(map(lambda j: design[j], seed)))
       #print("Barycenter: "+l_to_s(map(lambda j: design[j], seed)))
       Na = len(seed)
       group = list(seed)
    else:      
       return [], seed, report

# ______________________________CLASSIFICATION__________________________________

    i = 0
    Rg, Na_prior  = eye(M), Na
    
    #make_dir(pathjoin("TESTS",label,"plots","Clump"+str(Nc),""))

    #(round(asum(Rg),6) != round(asum(R_prior),6) and Na != Na_prior) (round(Rg,6) != round(R_prior,6) and Na != Na_prior)
    while (i == 0 or (round(asum(Rg),6) != round(asum(R_prior),6) and Na != Na_prior)) and i < 20 and Na > 2:

          R_prior = Rg

# *g --> group statistics

          ctg, devg, Sg, Rg = stats(data[group])

# Algorithm to replace zeroth deviations:

          if i == 0 : #and Na == Na_prior and any(zeroth): 
             for n, cov in ndenumerate(Sg):
                 if cov < minlim[n]:
                    Sg[n] = minlim[n]

          Na_prior = Na

          #plot_clump(i+1, [ctg*devt, devg*devt, Rg], elems[group], pathjoin(label,"plots","Clump"+str(Nc)), lim=[0.2,1.8], norm=[1,1e0], axis=[0.36,0.47,0.62,0.75,0.89])
             

# G hypothesis test:

          iSg, iRg = Invert(Sg), Invert(Rg)
          f = free(iRg)

          #f  = (M**2)/asum(Rg)
          #Rg = f/M

          group = filter(lambda x: x != None, \
                  imap(lambda ind, y: hyp_test(N,q1,f,ind,y,ctg,iSg), range(N), data))

          Na = len(group)

          report.append("Run "+str(i)+" Size: "+str(Na)+" A.D.: "+l_to_s(arround(devg, 3))+"\nf: "+str(f)+"\n")

          i += 1

    return group, seed, report

# End of the Central Method
