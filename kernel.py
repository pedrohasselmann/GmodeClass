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
from numpy import array, sqrt, amax, amin, eye, ndenumerate
from numpy import any as aany
from numpy import sum as asum
from numpy import round as arround
from plot_module import plot_clump
  
############ Neighboring Central Method: Recognition of Classes and Classification ###########
 
def classifying(q1, vlim, minlim, grid, design, data, devt, report):

    N = data.shape[0]    # Sample Number
    M = data.shape[1]

#_______________________________Whitenning the sample__________________________

    #elems = data
    data = data/devt

#________________________________Barycenter______________________________________

    seed = barycenter_density(data, grid, amax(data, axis=0), amin(data, axis=0), nmin=int(grid*30/M))

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

          if i == 0 and any([True for n, d in enumerate(devg) if d < sqrt(minlim[n,n])]) or aany(devg == 0e0):
             Sg = minlim

          Na_prior = Na

          #plot_clump(i+1, [ctg*devt, devg*devt, Rg], elems[group], pathjoin(label,"plots","Clump"+str(Nc)), lim=[0.2,1.8], norm=[1,1e0], axis=[0.36,0.47,0.62,0.75,0.89])
             

# G hypothesis test:

          iSg = Invert(Sg)
          f = free(Rg)
          #Rg = f/M

          group = filter(lambda x: x != None, \
                  imap(lambda ind, y: hyp_test(N,q1,f,ind,y,ctg,iSg), range(N), data))

          Na = len(group)

          report.append("Run "+str(i)+" Size: "+str(Na)+" A.D.: "+l_to_s(arround(devg, 3))+"\nf: "+str(f)+"\n")

          i += 1

    return group, seed, report

# End of the Central Method
