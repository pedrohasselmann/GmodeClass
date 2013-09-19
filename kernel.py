#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: central_method.py
'''Neighboring Point Central Method'''

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

from gmode_module import stats, hyp_test, Invert, free
from barycenter import barycenter_hist, barycenter_density
from file_module import l_to_s
from itertools import imap
from numpy import sqrt, amax, amin, eye, diagonal, ndenumerate, any, all
from numpy import sum as asum
from numpy import round as arround
#from plot_module import plot_clump

# Much faster look up than with lists, for larger lists

def filtering(X, criteria):
   @vectorize
   def verify(elem): return elem not in criteria
   return X[verify(X[:])]

def select(X, crt):
   @vectorize
   def verify(elem): return elem in criteria
   return X[verify(X[:])]
   
############ Neighboring Central Method: Recognition of Classes and Classification ###########
 
def clustering(q1, ulim, mlim, grid, design, data, devt, Rt, report):
    ''' Iterative procedure of clustering'''

    N = data.shape[0]    # Sample Size
    M = data.shape[1]    # Variable size

#_______________________________Whitenning the sample__________________________
    data = data/devt

#________________________________Barycenter______________________________________
    seed = barycenter_density(data, grid, amax(data, axis=0), amin(data, axis=0), sqrt(mlim))
    #seed = barycenter_hist(grid, design, data)
    
    if len(seed) > 2:
       report.append("Barycenter: "+l_to_s(design[seed])) # map(lambda j: design[j], seed)))
       #print("Barycenter: "+l_to_s(map(lambda j: design[j], seed)))
       Na = len(seed)
       cluster = list(seed)
    else:      
       return [], seed, report

# ______________________________CLASSIFICATION__________________________________

    i = 0
    Rc, Na_prior, Na_prior02  = eye(M), Na, Na
    index = range(N)

    while (i == 0 or (round(asum(Rc),6) != round(asum(R_prior),6) and Na != Na_prior and Na != Na_prior02)) and i < 40 and Na > 2:

          R_prior = Rc

# *c --> cluster statistics

          ctc, stdc, Sc, Rc = stats(data[cluster])

          Na_prior02 = Na_prior
          Na_prior   = Na
          
# Replace lower deviations than minimal limit:

          if i == 0 and all(Sc < mlim):
             Sc = mlim
                    
             stdc = sqrt(diagonal(Sc))

# Once upper limit is reached the iteration is haulted.
          if ulim < 1e0 and any(stdc >= ulim):
             #print(Na, seed, stdc)
             return cluster, seed, report

# G hypothesis test:

          iSc = Invert(Sc)
          f = free(Rc)
          #Rg = f/M

          if Na*f >= 30:
             cluster = filter(lambda x: x != None, \
                              imap(lambda ind, y: hyp_test(Na,q1,f,ind,y,ctc,iSc), index, data))
          
          else:
             return cluster, seed, report

          Na = len(cluster)
          
          report.append("Run "+str(i)+" Size: "+str(Na)+" S.D.: "+l_to_s(arround(stdc, 3))+"\nf: "+str(f)+"\n")
          
          i += 1

    return cluster, seed, report

# End of the Central Method
