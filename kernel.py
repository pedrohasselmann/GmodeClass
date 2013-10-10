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
from numpy import sqrt, amax, amin, diagonal, diagflat, ndenumerate, any, all
from numpy import sum as asum
from numpy import round as arround
from collections import deque
   
############ Neighboring Central Method: Recognition of Classes and Classification ###########
 
def clustering(q1, ulim, mlim, grid, design, data, report):
    ''' Iterative procedure of clustering'''

    N = data.shape[0]    # Sample Size
    M = data.shape[1]    # Variable size

#________________________________Barycenter______________________________________
    seed = barycenter_density(data, grid, amax(data, axis=0), amin(data, axis=0), sqrt(mlim))
    #seed = barycenter_hist(grid, design, data)
    
    if len(seed) > 2:
       report.append("Barycenter: "+l_to_s(design[seed])) # map(lambda j: design[j], seed)))
       #print("Barycenter: "+l_to_s(map(lambda j: design[j], seed)))
       Na = len(seed)
       cluster = list(seed)
    else:      
       return [], seed, report, 0.0

# ______________________________CLASSIFICATION__________________________________

    i = 0
    f_rec = deque()
    cluster = seed
    index = range(N)

    while (i == 0 or cluster[:] != cluster_prior[:]) and i < 40 and Na > 2: # (round(f,6) != round(f_prior,6) and Na != Na_prior and Na != Na_prior02)

# *c --> cluster statistics

          ct, std, S, r2 = stats(data[cluster])

          cluster_prior = cluster
          
# Replace lower deviations than minimal limit:
          if i == 0 and any(S < mlim):
             S = mlim                    
             std = sqrt(diagonal(S))

          iS = Invert(S)
          f = free(r2)

          f_rec.append(f)
          cluster_prior = cluster

          f_ = min(f_rec)
          
# Once upper limit is reached the iteration is haulted.
          if ulim < 1e0 and any(std >= ulim):
             return cluster, seed, report, Na*f

          report.append("Run "+str(i)+" Size: "+str(Na)+" S.D.: "+l_to_s(arround(std, 3))+"\nf: "+str(f)+"\n")
          
# G hypothesis test:

          if i == 0 or Na*f >= 30:
             cluster = filter(lambda x: x != None, \
                              imap(lambda ind, y: hyp_test(Na,q1,f_,ind,y,ct,iS), index, data))          
          else:
             return cluster, seed, report, Na*f

          Na = len(cluster)

# Discreetly increase std. until the seed start growing by itself.
          #if i == 0 and Na <= Na_prior:
          #   S = S + mlim
          #   std = sqrt(diagonal(S))
          #else:
          i += 1
             
    return cluster, seed, report, Na*f

# End of the Central Method
