#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: barycenter.py
# G-mode barycenter

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######
from numpy import Infinity, sqrt, median, array, amax, amin, linspace, fabs, argmax, argmin, where, histogramdd, unravel_index
from numpy import all as aall
from numpy import sum as asum
from scipy.spatial import cKDTree, KDTree
from itertools import combinations, tee, izip, imap, chain
from collections import Counter, deque

#  functions

boo = lambda x, xlim: x <= xlim[1] and x >= xlim[0]

# boolean function

def boolist(index, values, lim):
    if all([boo(item[0],item[1]) for item in izip(values,lim)]):
       return index

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def volume(lst):
    p = 1
    for i in lst: p *= i[1] - i[0]
    return p
#
# K-Nearest Neighbours
#

def barycenter_KNN(data, index, label):

    ''' KNN with closest pairs'''
    
    #print(len(index))
    tree = cKDTree(data[index],leafsize=int(data.shape[0]/2))
    
    closest = deque()
    for ind, item in enumerate(data[index]):
        i = int(tree.query(item, k=2)[1][1])
        if label[i] != label[ind]:
           closest.append(i)
        else:
           closest.append(None)
    
    score      = Counter(closest)
    popular    = score.most_common(2)
    common     = popular[0][0]
    neigh      = popular[0][1]
    if common == None: 
       common = popular[1][0]
       neigh  = popular[1][1]
    #print(popular)
    if neigh > 1:
       seed = [i for i, item in enumerate(closest) if item == common]   
       seed.append(common)
       return map(lambda j: index[j], set(seed))       
    else:
       return []

#
# Divide the sample in sectors, select the most crowded one,
# and find the most popular one and his friends to be the seed.
#

def barycenter_hist(grid, design, data):
    
    from numpy import histogramdd, unravel_index  

    upper = amax(data, axis=0)
    lower = amin(data, axis=0)
    rng   = range(data.shape[1])

    nbin = map(int,array([grid]*data.shape[1]))

    hist, edges = histogramdd(data,bins=nbin,range=tuple(zip(lower, upper)))

    ind = unravel_index(argmax(hist), hist.shape)

    if amax(hist) > 5:
       del hist
       limits = array([list(pairwise(edges[i])) for i in rng])  

       zone = [limits[i][j] for i, j in izip(rng,ind)]
    
       elems_ind = filter(lambda x: x != None, \
                          imap(lambda i, y: boolist(i,y,zone), range(data.shape[0]), data)) 

       return barycenter_KNN(data, elems_ind, design)
    else:
       grid = grid - 1
       if grid > 0:
          return barycenter_hist(grid, design, data)
       else:
          return []
    
#
# Initial seed by density. Recursive histogramdd.
#

def barycenter_density(data, grid, upper, lower, dens = 0e0, nmin = 6):

    rng   = range(data.shape[1])

    nbin = map(int,array([grid]*data.shape[1]))

    hist, edges = histogramdd(data,bins=nbin,range=tuple(zip(lower, upper)))

    limits = array([list(pairwise(edges[i])) for i in rng])

    ind = unravel_index(argmax(hist), hist.shape) 

    zone = array([limits[i,j] for i, j in izip(rng, ind)])

    density = amax(hist) / volume(zone)
    
    if density > dens and amax(hist) > nmin:
       zone = zone.T
       return barycenter_density(data, grid, zone[1], zone[0], density, nmin)

    else:
       #print(nmin, amax(hist), density)
       return filter(lambda x: x != None, \
                     imap(lambda i, y: boolist(i,y,zone), xrange(data.shape[0]), data))
      
