#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: barycenter.py
# G-mode barycenter

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######
from numpy import array, median, amax, amin, argmax, where, histogramdd, unravel_index, diagonal, fabs
from numpy import all as aall
from numpy import sum as asum
from scipy.spatial import cKDTree
from itertools import tee, izip, imap
from collections import Counter, deque

#  functions

# MAD
def mad(X, K=1.4826):
    return K*median(fabs(X - median(X, axis=0)), axis=0)

# boolean function

boo = lambda x, xlim: x <= xlim[1] and x >= xlim[0]

def mad(X, K=1.4826):
    return K*median(fabs(X - median(X, axis=0)), axis=0)

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

def barycenter_hist(grid, design, data, minsize=5): 
    ''' KNN with histogramdd.'''
    upper = amax(data, axis=0)
    lower = amin(data, axis=0)
    rng   = range(data.shape[1])

    nbin = map(int,array([grid]*data.shape[1]))

    hist, edges = histogramdd(data,bins=nbin,range=tuple(zip(lower, upper)))

    ind = unravel_index(argmax(hist), hist.shape)

    if amax(hist) > minsize:
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

def barycenter_density(data, grid, upper, lower, vmin, dens = 0e0):
    ''' Recursively search for the initial seed on the densest regions of the sample.
        Uses numpy.histogramdd to divide multidimensional samples on regions. '''
        
    rng   = range(data.shape[1])
    # Number of division on each variable
    nbin = map(int,array([grid]*data.shape[1]))
    # Define the number of points (hist) and the limits (edges) on each region.
    hist, edges = histogramdd(data,bins=nbin,range=tuple(zip(lower, upper)))
    # Transform the limits into a paired list.
    limits = array([list(pairwise(edges[i])) for i in rng])
    # Find the indexes of the region limits that contains the max number of points.
    ind = unravel_index(argmax(hist), hist.shape) 
    # Define the zone containing the maximum number of points.
    zone = array([limits[i,j] for i, j in izip(rng, ind)])
    # Calculate density:
    density = amax(hist) / volume(zone)
    # Elements inside the zone:
    inside = filter(lambda x: x != None, \
                    imap(lambda i, y: boolist(i,y,zone), xrange(data.shape[0]), data))
    # If density keep growing and vmin not reached, the function works recursively.
    # Stops when any of thoses conditions are not satisfied anymore.
    if density > dens and aall(mad(data[inside]) > diagonal(vmin)):
       zone = zone.T
       return barycenter_density(data, grid, zone[1], zone[0], vmin, density)

    else:
       return inside
      
