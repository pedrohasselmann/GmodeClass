#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: barycenter.py
# G-mode barycenter

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######
from numpy import Infinity, sqrt, median, array, amax, amin, linspace, fabs, argmax, argmin, where
from numpy import all as aall
from numpy import sum as asum
from scipy.spatial.distance import pdist, squareform as sqf
from scipy.spatial import cKDTree, KDTree
from itertools import combinations, product, tee, izip, imap, chain
from collections import Counter, deque
from gmode_module import stats

#  functions

boo = lambda x, xlim: x <= xlim[1] and x >= xlim[0]
G   = lambda z2, f: sqrt(2e0*z2) - sqrt(2e0*f - 1e0)
z2  = lambda G, f: ((G + sqrt(2e0*f - 1e0))**2)/2e0

# boolean function

def boolist(lim,values,index):
    if all([boo(item[0],item[1]) for item in izip(values,lim)]):
       return index

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

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

def barycenter_KNN2(data, index):

    ''' KNN with densest neighbourhood '''

    ct, dev, R = stats(data[index])
    f = R*data.shape[1]

    tree = KDTree(data[index])
    dist = median(pdist(data[index], 'seuclidean', V=dev))
    #print(G(dist,f))
    
    prior = 0
    for item in data[index]:
        neigh = tree.query_ball_point(item,r=dist)
        if len(neigh) > prior:
           seed = neigh
           prior = len(neigh)

    try:
       return map(lambda j: index[j], seed)
    except UnboundLocalError:
       return []

#
# find the closest pairs and its most common item.
#
def barycenter_pairs(data, index):

    closest = deque()
    pmatrix = pdist(data[index], 'seuclidean', V=stats(data)[1])

    i = 0   
    while i < len(index) and aall(pmatrix == Infinity) != True:
          minimum = amin(pmatrix)
          pair = where(sqf(pmatrix) == minimum)[0]
          closest.append(list(pair))
          pmatrix[argmin(pmatrix)] = Infinity
          i +=1
          
    score      = Counter(list(chain(*closest)))
    common     = score.most_common(1)
    seed       = deque()
    for pair in closest:
        if pair[0] == common[0][0] or pair[1] == common[0][0]:
           seed.extend(pair)

    return set(seed)

#
# Divide the sample in sectors, select the most crowded one,
# and find a clump of items to be seed.
#

def barycenter_hist(grid, design, data):
    
    from numpy import histogramdd, unravel_index
    
    # Best number of bin by Scott (1979).    

    upper = amax(data, axis=0)
    lower = amin(data, axis=0)

    nbin = map(int,array([grid]*data.shape[1]))

    hist, edges = histogramdd(data,bins=nbin,range=tuple(zip(lower, upper)))

    ind = unravel_index(argmax(hist), hist.shape)

    if amax(hist) > 5:
       del hist
       limits = array([list(pairwise(edges[i])) for i in xrange(data.shape[1])])  

       zone = [limits[i][j] for i, j in zip(range(data.shape[1]),ind)]
    
       elems_ind = filter(lambda x: x != None, \
                          imap(lambda i, y: boolist(zone,y,i), xrange(data.shape[0]), data)) 

       return barycenter_KNN(data, elems_ind, design)
    else:
       grid = grid - 1
       if grid > 0:
          return barycenter_hist(grid, design, data)
       else:
          return []

#
# Formal barycenter algorithm: Divide the sample in sector,
# calculate the moments and Gj of each sector and find
# the crowdest one, which will be the initial seed.
#

def barycenter(fa, q1, grid, data):
    from gmode_module import hyp_test, stats

    upper = amax(data, axis=0)
    lower = amin(data, axis=0)

    edges = [list(pairwise(linspace(lower[i],upper[i],grid))) for i in xrange(data.shape[1])] 

    seed, initial, prev_initial_len = list(), list(), 0

    for zones in product(*edges):

        elems_ind = filter(lambda x: x != None, \
                    imap(lambda i, y: boolist(zones,y,i), xrange(data.shape[0]), data)) 

        if len(elems_ind) > 5:
           
           elems = data[elems_ind]
           ct, dev, R = stats(elems)          
           initial = filter(lambda x: x != None, \
                     imap(lambda i, y: hyp_test(q1,fa,i,y,ct,dev,R), elems_ind, elems))
           
           if len(initial) > prev_initial_len:
              seed = initial
              prev_initial_len = len(seed)
    
    return seed
    
