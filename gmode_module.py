#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: gmode_func.py
# G-mode method functions

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######

from os import path
from numpy import sqrt, array, eye, matrix, dot, diagonal, diagflat, median, fabs, corrcoef, isnan, ravel
from numpy import all as aall
from numpy import sum as asum
from numpy.linalg.linalg import LinAlgError
from numpy.linalg import det 
from itertools import product, chain
from collections import deque
#import median_stats as md_sts
#from subroutines_gmode import pearson_r

from time import time
t0= time()

pathjoin = path.join
#median = md_sts.median
#mad    = md_sts.median_abs_dev

TINY = 1e-6

#
# Sample Statistics 
#

def free(R):
    return ravel( R.shape[1]/asum(R, axis=0) )

def pearson_R(X):
    r2 = corrcoef(zip(*X))**2

    if aall(isnan(r2)) : 
       r2 = eye(X.shape[1])
    else:
       whereNaN = isnan(r2)
       r2[whereNaN] = 1e0

    return r2

def Robust_R(X, ct, dev):
    ''' Shevlyakov 1997 - On a Robust Estimation of Correlation Coeficient'''
    X = (X - ct)/(dev + TINY)
    r2 = deque()
    
    for i in xrange(X.shape[1]):
        u  = median(fabs(X.T + X[:,i]), axis=1)**2
        v  = median(fabs(X.T - X[:,i]), axis=1)**2
        ri = (u - v)/(u + v)
        r2.append(ri**2)

    r2 = matrix(r2)

    if aall(isnan(r2)) : 
       r2 = matrix(eye(X.shape[1]))
    else:
       whereNaN = isnan(r2)
       r2[whereNaN] = 1e0

    return r2

def mad(X, ct, K=1.4826):
    return K*median(fabs(X - ct), axis=0)

def cov(X, ct, K=1.4826):
    X = X - ct
    return matrix( [median(X.T*X[:,i], axis=1)*K**2 for i in xrange(X.shape[1])] )
    

def stats(X):
               
    #X   = array(X)
    ct  = median(X, axis=0)      
    S   = cov(X, ct)
    dev = sqrt(diagonal(S))
    R   = Robust_R(X, ct, dev)

    return ct, dev, S, R
 

def Invert(A):

    try:
       iA = A.I
    except LinAlgError:
       iA = diagflat( 1e0/(diagonal(A) + TINY) )

    return iA

def shortcut(group, data): 
    return stats(map(lambda j: data[j], group))

#
# G estimator (Abramowitz 1962)
#

def G(N, f, X, ct, iS):

    ''' G parameter --> Transforms a X2 estimator to a Gaussian estimator '''

    #z2  = iR * asum( ( (X - ct) / (iS + TINY) )**2  )     # Z2 estimator
      
    X = X - ct
    z2 = fabs( X * ravel( dot(iS, X ) ) )

    if aall(N*f > 100e0):
       return  sqrt(2e0*z2) -  sqrt(2e0*f - 1e0)

    elif aall(N*f >= 30e0) and aall(N*f <= 100e0):
       return ((z2/f)**(1e0/3) - (1e0 - (2e0/9)*f))/sqrt((2e0/9)*f)
    
    elif aall(N*f < 30):
       return 9e9

#
# G Hypothesis test
#

def hyp_test(N, q1, f, index, x, ct, iS):
    if aall(G(N, f, x, ct, iS) < q1):
       #print(G(N, f, x, ct, iS))
       return index

def CollapseClassification(clusters, ID):

    ''' Returns a dictionary of classification per ID '''

    cat = dict()
    tax = 0
    for group in clusters:
        tax += 1
        t = 'T'+str(tax)
        for item in set(group):
            try:
               if cat.has_key(ID[item]):
                  cat[ID[item]].append(t)
               else:
                  cat[ID[item]] = [t,]
            except KeyError:
               print(item,ID.has_key(item),cat.has_key(item))
               break
    
    return cat



            
