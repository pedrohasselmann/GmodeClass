#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: support.py
# Shell and config manipulation
# Author: Pedro H. A. Hasselmann

from numpy import arange, concatenate, all, vectorize, vstack

# Much faster look up than with lists, for larger lists
def filtering(X, crt, i):
   criteria = set(crt)
   @vectorize
   def verify(elem): return elem not in criteria
   return X[verify(X[:,i])]

def select(X, crt, i):
   criteria = set(crt)
   @vectorize
   def verify(elem): return elem in criteria
   return X[verify(X[:,i])]

def group_by_node(D):
   ''' Binary Search Tree'''       
  
   import scipy.cluster.hierarchy as sch
   from numpy import arange, concatenate, all
    
   Y = sch.median(D)
   #Z = sch.to_tree(Y, rd=True)
   N = D.shape[0]
   #leaves = arange(0,N)
   sizes = iter(sorted(set(Y[:, -1])))
   order = {}
   Y = concatenate((Y, arange(N,2*N-1).reshape((Y.shape[0], 1))), 1)
   # Pairs' nodes:
   n = 0
   order[0] = Y[Y[:,-2] == 2e0][:, -1]
   s = sizes.next()
   # Nodes of higher order, where order is defined by total number of leaves:
   while True:
      try:
         print(n, s)
         s = sizes.next()
         n += 1
         order[n] = Y[Y[:,-2] == sizes][:,-1]
         Y = filtering(Y, order[n], -1)
         #print(order[n].size)
      except StopIteration:
         break
   
   return order
   
   
   
      