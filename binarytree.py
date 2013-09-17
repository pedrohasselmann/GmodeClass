#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: support.py
# Shell and config manipulation
# Author: Pedro H. A. Hasselmann

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

from numpy import array, arange, concatenate, all, insert, vectorize, vstack
from collections import deque
from itertools import izip

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

class BTree:
   ''' Convert scipy Linkage matrix to Binary Tree for easily cutting nodes'''
   
   def __init__(self, D):
      self.Y  = self.linkage(D)
      self.tree  = self.ordering_by_root()
      self.order = self.ordering_by_size()
      #print(self.Y)

   def linkage(self, D):
      ''' linkage matrix with derivative node n+i'''
  
      import scipy.cluster.hierarchy as sch
    
      Y = sch.median(D)
      #Z = sch.to_tree(Y, rd=True)
      N = D.shape[0]
      Y = concatenate((Y, arange(N,2*N-1).reshape((Y.shape[0], 1))), 1)
   
      return Y
   
   def ordering_by_size(self):      
      ''' Dictionary of nodes by order'''

      Y = self.Y
      sizes = iter(sorted(set(Y[:, -2])))
      order = dict()
      # Pairs' nodes:
      n = 0
      order[n] = list(Y[Y[:,-2] == 2e0][:, -1])
      s = sizes.next()
      # Nodes of higher order, where order is defined by total number of leaves:
      while True:
          try:
            #print(n, s)
            s = sizes.next()
            n += 1
            order[n] = list(Y[Y[:,-2] == s][:,-1])
            Y = filtering(Y, order[n], -1)
            #print(order[n].size)
          except StopIteration:
            break
   
      return order
   
   def ordering_by_root(self):
      ''' Get all leaves over a group of nodes.
          Nodes are dictionaries inside lists.'''

      Y = self.Y
      N = Y.shape[0]
      opp = {1:0,0:1}
      n = 0
      
      def get_sequence(node, bseq, order):
         ''' Recursively kernel.'''
         # b = {a:[{c:[]},{d:[]}], b:[{e:[]},{f:[]}]}
         
         if node > N:
           n += 1
           index = node - N -1
           
           order[n] = deque()
           bseq[Y[index][-1]] = deque()

           for branch in [0,1]:
              if Y[index][branch] > N:
                bseq[Y[index][-1]].append({Y[index][branch]:None})
                order[n].append(Y[index][branch])
              
                get_sequence(Y[index][branch], bseq[Y[index][-1]][branch], order)

              if Y[index][branch] < N:
                bseq[Y[index][-1]].append(Y[index][branch])
       
      tree = {}
      order = {0: Y[-1][-1]}
      get_sequence(Y[-1][-1], tree, order)
      return tree
             
   def cutting(self, branches):
      ''' Get all leaves over a group of nodes.'''

      Y = self.Y
      N = Y.shape[0]
      opp = {1:0,0:1}
      
      def get_leaves(node, leaves):
         ''' Recursively kernel.'''
      
         if node > N:
           index = node - N -1
    
           for branch in [0,1]:
              if Y[index][branch] > N:
                get_leaves(Y[index][branch], leaves)

              if Y[index][branch] < N:
                leaves.append(int(Y[index][branch]))
                #get_leaves(Y[index][opp[branch]], leaves)   
      
      ##################################################
      
      leaves = dict()

      for node in branches:
         leaves[node] = deque()
         get_leaves(node, leaves[node])
         leaves[node] = list(set(leaves[node]))
    
      return leaves

   def plot(self, group, stats, x = [0.36, 0.47, 0.62, 0.75, 0.89]):
      ''' Plot all leaves over a node'''
      
      import matplotlib.pyplot as plt
    
      Y = array(map(lambda x: x[0], stats))
      e = array(map(lambda x: x[1], stats))
      Y = insert(Y, 1, 1e0, axis=1)[group]
      e = insert(e, 1, 0e0, axis=1)[group]
      
      [plt.errorbar(x, y[0], yerr=y[1], color='k', fmt='o-') for y in izip(Y,e)]
      plt.show()
      
     
if __name__ == "__main__":

  from file_module import unpickle
  test = "q1.2_u0.7_m0.3_MOC3qnum"
  bt = BTree(unpickle(test,"D2"))
              