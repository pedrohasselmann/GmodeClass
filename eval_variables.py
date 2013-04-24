#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: eval_variable.py
# Evaluating variables;

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division


#__________________Part 2: Verifying the variables significance_________________________

def verifying(q2, all_groups, elems):

  Nc =len(all_groups)   # Class Number
  N  = len(elems)
  M  = len(elems[0])     # Variable number

  TINY = 1e-9

# Modules:
  from numpy import zeros, array, float32, sqrt, dot, ravel
  from numpy import sum as asum
  from gmode_module import stats, Invert


# Input:
  d2 = zeros((M,Nc,Nc), dtype=float32)
  Gc = zeros((M,Nc,Nc), dtype=float32)
  D2 = zeros((Nc,Nc),   dtype=float32)

# Hash elements of each group in a array:
  elems_group = list()
  for a in xrange(Nc): elems_group.append(elems[all_groups[a]])
  elems_group = array(elems_group)

# Calculate matrix terms:

  for a in range(0,Nc-1):

      Na = len(all_groups[a])
   
      # Statistics of group a:
      ct_a, dev_a, S_a, R_a = stats(elems_group[a])
      
      iR_a = Invert(R_a)

      for b in xrange(a+1,Nc):

          Nb = len(all_groups[b])
    
          # Statistics of group b:
          ct_b, dev_b, S_b, R_b = stats(elems_group[b])

          iR_b = Invert(R_b)          

          # Degrees of freedom:
          fab = (Nb - 1e0)*(M**2)/asum(R_a)
          fba = (Na - 1e0)*(M**2)/asum(R_b)
          
          # Calculating Z²i(a,b) e Z²i(b,a):
          Z2iab = asum( ( (elems_group[b] - ct_a)/(dev_a + TINY) )**2, axis=0 )
          #Z2iab = asum( (elems_group[b] - ct_a) * ravel( dot(iS_a , (elems_group[b] - ct_a).T) ), axis=0 )
          Z2iba = asum( ( (elems_group[a] - ct_b)/(dev_b + TINY) )**2, axis=0 )
          #Z2iba = asum( (elems_group[a] - ct_b) * ravel( dot(iS_b , (elems_group[a] - ct_b).T) ), axis=0 )

          # Calculating Z²(a,b) e Z2(b,a):

          Z2ab = asum( dot(iR_a, Z2iab) )
          Z2ba = asum( dot(iR_b, Z2iba) )


          for i in xrange(M):

              # Calculating d² e Gc :
              d2[i][a][b] = (Z2iab[i] + Z2iab[i])/(Na + Nb - 1e0)

              Gc[i][a][b] = sqrt(2e0*(Z2iab[i] + Z2iba[i])) - sqrt(2e0*(Na + Nb) - 1e0)

          D2[a][b] = (Z2ab + Z2ba)/(fab + fba - 1e0)
          D2[b][a] = (Z2ab + Z2ba)/(fab + fba - 1e0)


  return d2, Gc, D2

# End of the Part II
