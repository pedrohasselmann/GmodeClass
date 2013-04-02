#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: modules.py
# author: Pedro Henrique A. Hasselmann
# date: August 2012

# Escrevendo em python3 e usando python2.7:
from __future__ import print_function, unicode_literals, absolute_import, division
from os import path, mkdir
from itertools import izip

def l_to_s(list_of_numbers):
    ''' Return list in a likeable format '''
    return '  '.join(map(str, list(list_of_numbers)))

#
# Print pretty matrix format
#

def pretty_print(matrix):
    try:
       matrix = [[str(x) if x else "0" for x in row] for row in matrix]
       field_length = max(len(x) for row in matrix for x in row)
       output = "\n".join(" ".join("%%%ds" % field_length % x for x in row) for row in matrix)
    except TypeError:
       matrix = [str(x) if x else "0" for x in matrix]
       field_length = max(len(x) for x in matrix)
       output = " ".join("%%%ds" % field_length % x for x in matrix)
    return output

def FileQuery(query,filename):
  
   '''
    Write the ephemutil query in a file.
    WARNING: All MPCephem queries must be 3-level neasted dictionaries.
   '''

   log = open(path.join("Backup",str(filename)+".dat"),'w').write
  
   for item1 in dict.iterkeys(query):
      log("\n{0} \n\n".format(item1))
      for item2 in dict.iterkeys(query[item1]):
         log("{0:8}    ".format(item2)+"  ".join([str(item) for item in dict.itervalues(query[item1][item2])])+"\n")
  
   print(str(filename)+".dat has been generated.")

def FileKeys(dic,filename):

    log = open(path.join("Backup",str(filename)+".dat"),'w').write

    for item1 in dict.iterkeys(dic):
        log("{0} \n".format(item1))

    print(str(filename)+".dat has been generated.")

def WriteIt(pointer,text):

    pointer.write('\n'.join([line for line in text]))
    pointer.close()

def DumpQuery(query,filename):

   '''
      Dump file of a dictionary query produced by ephemutil package.
   '''
   try:
       mkdir(path.join("Backup",""))
   except OSError:
       pass

   import cPickle as pickle
   
   output = open(path.join("Backup","bkp_"+str(filename)+".pkl"),'wb')
   
   pickle.dump(query,output, 2)

   output.close()

def DumpRecover(filename):

   '''
      Recover a dump file produced by ephemutil package.
   '''
   
   import cPickle as pickle
   
   output = open(path.join("Backup",str(filename)+".pkl"),'rb')
   
   rec = pickle.load(output)
   
   output.close()
   
   return rec
