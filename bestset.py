#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: bestset.py
# Author: Pedro H. A. Hasselmann

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division
from os import path
from collections import deque

def fit(q1_range, u_range):
    import Gmode as Gm
    from numpy import arange
    
    gmode = Gm.Gmode()
    
    gmode.load_data(filename=path.join("SDSSMOC","lists","MOC4_3quartile_refl_num2.dat"))
    gmode.mlim = 0.3
    gmode.grid = 3
    
    form = "{0:.2f} {1:.1f}  {2:3}  {3:3}  {4:.4f}".format
    report = deque(["q1  ulim  Nc  E  R"])
    for ulim in arange(*u_range):
        for q1 in arange(*q1_range):
            gmode.q1 = q1
            gmode.ulim = ulim
            gmode.run(realtime_map='n', save='n')
            
            report.append(form(q1, ulim, len(gmode.cluster_stats), len(gmode.excluded), gmode.robust))

    out = open("tests.txt","w")

    text = '\n'.join(list(report))

    out.write(text)

def plot():
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    
    q1, ulim, Nc, excluded, robust = genfromtxt("tests.txt", unpack=True, dtype=None)
    
    name = ['{0:.2f} {1:.1f}'.format(*item) for item in zip(g1, ulim)]
    
    plt.figure(figsize=(10,10),dpi=70)
    
    plt.plot(robust, Nc, "k.", label="Gmode Tests")
    plt.text(robust, Nc, name)
    
    plt.xlabel("$N_{c}$")
    plt.ylabel("Robustness")
    plt.legend(loc=0)
    plt.show()
    plt.clf()
    
    

if __name__ == "__main__":
    
   fit([2.8,2.9,0.1],[0.7,0.8,0.1])
   plot()
    