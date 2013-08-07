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
    
    form = "{0:.2f} {1:.1f}  {2:3}  {3:4}  {4:.4f}".format
    report = deque(["q1  ulim  Nc  E  R"])
    for ulim in arange(*u_range):
        for q1 in arange(*q1_range):
            gmode.q1 = q1
            gmode.ulim = ulim
            gmode.run(realtime_map='n', save='n')
            
            report.append(form(q1, ulim, len(gmode.cluster_stats), len(gmode.excluded), gmode.robust))

    out = open(path.join("TESTS","tests.txt"),"w+")

    text = '\n'.join(list(report))

    out.write(text)

def plot(highlight=[]):
    import matplotlib.pyplot as plt
    from numpy import loadtxt
    
    q1, ulim, Nc, excluded, robust = loadtxt(path.join("TESTS","tests2.txt"), unpack=True, dtype=None, skiprows=1)
    
    name = ['({0:.2f},\n {1:.1f})'.format(*item) for item in zip(q1, ulim)]
    
    plt.figure(figsize=(10,10),dpi=100)
    
    plt.plot(Nc, robust, "k.", markersize=4, alpha=0.7, label="Gmode Tests ($G_{q_{1}}, upperlimit$)")

    if len(highlight) != 0:
       for n, texto in enumerate(name):
           if all([True if item[0] == q1[n] and item[1] == ulim[n] else False for item in highlight]):
              highl = plt.plot(Nc[n], robust[n], color="red", marker="o", markersize=4.5)
              plt.text(Nc[n], robust[n], texto, fontsize=8, fontweight='black', style='italic')

    highl[0].set_label("Highlighted tests")
    plt.xlabel("$N_{c}$")
    plt.ylabel("Robustness")
    plt.legend(loc=0, numpoints=1, scatterpoints=1)
    plt.show()
    plt.clf()
    
    

if __name__ == "__main__":
    
   #fit([1.9,2.0,0.1],[0.5,1.0,0.1])
   plot(highlight=[[2.5, 0.5]])
    