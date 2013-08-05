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

    out = open(path.join("TESTS","tests.txt"),"w")

    text = '\n'.join(list(report))

    out.write(text)

def plot(highlight=[]):
    import matplotlib.pyplot as plt
    from numpy import loadtxt
    
    q1, ulim, Nc, excluded, robust = loadtxt(path.join("TESTS","tests2.txt"), unpack=True, dtype=None, skiprows=1)
    
    name = ['({0:.2f},\n {1:.1f})'.format(*item) for item in zip(q1, ulim)]
    
    plt.figure(figsize=(10,10),dpi=100)
    
    plt.plot(Nc, robust, "k.", markersize=2, label="Gmode Tests ($G_{q_{1}}, upperlimit$)")
    [plt.text(Nc[n], robust[n], texto, fontsize=7, fontweight='black', style='italic') for n, texto in enumerate(name)]

    if len(highlight) != 0:
       for test in highlight:
           highl = plt.plot(test[0], test[1], color="red", marker="o", markersize=4.5)     

    highl[0].set_label("Highlighted tests")
    plt.xlabel("$N_{c}$")
    plt.ylabel("Robustness")
    plt.legend(loc=0, numpoints=1, scatterpoints=1)
    plt.show()
    plt.clf()
    
    

if __name__ == "__main__":
    
   #fit([1.6,1.9,0.1],[0.5,0.8,0.1])
   plot(highlight=[[44, 0.0469], [36, 0.0634]])
    