#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: bestset.py
# Author: Pedro H. A. Hasselmann

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division
from os import path
from collections import deque

def fit(where_data, q1_range, u_range):
    import Gmode as Gm
    from numpy import arange
    
    gmode = Gm.Gmode()
    
    gmode.load_data(filename=where_data)
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

def plot(highlight=None, **lim):
    import matplotlib.pyplot as plt
    from numpy import loadtxt
    
    q1, ulim, Nc, excluded, robust = loadtxt(path.join("TESTS","tests2.txt"), unpack=True, dtype=None, skiprows=1)
    
    name = ['({0:.2f}, {1:.1f})'.format(*item) for item in zip(q1, ulim)]

    plt.figure(figsize=(10,10),dpi=100)

    # Robustness vs Nc
    plt.subplot(210)
    plt.plot(Nc, robust, "ko", markersize=4, alpha=0.7, label="Gmode Tests ($G_{q_{1}}, upperlimit$)")


    for n, texto in enumerate(name):
        #if any([True if item[0] == q1[n] and item[1] == ulim[n] else False for item in highlight]):
        if Nc[n] >= lim["Nc"][0] and Nc[n] <= lim["Nc"][1] and robust[n] <= lim["robust"]:
           highl = plt.plot(Nc[n], robust[n], color="red", marker="o", markersize=7)
           plt.text(Nc[n], robust[n], texto, fontsize=8, fontweight='black', style='italic')  
           highl[0].set_label(texto)
    
    plt.xlabel("$N_{c}$")
    plt.ylabel("Robustness")
    plt.legend(loc=0, numpoints=1, scatterpoints=1, fontsize=9)

    # Excluded number vs Nc
    plt.subplot(211)
    plt.plot(Nc, excluded, "ko", markersize=4, alpha=0.7, label="Gmode Tests ($G_{q_{1}}, upperlimit$)")
    
    for n, texto in enumerate(name):
        #if any([True if item[0] == q1[n] and item[1] == ulim[n] else False for item in highlight]):
        if Nc[n] >= lim["Nc"][0] and Nc[n] <= lim["Nc"][1] and excluded[n] >= lim["excluded"][0] and excluded[n] <= lim["excluded"][1]:
           highl = plt.plot(Nc[n], excluded[n], color="red", marker="o", markersize=7)
           plt.text(Nc[n], excluded[n], texto, fontsize=8, fontweight='black', style='italic')  
           highl[0].set_label(texto)
       
    plt.xlabel("$N_{c}$")
    plt.ylabel("Number of Excluded")
    plt.legend(loc=0, numpoints=1, scatterpoints=1, fontsize=9)
    
    plt.show()
    plt.clf()
    
    

if __name__ == "__main__":
   
   data_path = path.join("SDSSMOC","lists","MOC4_3quartile_refl_num.dat")
   #fit(data_path, [1.2,1.8,0.1],[0.5,0.9,0.1])
   plot(Nc=[40,80],excluded=[1000,2000],robust=0.05)
    