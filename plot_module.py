#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: gmode_func.py
# G-mode plotting methods

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######

from os import path
#import scipy.stats as stats
from numpy import array, linspace, amin, amax
from numpy import insert as ainsert
from itertools import tee, izip, repeat, chain
import matplotlib.pyplot as plt
from collections import deque

pathjoin = path.join

#
# Plot Graphics of groups
#

def plot_distribution(data,*zones):

    plt.figure(figsize=(6,9),dpi=80)
    data = array(data)
    Na = data.shape[1]
    
    for n in xrange(Na):
        i = int(Na*100+10+n+1)
        plt.subplot(i)
        
        #x = linspace(0.5,amax(data, axis=0)[n],50)
        #dist = stats.gaussian_kde(data[:,n])
        
        #plt.plot(x, dist.evaluate(x))
        plt.hist(data[:,n],bins=50)
        if len(zones) != 0: plt.plot(zones[0][n],(1,1),'r*')
        plt.xlim(amin(data, axis=0)[n],amax(data, axis=0)[n])
        
    plt.show()

def plot_map(Nc, clump, seed, data, link, lim=[0.6, 1.8]):
    
    plt.figure(figsize=(6,10),dpi=60)
    
    data  = array(data)
    clump = array(clump)
    seed  = array(seed)
    Na = data.shape[1]    

    plt.title('Clump '+str(Nc))
    
    lbl = ['R_u','R_r','R_i','R_z']
    
    for n in xrange(1,Na):
        i = int((Na-1)*100+10+n)
        plt.subplot(i)

        plt.hexbin(data[:,n-1],data[:,n],gridsize=200,bins='log',mincnt=1)
        plt.plot(data[clump,n-1], data[clump,n], 'ro', data[seed,n-1], data[seed,n], 'go')
        plt.xlim(lim[0], lim[1]); plt.xlabel(lbl[n-1])
        plt.ylim(lim[0], lim[1]); plt.ylabel(lbl[n])

    cb=plt.colorbar(orientation='horizontal',fraction=0.10,pad=0.3,drawedges=False)
    cb.set_label('log10(N)')
    
    plt.savefig(pathjoin("TESTS",link,"maps",str(Nc)+'.png'),format='png')
    plt.clf()

def plot_clump(n, stats, data, label, lim, norm, axis):

    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

    ''' Plot the class unit with members and median values'''

    x, y = deque(), deque()
    
    plt.figure(figsize=(10,9),dpi=30)
    #plt.xlim(-1,data.shape[1])
    plt.ylim(lim[0],lim[1])
    plt.xlabel("variables")
    plt.title('Clump '+str(n)+' Na='+str(data.shape[0]))

    for item in data:
        if norm != None : item = ainsert(item,norm[0],norm[1])
        for xy in izip(pairwise(item),pairwise(axis)):
            y.extend(xy[0])
            y.append(None)
            x.extend(xy[1])
            x.append(None)

    plt.plot(x, y, 'k-')

    # Class Central Tendency e Total Deviation:
    if norm != None: 
       ct, dev = ainsert(stats[0],norm[0],norm[1]), ainsert(stats[1],norm[0],0e0)
    else:
       ct, dev = stats[0], stats[1]

    plt.errorbar(axis,ct,yerr = dev,fmt='o',color='r',ecolor='r',elinewidth=3)

    # Save figures:
    if n == 0:
       plt.show()
    else:
       try:
          plt.savefig(pathjoin("TESTS",label,"plots",str(len(data))+'_clump'+str(n)+'_'+label+'.png'),format='png')
          #plt.savefig(pathjoin("TESTS",label,'Graph'+str(n)+'.png'),format='png')
       except OverflowError:
          pass
    
    plt.clf()

#
# Dendograms
#

def dendrogram(D,label):
    ''' Plot dendrogram between classes'''
  
    import scipy.cluster.hierarchy as sch

    lbl = ['T'+str(n+1) for n in xrange(D.shape[0])]

    # Compute and plot dendrogram.
    fig = plt.figure(figsize=(8,14),dpi=100)
    ax1 = fig.add_axes([0.04,0.04,0.9,0.9])
    Y = sch.median(D)
    Z1 = sch.dendrogram(Y, orientation='right',labels=lbl)
    ax1.set_xticks([])


    plt.savefig(pathjoin("TESTS",label,"plots",'dendrogram_'+label+'.png'),format='png')

def dendrogram_cmap(D, label):
    ''' Plot dendrogram between classes and make a color map for distance matrix'''
  
    import scipy.cluster.hierarchy as sch

    lbl = ['T'+str(n+1) for n in xrange(D.shape[0])]

    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=(8,8),dpi=80)
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    Y = sch.complete(D)
    Z1 = sch.dendrogram(Y, orientation='right',labels=lbl)
    ax1.set_xticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.single(D)
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.Greys)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    plt.savefig(pathjoin("TESTS",label,"plots",'dendrogram_cmap_'+label+'.png'),format='png')

# END
