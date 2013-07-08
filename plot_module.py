#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: gmode_func.py
# G-mode plotting methods

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

##### IMPORT ######

from os import path
#import scipy.stats as stats
from numpy import array, amin, amax, degrees, sqrt, arctan2, zeros
from numpy import insert as ainsert
from numpy.linalg import eigh
from itertools import tee, izip
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from collections import deque
from support import config

pathjoin = path.join

option = config('config.cfg', 'PlotConfig')

lim   = map(float, option["lim"])

try:
   norm  = [int(option["norm"][0]), float(option["norm"][1])]
except ValueError:
   norm = None

axis   = map(float, option["axis"])
label  = option["label"]
xtitle = option["xtitle"]
ytitle = option["ytitle"]

#
# Plot clusters
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

def plot_cluster(cov, ct, q1, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        ct  : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
       ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = degrees(arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * q1 * sqrt(vals)
    ellip = Arc(xy=ct, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip

def plot_map(Nc, clump, seed, data, q1, ct, cov, link):
    ''' Density distribution of sample superposed with initial seed and identified cluster.'''    
    
    plt.figure(figsize=(6,10),dpi=60)
    
    data  = array(data)
    clump = array(clump)
    seed  = array(seed)
    M = data.shape[1]    

    plt.title('Cluster '+str(Nc))
    
    for n in xrange(1,M):
        i = int((M-1)*100+10+n)
        plt.subplot(i)

        plt.hexbin(data[:,n-1],data[:,n],gridsize=200,bins='log',mincnt=1)
        if seed.size != 0:
           plt.plot(data[clump,n-1], data[clump,n], 'ro', data[seed,n-1], data[seed,n], 'go', markersize=2.5)
           plot_cluster(cov[n-1:n+1, n-1:n+1], ct[n-1:n+1], 1e0, linestyle='solid', linewidth=1)          
           plot_cluster(cov[n-1:n+1, n-1:n+1], ct[n-1:n+1], q1, linestyle='solid', linewidth=2)
            
        plt.xlim(lim[0], lim[1])
        plt.ylim(lim[0], lim[1])
        
        if label != None:
           plt.ylabel(label[n])
           plt.xlabel(label[n-1])

    cb=plt.colorbar(orientation='horizontal',fraction=0.10,pad=0.3,drawedges=False)
    cb.set_label('log10(N)')
    
    plt.savefig(pathjoin("TESTS",link,"maps",str(Nc)+'.png'),format='png')
    plt.clf()

def plot_spectral(n, stats, data, link):
    ''' The spectral variable plot of central tendency and members of a given cluster.'''

    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

    x, y = deque(), deque()
    
    plt.figure(figsize=(10,9),dpi=30)
    #plt.xlim(-1,data.shape[1])
    plt.ylim(lim[0],lim[1])
    plt.xlabel(*xtitle)
    plt.ylabel(*ytitle)
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
          plt.savefig(pathjoin("TESTS",link,"plots",str(len(data))+'_clump'+str(n)+'_'+link+'.png'),format='png')
          #plt.savefig(pathjoin("TESTS",label,'Graph'+str(n)+'.png'),format='png')
       except OverflowError:
          pass
    
    plt.clf()

#
# Cluster Size Distribution
#

def histogram(Y, cluster_sizes, link):  
    ''' Plot histogram of cluster sizes and integrated standard deviation.'''
    
    from collections import OrderedDict as ordict
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    
    cl = ordict(sorted(cluster_sizes.iteritems(), key=lambda x: x[1]))
    #Y  = ordict(sorted(Y.iteritems(), key=lambda x: cl[x[0]]))
    #X  = ordict(zip( range(len(cluster_sizes)), cl.itervalues() ))

    fig = plt.figure(figsize=(15,8),dpi=80)
    #ax  = fig.add_subplot(1,1,1)

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)
    test = host.twinx()
    
    bottom = zeros(len(cl.keys()))

    size = host.bar(cl.keys(), cl.values(), align='center', width=0.4, color='black', label="Cluster Size", bottom=bottom)
    var  = test.plot(Y.keys(), Y.values(), color="red", marker='o', linestyle='-', linewidth=2, label="Cluster $\sigma$")

    host.set_xlabel("Identified Cluster")
    host.set_ylabel("Cluster Size")
    test.set_ylabel("Integrated Cluster $\sigma$")
    #test.set_ylim(0.0, 1.0)
    host.set_yscale('log') 
    
    #host.set_axisbelow(True)
    #host.set_xticks(X.keys(), cl.keys())

    host.legend(title=link)

    plt.savefig(pathjoin("TESTS",link,"plots",'hist_'+link+'.png'),format='png')
    plt.show()
    plt.close("all")    
        
#
# Dendograms
#

def dendrogram(D, link):
    ''' Plot dendrogram of cluster's parity.'''
  
    import scipy.cluster.hierarchy as sch

    lbl = [str(n+1) for n in xrange(D.shape[0])]

    # Compute and plot dendrogram.
    fig = plt.figure(figsize=(8,14),dpi=80)
    ax1 = fig.add_axes([0.04,0.04,0.9,0.9])
    Y = sch.median(D)
    Z1 = sch.dendrogram(Y, orientation='right',labels=lbl)
    ax1.set_xticks([])
    ax1.set_ylabel("Identified Cluster")
    ax1.set_xlabel("Cluster Parity")
    plt.title(link)

    #plt.show()
    plt.savefig(pathjoin("TESTS",link,"plots",'dendrogram_'+link+'.png'),format='png')
    plt.show()
    plt.close("all")


# END
