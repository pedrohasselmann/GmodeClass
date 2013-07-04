#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: G-mode.py
# G-mode Clustering method
# Author: Pedro H. A. Hasselmann

######################## GLOBAL IMPORT #######################################

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division


from os import path
from time import time
from support import make_dir
from scipy.stats import norm as normal
from collections import deque
from numpy import array, sum, sqrt, zeros, genfromtxt, float64, all
from gmode_module import stats, shortcut, Invert, hyp_test, free, cov
from file_module import l_to_s, pretty_print, WriteIt

pathjoin = path.join

# *s --> Interable with more than one element
# *t --> total sample statistic
# *g --> group/cluster statistic   

########################################################################
####################### Gmode Python Class #############################
########################################################################

class Gmode:   
     '''
     G-mode Multivariate Clustering Method
     -------------------------------------
     
     Developer : Pedro Henrique A Hasselman
     Method Developer : A. I. Gavrishin and A. Coradini
     
     WARNING: Minimum dependencies: Numpy 1.5, Scipy 0.9, matplotlib 1.0.1
     
     Parameters
     ----------
     q1 : Float. Input critical value.
     ulim : Float between 0.0 and 1.0. Upper Deviation Limit.
     mlim : Float between 0.0 and 1.0. Minimum Deviation Limit.
     grid : Integer. Number of folds in barycenter_density search.
     
     Input Format
     -----------
     The input file must be formatted as --> Designation / unique ID / variables / errors 

     Function
     --------
     Class Gmode
     Gmode.load_data()
     Gmode.run()
     Gmode.evaluate()
     Gmode.extension()
     Gmode.classification_per_id()
     Gmode.timeit()
     Gmode.classification()
     Gmode.writelog()
     Gmode.plot()
     Gmode.dendrogram()
     Gmode.histogram()
     '''

     __version__ = 1.0

     def __init__(self):

         if __name__ == '__main__':
            from support import main
            #print("main")
            par = main()

            self.filename = par.filename
            self.q1       = par.q1
            self.q2       = par.q2
            self.grid     = par.grid
            self.ulim     = par.ulim
            self.mlim     = par.mlim
            self.name     = par.name

            self.load()

         else:
            #print("imported")

            self.grid     = 2
            self.ulim     = 1e0
            self.mlim     = 1e0
            
         make_dir(pathjoin("TESTS",""))

     def load(self,**arg):

         if len(arg) == 0:
         
            q1    =  self.q1
            ulim  =  self.ulim or 1e0
            mlim  =  self.mlim or 1e0
            name  =  self.name

         else:

            q1    = arg['q1']
            ulim  = arg['ulim']   or 1e0
            mlim  = arg['mlim'] or 1e0
            name  = arg['name']

         if ulim != 1e0 and mlim == 1e0:
            self.label = 'q'+str(q1)+'_u'+str(ulim)+'_'+name
         elif mlim != 1e0 and ulim == 1e0:
            self.label = 'q'+str(q1)+'_m'+str(mlim)+'_'+name
         elif ulim != 1e0 and mlim != 1e0:
            self.label = 'q'+str(q1)+'_u'+str(ulim)+'_m'+str(mlim)+'_'+name
         else:
            self.label = 'q'+str(q1)+'_'+self.name

         mypath = pathjoin("TESTS",self.label)
         make_dir(mypath)

         self.clasf    = open(pathjoin(mypath,'gmode1_'+self.label+'.dat'),  'w')
         self.clasf2   = open(pathjoin(mypath,'gmode2_'+self.label+'.dat'),  'w')
         self.briefing = open(pathjoin(mypath,'clump_'+self.label+'.dat'),  'w')
         self.log      = open(pathjoin(mypath,   'log_'+self.label+'.dat'),  'w')

         make_dir(pathjoin(mypath,"plots",""))
         make_dir(pathjoin(mypath,"maps",""))

     ############################ Load Data ####################################

     def load_data(self,**arg):

         if len(arg) == 0:
            filename  = self.filename
         else:
            filename  = arg["fname"]

         from operator import itemgetter

         data = map(list,genfromtxt(filename, dtype=None))

         self.design   = map(itemgetter(0),data)
         self.uniq_id  = map(itemgetter(1),data)

         self.elems  = [array(item[2::2], dtype=float64) for item in data]
         self.errs   = [array(item[3::2], dtype=float64) for item in data]

         #self.elems  = [array(item[2:6], dtype=float64) for item in data]
         #self.errs   = [array(item[6:], dtype=float64) for item in data]

         self.indexs = range(len(self.design))

     ########################### START PROCEDURE #################################

     def run(self, realtime_map='y', **arg):

         from kernel import classifying
         from plot_module import plot_map

         if len(arg) == 0:

            if __name__ != '__main__': self.load()
            q1        =  self.q1
            grid      =  self.grid
            ulim      =  self.ulim
            mlim    =  self.mlim
            name      =  self.name

         else:

            q1      = arg['q1']
            grid    = arg['grid']
            ulim    = arg['ulim']
            mlim    = arg['mlim']
            name    = arg['name']

            self.load(**arg)

         #################################################   
         
         elems, design, indexs, excluded = list(), list(), list(), list()
         
         design.extend(self.design)
         elems.extend(self.elems)
         #errs.extend(self.errs)
         indexs.extend(self.indexs)

         t0 = time() # Start counting procedure time

         N=len(elems)    # Sample size
         M=len(elems[0]) # Variable size
         
         #print('grid: ',grid,"--> ", grid**(M))

         ##################################################
         
         ctt, devt, St, Rt = stats(elems)
         Se = cov(self.errs/devt, zeros(M), 1e0)
         
         mlim = mlim * Se

         #print(Se)

         ################# Write Into Log #################
         
         gmode_clusters =["Clump   N                median                      st. dev."]
         report = deque([" Sample size: "+str(N)+" Variable size: "+str(M)])
         report.append(" S.D.: "+str(devt))
         report.append("Upper Limit: "+str(ulim))
         report.append(" Minimum Deviation: "+str(mlim))
         report.append(" Confidence level q1: "+str(normal.cdf(q1) - normal.cdf(-q1)))
         report.append('grid: '+str(grid)+" --> "+str(grid**(M)))
         
         cluster_members   = deque()
         cluster_stats     = deque()

         #plot_map(0, [], [], elems, self.label)
         
         report.append('############################ Part I : Recognize Clusters and Classify ################################## \n ')

         ################### Cluster Recognition #################
               
         Nc = 0
         while Nc == 0 or N >= (M - 1):
               Nc+=1

               report.append('#################################### Clump '+str(Nc)+' ######################################### \n ')

               clump, seed, report = classifying(q1, ulim, mlim, grid, design, array(elems), devt, report)

               Na = len(clump)

               if Na > 3 or Na >= len(seed) and Na != 0:
                        #print("Barycenter size: ",len(seed))
                        #print(' N = ',N,'Nc = ',Nc,'Na = ',Na)
 
                        # Save cluster member indexes
                        cluster_members.append(map(lambda i: indexs[i], clump))

                        # save cluster statistics
                        cluster_stats.append(shortcut(clump, elems))

                        if realtime_map == 'y':
                           try:
                              plot_map(Nc, clump, seed, elems, self.label)
                           except IndexError:
                              pass
                        
                        # Exclude group members from the sample:
                        for i in clump:  elems[i], design[i], indexs[i] = None, None, None
                        elems  = filter(lambda x: x!=None, elems)
                        design = filter(lambda x: x!=None, design)
                        indexs = filter(lambda x: x!=None, indexs)
                        N = len(indexs)
  
                        #report.append("Statistical Significance")
                        report.append("\nC.T.: "+l_to_s(cluster_stats[-1][0])+"\nS.D.: "+l_to_s(cluster_stats[-1][1])+ \
                                        "\nSize: "+str(Na)+"       Left: "+str(N)+"\n")
                        
                        report.append("Cov. Matrix: \n"+str(cluster_stats[-1][2])+"\n")

                        gmode_clusters.append(str(Nc)+3*" "+str(Na)+3*" "+l_to_s(cluster_stats[-1][0])+3*" "+l_to_s(cluster_stats[-1][1]))
     
               else:
                        Nc -= 1
                        # Exclude clump members from the sample:
                        if len(seed) > 0 and Na > 0:
                           report.append("Failed Clump: "+l_to_s(map(lambda i: design[i], clump)))
                                        
                           excluded.extend(map(lambda i: indexs[i], clump))                         
                           
                           for i in clump:  elems[i], design[i], indexs[i] = None, None, None

                        elif len(seed) > 0 and Na == 0:
                           report.append("Failed Clump: "+l_to_s(map(lambda i: design[i], seed)))
                                        
                           excluded.extend(map(lambda i: indexs[i], seed))                         
                           
                           for i in seed:  elems[i], design[i], indexs[i] = None, None, None

                        elif len(seed) == 0:
                           break

                        elems  = filter(lambda x: x!=None, elems)
                        design = filter(lambda x: x!=None, design)
                        indexs = filter(lambda x: x!=None, indexs)
                        N = len(indexs)

         excluded.extend(indexs)

         report.append("######################### Excluded ###############################")
         report.append("Excluded Sample Size: "+str(len(excluded)))
         print("Number of Clusters: ", len(cluster_stats))
         print("Excluded Sample Size: ",len(excluded))
         
         # Setting in self
         self.t0              = t0
         self.t1              = time()
         self.report          = report
         self.gmode_clusters  = gmode_clusters
         self.cluster_members = cluster_members
         self.cluster_stats   = cluster_stats
         self.excluded        = excluded
         

     ################### Evaluate Variables and discriminate it #####################

     def evaluate(self,**arg):

         if len(arg) == 0:
            q2    = self.q2
         else:
            q2    = arg["q2"]

         if len(self.cluster_members) > 1:

            from eval_variables import verifying

            cluster_members = self.cluster_members
            elems      = self.elems
            #errs      = self.errs
            report     = self.report

            report.append('\n############################## Part II : Verifying the variable significance ###############################\n')

            report.append("Confidence level q2: "+str(normal.cdf(q2) - normal.cdf(-q2)))

            d2, Gc, D2 = verifying(q2, cluster_members, array(elems)/(stats(elems)[1])) 

            report.append('Matrix D2: \n'+pretty_print(D2))

            j = 0
            for i in range(len(elems[0])):
                report.append('\nMatrix Gc for variable '+str(i+1)+10*" "+' Weight: '+str(d2[i].sum()/d2.sum())+'\n'+pretty_print(Gc[i]))

                if all(Gc[i] < q2):
                   report.append('\n Variable '+str(i+1)+' is statistically redundant.')
                   print('Variable '+str(i+1)+' is statistically redundant.')
                   j += 1
            
            # setting in self

            self.Gc = Gc
            self.D2 = D2

     ###### Fulchignoni et al. (2000) Extension ######

     def extension(self,**arg):
         from itertools import  imap
         from plot_module import plot_map

         if len(arg) == 0:
            q1   = self.q1
         else:
            q1   = arg['q1']

         cluster_members = self.cluster_members
         excluded   = self.excluded
         elems      = self.elems
         sample = map(lambda j: elems[j], excluded)
         
         N = deque()

         for n, st in enumerate(self.cluster_stats):
             iS = Invert(st[2]) #, Invert(st[3])
             f = free(st[3])
             #iR = f/st[1].size
             #iS = st[1]
             selected = filter(lambda x: x != None, \
                               imap(lambda ind, y: hyp_test(len(cluster_members[n]),q1,f,ind,y,st[0],iS), excluded, sample))

             if len(selected) != 0: 
                plot_map(1001+n, excluded, selected, elems, self.label)
                cluster_members[n].extend(selected)
                N.extend(selected)
         
         N = set(N)

         self.report.append("\n Reclassified Excluded Sample Size: "+str(len(N)))
         print("Excluded : ", len(sample) - len(N))
         self.report.append("\n Totally Excluded: "+str(len(sample) - len(N)))


     ################# OUTPUT #####################

     # Write classifications into a file:
     def classification(self):
         cluster_members = self.cluster_members
         writing         = self.clasf.write
         design          = self.design
         uniq_id         = self.uniq_id
         
         [[writing(str('{0:7} {1:>10} {2:7} '+str(n+1)+'\n').format(ind,design[ind],uniq_id[ind])) for ind in cluster_members[n]] for n in range(len(cluster_members))]
         self.clasf.close()   

     def classification_per_id(self):
         from gmode_module import collapse_classification

         text = deque()
         
         catalogue = collapse_classification(self.cluster_members,self.design)

         form = "{0:>10} {1}".format
         [text.append(form(each,l_to_s(catalogue[each]))) for each in catalogue.keys()]
         
         WriteIt(self.clasf2,text)

     ################### Log #######################

     # Write into a file:
     def writelog(self):
         
         WriteIt(self.log,      self.report)
         WriteIt(self.briefing, self.gmode_clusters)

     ################### Plot #######################
     
     def plot(self):
         from plot_module import plot_clump
         from matplotlib.pyplot import close
         
         for n in xrange(len(self.cluster_members)):
             elems_group = array(map(lambda j: self.elems[j], self.cluster_members[n]))
             
             plot_clump(n+1, self.cluster_stats[n], elems_group, self.label)
         
         close("all")
         
     def dendrogram(self):
         from plot_module import dendrogram
          
         dendrogram(self.D2, self.label)
     
     def histogram(self):
         from plot_module import histogram
         
         cluster_sizes, cluster_stats = dict(), dict()
         for n, cluster in enumerate(self.cluster_members): cluster_sizes[n+1] = len(cluster)
         for n, cluster in enumerate(self.cluster_stats): cluster_stats[n+1] = sqrt(sum(cluster[1]**2))
         
         histogram(cluster_stats, cluster_sizes, self.label)

     def timeit(self):
         # Total processing time:
         t = (time() - self.t0)/60e0
         self.report.append('total processing time: '+str(t)+' min')
         print('total processing time: '+str(t)+' min')

# END of the method

if __name__ == '__main__':
  
   gmode  = Gmode()
   gmode.load_data()
   gmode.run()
   gmode.evaluate()
   #gmode.extension()
   gmode.classification_per_id()
   gmode.classification()
   gmode.plot()
   gmode.timeit()
   gmode.writelog()
   gmode.dendrogram()
   gmode.histogram()