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
from numpy import median, matrix, array, sum, sqrt, diagflat, diagonal, ones, zeros, genfromtxt, float64, all, copy, delete
from file_module import l_to_s, pretty_print

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
     
     Algorithm Developer : Pedro Henrique A Hasselman (Hasselmann et al. 2013)
     Original Method Developers : A. I. Gavrishin and A. Coradini (Coradini et al. 1977)
     
     WARNING: Minimal dependencies: Numpy 1.5, Scipy 0.9, matplotlib 1.0.1
     
     Parameters
     ----------
     q1 : Float. Input critical value.
     ulim : Float between 0.0 and 1.0. Upper Std. Deviation Limit.
     mlim : Float between 0.0 and 1.0. Minimum Std. Deviation Limit.
     grid : Integer. Number of folds in barycenter_density search.
     name : Test label.
     
     Input Format
     -----------
     The input file must be formatted as --> Designation / unique ID / variables / errors
     
     Estimator
     ---------
     
     Mahalonobis distance --> (X - b)' * S^-1 * (X - b)
     S --> median covariance distance

     Function
     --------
     Class Gmode
     Gmode.load_data(filename)
     Gmode.run(realtime_map='n', save='y', **arg)
     Gmode.evaluate(q2=None)
     Gmode.extension(q1=None)
     Gmode.correspondence(templates, template_name, q1=None, artifact=None)
     Gmode.robustness()
     Gmode.classification_per_id()
     Gmode.timeit()
     Gmode.classification()
     Gmode.writelog()
     Gmode.plot()
     Gmode.dendrogram()
     Gmode.histogram()
     '''

     __version__ = 1.4

     def __init__(self):

         if __name__ == '__main__':
            from support import main
            print("main")
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
            print("imported")

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
            ulim  = arg['ulim']  or 1e0
            mlim  = arg['mlim']  or 1e0
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
         make_dir(pathjoin(mypath,"plots",""))
         make_dir(pathjoin(mypath,"maps",""))

     ############################ Load Data ####################################

     def load_data(self,**arg):

         if len(arg) == 0:
            filename  = self.filename
         else:
            filename  = arg["filename"]

         from operator import itemgetter

         data = map(list,genfromtxt(filename, dtype=None))

         self.design   = array(map(itemgetter(0),data))
         self.uniq_id  = array(map(itemgetter(1),data))

         self.elems  = array([array(item[2::2], dtype=float64) for item in data])
         self.errs   = array([array(item[3::2], dtype=float64) for item in data])

         self.indexs = range(len(self.design))

     ########################### START PROCEDURE #################################

     def run(self, realtime_map='n', save='y', **arg):

         from kernel import clustering
         from plot_module import plot_map
         from gmode_module import stats, cov, free

         if len(arg) == 0:

            q1     =  self.q1
            grid   =  self.grid
            ulim   =  self.ulim
            mlim   =  self.mlim

         else:

            q1    = arg['q1']
            grid  = arg['grid']
            ulim  = arg['ulim']
            mlim  = arg['mlim']

         if save == 'y': self.load(**arg)

         print(q1, mlim, ulim) #30/free(Rt))
         #################################################   
         
         design = copy(self.design)
         elems  = copy(self.elems)
         #errs.extend(self.errs)
         indexs = copy(self.indexs)

         t0 = time() # Start counting procedure time

         N=len(elems)    # Sample size
         M=len(elems[0]) # Variable size

         ##################################################
         
         ctt, devt, St, r2t = stats(elems)
         #Se = cov(self.errs/devt, zeros(M), 1e0)
         
         #mlim = (mlim**2) * Se
         mlim = matrix(diagflat(mlim*ones(ctt.size)))

         ################# START REPORT #################

         #print('mlim: ',sqrt(diagonal(mlim)))
         #print('Se: ',sqrt(diagonal(Se)))
         clusters_report =["Clump   N                median                      st. dev."]
         report = deque([" Sample size: "+str(N)+" Variable size: "+str(M)])
         report.append(" S.D.: "+str(devt))
         report.append("Upper Limit: "+str(ulim))
         report.append(" Minimum Deviation: "+str(diagonal(mlim)))
         report.append(" Confidence level q1: "+str(normal.cdf(q1) - normal.cdf(-q1)))
         report.append('grid: '+str(grid)+" --> "+str(grid**(M)))

         excluded          = deque()
         failed_seed       = deque()
         cluster_members   = deque()
         cluster_stats     = deque()
         n_failedseeds     = 0

         if realtime_map == 'y': plot_map(0, [], [], elems, q1, [], [], self.label)
         
         report.append('############################ Part I : Recognize Clusters and Classify ################################## \n ')

         ################### Cluster Recognition #################
               
         Nc = 0
         while Nc == 0 or N >= (M - 1):
               Nc+=1
               report.append('#################################### Clump '+str(Nc)+' ######################################### \n ')
               cluster, seed, report, freedom = clustering(q1, ulim, mlim, grid, design, elems/devt, report) # whitten happenning here

               Na = len(cluster)

               if Na > 2 and freedom >= 30:

                     #press = raw_input("press enter") 
                     # Save cluster member indexes
                     cluster_members.append(indexs[cluster])

                     # save cluster statistics
                     cluster_stats.append(stats(elems[cluster]))
                        
                     # Exclude group members from the sample:
                     elems  = delete(elems, cluster, 0)
                     design = delete(design, cluster, 0)
                     indexs = delete(indexs, cluster, 0)

                     if realtime_map == 'y':
                       print(Nc, "Seed size: ",len(seed),'Na= ',Na,' N= ',N,' f= ',freedom)
                       try:
                         plot_map(Nc, cluster, seed, elems, q1, cluster_stats[-1][0], cluster_stats[-1][2], self.label)
                       except IndexError:
                         pass

                     # appending into logs
                     report.append("\nC.T.: "+l_to_s(cluster_stats[-1][0])+"\nS.D.: "+l_to_s(cluster_stats[-1][1])+ \
                        "\nSize: "+str(Na)+"       Left: "+str(N)+"\nCov. Matrix: \n"+str(cluster_stats[-1][2])+"\n")

                     clusters_report.append(str(Nc)+3*" "+str(Na)+3*" "+l_to_s(cluster_stats[-1][0])+3*" "+l_to_s(cluster_stats[-1][1]))
               else:

                     Nc-=1
                     # Exclude clump members from the sample:
                     if len(seed) > 2 and Na > 2: # Has initial seed and members.

                       report.append("Failed Clump: "+l_to_s(design[cluster])) #map(lambda i: design[i], cluster)))
                                        
                       failed_seed.append(set(indexs[cluster])) #map(lambda i: indexs[i], cluster))
                       n_failedseeds = n_failedseeds + len(cluster)
                           
                       elems  = delete(elems, cluster, 0)
                       design = delete(design, cluster, 0)
                       indexs = delete(indexs, cluster, 0)

                     elif len(seed) > 2 and Na < 3: # Has initial seed but no members.

                       report.append("Failed Seed: "+l_to_s(design[seed]))
                                        
                       failed_seed.append(set(indexs[seed]))                        
                       n_failedseeds = n_failedseeds + len(cluster)

                       elems  = delete(elems, seed, 0)
                       design = delete(design, seed, 0)
                       indexs = delete(indexs, seed, 0)

                     elif len(seed) < 3: # It does not have initial seed.
                       break

               N = len(indexs)


         excluded.extend(indexs)

         report.append("######################### Excluded ###############################")
         report.append("Excluded Sample Size: "+str(len(excluded)))
         report.append("Failed Seeds: "+str(n_failedseeds))
         print("Number of Clusters: ", len(cluster_stats))
         print("Excluded Sample Size: ",len(excluded))
         print("Failed Seeds: ",len(failed_seed))
     
         # Setting in self
         self.t0              = t0
         # logs
         self.report          = report
         self.clusters_report = clusters_report
         # python objects
         self.cluster_members = cluster_members
         self.cluster_stats   = cluster_stats
         self.excluded        = excluded
         self.failed_seed     = failed_seed
         
         # Robustness        
         self.robustness()
         print("Robustness: ", self.robust)
         report.append("Robustness: "+str(self.robust))
         

     ################### Evaluate Variables and discriminate them #####################

     def evaluate(self, q2=None):

         if q2 == None: q2 = self.q1

         if len(self.cluster_members) > 1:

            from eval_variables import distance
            from gmode_module import mad
            from file_module import pickle

            elems      = copy(self.elems)
            dev        = mad(elems, median(elems, axis=0))
            #errs      = self.errs

            self.report.append('\n############################## Part II : Verifying the variable significance ###############################\n')

            self.report.append("Confidence level q2: "+str(normal.cdf(q2) - normal.cdf(-q2)))

            d2, Gc, D2 = distance(self.cluster_members, self.cluster_stats, elems/dev) 

            j = 0
            for i in range(len(elems[0])):
                self.report.append('\nMatrix Gc for variable '+str(i+1)+10*" "+' Weight: '+str(d2[i].sum()/d2.sum())) #+pretty_print(Gc[i]))
                

                if all(Gc[i] < q2):
                   self.report.append('\n Variable '+str(i+1)+' is statistically redundant.')
                   print('Variable '+str(i+1)+' is statistically redundant.')
                   j += 1

            pickle(D2, self.label, "D2")
            pickle(Gc, self.label, "Gc")
                   
     ###### Fulchignoni et al. (2000) Extension ######
 
     def extension(self,q1=None):      
         from itertools import  imap, chain
         from gmode_module import Invert, free, hyp_test
 
         if q1 == None:  q1  = self.q1
           
         cluster_members = self.cluster_members
         sample = self.elems[self.excluded]
          
         self.reclass = deque()
         for n, st in enumerate(self.cluster_stats):
             self.reclass.append(list())
             iS = Invert(st[2]) #, Invert(st[3])
             f = free(st[3])
             size = len(cluster_members[n])
             
             selected = filter(lambda x: x != None, \
                               imap(lambda ind, y: hyp_test(size,q1,f,ind,y,st[0],iS), self.excluded, sample))

 
             if len(selected) != 0: 
                self.reclass[n].extend(selected)
         
         N = set(chain.from_iterable(reclass))
 
         self.report.append("\n Reclassified Excluded Sample Size: "+str(len(N)))
         print("Reclassified : ",len(N))
         self.report.append("\n Totally Excluded: "+str(len(sample) - len(N)))

     ###### GIVE MEANING TO YOUR CLUSTERS ######

     def correspondence(self, templates, template_name, q1=None, artifact=None, var=None):
         ''' Fulchignoni et al. (2000) extension used to give a correspondence to clusters'''
       
         from itertools import  imap
         from gmode_module import Invert, free, hyp_test
         from file_module import pickle, writedict
         import cPickle as pkl

         if q1 == None:  q1  = self.q1

         cluster_members = self.cluster_members
         cluster_stats   = self.cluster_stats
         
         templ = pkl.load(open(templates,'rb')) 
         
         interpretation = dict()

         for n, stat in enumerate(cluster_stats):
             iS = Invert(stat[2][var, :][:, var])
             f = free(stat[3][var, :][:, var])
             size =  len(cluster_members[n])

             selected = filter(lambda x: x != None, \
                               imap(lambda key, y: hyp_test(size, q1, f, key, y[var], stat[0][var], iS), templ.keys(), templ.values()))
           
             interpretation[n+1] = selected

         writedict(interpretation,open(pathjoin("TESTS",self.label,'correspondence_q'+str(q1)+'_'+template_name+'.dat'),  'w'))
         pickle(interpretation, self.label, "correspondence_q"+str(q1)+'_'+template_name)
             
     ############### ROBUSTNESS PARAMETER ##################
     
     def robustness(self):
         ''' Measuring robustness. '''
       
         from gmode_module import robust_parameter as rpar
         
         self.robust = rpar(self.cluster_members, self.cluster_stats, array(self.elems))
         
     ################# OUTPUT #####################

     def classification(self):
         ''' Write classifications into a file '''
         
         cluster_members = self.cluster_members
         f               = open(pathjoin("TESTS",self.label,'gmode1_'+self.label+'.dat'),  'w')
         design          = self.design
         uniq_id         = self.uniq_id
         
         [[f.write(str('{0:7} {1:>10} {2:7} '+str(n+1)+'\n').format(ind,design[ind],uniq_id[ind])) for ind in cluster_members[n]] for n in range(len(cluster_members))]
         f.close()   

     def classification_per_id(self):
         from gmode_module import collapse_classification
         from file_module import writeit

         text = deque()
         
         catalogue = collapse_classification(self.cluster_members,self.design)

         form = "{0:>10} {1}".format
         [text.append(form(each,l_to_s(catalogue[each]))) for each in catalogue.keys()]
         
         writeit(text, open(pathjoin("TESTS",self.label,'gmode2_'+self.label+'.dat'),  'w'))

     ################### Log #######################

     # Write into a file:
     def writelog(self):
         from file_module import writeit, writedict
         from file_module import pickle

         mypath = pathjoin("TESTS",self.label)
             
         writeit(self.report,          open(pathjoin(mypath, 'log_'+self.label+'.dat'),      'w'))
         writeit(self.clusters_report, open(pathjoin(mypath, 'cluster_'+self.label+'.dat'),  'w'))
         
         pickle(self.cluster_stats,   self.label, "cluster_stats")
         pickle(self.cluster_members, self.label, "cluster_members")
         pickle(self.excluded,        self.label, "excluded")

     ################### Plot #######################
     
     def plot(self):
         from plot_module import plot_spectral, mosaic
         from matplotlib.pyplot import close
         
         for n, cl in enumerate(self.cluster_members):
             elems_group = array(map(lambda j: self.elems[j], cl))
             
             plot_spectral(n+1, self.cluster_stats[n], elems_group, self.label)
         
         mosaic(self.cluster_members, self.elems, self.label)
         close("all")
         
     def dendrogram(self):
         from plot_module import dendrogram
         from file_module import unpickle
         
         dendrogram(unpickle(self.label,"D2"), self.label)
     
     def histogram(self):
         from plot_module import histogram
         
         cluster_sizes, cluster_stats = dict(), dict()
         # cluster size dictionary:
         for n, cluster in enumerate(self.cluster_members): cluster_sizes[n+1] = len(cluster)
         # cluster variance dictionary:
         for n, cluster in enumerate(self.cluster_stats): cluster_stats[n+1] = sqrt(sum(cluster[1]**2))
         
         histogram(cluster_stats, cluster_sizes, self.label)

     def timeit(self):
         # Total processing time:
         t = (time() - self.t0)/60e0
         self.report.append('total processing time: '+str(t)+' min')
         print('total processing time: '+str(t)+' min')

# END

if __name__ == '__main__':
   gmode  = Gmode()
   gmode.load_data()
   gmode.run(realtime_map="y", save="y")
   gmode.evaluate()
   gmode.classification_per_id()
   gmode.classification()
   gmode.timeit()
   gmode.writelog()
   gmode.plot()
   gmode.dendrogram()
   gmode.histogram()
