#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: G-mode.py
# Statistical Classification method G-mode

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

######################## GLOBAL IMPORT #######################################

# Escrevendo em python3 e usando python2.6:
from __future__ import print_function, unicode_literals, absolute_import, division

import os
from time import time
from scipy import stats as scp_sts
from collections import deque
from numpy import array, zeros, median, diagonal, sqrt, genfromtxt, loadtxt, float64, all
from numpy import sum as asum
from gmode_module import stats, shortcut, Invert, G, hyp_test, free, cov
from file_module import l_to_s, pretty_print, WriteIt 
from plot_module import plot_map

# Global variables:
pathjoin = os.path.join

# *s --> Interable with more than one element
# *t --> total sample statistic

########################## Shell Parameters ##################################
def main():
    import optparse

    parser = optparse.OptionParser()
    parser.add_option('-i','--in', dest="arq", default=pathjoin('Barucci_test','Birlan_sample.dat'), type="str", help="File name and path")
    parser.add_option('--q1', action="store", dest="q1", default=2.5, type="float", help="Confidence level parameter 1")
    parser.add_option('--q2', action="store", dest="q2", default=2.5, type="float", help="Confidence level parameter 2")
    parser.add_option('-z','--zones', action="store", dest="zones", default=1, type="int", help="Number of zones")
    parser.add_option('-v','--vlim', action="store", dest="vlim", default=1e0, type="float", help="Limited Deviation")
    parser.add_option('-m','--minlim', action="store", dest="minlim", default=1e0, type="float", help="Minimum Deviation")
    parser.add_option('-n','--name', action="store", dest="name", default='birlan', type="str", help="Label for output customization")

    parser.add_option('-f','--from', action="store", dest="test", default=None, type="str")
    parser.add_option('-c','--clump', action="store", dest="clump", default=None, type="str")

    par, remain = parser.parse_args()

    return par

def make_dir(pth):
    try:
        os.mkdir(pth)
    except OSError:
        pass    

########################################################################
####################### Gmode Python Class #############################
########################################################################

class Gmode:

     __version__ = 3.5

     def __str__(self):
         return ' #########################################################################################################\n   \
                  ################ Statistical Classification Method G-mode - version 3.5 for Python 2.7.2 ################\n   \
                  ################         Program Developer: Pedro Henrique A Hasselmann                  ################\n   \
                  ################       Method Developer: A. I. Gavrishin and A. Coradini                 ################\n   \
                  #########################################################################################################\n   \
                  WARNING: Minimal python packages dependencies: Numpy 1.5, Scipy 0.9, matplotlib 1.0.1\n                       \
                  The input file must be formatted as --> Designation / unique ID / variables / errors \n '

     def __init__(self):

         if __name__ == '__main__': self.Load()

         make_dir(pathjoin("TESTS",""))

     def Load(self,**arg):

         from math import log10

         self.unmeaning_var = list()

         if __name__ == '__main__':
            print("main")

            self.par = main()

            q1        =  self.par.q1
            zones     =  self.par.zones
            vlim      =  self.par.vlim
            minlim    =  self.par.minlim
            self.name =  self.par.name

         else:
            print("imported")

            q1         = arg["q1"]
            zones      = arg["zones"]
            vlim       = arg["vlim"]
            minlim     = arg["minlim"]
            self.name  = arg["name"]

         if vlim != 1e0:
             self.label = 'q'+str(q1)+'_z'+str(round(log10(zones),2))+'_v'+str(vlim)+'_'+self.name
         elif minlim != 1e0:
             self.label = 'q'+str(q1)+'_z'+str(round(log10(zones),2))+'_m'+str(minlim)+'_'+self.name
         else:
             self.label = 'q'+str(q1)+'_z'+str(round(log10(zones),2))+'_'+self.name

         mypath = pathjoin("TESTS",self.label)
         make_dir(mypath)

         self.clasf    = open(pathjoin(mypath,'gmode1_'+self.label+'.dat'),  'w')
         self.clasf2   = open(pathjoin(mypath,'gmode2_'+self.label+'.dat'),  'w')
         self.briefing = open(pathjoin(mypath,'clump_'+self.label+'.dat'),  'w')
         self.log      = open(pathjoin(mypath,   'log_'+self.label+'.dat'),  'w')

         make_dir(pathjoin(mypath,"plots",""))
         make_dir(pathjoin(mypath,"maps",""))

     ############################ Load Data ####################################

     def LoadData(self,**arg):

         if len(arg) == 0:
            filename  = self.par.arq
         else:
            filename  = arg["file"]

         from time import time
         from operator import getitem, itemgetter

         data = map(list,genfromtxt(filename, dtype=None))

         self.design   = map(itemgetter(0),data)
         self.uniq_id  = map(itemgetter(1),data)

         self.elems  = [array(item[2::2], dtype=float64) for item in data]
         self.errs   = [array(item[3::2], dtype=float64) for item in data]

         #self.elems  = [array(item[2:6], dtype=float64) for item in data]
         #self.errs   = [array(item[6:], dtype=float64) for item in data]

         self.indexs = range(len(self.design))

     ############################ Load Sample ####################################

     def LoadClump(self,**arg):

         if len(arg) == 0:
            T       = self.par.clump
            test    = self.par.test
         else:
            T       = arg["clump"]
            test    = arg["test"]

         if T is not None and test is not None:
            mypath = pathjoin("TESTS",test)
            ind, tax = loadtxt(pathjoin(mypath,'gmode1_'+test+'.dat'), unpack=True, dtype=str, usecols=[0,-1])

            uniq_id, design, elems, indexs = list(), list(), list(), list()
            for i, classify in zip(ind, tax):
                i = int(i)
                if classify == T:
                   uniq_id.append(self.uniq_id[i])
                   design.append(self.design[i])
                   elems.append(self.elems[i])

            # setting in self
            self.uniq_id = uniq_id
            self.design  = design
            self.elems   = elems
            self.indexs  = range(len(self.design))

     ########################### START PROCEDURE #################################

     def Run(self,**arg):

         from kernel import classifying 

         if len(arg) == 0:

            q1      = self.par.q1
            zones   = self.par.zones
            vlim    = self.par.vlim
            minlim  = self.par.minlim

         else:

            q1      = arg['q1']
            zones   = arg['zones']
            vlim    = arg['vlim']
            minlim  = arg['minlim']
            self.Load(**arg)
         
         elems, design, indexs, excluded = list(), list(), list(), list()
         
         design.extend(self.design)
         elems.extend(self.elems)
         #errs.extend(self.errs)
         indexs.extend(self.indexs)

         t0 = time() # Start counting procedure time

         N=len(elems)    # Sample size
         M=len(elems[0]) # Variable size

         # Find best number of grid:
         grid = 2
         while grid**(M) < zones:  grid += 1
         grid -= 1
         
         print('grid: ',grid,"--> ", grid**(M),' < ',zones)

         ########## Sample statistical moments ##############

         ctt, devt, St, Rt = stats(elems)
         #Se = array(cov(self.errs,zeros(M), 1e0))
         Se = median( self.errs, axis=0)
         #Se_St = sqrt(Se)/sqrt(array(St))
         Se_St = Se/devt

         vlim   = vlim   * Se_St
         minlim = minlim * Se_St
         #print(median(self.errs, axis=0)/ devt)
         #print(vlim)
         print(minlim)

         #fa = N * free(Rt)         # Degrees of Freedom of all sample

         groups_syntesis =["Clump   N                mean                      dev"]
         report = deque([" Sample size: "+str(N)+" Variable size: "+str(M)])
         report.append(" A.D.: "+str(devt))
         report.append("Variance Limit: "+str(vlim))
         report.append(" Minimum Deviation: "+str(minlim))
         report.append(" Confidence level q1: "+str(scp_sts.norm.cdf(q1) - scp_sts.norm.cdf(-q1)))
         report.append('grid: '+str(grid)+"-->"+str(grid**(M))+'<'+str(zones))
         #print("Maximum number of excluded: ",int(N*(1e0 - scp_sts.norm.cdf(q1))))

         ######################################################
         
         all_groups   = deque()
         all_stats_groups = deque()

         report.append('############################ Part I : Indentifying Classes and Classifying ################################## \n ')

         ################### Indentify Groups #################
               
         Nc = 0
         while Nc == 0 or N >= (M - 1):
               Nc+=1

               report.append('#################################### Clump '+str(Nc)+' ######################################### \n ')

               group, seed, report = classifying(q1, vlim, minlim, grid, design, array(elems), devt, report)

               Na = len(group)

               if Na > 3 or Na >= len(seed):
                        print("Barycenter size: ",len(seed))
                        print(' N = ',N,'Nc = ',Nc,'Na = ',Na)
 
                        '''try:
                          plot_map(Nc, group, seed, elems, self.label)
                        except IndexError:
                          pass '''
 
                        # Save group member indexes
                        all_groups.append(map(lambda i: indexs[i], group))

                        # save group statistics
                        all_stats_groups.append(shortcut(group, elems))

                        # Exclude group members from the sample:
                        for i in group:  elems[i], design[i], indexs[i] = None, None, None
                        elems  = filter(lambda x: x!=None, elems)
                        design = filter(lambda x: x!=None, design)
                        indexs = filter(lambda x: x!=None, indexs)
                        N = len(indexs)
  
                        #report.append("Statistical Significance")
                        report.append("\nC.T.: "+l_to_s(all_stats_groups[-1][0])+"\nA.D.: "+l_to_s(all_stats_groups[-1][1])+ \
                                        "\nSize: "+str(Na)+"       Left: "+str(N)+"\n")

                        groups_syntesis.append("T"+str(Nc)+3*" "+str(Na)+3*" "+l_to_s(all_stats_groups[-1][0])+3*" "+l_to_s(all_stats_groups[-1][1]))
     
               else:
                        Nc -= 1
                        # Exclude group members from the sample:
                        if len(seed) > 0 and Na > 0:
                           report.append("Failed Clump: "+l_to_s(map(lambda i: design[i], group)))
                                        
                           excluded.extend(map(lambda i: indexs[i], group))                         
                           
                           for i in group:  elems[i], design[i], indexs[i] = None, None, None

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
         report.append(" Excluded Sample Size: "+str(len(excluded)))

         # Setting in self

         self.t0 = t0
         self.t1  = time()
         self.report = report
         self.groups_syntesis = groups_syntesis
         self.all_groups = all_groups
         self.all_stats_groups = all_stats_groups
         self.excluded  = excluded
         

     ################### Evaluate Variables and discriminate it #####################
     def Evaluate(self,**arg):

         if len(arg) == 0:
            q2    =  self.par.q1
         else:
            q2    = arg["q2"]

         if len(self.all_groups) > 1:

            from eval_variables import verifying

            all_groups = self.all_groups
            elems      = self.elems
            #errs      = self.errs
            report     = self.report

            report.append('\n############################## Part II : Verifying the variable significance ###############################\n')

            report.append("Confidence level q2: "+str(scp_sts.norm.cdf(q2) - scp_sts.norm.cdf(-q2)))

            d2, Gc, D2 = verifying(q2, all_groups, array(elems)/(stats(elems)[1])) 

            report.append('Matrix D2: \n'+pretty_print(D2))

            j = 0
            for i in range(len(elems[0])):
                report.append('\nMatrix Gc for variable '+str(i+1)+10*" "+' Weight: '+str(d2[i].sum()/d2.sum())+'\n'+pretty_print(Gc[i]))

                if all(Gc[i] < q2):
                   self.unmeaning_var.append(i)
                   report.append('\n Variable '+str(i+1)+' is statistically redundant.')
                   print('Variable '+str(i+1)+' is statistically redundant.')
                   j += 1
            
            # setting in self
            
            self.d2 = d2
            self.Gc = Gc
            self.D2 = D2

     ################################ Fulchignoni et al. (2000) Extension ############################

     def Extension(self,**arg):
         from itertools import  imap

         if len(arg) == 0:
            q1      = self.par.q1
         else:
            q1      = arg['q1']

         sample = map(lambda j: self.elems[j], self.excluded)
         
         N = deque()

         for n, st in enumerate(self.all_stats_groups):
             #iS, iR = Invert(st[2]), Invert(st[3])
             #f = free(st[3])
             f  = (st[1].size**2)/asum(st[3])
             iR = f/st[1].size
             iS = st[1]
             selected = filter(lambda x: x != None, \
                               imap(lambda ind, y: hyp_test(q1,f,ind,y,st[0],iS,iR), self.excluded, sample))

             if len(selected) != 0: 
                plot_map(1001+n,self.excluded,selected,self.elems, self.label)
                self.all_groups[n].extend(selected)
                N.extend(selected)
         
         N = set(N)

         self.report.append("\n Reclassified Excluded Sample Size: "+str(len(N)))
         print("Excluded : ", len(sample) - len(N))
         self.report.append("\n Totally Excluded: "+str(len(sample) - len(N)))


     ################# OUTPUT #####################

     # Write classifications into a file:
     def Classification(self):
         all_groups = self.all_groups
         writing    = self.clasf.write
         design     = self.design
         uniq_id    = self.uniq_id
         
         [[writing(str('{0:7} {1:>10} {2:7} T'+str(n+1)+'\n').format(ind,design[ind],uniq_id[ind])) for ind in self.all_groups[n]] for n in range(len(self.all_groups))]
         self.clasf.close()   

     def ClassificationPerID(self):
         from gmode_module import CollapseClassification

         text = deque()
         
         catalogue = CollapseClassification(self.all_groups,self.design)

         form = "{0:>10} {1}".format
         [text.append(form(each,l_to_s(catalogue[each]))) for each in catalogue.keys()]
         
         WriteIt(self.clasf2,text)

     ################### Log #######################

     # Write into a file:
     def WriteLog(self):
         
         WriteIt(self.log,self.report)
         WriteIt(self.briefing,self.groups_syntesis)

     ################### Plot #######################
     
     def Plot(self, lim=[0.2,1.8], norm=[1, 1e0], axis=[0.36,0.47,0.62,0.75,0.89]):
         from plot_module import plot_clump
         
         for n in xrange(len(self.all_groups)):
             elems_group = array(map(lambda j: self.elems[j], self.all_groups[n]))
             
             plot_clump(n+1, self.all_stats_groups[n], elems_group, self.label, lim, norm, axis)
         
     def Dendrogram(self):
         from plot_module import dendrogram
          
         dendrogram(self.D2,self.label)

     def TimeIt(self):
         # Total processing time:
         t = (time() - self.t0)/60e0
         self.report.append('general processing time: '+str(t)+' min')
         print('general processing time: '+str(t)+' min')

# END of the method

if __name__ == '__main__':
  
   gmode  = Gmode()
   load   = gmode.LoadData()
   samp   = gmode.LoadClump()
   run    = gmode.Run()
   ex     = gmode.Extension()
   ev     = gmode.Evaluate()
   col    = gmode.ClassificationPerID()
   end    = gmode.TimeIt()
   classf = gmode.Classification()
   log    = gmode.WriteLog()
   plot   = gmode.Plot()
   dendo  = gmode.Dendrogram()
