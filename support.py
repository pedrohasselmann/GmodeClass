#!/usr/bin/python
# -*- coding: utf-8 -*-
# Filename: support.py
# Shell and config manipulation
# Author: Pedro H. A. Hasselmann

import os

########################## Shell Parameters ##################################

def main():
    import optparse

    pathjoin = os.path.join
    
    parser = optparse.OptionParser()
    parser.add_option('-i','--in', dest="filename", default=pathjoin('Barucci_test','Birlan_sample.dat'), type="str", help="path/to/file")
    parser.add_option('--q1', action="store", dest="q1", default=2.5, type="float", help="Confidence level parameter 1")
    parser.add_option('--q2', action="store", dest="q2", default=2.5, type="float", help="Confidence level parameter 2")
    parser.add_option('-g','--grid', action="store", dest="grid", default=3, type="int", help="Grid")
    parser.add_option('-u','--ulim', action="store", dest="ulim", default=1e0, type="float", help="Limited Deviation")
    parser.add_option('-m','--mlim', action="store", dest="mlim", default=1e0, type="float", help="Minimum Deviation")
    parser.add_option('-n','--name', action="store", dest="name", default='birlan', type="str", help="Label for output customization")

    par, remain = parser.parse_args()
         
    return par

def config(filename, section):
    import ConfigParser
    
    Config = ConfigParser.ConfigParser()
    Config.read(filename)
    
    option = dict(Config.items(section))
    
    for key in option.keys():
        option[key] = option[key].split(",")
    
    return option
    
def make_dir(pth):
    try:
        os.mkdir(pth)
    except OSError:
        pass 
             
# END