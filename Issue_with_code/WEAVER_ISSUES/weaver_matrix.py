# -*- coding: utf-8
import os
import array
from array import *
import numpy as np
import glob
import math
from math import pi, cos, asinh, sqrt
import random
import ROOT
#from ROOT import *
import sys
import itertools
from itertools import *
from optparse import OptionParser
#import cppyy.ll
# Obviously can only be run in a CMSSW Framework
from DataFormats.FWLite import * 
from HLTrigger import *


files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")

dir = file1.Get("tree")
tree = dir.Get("Events")

i = 0
not_indexed =0
b_tagging_matrix = np.zeros((15,10))
not_b_tagging_matrix = np.zeros((15,10))
pt_cuts = np.array([15, 20, 26, 35, 46, 61, 80, 106, 141, 186, 247, 326, 432, 571, 756, 1000])
eta_cuts = np.array([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    i +=1
    pt_value = tree.j_pt
    eta_value = tree.j_eta
    btagging_value = tree.j_nBHadrons

    pt_index_array = np.where(pt_cuts < pt_value)
    eta_index_array = np.where(eta_cuts < eta_value)

    pt_index = len(pt_index_array[0])
    eta_index = len(eta_index_array[0])

    if pt_index < 16 and eta_index <11:
        if btagging_value == 1:
            b_tagging_matrix[pt_index-1][eta_index-1] +=1
        else:
            not_b_tagging_matrix[pt_index-1][eta_index-1] +=1
    else:
        not_indexed +=1



print("The Matrix for jets that are btagged")
print(b_tagging_matrix)
print("the Matrix for jets that are not btagged")
print(not_b_tagging_matrix)
print("Number of jets that are not indexed")
print(not_indexed)
print(i)      
