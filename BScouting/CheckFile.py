# -*- coding: utf-8 -*-
# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
import ROOT
from ROOT import *
import sys
import itertools
from itertools import *
from optparse import OptionParser
# Obviously can only be run in a CMSSW Framework
from DataFormats.FWLite import * 
from HLTrigger import *
​
# Tools:-----------------------------------------------------#
def HardGet(e, L, H): # shorthand def for getting collection from event
        e.getByLabel(L, H)
        if H.isValid() and len(H.product()) > 0: 
                return H.product()
        return False
​
# Actual:-----------------------------------------------------#
files = glob.glob( sys.argv[1]+"*root" )
events = Events (files)
print(files)
GH = Handle("vector<reco::GenParticle>")
GL = ("GenParticles", "", "HLT")
​
for event in events:
    print(10*"=")
    GP = HardGet(event, GL, GH)
    for ig in GP:
            print(ig.pdgId())
