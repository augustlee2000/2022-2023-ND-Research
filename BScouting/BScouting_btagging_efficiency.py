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


# Open the events tree in the the file, python CheckFile.py output.root is the correct file has to be done on the cms
files = glob.glob( sys.argv[1]+"*root" )
events = Events (files)

#giving these more significant variable names could be useful
#The Handle and legend? This allows us to get a certain branch (gen)
GH = Handle("vector<reco::GenParticle>")
GL = ("GenParticles", "", "HLT")#This will have to be changed to fit our needs

#The Handle and legend? This allows us to get a certain branch (sim)
SH = Handle("vector<reco::GenParticle>")
SL = ("GenParticles", "", "HLT")#This will have to be changed to fit our needs

correct_btagged_jets = 0

total_btagged_jets = 0

for event in events:
    print(10*"=")
    GP = HardGet(event, GL, GH)
    SP = HardGet(event, SH, SL)
    for ig in GP: # this is going to get what what type of particles are in each event
        gen_id = ig.pdgId()
    for iS in SP: #This line might be slightly more complicated, I probably need to to get the btagged information 
        #btag_array_i = tree.Jet_btagCSVV2
        #btag_array = np.frombuffer(btag_array_i, dtype=np.float32)
        #btag_sorted = np.where(btag_array > 0.7527)#tight btagging on Deep CSV
        #btag_sorted = btag_sorted[0]

        #the code below will give us an array that tells us which indencies are correctly btagged (this code is going to have to be changed to fit our needs)

        sim_id = iS.pdgId()

    #Now that we have our gen particles and our sim btagged we want to check our efficency rate

    for j in range(btag_sorted.size):
        total_btagged_jets += 1 
        if gen_id[j] == 5:
            correct_btagged_jets += 1

Print("The total correct btagged jets: ", correct_btagged_jets, "Out of: ", total_btagged_jets, "total jets")
Print("This give an efficency of: ", (correct_btagged_jets / total_btagged_jets)*100)





