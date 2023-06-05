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

# Tools:-----------------------------------------------------#

def MakeHist(N,colors,names,title,x_axis,y_axis,low,high,bins):
    hist = []
    for i in range(N):
        hists = ROOT.TH1F(names[i],title,bins,low,high)
        hists.SetLineColor(colors[i])
        #hists.SetFillColor(colors[i])
        hists.GetXaxis().SetTitle(x_axis)
        hists.GetYaxis().SetTitle(y_axis)
        hist.append(hists)
    return hist

def HardGet(e, L, H): # shorthand def for getting collection from event
        e.getByLabel(L, H)
        if H.isValid() and len(H.product()) > 0: 
                return H.product()
        return False

def GetDir(events,GL,GH):
    for event in events:
     #if statement is only here because I am being super lazy and don't feel like remaking this event loop to go over only the first event
        event.getByLabel(GL,GH)
        GP = HardGet(event, GL, GH)
        return dir(GP[0])
        break

def Hist(events,GL,GH, variable,hist):
      for event in events:
            GP = HardGet(event, GL,GH)
            for ig in GP:
                  hist.Fill(ig.variable())

# Actual:-----------------------------------------------------#
files = glob.glob( sys.argv[1] ) #+"*root"
events = Events(files)
print(files)
Handle_string = input("Please input the Type: ")
label_string_1 = input("Please input the Module: ")
label_string_2 = input("Please input the Label: ")
label_string_3 = input("Please input the Process: ")

GH = Handle(Handle_string) #vector<Run3ScoutingPFJet>
GL = (label_string_1, label_string_2, label_string_3) 

yes_no = input("Whould you like the directory (y for yes; n for no): ")
#This portion will get you all the variables that can be called!
if yes_no =="y":
    print(GetDir(events, GL, GH))


histogram = input("whould you want a histogram? (y for yes; n for no)")

if histogram == "y":
    variable = input("What variable would you like the histogram of? ") 
    #This is a very fast and dirty way to get a histogram from the file if you don't have access to a TBrowswer or justing being lazy
    events = Events(files)

    colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
    names = [variable]
    x_starting = float(input("What is your minimium x value: "))
    x_ending = float(input("What is your maximum x value: "))
    bins = 500

    # (number of histograms, colors array, name array, title, x axis, y axis, starting value, ending value, bins)
    hist_array = MakeHist(1,colors,names,variable ,variable,"Instances",x_starting,x_ending,bins) 


    #Event loop
    for event in events:
        GP = HardGet(event, GL,GH)
        test = GetDir(events, GL, GH)

        for ig in GP:
                value = getattr(ig,variable)()
                hist_array[0].Fill(value)

    #Drawing the histogram
    c1=ROOT.TCanvas()
    c1.Draw()
    c1.cd()
    hist_array[0].Draw("hist_same")
    c1.Print("AOD_Variable_Output.png") #remember to change this if you are doing multiple runs
    c1.cd()

