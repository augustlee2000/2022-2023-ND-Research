# -*- coding: utf-8
import os
import array
from array import *
import numpy as np
import glob
import math
from math import pi
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


def CMSLabel(x,y,color,name):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(ROOT.kBlack)
    l.DrawLatex(x,y,"CMS")
    
    
    delx = (0.115*496*ROOT.gPad.GetWh())/(472*ROOT.gPad.GetWw()) #696
    
    t =ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(42)
    t.SetTextColor(ROOT.kBlack)
    t.DrawLatex(x+delx,y,name)

def Label(x,y,name):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(42)
    l.SetTextSize(0.03)
    l.SetTextColor(ROOT.kBlack)
    l.DrawLatex(x,y,name)

def legends(hists):
    print("name = ",hists.GetName())
    hists_legend.AddEntry(hists,hists.GetName(),"lp") 

def MakeHist(N,colors,names,title,x_axis,y_axis,high,low,bins):
    hist = []
    for i in range(N):
        hists = ROOT.TH1F(names[i],title,bins,low,high)
        hists.SetLineColor(colors[i])
        #hists.SetFillColor(colors[i])
        hists.GetXaxis().SetTitle(x_axis)
        hists.GetYaxis().SetTitle(y_axis)
        hist.append(hists)
    return hist

def vec_creator(phi,eta,mass,pt):
    vec_temp = ROOT.TLorentzVector()
    vec_temp.SetPtEtaPhiM(pt,eta,phi,mass)
    
    return vec_temp

files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")
tree = file1.Get("Events")
#tree.Print("GenJet*")
#tree.Print()


hists_legend = ROOT.TLegend(.65,.70,.90,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)

colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["Parton_vs_Gen_angle", "Parton_vs_Scouting_angle"]  

hist_array = MakeHist(2,colors,names,"Angle Btween Parton and associated jet" ,"Angle (Rad)","Normalized Instances",.1,0,50)
legends(hist_array[0])
legends(hist_array[1])



for ientry in range(tree.GetEntries()):
    #print("="*10)
    tree.GetEntry(ientry)

    MyParticles = []

    for gp in range(tree.nGenPart):
        #if math.fabs(tree.GenPart_pdgId[gp]) in (1,2,3,4,5):# and tree.GenPart_status[gp] in (62,52):
        if tree.GenPart_statusFlags[gp] == 4481 and math.fabs(tree.GenPart_pdgId[gp]) == 5 and math.fabs(tree.GenPart_pdgId[tree.GenPart_genPartIdxMother[gp]]) == 1000024:
            TVc = ROOT.TLorentzVector(tree.GenPart_pt[gp], tree.GenPart_eta[gp], tree.GenPart_phi[gp], tree.GenPart_mass[gp])
            MyParticles.append(TVc)
    
    MyJets = []
    for j in range(tree.nScoutingJet):
        TVj = ROOT.TLorentzVector(tree.ScoutingJet_pt[j],tree.ScoutingJet_eta[j],tree.ScoutingJet_phi[j],tree.ScoutingJet_mass[j])
        MyJets.append(TVj)
    
    MyGen = []
    for l in range(tree.nGenJet):
        TVl = ROOT.TLorentzVector(tree.GenJet_pt[l],tree.GenJet_eta[l],tree.GenJet_phi[l],tree.GenJet_mass[l])
        MyGen.append(TVl)

    MyParticles_pt = []
    for m in range(len(MyParticles)):
        MyParticles_pt.append(MyParticles[m].Pt())
    
    MyJets_pt = []
    for m in range(len(MyJets)):
        MyJets_pt.append(MyJets[m].Pt())

    MyGen_pt = []
    for m in range(len(MyGen)):
        MyGen_pt.append(MyGen[m].Pt())
    
    mp = len(MyParticles)
    mj = len(MyJets)
    mg = len(MyGen)

    if mp  > 0:

        #looking for best pairs for parton and gen particle
        ind = list(product(range(mp),range(mg)))

        pt_diff = []
        for i in range(len(ind)):
            tem = ind[i]
            pt = abs((MyParticles_pt[tem[0]] - MyGen_pt[tem[1]])/MyGen_pt[tem[1]])
            pt_diff.append(pt)   
        min_arg=pt_diff.index(min(pt_diff))
        min_ind=ind[min_arg]
        hist_array[0].Fill(MyParticles[min_ind[0]].DeltaR(MyGen[min_ind[1]]))

        

        #looking for bet fit for parton and scouting particle
        #print("*"*10)
        ind1 = list(product(range(mp),range(mj)))

        pt_diff1 = []
        for i in range(len(ind1)):
            tem = ind1[i]
            pt = abs((MyParticles_pt[tem[0]] - MyJets_pt[tem[1]])/MyJets_pt[tem[1]])
            pt_diff1.append(pt)    
        min_arg1=pt_diff1.index(min(pt_diff1))
        min_ind1=ind1[min_arg1]
        hist_array[1].Fill(MyParticles[min_ind1[0]].DeltaR(MyJets[min_ind1[1]]))

ROOT.gStyle.SetOptStat(0000000)
c1=ROOT.TCanvas("","",800,600)
c1.Draw()
c1.cd()
hist_array[0].Draw("same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c1.Print("Part_Gen_angle.png")
c1.cd()    

ROOT.gStyle.SetOptStat(0000000)
c2=ROOT.TCanvas("","",800,600)
c2.Draw()
c2.cd()
hist_array[1].Draw("same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c2.Print("Part_Scouting_angle.png")
c2.cd()    
        

   