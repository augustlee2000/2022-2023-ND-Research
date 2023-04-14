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

def HistNorm(hist):
    integral = hist.Integral()
    hist.Scale(1./integral)

files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")
tree = file1.Get("Events")
#tree.Print("GenJet*")
#tree.Print()

#h1 = ROOT.TH2F("Parton vs Gen", "Parton vs Gen", 50, 0, 1000, 50, 0, 1000)
#h2 = ROOT.TH2F("Parton vs Scouting", "Parton vs Scouting", 50, 0, 1000, 50, 0, 1000)

hists_legend = ROOT.TLegend(.65,.70,.85,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)




colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["Parton vs Gen", "Parton vs Scouting", "Parton vs AK8","Gen vs Scouting" ]  

hist_array_pt = MakeHist(4,colors,names,"Pt Difference" ,"Pt (GeV)","Normalized Instances",1.5,-1.5,50)
hist_array_mass = MakeHist(4,colors,names,"Mass Difference" ,"Mass (GeV)","Normalized Instances",1.5,-1.5,50)
hist_array_eta = MakeHist(4,colors,names,"Eta Difference" ,"Eta","Normalized Instances",1.5,-1.5,50)
hist_array_angle = MakeHist(3,colors,names,"Angle between Parton and Jet" ,"Angle (rad)","Normalized Instances",.1,0,50)

legends(hist_array_pt[0])
legends(hist_array_pt[1])
legends(hist_array_pt[2])
legends(hist_array_pt[3])


for ientry in range(tree.GetEntries()):
    print("="*10)
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

    MyFatJets = []
    for l in range(tree.nGenJetAK8):
        TVf = ROOT.TLorentzVector(tree.GenJetAK8_pt[l],tree.GenJetAK8_eta[l],tree.GenJetAK8_phi[l],tree.GenJetAK8_mass[l])
        MyFatJets.append(TVf)    
    
    mp = len(MyParticles)
    mj = len(MyJets)
    mg = len(MyGen)
    mf = len(MyFatJets)

    if mp  > 0:

       #looking for best pairs for parton and gen particle
        ind_g = list(map(list,product(range(mp),range(mg))))
        for i in range(len(ind_g)):
            ind_g[i].append(MyParticles[ind_g[i][0]].DeltaR(MyGen[ind_g[i][1]]))  
        ind_g.sort(key=lambda x:x[2])

        

        #looking for bet fit for parton and scouting particle
        ind_s = list(map(list,product(range(mp),range(mj))))
        for i in range(len(ind_s)):
            ind_s[i].append(MyParticles[ind_s[i][0]].DeltaR(MyJets[ind_s[i][1]]))  
        ind_s.sort(key=lambda x:x[2])

        #looking for best fit between parton and fat jet
        ind_f = list(map(list,product(range(mp),range(mf))))
        for i in range(len(ind_f)):
            ind_f[i].append(MyParticles[ind_f[i][0]].DeltaR(MyJets[ind_f[i][1]]))  
        ind_f.sort(key=lambda x:x[2])



        #Look at the pt difference
        
        #parton vs gen
        
        pt_diff = (MyParticles[ind_g[0][0]].Pt() - MyGen[ind_g[0][1]].Pt()) / MyGen[ind_g[0][1]].Pt()

        hist_array_pt[0].Fill(pt_diff)
        
        #Parton vs Scouting
        
        pt_diff = (MyParticles[ind_s[0][0]].Pt() - MyJets[ind_s[0][1]].Pt()) / MyJets[ind_s[0][1]].Pt()

        hist_array_pt[1].Fill(pt_diff)

        #Parton vs fat jet

        pt_diff = (MyParticles[ind_f[0][0]].Pt() - MyFatJets[ind_f[0][1]].Pt()) / MyFatJets[ind_f[0][1]].Pt()

        hist_array_pt[2].Fill(pt_diff)
        
        #Gen vs Scouting 
        
        pt_diff = (MyGen[ind_g[0][1]].Pt() - MyJets[ind_s[0][1]].Pt()) / MyJets[ind_s[0][1]].Pt()
        
        hist_array_pt[3].Fill(pt_diff)

        #Looking at the mass difference
        
        #parton vs gen
        
        mass_diff = (MyParticles[ind_g[0][0]].M() - MyGen[ind_g[0][1]].M()) / MyGen[ind_g[0][1]].M()

        hist_array_mass[0].Fill(mass_diff)
        
        #Parton vs Scouting
        
        mass_diff = (MyParticles[ind_s[0][0]].M() - MyJets[ind_s[0][1]].M()) / MyJets[ind_s[0][1]].M()

        hist_array_mass[1].Fill(mass_diff)

        #Parton vs fat jet

        mass_diff = (MyParticles[ind_f[0][0]].M() - MyFatJets[ind_f[0][1]].M()) / MyFatJets[ind_f[0][1]].M()

        hist_array_mass[2].Fill(mass_diff)
        
        #Gen vs Scouting 
        
        mass_diff = (MyGen[ind_g[0][1]].M() - MyJets[ind_s[0][1]].M()) / MyJets[ind_s[0][1]].M()

        hist_array_mass[3].Fill(mass_diff)
        
        
        #Looking at the Eta difference
        
        #parton vs gen
        
        eta_diff = (MyParticles[ind_g[0][0]].Eta() - MyGen[ind_g[0][1]].Eta()) / MyGen[ind_g[0][1]].Eta()

        hist_array_eta[0].Fill(eta_diff)
        
        #Parton vs Scouting
        
        eta_diff = (MyParticles[ind_s[0][0]].Eta() - MyJets[ind_s[0][1]].Eta()) / MyJets[ind_s[0][1]].Eta()

        hist_array_eta[1].Fill(eta_diff)

        #Parton vs fat jet

        eta_diff = (MyParticles[ind_f[0][0]].Eta() - MyFatJets[ind_f[0][1]].Eta()) / MyFatJets[ind_f[0][1]].Eta()

        hist_array_eta[2].Fill(eta_diff)
        
        #Gen vs Scouting 
        
        eta_diff = (MyGen[ind_g[0][1]].Eta() - MyJets[ind_s[0][1]].Eta()) / MyJets[ind_s[0][1]].Eta()

        hist_array_eta[3].Fill(eta_diff)

        hist_array_angle[0].Fill(ind_g[0][2])
        hist_array_angle[1].Fill(ind_s[0][2])
        hist_array_angle[2].Fill(ind_f[0][2])











HistNorm(hist_array_pt[0])
HistNorm(hist_array_pt[1])
HistNorm(hist_array_pt[2])
HistNorm(hist_array_pt[3])

HistNorm(hist_array_mass[0])
HistNorm(hist_array_mass[1])
HistNorm(hist_array_mass[2])
HistNorm(hist_array_mass[3])

HistNorm(hist_array_eta[0])
HistNorm(hist_array_eta[1])
HistNorm(hist_array_eta[2])
HistNorm(hist_array_eta[3])

HistNorm(hist_array_angle[0])
HistNorm(hist_array_angle[1])
HistNorm(hist_array_angle[2])









ROOT.gStyle.SetOptStat(0000000)

c1=ROOT.TCanvas()
c1.Draw()
c1.cd()
hist_array_pt[3].Draw("hist_same")
hist_array_pt[1].Draw("hist_same")
hist_array_pt[2].Draw("hist_same")
hist_array_pt[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c1.Print("pt_diff.png")
c1.cd()


c2=ROOT.TCanvas()
c2.Draw()
c2.cd()
hist_array_mass[3].Draw("hist_same")
hist_array_mass[1].Draw("hist_same")
hist_array_mass[2].Draw("hist_same")
hist_array_mass[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c2.Print("mass_diff.png")
c2.cd()

c3=ROOT.TCanvas()
c3.Draw()
c3.cd()
hist_array_eta[3].Draw("hist_same")
hist_array_eta[1].Draw("hist_same")
hist_array_eta[2].Draw("hist_same")
hist_array_eta[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c3.Print("eta_diff.png")
c3.cd()    

c4=ROOT.TCanvas()
c4.Draw()
c4.cd()
hist_array_angle[2].Draw("hist_same")
hist_array_angle[1].Draw("hist_same")
hist_array_angle[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c4.Print("angle.png")
c4.cd()  

# ROOT.gStyle.SetOptStat(0000000)
# c4=ROOT.TCanvas("","",800,600)
# c4.Draw()
# c4.cd()
# hist_array[1].Draw("same")
# #hist_array[3].Draw("same")
# #hists_legend.Draw()
# CMSLabel(.65,.85,ROOT.kBlack,"Internal")
# c4.Print("Parton_vs_Scouting_pt.png")
# c4.cd()   
    
    
    
    
#     #double-for-loop(MyJets, MyParticles):
#     #FindClosests: J.DeltaR(P)
   