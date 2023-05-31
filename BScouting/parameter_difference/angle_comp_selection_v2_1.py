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

def vec_creator(pt,eta,phi,mass):
    vec_temp = ROOT.TLorentzVector()
    vec_temp.SetPtEtaPhiM(pt,eta,phi,mass)
    
    return vec_temp

def HistNorm(hist):
    integral = hist.Integral()
    hist.Scale(1./integral)

def angle_sort(n1,n2,v1,v2):
    ind = list(map(list,product(range(n1),range(n2))))
    for i in range(len(ind)):
        ind[i].append(v1[ind[i][0]].DeltaR(v2[ind[i][1]]))  
    ind.sort(key=lambda x:x[2])
    return ind

def difference(v1, v2, command):
    diff = (command(v1) - command(v2)) / command(v2)
    return diff

def Pt(v1):
    return v1.Pt()

def Mass(v1):
    return v1.M()

def Eta(v1):
    return v1.Eta()

files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")
dir = file1.Get("mmtree")
tree = dir.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.85,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)


colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["Parton vs Gen", "Parton vs Scouting", "Parton vs AK8","Gen vs Scouting" ]  

hist_array_pt = MakeHist(4,colors,names,"Pt Difference" ,"Pt Difference","Normalized Instances",1.5,-1.5,50)
hist_array_mass = MakeHist(4,colors,names,"Mass Difference" ,"Mass Difference","Normalized Instances",5,-5,50)
hist_array_eta = MakeHist(4,colors,names,"Eta Difference" ,"Eta Difference","Normalized Instances",1.5,-1.5,50)
hist_array_angle = MakeHist(3,colors,names,"Angle between Parton and Jet" ,"Angle (rad)","Normalized Instances",.4,0,25)

bins = 25
h1 = ROOT.TH2F("Pt Difference Parton vs Gen", "Pt Difference Parton vs Gen", bins, 0, 600, bins, 0, 600)
h2 = ROOT.TH2F("Pt Difference Parton vs Scouting", "Pt Difference Parton vs Scouting", bins, 0, 600, bins, 0, 600)
h3 = ROOT.TH2F("Pt Difference Parton vs AK8", "Pt Difference Parton vs AK8", bins, 0, 600, bins, 0, 600)
h4 = ROOT.TH2F("Pt Difference Gen vs Scouting", "Pt Difference Gen vs Scouting", bins, 0, 600, bins, 0, 600)

legends(hist_array_pt[0])
legends(hist_array_pt[1])
legends(hist_array_pt[2])
legends(hist_array_pt[3])


for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    #print("="*10)

    #We can make a small optimization here with vec_creator function
    MyParticles = []
    for gp in range(tree.n_GenPart):
        #if math.fabs(tree.GenPart_pdgId[gp]) in (1,2,3,4,5):# and tree.GenPart_status[gp] in (62,52):
        if tree.GenPart_StatusFlags[gp] ==23 and math.fabs(tree.GenPart_pdgId[gp]) == 5 and math.fabs(tree.GenPart_genPartIdxMother[gp]) == 6:
            MyParticles.append(vec_creator(tree.GenPart_pt[gp], tree.GenPart_eta[gp], tree.GenPart_phi[gp], tree.GenPart_m[gp]))
    
    MyJets = []
    for j in range(tree.n_jet):
        MyJets.append(vec_creator(tree.Jet_pt[j],tree.Jet_eta[j],tree.Jet_phi[j],tree.Jet_m[j]))
    
    MyGen = []
    for l in range(tree.n_GenJetAK4):
        MyGen.append(vec_creator(tree.GenJetAK4_pt[l],tree.GenJetAK4_eta[l],tree.GenJetAK4_phi[l],tree.GenJetAK4_m[l]))

    MyFatJets = []
    for l in range(tree.n_GenJetAK8):
        MyFatJets.append(vec_creator(tree.GenJetAK8_pt[l],tree.GenJetAK8_eta[l],tree.GenJetAK8_phi[l],tree.GenJetAK8_m[l]))    
    
    mp = len(MyParticles)
    mj = len(MyJets)
    mg = len(MyGen)
    mf = len(MyFatJets)

    if mp  > 0:

        ind_g = angle_sort(mp,mg,MyParticles,MyGen)
        ind_s = angle_sort(mp,mj,MyParticles,MyJets)
        ind_f = angle_sort(mp,mf,MyParticles,MyFatJets)

        #Looking at the pt difference
    
        hist_array_pt[0].Fill(difference(MyParticles[ind_g[0][0]],MyGen[ind_g[0][1]],Pt))
        hist_array_pt[1].Fill(difference(MyParticles[ind_s[0][0]],MyJets[ind_s[0][1]],Pt))
        hist_array_pt[2].Fill(difference(MyParticles[ind_f[0][0]],MyFatJets[ind_f[0][1]],Pt))        
        hist_array_pt[3].Fill(difference(MyGen[ind_g[0][1]],MyJets[ind_s[0][1]],Pt))

        #makeing 2-d Plots

        h1.Fill(MyParticles[ind_g[0][0]].Pt(), MyGen[ind_g[0][1]].Pt())
        h2.Fill(MyParticles[ind_s[0][0]].Pt(),MyJets[ind_s[0][1]].Pt())
        h3.Fill(MyParticles[ind_f[0][0]].Pt(),MyFatJets[ind_f[0][1]].Pt())
        h4.Fill(MyGen[ind_g[0][1]].Pt(), MyJets[ind_s[0][1]].Pt())

        #Looking at the mass difference
        
        hist_array_mass[0].Fill(difference(MyParticles[ind_g[0][0]],MyGen[ind_g[0][1]],Mass))
        hist_array_mass[1].Fill(difference(MyParticles[ind_s[0][0]],MyJets[ind_s[0][1]],Mass))
        hist_array_mass[2].Fill(difference(MyParticles[ind_f[0][0]],MyFatJets[ind_f[0][1]],Mass))        
        hist_array_mass[3].Fill(difference(MyGen[ind_g[0][1]],MyJets[ind_s[0][1]],Mass))

        #Looking at the Eta difference

        hist_array_eta[0].Fill(difference(MyParticles[ind_g[0][0]],MyGen[ind_g[0][1]],Eta))
        hist_array_eta[1].Fill(difference(MyParticles[ind_s[0][0]],MyJets[ind_s[0][1]],Eta))
        hist_array_eta[2].Fill(difference(MyParticles[ind_f[0][0]],MyFatJets[ind_f[0][1]],Eta))        
        hist_array_eta[3].Fill(difference(MyGen[ind_g[0][1]],MyJets[ind_s[0][1]],Eta))

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
hist_array_pt[0].Draw("hist_same")
hist_array_pt[1].Draw("hist_same")
hist_array_pt[3].Draw("hist_same")
hist_array_pt[2].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c1.Print("pt_diff.png")
c1.cd()

c2=ROOT.TCanvas()
c2.Draw()
c2.cd()
hist_array_mass[2].Draw("hist_same")
hist_array_mass[1].Draw("hist_same")
hist_array_mass[0].Draw("hist_same")
hist_array_mass[3].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c2.Print("mass_diff.png")
c2.cd()

c3=ROOT.TCanvas()
c3.Draw()
c3.cd()
hist_array_eta[0].Draw("hist_same")
hist_array_eta[1].Draw("hist_same")
hist_array_eta[3].Draw("hist_same")
hist_array_eta[2].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c3.Print("eta_diff.png")
c3.cd()    

c4=ROOT.TCanvas()
c4.Draw()
c4.cd()
hist_array_angle[0].Draw("hist_same")
hist_array_angle[1].Draw("hist_same")
hist_array_angle[2].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c4.Print("angle.png")
c4.cd()

h1.GetXaxis().SetTitle("Pt (GeV) (Parton)")
h1.GetYaxis().SetTitle("Pt (GeV) (Gen)")
c5=ROOT.TCanvas()
c5.Draw()
c5.cd()
h1.Draw("COLZ")
c5.Print("2d_Pt_parton_gen.png")
c5.cd()

h2.GetXaxis().SetTitle("Pt (GeV) (Parton)")
h2.GetYaxis().SetTitle("Pt (GeV) (Scouting)")
c6=ROOT.TCanvas()
c6.Draw()
c6.cd()
h2.Draw("COLZ")
c6.Print("2d_Pt_parton_scouting.png")
c6.cd()

h3.GetXaxis().SetTitle("Pt (GeV) (Parton)")
h3.GetYaxis().SetTitle("Pt (GeV) (AK8)")
c7=ROOT.TCanvas()
c7.Draw()
c7.cd()
h3.Draw("COLZ")
c7.Print("2d_Pt_parton_ak8.png")
c7.cd()

h4.GetXaxis().SetTitle("Pt (GeV) (Gen)")
h4.GetYaxis().SetTitle("Pt (GeV) (Scouting)")
c8=ROOT.TCanvas()
c8.Draw()
c8.cd()
h4.Draw("COLZ")
c8.Print("2d_Pt_gen_scouting.png")
c8.cd()
