# -*- coding: utf-8 -*-
# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, cos, sin
import random
import ROOT
from ROOT import *
import sys
import itertools
from itertools import *
from optparse import OptionParser
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
    for i in range(len(hists)):
        print("name = ",hists[i].GetName())
        hists_legend.AddEntry(hists[i],hists[i].GetName(),"lp") 

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
    for i in range(len(hist)):
        integral = hist[i].Integral()
        hist[i].Scale(1./integral)

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

def Btagging(tree,fitting):
    csv = np.array([])
    for i in range(tree.n_jet):
        np.append(csv, tree.Jet_csv[i])
    csv_f = np.where(csv > fitting)
    return csv_f[0]

def mu(W_mass,pt_l,pt_v, delta_phi):
    return ((W_mass**2)/2 + pt_l * pt_v * cos(delta_phi))

def MET_pz(W_mass, pt_l , pt_v, pz_l, E_l, l_phi, v_phi): # l =  muon, while v =  neutrino or MET
    delta_phi = l_phi - v_phi
    mu_f = mu(W_mass,pt_l,pt_v,delta_phi)
    square_root = ((((mu_f**2) * (pz_l**2))/(pt_l**4))  - (((E_l**2)*(pt_v**2) - (mu_f**2))/(pt_l**2)))
    part_one = ((mu_f)*(pz_l)/(pt_l**2))

    if square_root < 0:
        return part_one
    else:
        random_number = random.randint(0,1) # can be improved 
        if random_number  == 0:
            return part_one + sqrt(square_root)
        else:
            return part_one - sqrt(square_root)

#basic trapizoid rule
def trapizoid_rule(x,f,N):
    integral  = 0.0
    for i in range(N):
        A_k = 0.5*(x[i+1] - x[i])*(f[i] + f[i+1])
        integral +=A_k
    return integral     

def MET_eta(p_z,p_t):
    return asinh((p_z)/(p_t))

def MET_mass(p_t,phi): #this is assuming m = 0 for a neutrino
    return (p_t*cos(phi))**2 + (p_t*sin(phi))**2


files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")

dir = file1.Get("mmtree")
tree = dir.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.85,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)

colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["Top Mass" ]  

hist_array_top_mass = MakeHist(1,colors,names,"Top Mass reconstructed from the Muon, Met, associated Jet" ,"Mass (GeV)","Normalized Instances",300,0,50)



for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    if tree.n_mu > 0:
        if tree.n_jet > 0:
            muon_vector = []
            muon_vector.append(vec_creator(tree.Muon_pt[0],tree.Muon_eta[0],tree.Muon_phi[0],tree.Muon_m[0]))
            E_l = muon_vector[0].Energy()
            Pz_l = muon_vector[0].Pz()
            METpz = MET_pz(80.377, tree.Muon_pt[0],tree.PFMet_Pt,Pz_l,E_l,tree.Muon_phi[0],tree.PFMet_Phi)
            Met_et =MET_eta(METpz,tree.PFMet_Pt)
            Met_m = MET_mass(tree.PFMet_Pt, tree.PFMet_Phi) 

            Met_vector = vec_creator(tree.PFMet_Pt, Met_et,tree.PFMet_Phi, Met_m)

            MyJets = [] #jets
            for j in range(tree.n_jet):
                MyJets.append(vec_creator(tree.Jet_pt[j],tree.Jet_eta[j],tree.Jet_phi[j],tree.Jet_m[j]))
            
            mj = len(MyJets)
            ind_s = angle_sort(1,mj,muon_vector,MyJets)

            jet_vector = MyJets[ind_s[0][1]]

            Top_vector = muon_vector[0] + Met_vector + jet_vector

            top_mass  = Top_vector.M()

            hist_array_top_mass[0].Fill(jet_vector.M())


HistNorm(hist_array_top_mass)
legends(hist_array_top_mass)

c4=ROOT.TCanvas()
c4.Draw()
c4.cd()
hist_array_top_mass[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c4.Print("Top_mass.png")
c4.cd()
