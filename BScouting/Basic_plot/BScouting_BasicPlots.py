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
    

def MET_eta(p_z,p_t):
    return asinh((p_z)/(p_t))


files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")

dir = file1.Get("mmtree")
tree = dir.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.85,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)


colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["Jets Pt", "B-tagged Pt","leading Muon PT"] 
names2 = ["Number of Muons"]
names3 = ["Missing Energy"]
names4 = ["Delta Eta between Met and leading muon"]

hist_array_pt = MakeHist(3,colors,names,"Pt" ,"Pt (GeV)","Normalized Instances",100,0,50)
hist_array_n_muons = MakeHist(1,colors,names2,"Number of Muons","Muons", "Normalized Instatnces", 5,0,5)
hist_array_energy  = MakeHist(1,colors,names3,"Missing Energy in Event","Energy (GeV)", "Normalized Instatnces", 300,0,50)
hist_array_delta_eta  = MakeHist(1,colors,names3,"Delta Eta between MET and Leading Muon","Delta Eta ", "Normalized Instatnces", 10,0,50)

legends(hist_array_pt[0])
legends(hist_array_pt[1])


for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    btagging_ind = Btagging(tree, 0.7527 ) #tight BTagging
    
    if btagging_ind.size > 0:
        hist_array_pt[1].Fill(tree.Jet_pt[btagging_ind[0]])
    
    for j in range(tree.n_jet):
        hist_array_pt[0].Fill(tree.Jet_pt[j])

    if tree.n_mu > 0:
        hist_array_pt[2].Fill(tree.Muon_pt[0])

    hist_array_n_muons[0].Fill(tree.n_mu)


    hist_array_energy[0].Fill(tree.PFMet_Pt)

    if tree.n_mu > 0:
        muon_vector = vec_creator(tree.Muon_pt[0],tree.Muon_eta[0],tree.Muon_phi[0], tree.Muon_m[0])
        E_l = muon_vector.Energy()
        Pz_l = muon_vector.Pz()
        METpz = MET_pz(80.377, tree.Muon_pt[0],tree.PFMet_Pt,Pz_l,E_l,tree.Muon_phi[0],tree.PFMet_Phi)

        hist_array_delta_eta[0].Fill(abs(MET_eta(METpz,tree.PFMet_Pt) - tree.Muon_eta[0]))


HistNorm(hist_array_pt[0])
#HistNorm(hist_array_pt[1])
HistNorm(hist_array_pt[2])

HistNorm(hist_array_n_muons[0])
HistNorm(hist_array_energy[0])
# HistNorm(hist_array_delta_eta[0])

ROOT.gStyle.SetOptStat(0000000)

c1=ROOT.TCanvas()
c1.Draw()
c1.cd()
hist_array_pt[1].Draw("hist_same")
hist_array_pt[0].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c1.Print("pt_Btagging_Jets.png")
c1.cd()

c2=ROOT.TCanvas()
c2.Draw()
c2.cd()
hist_array_n_muons[0].Draw("hist_same")
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c2.Print("Number_of_Muons.png")
c2.cd()

c3=ROOT.TCanvas()
c3.Draw()
c3.cd()
hist_array_energy[0].Draw("hist_same")
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c3.Print("MET_ttbar.png")
c3.cd() 

c4=ROOT.TCanvas()
c4.Draw()
c4.cd()
hist_array_delta_eta[0].Draw("hist_same")
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c4.Print("Delta_Eta_MET_Muon.png")
c4.cd()

c5=ROOT.TCanvas()
c5.Draw()
c5.cd()
hist_array_pt[2].Draw("hist_same")
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c5.Print("Leading_Muon_Pt.png")
c5.cd()
