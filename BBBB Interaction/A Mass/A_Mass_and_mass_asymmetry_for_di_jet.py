#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
from ROOT import TMVA,TLatex
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import cppyy.ll

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
    

def PhiAngle(phi):
    while phi>pi:
        phi = phi-pi
    return phi

def b_tagging(jet_1,jet_2,BEff_77):
    if BEff_77[jet_1] =="\x00" or BEff_77[jet_2] =="\x00":
        return 0
    else:
        return 1

def HistNorm(hist):
    integral = hist.Integral()
    hist.Scale(1./integral)
    

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

def AddCMSLumi(pad, fb, extra):
    cmsText     = "CMS " + extra
    cmsTextFont   = 61  
    lumiTextSize     = 0.45
    lumiTextOffset   = 0.15
    cmsTextSize      = 0.5
    cmsTextOffset    = 0.15
    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    e = 0.025
    pad.cd()
    lumiText = str(fb)+" fb^{-1} (13 TeV)"
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)	
    extraTextSize = 0.76*cmsTextSize
    latex.SetTextFont(42)
    latex.SetTextAlign(31) 
    latex.SetTextSize(lumiTextSize*t)
    latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)
    pad.cd()
    latex.SetTextFont(cmsTextFont)
    latex.SetTextSize(cmsTextSize*t)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.1265, 0.825, cmsText)
    pad.Update()


def vec_creator(phi,eta,mass,pt):
    vec_temp = ROOT.TLorentzVector()
    vec_temp.SetPtEtaPhiM(pt,eta,phi,mass)
    
    return vec_temp
    

file1 = ROOT.TFile("XtoAAto4b_X3000_A50.root","READ")
tree = file1.Get("Events")


h1 = ROOT.TH2F("Di-Jet mass vs mass asymmetry", "XtoAAto4b_X3000_A50: Di-Jet mass vs mass asymmetry ", 200, 0, 800, 50, 0, 1)

h1.GetXaxis().SetTitle("di-jet mass (Gev) ")
h1.GetYaxis().SetTitle("mass asymmetry")

for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    
    if tree.nJet > 3:
        
        btag_array_i = tree.Jet_btagCSVV2
        btag_array = np.frombuffer(btag_array_i, dtype=np.float32)
        
        btag_sorted = np.where(btag_array > 0.7527)#tight btagging on Deep CSV
        btag_sorted = btag_sorted[0]
        
        if btag_sorted.size > 3:
            
            eta_array_intital = tree.Jet_eta
            phi_array_intital = tree.Jet_phi
            pt_array_intital = tree.Jet_pt
            mass_array_intital = tree.Jet_mass
    
            eta_array = np.frombuffer(eta_array_intital, dtype=np.float32)
            phi_array = np.frombuffer(phi_array_intital, dtype=np.float32)
            pt_array = np.frombuffer(pt_array_intital, dtype=np.float32)
            mass_array = np.frombuffer(mass_array_intital, dtype=np.float32)
            
            vector_array = np.array([])
            
            for q in range(btag_sorted.size):
                vector  = vec_creator(phi_array[q],eta_array[q],mass_array[q],pt_array[q])
                
                vector_array = np.append(vector_array, vector)
                
            pair_array = np.array([])
            M = np.size(vector_array)
            for i in range(M):
                j=i+1
                while j <M:
                    pair_array = np.append(pair_array,[i,j])
                    j+=1
            rows = int(len(pair_array)/2)    
            pair_array = np.resize(pair_array, (rows,2))
            
            for x in range(rows):
                for y in range(rows):
                    a = int(pair_array[x,0])
                    b = int(pair_array[x,1])
                    c = int(pair_array[y,0])
                    d = int(pair_array[y,1])
                    if ((a == c) or (a == d) or (b == c) or (b == d) or (c < a) ):
                        a=a
                    else:
                        vec_final_1 = vector_array[a] + vector_array[b]
                        vec_final_2 = vector_array[c] + vector_array[d]
                        
                        vec_final_1_mass = vec_final_1.M()
                        vec_final_2_mass = vec_final_2.M()
                        
                        mass_asymmetry =  abs((vec_final_1_mass - vec_final_2_mass)/(vec_final_1_mass + vec_final_2_mass))
                        
                        h1.Fill(vec_final_1_mass, mass_asymmetry)
                        h1.Fill(vec_final_2_mass, mass_asymmetry)
                        
                        
                        
ROOT.gStyle.SetOptStat(0000000)            
c1=ROOT.TCanvas()
c1.Draw()
c1.cd()
h1.Draw("COLZ")
c1.cd()


# In[ ]:




