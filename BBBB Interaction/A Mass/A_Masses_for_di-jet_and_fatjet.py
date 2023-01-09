#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# In[17]:


file1 = ROOT.TFile("XtoAAto4b_X3000_A50.root","READ")
tree = file1.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.90,.89)
hists_legend.Clear()
hists_legend.SetLineWidth(0)

colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["A Mass from Fat Jets" ]  

hist_array = MakeHist(1,colors,names,"XtoAAto4b_X3000_A50: A Mass" ,"Mass (GeV)","Normalized Instances",100,0,75)
legends(hist_array[0])

i = 0
j = 1
k = 2
l = 3


for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    
    if tree.nFatJet > 1:
        
        btag_array_i = tree.FatJet_btagHbb
        btag_array = np.frombuffer(btag_array_i, dtype=np.float32)
        
        btag_sorted = np.where(btag_array > 0.7527) #tight btagging on Deep CSV 0.7527
        btag_sorted = btag_sorted[0]
        
        if btag_sorted.size > 1:
            
            FatJet_msoftdrop = np.frombuffer(tree.FatJet_msoftdrop, dtype=np.float32)
            
            hist_array[0].Fill(FatJet_msoftdrop[btag_sorted[0]])

            
            
HistNorm(hist_array[0])

ROOT.gStyle.SetOptStat(0000000)
c1=ROOT.TCanvas("","",800,600)
c1.Draw()
c1.cd()
hist_array[0].Draw("hist_same")
hists_legend.Draw()

AddCMSLumi(c1, 138, "Internal")
c1.cd()    


# In[33]:


file1 = ROOT.TFile("XtoAAto4b_X800_A50.root","READ")
tree = file1.Get("Events")


h1 = ROOT.TH2F("Di-Jet mass vs mass asymmetry", "XtoAAto4b_X1000_A25: Di-Jet mass vs mass asymmetry ", 25, 0, 100, 50, 0, 1)

h1.GetXaxis().SetTitle("di-jet mass (Gev) ")
h1.GetYaxis().SetTitle("mass asymmetry")


i = 0
j = 3

for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    
    if tree.nJet > 3:
        
        btag_array_i = tree.Jet_btagCSVV2
        btag_array = np.frombuffer(btag_array_i, dtype=np.float32)
        
        btag_sorted = np.where(btag_array > 0.7527)#tight btagging on Deep CSV
        btag_sorted = btag_sorted[0]
        
        if btag_sorted.size > 3:
            
            vec_1 = ROOT.TLorentzVector()
            vec_2 = ROOT.TLorentzVector()
            
            eta_array_intital = tree.Jet_eta
            phi_array_intital = tree.Jet_phi
            pt_array_intital = tree.Jet_pt
            mass_array_intital = tree.Jet_mass
    
            eta_array = np.frombuffer(eta_array_intital, dtype=np.float32)
            phi_array = np.frombuffer(phi_array_intital, dtype=np.float32)
            pt_array = np.frombuffer(pt_array_intital, dtype=np.float32)
            mass_array = np.frombuffer(mass_array_intital, dtype=np.float32)
            
            vec_1.SetPtEtaPhiM(pt_array[btag_sorted[i]],eta_array[btag_sorted[i]],phi_array[btag_sorted[i]],mass_array[btag_sorted[i]])
            vec_2.SetPtEtaPhiM(pt_array[btag_sorted[j]],eta_array[btag_sorted[j]],phi_array[btag_sorted[j]],mass_array[btag_sorted[j]])
            
            di_jet_mass = vec_1 + vec_2
            
            mass_asymmetry = abs(mass_array[btag_sorted[i]] - mass_array[btag_sorted[j]])/(mass_array[btag_sorted[i]] + mass_array[btag_sorted[j]])
            
            #print(di_jet_mass.M())
            
            h1.Fill(di_jet_mass.M(), mass_asymmetry)
            
            
            
            
            
c1=ROOT.TCanvas()
c1.Draw()
c1.cd()
h1.Draw("COLZ")
c1.cd()


# In[ ]:




