#!/usr/bin/env python
# coding: utf-8

# In[15]:


import ROOT
from ROOT import TMVA
import numpy as np
from math import pi
import matplotlib.pyplot as plt
import cppyy.ll

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


# In[22]:


file1 = ROOT.TFile("X_bbbb_v2.root","READ")
tree1 = file1.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.85,.85)
hists_legend.Clear()
hists_legend.SetLineWidth(0)

colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack]
names = ["BBBB Mass"]  

hist_array = MakeHist(1,colors,names,"X-> AA ->BBBB: 4B Mass" ,"Mass (GeV)","Normalized Instances",4000,0,50)
legends(hist_array[0])


i = 0
j = 1
k = 2
l = 3

for ientry in range(tree.GetEntries()):
    tree.GetEntry(ientry)
    vec_1 = ROOT.TLorentzVector()
    vec_2 = ROOT.TLorentzVector()
    vec_3 = ROOT.TLorentzVector()
    vec_4 = ROOT.TLorentzVector()
    
    eta_array_intital = tree.Jet_eta
    phi_array_intital = tree.Jet_phi
    pt_array_intital = tree.Jet_pt
    mass_array_intital = tree.Jet_mass
    
    eta_array = np.frombuffer(eta_array_intital, dtype=np.float32)
    phi_array = np.frombuffer(phi_array_intital, dtype=np.float32)
    pt_array = np.frombuffer(pt_array_intital, dtype=np.float32)
    mass_array = np.frombuffer(mass_array_intital, dtype=np.float32)
    
    if eta_array.size > 3: 
        vec_1.SetPtEtaPhiM(pt_array[i],eta_array[i],phi_array[i],mass_array[i])
        vec_2.SetPtEtaPhiM(pt_array[j],eta_array[j],phi_array[j],mass_array[j])
        vec_3.SetPtEtaPhiM(pt_array[k],eta_array[k],phi_array[k],mass_array[k])
        vec_4.SetPtEtaPhiM(pt_array[l],eta_array[l],phi_array[l],mass_array[l])
    
        vec_final = vec_1 + vec_2 +vec_3 + vec_4
    
        hist_array[0].Fill(vec_final.M())
    
HistNorm(hist_array[0])

ROOT.gStyle.SetOptStat(0000000)
c1=ROOT.TCanvas("","",800,600)
c1.Draw()
c1.cd()
hist_array[0].Draw("same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c1.cd()    


# In[ ]:




