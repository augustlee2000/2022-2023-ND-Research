# -*- coding: utf-8 -*-
# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt, cos
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
    ind.sort(key=lambda x:x[0])
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


files = sys.argv[1]
file1 = ROOT.TFile(files,"READ")
dir = file1.Get("mmtree")
tree = dir.Get("Events")

hists_legend = ROOT.TLegend(.65,.70,.85,.80)
hists_legend.Clear()
hists_legend.SetLineWidth(0)

colors=[ROOT.kBlue,ROOT.kRed, ROOT.kBlack, ROOT.kGreen]
names = ["True Positive", "True Negative", "False Positive","False Negative" ]
hist_array_angle = MakeHist(4,colors,names,"Angle between Parton and Jet" ,"Angle (rad)","Normalized Instances",1,0,25)



CSV= np.arange(0,1.01,.01)

True_positive_rate_array  = np.array([])
False_positive_rate_array = np.array([])

for i in range(len(CSV)):
    false_positive = 0
    false_negative = 0
    true_positive = 0
    true_negative = 0

    for ientry in range(tree.GetEntries()):
        tree.GetEntry(ientry)
        MyParticles = []
        for gp in range(tree.n_GenPart):
            #if math.fabs(tree.GenPart_pdgId[gp]) in (1,2,3,4,5):# and tree.GenPart_status[gp] in (62,52):
            if tree.GenPart_StatusFlags[gp] ==23 and math.fabs(tree.GenPart_pdgId[gp]) == 5 and math.fabs(tree.GenPart_genPartIdxMother[gp]) == 6:
                MyParticles.append(vec_creator(tree.GenPart_pt[gp], tree.GenPart_eta[gp], tree.GenPart_phi[gp], tree.GenPart_m[gp]))

        MyJets = [] #jets
        for j in range(tree.n_jet):
            MyJets.append(vec_creator(tree.Jet_pt[j],tree.Jet_eta[j],tree.Jet_phi[j],tree.Jet_m[j]))
        

        mp = len(MyParticles)
        mj = len(MyJets)
        if mp  > 0:
            ind_s = angle_sort(mp,mj,MyParticles,MyJets)
            for l in range(mp):
                for j in range(tree.n_jet):
                    j += tree.n_jet * l
                    if j == tree.n_jet * l:
                        if tree.Jet_csv[ind_s[j][1]] > CSV[i]:
                            true_positive +=1
                            hist_array_angle[0].Fill(ind_s[j][2])
                        else:
                            false_negative +=1
                            hist_array_angle[3].Fill(ind_s[j][2])
                    else:
                        if tree.Jet_csv[ind_s[j][1]] > CSV[i]:
                            false_positive +=1
                            hist_array_angle[2].Fill(ind_s[j][2])
                        else:
                            true_negative +=1
                            hist_array_angle[1].Fill(ind_s[j][2])
        
    True_positive_rate = true_positive.sum() / (true_positive.sum() + false_negative.sum())
    False_positive_rate = false_positive.sum() / (false_positive.sum() + true_negative.sum())

    True_positive_rate_array = np.append(True_positive_rate_array,True_positive_rate)
    False_positive_rate_array = np.append(False_positive_rate_array, False_positive_rate)

False_positive_rate_array = 1-False_positive_rate_array
N=100

roc_integral = trapizoid_rule(False_positive_rate_array,True_positive_rate_array, N )

print("The value of the ROC Integral is:  ",roc_integral)

plt.scatter(False_positive_rate_array,True_positive_rate_array)  # I need to see if I can actually make a matplotlib graph on a ssh
plt.title("ROC Curve for Signal vs vbfH")
plt.xlabel("1-FPR")
plt.ylabel("TPR")
plt.show()


HistNorm(hist_array_angle)
legends(hist_array_angle)

c4=ROOT.TCanvas()
c4.Draw()
c4.cd()
hist_array_angle[0].Draw("hist_same")
hist_array_angle[1].Draw("hist_same")
hist_array_angle[2].Draw("hist_same")
hist_array_angle[3].Draw("hist_same")
hists_legend.Draw()
CMSLabel(.65,.85,ROOT.kBlack,"Internal")
c4.Print("Angle_multiHist_tpr.png")
c4.cd()
