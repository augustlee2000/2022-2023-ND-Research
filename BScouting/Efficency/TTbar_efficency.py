import numpy as np
import pandas as pd
import uproot
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, confusion_matrix
from math import pi, sin, cos, sqrt, sinh, log
import glob
import hist
from itertools import combinations, combinations_with_replacement, permutations
import ROOT




def chi2_function(top_l, top_h):
    chi2 = (172.76 - top_l)**2 + (172.76 - top_h)**2
    return chi2/172.76
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

def legends(hists):
    print("name = ",hists.GetName())
    hists_legend.AddEntry(hists,hists.GetName(),"lp") 


def neutrino(jet,muon, pt_v, phi_v):
    pt_l = muon.Pt()
    phi_l = muon.Phi()
    pz_l = muon.Pz()
    E_l = muon.E()
    mu = mu_calc(pt_l, pt_v, phi_l, phi_v)
    square_root = (((mu**2 * pz_l**2)/(pt_l**4)) - (((E_l**2 * pt_v**2) - mu**2) /(pt_l**2)))
    real_part = ((mu*pz_l)/(pt_l**2))
    if square_root < 0:
        pz_v = real_part
        neutrino = neutrino_vec(pt_v*cos(phi_v), pt_v*sin(phi_v), pz_v)

        top = jet + muon + neutrino

        chi2 = (top.M() - 172.5)**2 / (25)

        return chi2
        
    else:
        chi2 = []
        pz_v = [real_part + sqrt(square_root), real_part - sqrt(square_root) ]

        neutrino_0 = neutrino_vec(pt_v*cos(phi_v), pt_v*sin(phi_v), pz_v[0])
        neutrino_1 = neutrino_vec(pt_v*cos(phi_v), pt_v*sin(phi_v), pz_v[1])

        top_0 = jet + muon + neutrino_0
        top_1 = jet + muon + neutrino_1

        chi2.append((top_0.M() - 172.5)**2 / (25))
        chi2.append((top_1.M() - 172.5)**2 / (25))

        return min(chi2)

                  
def mu_calc(pt_l, pt_v, phi_l, phi_v):
    return ( 3230.23 + (pt_l * pt_v* cos(phi_v - phi_l)))

def neutrino_vec(px,py,pz):
    v_E = sqrt(px**2 + py**2 + pz**2)
    v = ROOT.TLorentzVector()
    v.SetPxPyPzE(px,py,pz,v_E)
    return v

def btagging_efficenty_v2(file):  #, hist_0, hist_1, hist_2):
    infile_name = file #takes in the current file
    infile = uproot.open(infile_name)
    tree = infile['Events'] #Get the event tree

    #Get the ML scores if it greater than .3 then it is 1 if it less then is 0
    b_jets = tree["ScoutingJet_particleNet_prob_b"].array()
    b_jets_btagged  = np.where(b_jets > .5 ,1, 0)


    #applying muon and jet cuts on the data
    muons = tree["nScoutingMuon"].array()
    B_jet_mask = np.sum(b_jets_btagged == 1, axis=1)
    b_jets_cut = B_jet_mask > -1
    number_of_events_cut = np.array([len(jets) for jets in b_jets]) == 4
    muons_cut = muons == 1

    jet_pt = tree["ScoutingJet_pt"].array()
    pt_cut = jet_pt >= 30
    filtered_jet_pt = jet_pt[pt_cut]
    num_jets_after_cut = ak.num(filtered_jet_pt)
    number_of_events_cut = num_jets_after_cut == 4
    

    met_pt_cut = tree["ScoutingMET_pt"].array() > 30

    event_mask = b_jets_cut & number_of_events_cut & muons_cut & met_pt_cut

    #making a cut on individual muons

    #making the cut on all the branches that I need
    cut_b_jets = b_jets_btagged[event_mask]
    cut_phi = tree["ScoutingJet_phi"].array()[event_mask]
    cut_eta = tree["ScoutingJet_eta"].array()[event_mask]
    cut_mass = tree["ScoutingJet_mass"].array()[event_mask]
    cut_pt = tree["ScoutingJet_pt"].array()[event_mask]

    n_muon = tree["nScoutingMuon"].array()[event_mask]
    muon_pt = tree["ScoutingMuon_pt"].array()[event_mask]
    muon_eta = tree["ScoutingMuon_eta"].array()[event_mask]
    muon_phi = tree["ScoutingMuon_phi"].array()[event_mask]
    muon_mass = tree["ScoutingMuon_m"].array()[event_mask]
    met_pt = tree["ScoutingMET_pt"].array()[event_mask]
    met_phi = tree["ScoutingMET_phi"].array()[event_mask]
    muon_ecal_iso = tree["ScoutingMuon_ecalIso"].array()[event_mask]
    muon_hcal_iso =tree["ScoutingMuon_hcalIso"].array()[event_mask]
    muon_track_iso = tree["ScoutingMuon_trackIso"].array()[event_mask]

    muon_iso = (muon_ecal_iso + muon_hcal_iso + muon_track_iso)/muon_pt

    sigma2 = np.array([[25, 1], [1, 25]])
    idx_combo = [[0,1,2,3],[0,2,1,3],[0,3,1,2], [1,0,2,3], [1,2,0,3], [1,3,0,2], [2,0,1,3], [2,1,0,3], [2,3,0,1], [3,0,1,2], [3,1,0,2], [3,2,0,1]]
    tag = 0
    probe = 0
    veto = 0
    #loop over all events that are left
    for i in range(100): #len(cut_b_jets)
        if muon_pt[i][0] > 30 and muon_iso[i] < .3: #and ((muon_pt[i][0]/ cut_pt[i][0]) < .98 and (muon_pt[i][0]/ cut_pt[i][0]) > 1.02) and ((muon_pt[i][0]/ cut_pt[i][1]) < .98 and (muon_pt[i][0]/ cut_pt[i][1]) > 1.02) and ((muon_pt[i][0]/ cut_pt[i][2]) < .98 and (muon_pt[i][0]/ cut_pt[i][2]) > 1.02) and ((muon_pt[i][0]/ cut_pt[i][3]) < .98 and (muon_pt[i][0]/ cut_pt[i][3]) > 1.02):
            chi2 = []
            
            jet_0 = ROOT.TLorentzVector()
            jet_1 = ROOT.TLorentzVector()
            jet_2 = ROOT.TLorentzVector()
            jet_3 = ROOT.TLorentzVector()
            muon = ROOT.TLorentzVector()
            jet_0.SetPtEtaPhiM(cut_pt[i][0], cut_eta[i][0], cut_phi[i][0], cut_mass[i][0])
            jet_1.SetPtEtaPhiM(cut_pt[i][1], cut_eta[i][1], cut_phi[i][1], cut_mass[i][1])
            jet_2.SetPtEtaPhiM(cut_pt[i][2], cut_eta[i][2], cut_phi[i][2], cut_mass[i][2])
            jet_3.SetPtEtaPhiM(cut_pt[i][3], cut_eta[i][3], cut_phi[i][3], cut_mass[i][3])
            muon.SetPtEtaPhiM(muon_pt[i][0], muon_eta[i][0], muon_phi[i][0], muon_mass[i][0])
            jet_array = [jet_0, jet_1, jet_2, jet_3]
            lep_top_0 = neutrino(jet_0, muon, met_pt[i], met_phi[i])
            lep_top_1 = neutrino(jet_1, muon, met_pt[i], met_phi[i])
            lep_top_2 = neutrino(jet_2, muon, met_pt[i], met_phi[i])
            lep_top_3 = neutrino(jet_3, muon, met_pt[i], met_phi[i])
            lep_top_chi2 = [lep_top_0, lep_top_1, lep_top_2, lep_top_3]

            for j in range(len(idx_combo)): # 0: leptonic b, 1: hadronic b, 2: w, 3: w index
                w_quark = jet_array[idx_combo[j][2]] + jet_array[idx_combo[j][3]]
                w_mass = w_quark.M()
                top_quark = w_quark + jet_array[idx_combo[j][1]]
                top_mass = top_quark.M()

                had_top_chi2 = (top_mass - 172.5)**2/172.5 + (w_mass - 80.385)**2/80.385

                chi2_final = -1*log(had_top_chi2) - log(lep_top_chi2[idx_combo[j][0]]) 

                chi2.append(chi2_final * -1)
            min_idx = np.argmin(chi2)
            min_jet_idx = idx_combo[min_idx]
            if cut_b_jets[i][min_jet_idx[1]] == 1:
                tag +=1
                if cut_b_jets[i][min_jet_idx[0]] == 1:
                    probe +=1
                else:
                    veto +=1
    print(probe)
    print(probe / (tag))

btagging_efficenty_v2('/hadoop/store/user/aulee/TT_TuneCP5_13p6TeV_powheg-pythia8/btagging_V1_TTbar_1/240430_194500/0000/step1_1.root')




