import ROOT
import json
from ROOT import TMVA,TLatex,Math
import numpy as np
from math import pi, sin, cos, sqrt, sinh
import matplotlib.pyplot as plt
from array import array

import onnxruntime

# Load preprocess.json
with open("preprocess.json", "r") as f:
    preprocess_info = json.load(f)



file1 = ROOT.TFile("weirdtracks.root","READ")
dir = file1.Get("tree")
tree = dir.Get("Events")

file_temp = ROOT.TFile("output.root", 'recreate')
newtree = tree.CloneTree(0)

btagger_score = ROOT.std.vector('float')()
newtree.Branch("btagger_score", btagger_score)



onnx_model_path = "very_bad_btagger.onnx"
session = onnxruntime.InferenceSession(onnx_model_path)

label_name = session.get_outputs()[0].name

print(label_name)

for ientry in range(10): #tree.GetEntries()
    tree.GetEntry(ientry)
    for i in range(tree.j_no):
        jet_px = tree.j_pt[i] * cos(tree.j_phi[i])
        jet_py = tree.j_pt[i] * sin(tree.j_phi[i])
        jet_pz = tree.j_pt[i] * sinh(tree.j_eta[i])
        
        jet_dir_temp = ROOT.Math.XYZVector(jet_px, jet_py, jet_pz)
        jet_dir = jet_dir_temp.Unit()
        jet_dir3 = ROOT.TVector3(jet_px, jet_py, jet_pz)

        
        pfcand_pt_log_nopuppi_array = np.array([])
        pfcand_e_log_nopuppi_array = np.array([])
        pfcand_etarel_array = np.array([])
        pfcand_phirel_array = np.array([])
        pfcand_abseta_array = np.array([])
        pfcand_charge_array = np.array([],int)
        pfcand_lostInnerHits_array = np.array([])
        pfcand_normchi2_array = np.array([])
        pfcand_quality_array = np.array([])
        pfcand_dz_array = np.array([])
        pfcand_dzsig_array = np.array([])
        pfcand_dxy_array = np.array([])
        pfcand_dxysig_array = np.array([])
        pfcand_btagEtoRel_array = np.array([])
        pfcand_btagPtRatio_array = np.array([])
        pfcand_btagPParRatio_array = np.array([])
        
        #looping over all the pfcand
        
        pfcand_mask = np.array([],int)
        for j in range(len(tree.pfcand_index[i])): 
            pfcand_spot = int(tree.pfcand_index[i][j])
            
            pfcand_mask = np.append(pfcand_mask, 1)

            pfcand_pt_log_nopuppi_array = np.append(pfcand_pt_log_nopuppi_array,tree.pfcand_pt_log_nopuppi[pfcand_spot])
            pfcand_e_log_nopuppi_array = np.append(pfcand_e_log_nopuppi_array, tree.pfcand_e_log_nopuppi[pfcand_spot])
            pfcand_etarel_array = np.append(pfcand_etarel_array, (tree.pfcand_abseta[pfcand_spot]-tree.j_eta[i])*np.sign(tree.j_eta[i]))
            pfcand_phirel_array = np.append(pfcand_phirel_array, ((tree.pfcand_absphi[pfcand_spot]*np.sign(tree.j_eta[i]))-tree.j_phi[i]))
            pfcand_abseta_array = np.append(pfcand_abseta_array, tree.pfcand_abseta[pfcand_spot])
            pfcand_charge_array = np.append(pfcand_charge_array, tree.pfcand_charge[pfcand_spot])

            if tree.pfcand_normchi2[pfcand_spot] > 900:
                pfcand_lostInnerHits_array = np.append(pfcand_lostInnerHits_array, 0)
                pfcand_normchi2_array = np.append(pfcand_normchi2_array, 0)
                pfcand_quality_array = np.append(pfcand_quality_array, 0)
                pfcand_dz_array = np.append(pfcand_dz_array, 0)
                pfcand_dzsig_array = np.append(pfcand_dzsig_array, 0)
                pfcand_dxy_array = np.append(pfcand_dxy_array, 0)
                pfcand_dxysig_array = np.append(pfcand_dxysig_array, 0)
                pfcand_btagEtoRel_array = np.append(pfcand_btagEtoRel_array, 0)
                pfcand_btagPtRatio_array = np.append(pfcand_btagPtRatio_array, 0)
                pfcand_btagPParRatio_array = np.append(pfcand_btagPParRatio_array , 0)
                    

            else:            
                trk_eta = tree.pfcand_trk_eta[pfcand_spot]
                trk_phi = tree.pfcand_trk_phi[pfcand_spot]
                trk_pt = tree.pfcand_trk_pt[pfcand_spot]

                trk_px = trk_pt * cos(trk_phi)
                trk_py = trk_pt * sin(trk_phi)
                trk_pz = trk_pt * sinh(trk_eta)

                trk_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz)
                mom = ROOT.Math.XYZVector(trk_px,trk_py, trk_pz)
                    
                pfcand_lostInnerHits_array = np.append(pfcand_lostInnerHits_array, tree.pfcand_lostInnerHits[pfcand_spot])
                pfcand_normchi2_array = np.append(pfcand_normchi2_array, tree.pfcand_normchi2[pfcand_spot])
                pfcand_quality_array = np.append(pfcand_quality_array, tree.pfcand_quality[pfcand_spot])
                pfcand_dz_array = np.append(pfcand_dz_array, tree.pfcand_dz[pfcand_spot])
                pfcand_dzsig_array = np.append(pfcand_dzsig_array, tree.pfcand_dzsig[pfcand_spot])
                pfcand_dxy_array = np.append(pfcand_dxy_array, tree.pfcand_dxy[pfcand_spot])
                pfcand_dxysig_array = np.append(pfcand_dxysig_array, tree.pfcand_dxysig[pfcand_spot])
                pfcand_btagEtoRel_array = np.append(pfcand_btagEtoRel_array, jet_dir.Eta() - mom.Eta())


                if trk_mag == 0:
                    pfcand_btagPtRatio_array = np.append(pfcand_btagPtRatio_array, 0)
                    pfcand_btagPParRatio_array = np.append(pfcand_btagPParRatio_array , 0)
                else:     
                    mom3 = ROOT.TVector3(trk_px,trk_py,trk_pz)
                    pfcand_btagPtRatio_array = np.append(pfcand_btagPtRatio_array, mom3.Perp(jet_dir3)/trk_mag)
                    pfcand_btagPParRatio_array = np.append(pfcand_btagPParRatio_array , (trk_px*jet_px + trk_py*jet_py + trk_pz*jet_pz) / trk_mag)
                        
       
        pf_features = {
            'pfcand_pt_log_nopuppi': pfcand_pt_log_nopuppi_array,
            'pfcand_e_log_nopuppi': pfcand_e_log_nopuppi_array,
            'pfcand_etarel': pfcand_etarel_array,
            'pfcand_phirel': pfcand_phirel_array,
            'pfcand_abseta': pfcand_abseta_array,
            'pfcand_charge': pfcand_charge_array,
            'pfcand_lostInnerHits':  pfcand_lostInnerHits_array,
            'pfcand_normchi2': pfcand_normchi2_array,
            'pfcand_quality': pfcand_quality_array,
            'pfcand_dz': pfcand_dz_array,
            'pfcand_dzsig': pfcand_dzsig_array,
            'pfcand_dxy': pfcand_dxy_array,
            'pfcand_dxysig': pfcand_dxysig_array,
            'pfcand_btagEtoRel': pfcand_btagEtoRel_array,
            'pfcand_btagPtRatio': pfcand_btagPtRatio_array,
            'pfcand_btagPParRatio': pfcand_btagPParRatio_array
        }
        
        pf_points = {
            'pfcand_etarel': pfcand_etarel_array,
            'pfcand_phirel': pfcand_phirel_array
        }
        
        pf_mask = {
            'pfcand_mask': pfcand_mask
        }
        
        inputs = {
            'pf_points': pf_points,
            'pf_features': pf_features,
            'pf_mask': pf_mask
        }

        
        # Run inference
        output = session.run(["label_b","label_bb","label_c","label_cc","label_uds","label_g","label_undef"],inputs)

        
    
    
    
    
    
    btagger_score.push_back(ientry)
    
    newtree.Fill()
    
    temp_vector.clear()
    
file_temp.Write()
file_temp.Close()
