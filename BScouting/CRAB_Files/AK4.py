import FWCore.ParameterSet.Config as cms
#from PhysicsTools.NanoAOD.run3scouting_cff import *

from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register('inputDataset',
    '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Input dataset"
)

process = cms.Process("LL")
params.parseArguments()

#process.load("PhysicsTools.NanoAOD.run3scouting_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000)) #47168
# process.source = cms.Source("PoolSource",
# 	fileNames = cms.untracked.vstring("/store/mc/Run3Summer22MiniAODv3/TT_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/0f0274b6-b4cd-4ed3-8bcf-3744e0704bcc.root"),
#         secondaryFileNames = cms.untracked.vstring(
#          	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/00c0a31b-fd93-4bb0-87e8-02f6921a2607.root",
#          	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/36607336-57c8-4c1a-be5d-f4f978478ec8.root",
#          	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/9d61998c-e9dc-4952-9674-76171199ee85.root",
#          	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/d2419612-82cb-4035-9802-fb31ab4a372b.root",
#             "/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/e792db21-006d-4ec8-a8e9-e58584bfb357.root"
#         )
# )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000)) #47168
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring("/store/mc/Run3Summer22MiniAODv3/BulkGravitonToHH_MX600_MH22_TuneCP5_13p6TeV_madgraph-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12-v4/30000/0ea2defc-f6da-4fb1-95a8-280555a62b40.root"),
        secondaryFileNames = cms.untracked.vstring(
         	"/store/mc/Run3Summer22DRPremix/BulkGravitonToHH_MX600_MH22_TuneCP5_13p6TeV_madgraph-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12-v4/30000/07fa7492-e0a3-4087-b0a5-3d46b8344ecb.root",
         	"/store/mc/Run3Summer22DRPremix/BulkGravitonToHH_MX600_MH22_TuneCP5_13p6TeV_madgraph-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12-v4/30000/412f0428-d625-4519-b5cd-a91ac401f29f.root"

             
        )
)


process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test.root")
)


# # # Create SoftDrop pruned GEN jets
# process.load('PhysicsTools.NanoAOD.jetMC_cff')
# process.genJetSequence = cms.Sequence(
#    process.patJetPartonsNano+
#    process.genJetFlavourAssociation,
#    process.ak8GenJetsWithNu+
#    process.ak8GenJetsWithNuSoftDrop
# )

# Create SoftDrop pruned GEN jets
from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets
process.ak8GenJetsWithNu = ak8GenJets.clone(
    src='packedGenParticles',
    rParam=cms.double(0.8),
    jetPtMin=100.0
)

process.ak8GenJetsWithNuSoftDrop = process.ak8GenJetsWithNu.clone(
    useSoftDrop=cms.bool(True),
    zcut=cms.double(0.1),
    beta=cms.double(0.0),
    R0=cms.double(0.8),
    useExplicitGhosts=cms.bool(True)
)

# # Create SoftDrop pruned GEN jets
process.load('PhysicsTools.NanoAOD.jetMC_cff')
process.genJetSequence = cms.Sequence(
   process.patJetPartonsNano+
   process.genJetFlavourAssociation+
   process.ak8GenJetsWithNu+
   process.ak8GenJetsWithNuSoftDrop
)

# process.genJetSequence = cms.Sequence(
#     process.ak8GenJetsWithNu+
#     process.ak8GenJetsWithNuSoftDrop
# )

# Create ParticleNet ntuple
process.tree = cms.EDAnalyzer("AK4JetNtupleProducer",
      isQCD = cms.untracked.bool( '/QCD_' in params.inputDataset ),
      gen_jets = cms.InputTag( "genJetFlavourAssociation" ),
      pf_candidates = cms.InputTag( "hltScoutingPFPacker" ),
      fgen_jets = cms.InputTag( "ak8GenJetsWithNuSoftDrop" ),
      gen_candidates = cms.InputTag( "prunedGenParticles" ),
      muons = cms.InputTag( "hltScoutingMuonPacker" ),
      pfMet            = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
      pfMetPhi         = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
      gen_jet_data     = cms.InputTag("ak4GenJetsNoNu"),

)

process.p = cms.Path(process.genJetSequence*process.tree)