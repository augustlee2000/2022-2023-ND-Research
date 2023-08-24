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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10)) #47168
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring("/store/mc/Run3Summer22MiniAODv3/TT_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/0f0274b6-b4cd-4ed3-8bcf-3744e0704bcc.root"),
        secondaryFileNames = cms.untracked.vstring(
         	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/00c0a31b-fd93-4bb0-87e8-02f6921a2607.root",
         	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/36607336-57c8-4c1a-be5d-f4f978478ec8.root",
         	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/9d61998c-e9dc-4952-9674-76171199ee85.root",
         	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/d2419612-82cb-4035-9802-fb31ab4a372b.root",
            "/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/e792db21-006d-4ec8-a8e9-e58584bfb357.root"
        )
)



# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))
# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring("/store/mc/Run3Summer22EEMiniAODv3/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30000/b5d16dec-30fc-4a1e-a287-fd2b57c6c2bc.root"),
#         secondaryFileNames = cms.untracked.vstring(
#         "/store/mc/Run3Summer22EEDRPremix/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30003/4c74832b-a3f0-48fb-9a76-e871fed3c151.root",
#         "/store/mc/Run3Summer22EEDRPremix/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30003/62d3132b-23bb-4001-9e67-27641ec3d098.root",
#         "/store/mc/Run3Summer22EEDRPremix/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30003/bd759eb2-5a6f-4582-ab32-fe5463891960.root",
#         "/store/mc/Run3Summer22EEDRPremix/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30003/c1bcdff5-73a5-460a-931c-b51fca6fc49b.root",
#         "/store/mc/Run3Summer22EEDRPremix/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_postEE_v1-v2/30003/f45d6172-aab1-45f4-b07a-55a561225b3b.root"
#         )
# )



process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test.root")
)


# Create SoftDrop pruned GEN jets
process.load('PhysicsTools.NanoAOD.jetMC_cff')
process.genJetSequence = cms.Sequence(
   process.patJetPartonsNano+
   process.genJetFlavourAssociation
)

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

process.genJetSequence = cms.Sequence(
    process.ak8GenJetsWithNu+
    process.ak8GenJetsWithNuSoftDrop
)

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

)

process.p = cms.Path(process.genJetSequence*process.tree)
