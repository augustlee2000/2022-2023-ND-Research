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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100)) #47168
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring("/store/mc/Run3Summer22MiniAODv4/TT_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/00c263d2-32c5-41f0-8f94-5cda0a135237.root"),
        secondaryFileNames = cms.untracked.vstring(
         	"/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12-v3/70000/050baf38-d1d8-43b5-93fc-6b206dcef3d3.root"
        )
)

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000)) #47168
# process.source = cms.Source("PoolSource",
# 	fileNames = cms.untracked.vstring(),
#         secondaryFileNames = cms.untracked.vstring()
# )


process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string("test.root")
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

# # Create SoftDrop pruned GEN jets
process.load('PhysicsTools.NanoAOD.jetMC_cff')
process.genJetSequence = cms.Sequence(
   process.patJetPartonsNano+
   process.genJetFlavourAssociation+
   process.ak8GenJetsWithNu+
   process.ak8GenJetsWithNuSoftDrop
)

# Create ParticleNet ntuple
process.tree = cms.EDAnalyzer("AK4JetNtupleProducer",
      isQCD = cms.untracked.bool( False ),
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
