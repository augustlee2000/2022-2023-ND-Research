import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.nanogen_cff import *
# Set parameters externally
from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
    'isMC',
    True,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether the sample is simulation or data'
)

params.register(
    'useWeights',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to use the events weights from a Monte Carlo generator'
)

params.register(
    'filterTrigger',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to fire a trigger used in the analysis'
)

params.register(
    'filterMuons',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to ask the event to contain at least two muons'
)

params.register(
    'reducedInfo',
    False,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Flag to indicate whether or not to store just the reduced information'
)

params.register(
    'trigProcess',
    'HLT',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagData',
    '<not available>',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'GlobalTagMC',
    '112X_mcRun3_2021_realistic_v16',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'Process name for the HLT paths'
)

params.register(
    'xsec',
    0.001,
    VarParsing.multiplicity.singleton,VarParsing.varType.float,
    'Cross-section for a Monte Carlo Sample'
)

# Define the process
process = cms.Process("LL")

# Parse command line arguments
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 5

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary      = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound')
)


# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(params.maxEvents) ) #params.maxEvents

# Input EDM files
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring([
	'/store/mc/Run3Summer22DRPremix/TT_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/124X_mcRun3_2022_realistic_v12_ext1-v3/2810000/00c0a31b-fd93-4bb0-87e8-02f6921a2607.root'
	])
)

# Dataset: /TT_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22DRPremix-124X_mcRun3_2022_realistic_v12_ext1-v3/AODSIM

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')

# Load the global tag
from Configuration.AlCa.GlobalTag import GlobalTag
if params.isMC :
    process.GlobalTag.globaltag = params.GlobalTagMC
else :
    process.GlobalTag.globaltag = params.GlobalTagData

# Define the services needed for the treemaker
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(params.outputFile)
)

# Tree for the generator weights
'''
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheInfo = cms.InputTag("externalLHEProducer"),
    genInfo = cms.InputTag("generator"),
    useLHEWeights = cms.bool(params.useWeights)
)
'''

#from DarkPhotonAnalysis.DimuonAnalysis2018.TriggerPaths_cfi import getL1Conf
#L1Info = ['L1_DoubleMu_12_5', 'L1_DoubleMu_15_7', 'L1_HTT200er', 'L1_HTT255er', 'L1_HTT280er', 'L1_HTT320er', 'L1_HTT360er', 'L1_ETT2000', 'L1_HTT400er', 'L1_HTT450er', 'L1_SingleJet180', 'L1_SingleJet200', 'L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5', 'L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5', 'L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7', 'L1_DoubleMu4_SQ_OS_dR_Max1p2', 'L1_SingleEG36er2p5', 'L1_SingleLooseIsoEG28er2p1', 'L1_SingleEG8er2p5', 'L1_SingleEG10er2p5', 'L1_SingleEG15er2p5', 'L1_SingleEG26er2p5', 'L1_SingleEG28_FWD2p5', 'L1_DoubleEG4_er1p2_dR_Max0p9', 'L1_DoubleEG4p5_er1p2_dR_Max0p9', 'L1_DoubleEG5_er1p2_dR_Max0p9', 'L1_DoubleEG5p5_er1p2_dR_Max0p8', 'L1_DoubleEG7_er1p2_dR_Max0p8', 'L1_DoubleEG7p5_er1p2_dR_Max0p7', 'L1_DoubleEG_15_10_er2p5', 'L1_DoubleEG_20_10_er2p5', 'L1_DoubleEG_22_10_er2p5', 'L1_DoubleEG_25_12_er2p5', 'L1_DoubleEG_25_14_er2p5', 'L1_DoubleEG_27_14_er2p5', 'L1_DoubleEG_LooseIso22_12_er2p5', 'L1_DoubleEG_LooseIso25_12_er2p5', 'L1_TripleEG_18_17_8_er2p5', 'L1_TripleEG_18_18_12_er2p5', 'L1_TripleEG16er2p5', 'L1_DoubleEG8er2p5_HTT300er']
L1Info = ['L1_DoubleMu_12_5','L1_DoubleMu_15_7','L1_HTT200er','L1_HTT255er','L1_HTT280er','L1_HTT320er','L1_HTT360er','L1_ETT2000','L1_HTT400er','L1_HTT450er','L1_SingleJet180','L1_SingleJet200','L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5','L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5','L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5','L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7','L1_DoubleMu4_SQ_OS_dR_Max1p2','L1_SingleEG36er2p5','L1_SingleLooseIsoEG28er2p1']

# process.genParticlesMerged = cms.EDProducer("MergedGenParticleProducer",   I'm not sure why or how I added this line 
#     inputPruned = cms.InputTag("prunedGenParticles"),
#     inputPacked = cms.InputTag("packedGenParticles")
# )

# Make tree
process.mmtree = cms.EDAnalyzer('ScoutingNanoAOD',

    	triggerresults   = cms.InputTag("TriggerResults", "", params.trigProcess),
        doL1 = cms.bool(True),
        triggerConfiguration = cms.PSet(
    		hltResults            = cms.InputTag('TriggerResults','', params.trigProcess),
    		l1tResults            = cms.InputTag(''),
    		daqPartitions         = cms.uint32(1),
    		l1tIgnoreMaskAndPrescale = cms.bool(False),
    		throw                 = cms.bool(False)
  	),
	ReadPrescalesFromFile = cms.bool( False ),
        AlgInputTag       = cms.InputTag("gtStage2Digis"),
        l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
        l1Seeds           = cms.vstring(L1Info),
	#vertices         = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
	muons            = cms.InputTag("muons"),
	electrons        = cms.InputTag("gedGsfElectrons"),
        photons          = cms.InputTag("gedPhotons"),
	pfcands          = cms.InputTag("particleFlow"),
	pfjets           = cms.InputTag("ak4PFJets"),
        tracks           = cms.InputTag("generalTracks"),
        primaryVertices  = cms.InputTag("offlinePrimaryVertices"),
        displacedVertices  = cms.InputTag("inclusiveSecondaryVertices"),
        pfMet            = cms.InputTag("hltScoutingPFPacker","pfMetPt"),
        pfMetPhi         = cms.InputTag("hltScoutingPFPacker","pfMetPhi"),
        rho         = cms.InputTag("hltScoutingPFPacker","rho"),    
    GenParticle         = cms.InputTag("genParticles"),
    GenJetAK4         = cms.InputTag("ak4GenJets"),
    GenJetAK4NoNu         = cms.InputTag("ak4GenJetsNoNu"),
    GenJetAK8         = cms.InputTag("ak8GenJets"),
    	#pileupinfo       = cms.InputTag("addPileupInfo"),
    	#geneventinfo     = cms.InputTag("generator"),

)
#process  = customizeNanoGEN(process) #this seems to enter the code but never wrote a new branch so I don't know what that is all about


process.p = cms.Path(process.gtStage2Digis*process.mmtree)
