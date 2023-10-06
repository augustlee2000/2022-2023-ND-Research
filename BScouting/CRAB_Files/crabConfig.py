import CRABClient
from CRABClient.UserUtilities import config
from WMCore.Configuration import Configuration

config = config()

config.General.requestName = 'bscouting_ntuple_analysis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Run3ScoutingJetTagging/Analysis/test/AK4.py'

config.Data.inputDataset = '/BulkGravitonToHH_MX700_MH91_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v3/MINIAODSIM' #input data set
config.Data.useParent = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
config.Data.outputDatasetTag = 'bscouting_ntuple_analysis'

config.Site.storageSite = 'T3_US_NotreDame'


