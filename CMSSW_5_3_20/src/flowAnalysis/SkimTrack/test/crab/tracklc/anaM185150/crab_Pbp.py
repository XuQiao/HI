from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'FlowLYZHMPbp'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'skimtrack_185150_cfi.py'
config.JobType.outputFiles = ['skimTreeTrack.root']

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/PAHighPt/davidlw-PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_ReTracking_v18-28b2b9cce04ec3f20baeb96fbd2295a8/USER'
#config.Data.lumiMask = ''
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 50
#config.Data.runRange = '193093-193999'
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.publication = False
#config.Data.publishDataName = ''
config.Data.outLFNDirBase = '/store/user/qixu/flow/PACorrHM/skim/tracklc/multiM185150/FlowLYZHMPbp'

config.section_('Site')
config.Site.storageSite = 'T2_US_Vanderbilt'
