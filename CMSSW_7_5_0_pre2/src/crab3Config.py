from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MultiStrangePYTHIA_step1'
config.General.workArea = 'crab_projects_nofilter'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'PYTHIA8_TuneCUETP8M1_5TeV_XiFilter_cff_GEN_SIM.py'

config.Data.outputPrimaryDataset = 'MinBias'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 20000
NJOBS = 20  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/honi/RpPbInput/MultiStrangeFilteredPYTHIA'
config.Data.publication = True
config.Data.outputDatasetTag = 'MultistrangeFilteredPYTHIA_GENSIM_755patch3'

config.Site.storageSite = 'T2_US_Vanderbilt'


#crab submit -c crab3Config.py
#crab status -d crab_projects/xxxx
