import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')

process.GlobalTag.globaltag = 'STARTHI53_V25::All'

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltZeroBiasSingleTrack = process.hltHighLevel.clone()
process.hltZeroBiasSingleTrack.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]

from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)
process.HeavyIonGlobalParameters = cms.PSet(
  centralityVariable = cms.string("HFtowersPlusTrunc"),
  nonDefaultGlauberModel = cms.string(""),
  centralitySrc = cms.InputTag("pACentrality")
  )

process.options = cms.untracked.PSet(
# after crash
SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 20

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/user/tuos/MC_HIJING_pPb/1EBA3EC5-D375-E211-B9BA-00266CFAE228.root'
    )
)

# histogram output
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("DS_Epos.root"))


process.demo = cms.EDAnalyzer('AnaEff'
)


process.p = cms.Path(process.demo)
