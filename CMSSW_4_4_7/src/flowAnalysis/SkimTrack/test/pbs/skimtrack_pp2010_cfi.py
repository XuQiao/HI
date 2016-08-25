import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'GR_E_V33A::All'
#process.GlobalTag.globaltag = 'GR_P_V41::All'
process.GlobalTag.globaltag = 'GR_R_44_V15::All'

process.options = cms.untracked.PSet(
# after crash
SkipEvent = cms.untracked.vstring('ProductNotFound')	
)

# Input Source
process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(
#'/store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_v12/331639eac1650e369e7dacf443369878/pPb_HM_100_1_tMI.root'
#'file:/scratch/tuos/tmp/ppGoodColl_661_1_Miw.root'
),
                            #duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )
#mylist = FileUtils.loadListFromFile ('filelistAll.dat')
#for fname in mylist:
# process.source.fileNames.append('file:%s' % (fname))
#process.source.fileNames.append('%s' % (fname))
filename = "filelistAll.dat"
mylist = FileUtils.loadListFromFile (filename)
for fname in mylist[int(sys.argv[3]):int(sys.argv[4])]:
 process.source.fileNames.append('file:%s' % (fname))

# histogram output
filenameout="/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PP2010Data/pbs/pPbskimAnaTree_"+sys.argv[2]+".root"
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string(filenameout))

process.demo = cms.EDAnalyzer('SkimTrack',
   srcTrack = cms.untracked.InputTag("generalTracks"),
   vzcutmax = cms.untracked.double(15.0),
   vzcutmin = cms.untracked.double(-0.1),
   etaMin = cms.untracked.double(-2.4),
   etaMax = cms.untracked.double(2.4),
   ptMin = cms.untracked.double(0.3),
   ptMax = cms.untracked.double(6.0),
   quality = cms.untracked.int32(0),
   ptEcut = cms.untracked.double(0.1),
   dzcut = cms.untracked.double(3.0),
   dxycut = cms.untracked.double(3.0),
   trkbin = cms.untracked.vint32(150,100)
#   d0cut = cms.untracked.double(3.0),
#   chi2cut = cms.untracked.double(40.0),
#   jetptcut = cms.untracked.double(100.0),
#   minMult = cms.untracked.int32(220),
#   maxMult = cms.untracked.int32(240)
)

#process.p = cms.Path(process.PAcollisionEventSelection * process.demo)
#process.p = cms.Path(process.pACentrality * process.demo)
#process.p = cms.Path(process.demo)
process.p = cms.Path(process.demo)
#process.p = cms.Path(process.hltHM * process.pACentrality * process.demo)
#process.p = cms.Path(process.PAcollisionEventSelection * process.pACentrality * process.demo)
#process.p = cms.Path(process.hltZeroBiasSingleTrack * process.PAcollisionEventSelection * process.pACentrality * process.demo)

