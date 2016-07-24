import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import sys
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#  '/store/hidata/HIRun2013A/PAHighPt/RECO/PromptReco-v1/000/210/634/FA4E6B7E-7366-E211-8DD0-0019B9F581C9.root'
#  'file:/cms/store/hidata/HIRun2013/PAHighPt/RECO/PromptReco-v1/000/210/498/00000/1E825832-FA64-E211-8F9C-003048CF9B28.root'
),
				)
filename = "filelist"+sys.argv[3]+".dat"
mylist = FileUtils.loadListFromFile (filename)
for fname in mylist[int(sys.argv[4]):int(sys.argv[5])]:
 process.source.fileNames.append('file:%s' % (fname))

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V43F::All'
if sys.argv[3] == "pPbReReco":
    process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_210498-210658_HI_PromptReco_Collisions13_JSON_v2.txt').getVLuminosityBlockRange()
elif sys.argv[3] == "pPb":
    process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_210676-211256_HI_PromptReco_Collisions13_JSON_v2.txt').getVLuminosityBlockRange()
elif sys.argv[3] == "Pbp":
    process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_211313-211631_HI_PromptReco_Collisions13_JSON_v2.txt').getVLuminosityBlockRange()

# =============== Import Sequences =====================
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = [
#		"HLT_*PAPixelTracks_Multiplicity*_v*"
               "HLT_PAPixelTracks_Multiplicity100_v*",  #120-150
               "HLT_PAPixelTracks_Multiplicity130_v*",  #150-185
               "HLT_PAPixelTracks_Multiplicity160_v*",  #185-220
               "HLT_PAPixelTracks_Multiplicity190_v*"  #220-260
#               "HLT_PAPixelTracks_Multiplicity220_v*"   #260-300
              ]
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.PAprimaryVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),   # otherwise it won't filter the events
)

process.PAcollisionEventSelection = cms.Sequence(process.hfCoincFilter *
                                         process.PAprimaryVertexFilter *
                                         process.NoScraping
                                         )

process.eventFilter_HM = cms.Sequence( 
    process.hltHM *
    process.PAcollisionEventSelection *
    process.pileupVertexFilterCutGplus
)

#process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#process.load("Analysis.Skim2013.ppanalysisSkimContentFull_cff")
#process.output_HM = cms.OutputModule("PoolOutputModule",
#    outputCommands = process.analysisSkimContent.outputCommands,
#    fileName = cms.untracked.string('pPb_HM.root'),
#    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
#    dataset = cms.untracked.PSet(
#      dataTier = cms.untracked.string('AOD'),
#      filterName = cms.untracked.string('pPb_HM'))
#)
#####combine###############
# histogram output
filenameout = "/lio/lfs/cms/store/user/qixu/flow/PACorrHM/skim/PFcandi/tracktest/pbs/M260220/FlowHM"+sys.argv[3]+"/pPbskimAnaTree_"+sys.argv[2]+".root"
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string(filenameout))

process.demo = cms.EDAnalyzer('SkimPFCand',
   srcTrack = cms.untracked.InputTag("generalTracks"),
   pfCandidateLabel = cms.untracked.InputTag("particleFlow"),
   pfetaMin = cms.untracked.double(-5.0),
   pfetaMax = cms.untracked.double(5.0),
   pfptMin = cms.untracked.double(0.1),
   pfptMax = cms.untracked.double(20.0),
   trkbin = cms.untracked.vint32(260,220),
   doCharge = cms.untracked.int32(0)
)
###########################

#process.output_HM_step = cms.EndPath(process.output_HM)
#process.schedule = cms.Schedule(
#    process.eventFilter_HM_step,
#    process.output_HM_step
#)
process.p = cms.Path(process.eventFilter_HM * process.demo)
