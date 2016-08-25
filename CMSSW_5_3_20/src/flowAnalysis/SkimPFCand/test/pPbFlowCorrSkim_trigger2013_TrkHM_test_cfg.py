import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
process = cms.Process("ANASKIM")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#  '/store/hidata/HIRun2013A/PAHighPt/RECO/PromptReco-v1/000/210/634/FA4E6B7E-7366-E211-8DD0-0019B9F581C9.root'
  'file:/cms/store/hidata/HIRun2013/PAHighPt/RECO/PromptReco-v1/000/210/498/00000/1E825832-FA64-E211-8F9C-003048CF9B28.root'
)
)
#mylist = FileUtils.loadListFromFile ('filelistAll.dat')
#for fname in mylist:
# process.source.fileNames.append('file:%s' % (fname))

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(3910))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V43F::All'

# =============== Import Sequences =====================
#process.load("Analysis.SkimPFCand.ppExtraReco_cff")
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = [
#		"HLT_*PAPixelTracks_Multiplicity*_v*"
               "HLT_PAPixelTracks_Multiplicity100_v*"  #120-150
#               "HLT_PAPixelTracks_Multiplicity130_v*",  #150-185
#               "HLT_PAPixelTracks_Multiplicity160_v*",  #185-220
#               "HLT_PAPixelTracks_Multiplicity190_v*",  #220-260
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

process.pACentrality_step = cms.Path( process.eventFilter_HM * process.pACentrality)
process.pACentrality.producePixelhits = False

########## V0 candidate rereco ###############################################################
process.generalV0CandidatesNew = process.generalV0Candidates.clone (
    tkNhitsCut = cms.int32(0),
    tkChi2Cut = cms.double(7.0),
    dauTransImpactSigCut = cms.double(0.5),
    dauLongImpactSigCut = cms.double(0.5),
    xiVtxSignificance3DCut = cms.double(0.0),
    xiVtxSignificance2DCut = cms.double(0.0),
    vtxSignificance2DCut = cms.double(0.0),
    vtxSignificance3DCut = cms.double(4.0)
)
process.v0rereco_step = cms.Path( process.eventFilter_HM * process.generalV0CandidatesNew )

########## ReTracking #########################################################################
process.generalTracksLowPt = process.generalTracks.clone()
process.iterTracking.replace(process.generalTracks,process.generalTracksLowPt)
process.offlinePrimaryVerticesLowPt = process.offlinePrimaryVertices.clone( TrackLabel = cms.InputTag("generalTracksLowPt") )
process.generalV0CandidatesLowPt = process.generalV0CandidatesNew.clone(
  trackRecoAlgorithm = cms.InputTag('generalTracksLowPt'),
  vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVerticesLowPt')
)

process.lowPtTripletStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.075
process.detachedTripletStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.075

process.dedxTruncated40LowPt = process.dedxTruncated40.clone(
    tracks                     = cms.InputTag("generalTracksLowPt"),
    trajectoryTrackAssociation = cms.InputTag("generalTracksLowPt")
)
process.dedxHarmonic2LowPt = process.dedxHarmonic2.clone(
    tracks                     = cms.InputTag("generalTracksLowPt"),
    trajectoryTrackAssociation = cms.InputTag("generalTracksLowPt")
)
process.dedxDiscrimASmiLowPt = process.dedxDiscrimASmi.clone(
    tracks                     = cms.InputTag("generalTracksLowPt"),
    trajectoryTrackAssociation = cms.InputTag("generalTracksLowPt")
)
process.trackingGlobalReco.replace(process.dedxTruncated40,process.dedxTruncated40LowPt)
process.trackingGlobalReco.replace(process.dedxHarmonic2,process.dedxHarmonic2LowPt)
process.trackingGlobalReco.replace(process.dedxDiscrimASmi,process.dedxDiscrimASmiLowPt)

process.reTracking = cms.Sequence(
   process.siPixelRecHits *
   process.siStripMatchedRecHits *
   process.recopixelvertexing *
   process.trackingGlobalReco *
   process.offlinePrimaryVerticesLowPt *
   process.generalV0CandidatesLowPt
)

process.reTracking_step = cms.Path( process.eventFilter_HM * process.reTracking )

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
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("pPbskimAnaTree_no.root"))

process.demo = cms.EDAnalyzer('SkimPFCand',
   srcTrack = cms.untracked.InputTag("generalTracks"),
   pfCandidateLabel = cms.untracked.InputTag("particleFlow"),
   pfetaMin = cms.untracked.double(-5.0),
   pfetaMax = cms.untracked.double(5.0),
   pfptMin = cms.untracked.double(0.1),
   pfptMax = cms.untracked.double(10.0),
   trkbin = cms.untracked.vint32(150,120),
   doCharge = cms.untracked.int32(0)
)
###########################

#process.output_HM_step = cms.EndPath(process.output_HM)
#process.schedule = cms.Schedule(
#    process.eventFilter_HM_step,
#    process.output_HM_step
#)
process.p = cms.Path(process.eventFilter_HM * process.demo)
process.schedule = cms.Schedule(
#	process.reTracking_step,
#        process.v0rereco_step,
	process.p
)

