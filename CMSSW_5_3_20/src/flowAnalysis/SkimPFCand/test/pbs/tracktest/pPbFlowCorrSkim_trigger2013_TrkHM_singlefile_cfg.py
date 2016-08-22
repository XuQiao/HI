import FWCore.ParameterSet.Config as cms
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
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 200

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
#  '/store/user/davidw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_v10/75efa6cce8565664c7fea989d6c24a72/pPb_HM_285_1_lwi.root'
#'/store/group/phys_heavyions/icali/PAPhysics/pAPilotRun_Run202792GoodLumis_RAWRECO_L1Em_PrescaleActiveBitsSkimNoZB_CMSSW528_V94_FinalWorkflow_200kHz_v2_v1_v2/f3394926c5028783289fd2cd57b36909/PAPhysics_RAWRECO_inRECO_9_1_Oxr.root',
  '/store/hidata/HIRun2013A/PAHighPt/RECO/PromptReco-v1/000/210/634/FA4E6B7E-7366-E211-8DD0-0019B9F581C9.root'
)
)

# =============== Other Statements =====================
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.GlobalTag.globaltag = 'GR_P_V43F::All'

# =============== Import Sequences =====================
#process.load("davidlw.HighPtFlow2011.ppExtraReco_cff")
process.load("Analysis.SkimPFCand.ppExtraReco_cff")
#process.load("davidlw.HighPtFlow2011.PAHighMultiplicityPileUpFilter_cff")
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')
process.load("UserCode.FerencSiklerEnergyLoss.EnergyLossProducer_cff")

#Trigger Selection
### Comment out for the timing being assuming running on secondary dataset with trigger bit selected already
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM.HLTPaths = ['HLT_*PixelTracks_Multiplicity*_v*']
process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.eventFilter_HM = cms.Sequence( 
    process.hltHM *
    process.PAcollisionEventSelection *
    process.pileupVertexFilterCutGplus
)

process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

process.pACentrality_step = cms.Path( process.eventFilter_HM * process.pACentrality)
process.pACentrality.producePixelhits = False

process.dEdx_step = cms.Path( process.eventFilter_HM * process.dEdx )

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

########## MinBias tracking ###################################################################
# Tracker local reco
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

# Beamspot
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoPixelVertexing.PixelLowPtUtilities.MinBiasTracking_cff")

# Agglomerative vertexing
import UserCode.FerencSiklerVertexing.NewVertexProducer_cfi

process.allVertices = UserCode.FerencSiklerVertexing.NewVertexProducer_cfi.newVertices.clone()
process.allVertices.TrackCollection = 'allTracks'
process.allVertices.PtMin  = cms.double(0.1)

# V0
process.pixelVZeros.trackCollection  = cms.InputTag("allTracks")
process.pixelVZeros.vertexCollection = cms.InputTag("allVertices")
process.pixelVZeros.minImpactPositiveDaughter = cms.double(0.2)
process.pixelVZeros.minImpactNegativeDaughter = cms.double(0.2)

# Cluster-shape filter re-run offline
import HLTrigger.special.hltPixelClusterShapeFilter_cfi

process.hltPixelClusterShapeFilter = HLTrigger.special.hltPixelClusterShapeFilter_cfi.hltPixelClusterShapeFilter.clone()
process.hltPixelClusterShapeFilter.inputTag = "siPixelRecHits"

process.generalV0CandidatesMB = process.generalV0Candidates.clone( trackRecoAlgorithm = cms.InputTag('allTracks') )

process.MBtracking = cms.Sequence(
                        process.siPixelRecHits
                      * process.hltPixelClusterShapeFilter
                      * process.siStripMatchedRecHits
                      * process.minBiasTracking
                      * process.allVertices
                      * process.generalV0CandidatesMB
                       )

process.MBtracking_step = cms.Path( process.eventFilter_HM * process.MBtracking )

###############################################################################################

process.load("davidlw.HighPtFlow2011.ppanalysisSkimContentFull_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('pPb_HM.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('AOD'),
      filterName = cms.untracked.string('pPb_HM'))
)

process.output_HM_step = cms.EndPath(process.output_HM)

process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
#    process.extraTrks_HM_step,
    process.pACentrality_step,
#    process.dEdx_step,
#    process.MBtracking_step, 
    process.reTracking_step,
    process.v0rereco_step,
#    process.v0rereco_refitted_step,
    process.output_HM_step
)
