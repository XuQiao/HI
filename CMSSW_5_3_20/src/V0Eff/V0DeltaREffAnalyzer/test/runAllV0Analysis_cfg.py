import FWCore.ParameterSet.Config as cms

process = cms.Process('V0COLL')

process.load('CondCore.DBCommon.CondDBSetup_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
#process.load('Configuration/TrackingAnalysis/trackingParticles_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.load("myAnalyzers/V0RecoAnalyzer/runV0Analysis_cff")
process.analyzeKshort.writeTree = cms.bool(True)
process.analyzeLambda.writeTree = cms.bool(True)
process.localV0Candidates.vtxSignificance2DCut = cms.double(10.0)

process.ntupleXi = cms.EDAnalyzer("PiLambda",
 HepMC       = cms.untracked.string('source'),
 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
 TrackHypo   = cms.untracked.string('pion'),
 V0Collection = cms.untracked.string('localV0Candidates')
)

process.ntupleOm = cms.EDAnalyzer("PiLambda",
 HepMC       = cms.untracked.string('source'),
 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
 TrackHypo   = cms.untracked.string('kaon'),
 V0Collection = cms.untracked.string('localV0Candidates')
)

process.genParticlesPlusSim = cms.EDProducer("GenPlusSimParticleProducer",
        src           = cms.InputTag("g4SimHits"), # use "famosSimHits" for FAMOS
        setStatus     = cms.int32(8),
#        particleTypes = cms.vstring(""),
        filter        = cms.vstring("pt > 0.0"),  # just for testing
        genParticles   = cms.InputTag("genParticles")
)

process.load('V0Eff.V0DeltaREffAnalyzer.V0DeltaREffAnalyzer_cfi')

# Hits matching, for comparison, if the MC sample has TrackingParticles (ONLY!)
#process.load('SimTracker.TrackAssociation.TrackAssociatorByHits_cfi')
#process.load('SimTracker.TrackAssociation.trackMCMatch_cfi')
#process.load('SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cff')

#process.matchVeesByHits = cms.EDAnalyzer('EffAnalyzer',
#    associationModule = cms.string('TrackAssociatorByHits'),
#    trackingParticleModule = cms.string("mergedtruth"),
#    trackingParticleProduct = cms.string('MergedTrackTruth'),
#    bestMatchByMaxValue = cms.bool(True),
#    recoTrackModule = cms.string('generalTracks'),
#    v0Collection = cms.InputTag('generalV0Candidates:Kshort')
#)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/nfs/data35/cms/drell/7TeVMC/54E25B24-F33D-DF11-BAD0-002618943836.root'
    )#,
#   inputCommands = cms.untracked.vstring("keep *",
#                                         "drop *_gctDigis_*_*",
#                                         "drop *_gtDigis_*_*",
#					 "drop *_*_*_HLT")
)


process.GlobalTag.globaltag = cms.string('START3X_V26A::All')
#process.GlobalTag.globaltag = cms.string('START3X_V22A:All')

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *',
					   'keep recoBeamSpot_*_*_*',
					   'keep *_generalTracks_*_*',
					   'keep recoVertexs_*_*_*',
					   'keep recoVertexCompositeCandidates_*_*_*',
                                           'keep patCompositeCandidates_*_*_*'),
    fileName = cms.untracked.string('v0Collection.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('vee_step')
    )
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('analysis.root')
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
#process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
# Bit zero doesn't work for MC
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.2)
)

process.vee_step = cms.Path(process.hltLevel1GTSeed * process.noScraping * process.v0analysis * process.ntupleXi * process.ntupleOm)
process.match_step = cms.Path(process.hltLevel1GTSeed * process.noScraping * process.genParticlesPlusSim * process.matchVeesByDeltaR)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.vee_step, process.match_step, process.endjob_step, process.out_step)
#process.schedule = cms.Schedule(process.vee_step, process.match_step, process.endjob_step)
