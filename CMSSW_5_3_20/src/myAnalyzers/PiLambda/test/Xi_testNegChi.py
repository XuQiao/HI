import FWCore.ParameterSet.Config as cms

process = cms.Process("XiRECO")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(10000) )
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/F6361662-B9ED-DE11-877E-002618943950.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/F2768061-B9ED-DE11-A4E6-002618943864.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/EC82CF5C-B9ED-DE11-949C-002618943900.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/EA208C63-B9ED-DE11-BA1F-00261894386C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/E8416B42-B9ED-DE11-867B-00261894397F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/E266AC46-B9ED-DE11-91C3-002618943900.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/E261AF48-B9ED-DE11-B614-002618FDA208.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/E2383F82-B9ED-DE11-88A2-002618943944.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/DA317A74-B9ED-DE11-ABC1-00248C0BE014.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/D27C3292-B9ED-DE11-B457-002618943868.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/D25ED03F-B9ED-DE11-AF35-0026189438BA.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/CE144D48-B9ED-DE11-9AAF-002618943884.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/CC7F9E8D-B9ED-DE11-9356-00261894396A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/C8AA8F3E-B9ED-DE11-95F5-002618943944.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/C8A8AC5F-B9ED-DE11-8A5B-002618943864.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/BA359B47-B9ED-DE11-A00F-0026189438DE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/B870B05D-B9ED-DE11-9140-002618943932.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/B2C21F52-B9ED-DE11-9FA2-002618943930.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/AEF26C42-B9ED-DE11-B27B-002618943907.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/AE3BCF55-B9ED-DE11-9E83-00261894388D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/AA55D651-B9ED-DE11-9E5A-00261894382D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/A64B795D-B9ED-DE11-AE61-0026189438DE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/A2767553-B9ED-DE11-9381-00248C0BE014.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/9E23105B-B9ED-DE11-AFAD-002618943907.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/9C264FD4-2CEE-DE11-9CB1-002618943945.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/98DB7543-B9ED-DE11-B036-002618943867.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/9457AE65-B9ED-DE11-B87D-002618943884.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/92680451-B9ED-DE11-B35D-002618943961.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/8AB0154B-B9ED-DE11-AD43-002618943900.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec19thSkim_341_v2/0006/8645194D-B9ED-DE11-92A3-00261894396C.root'
)
)

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange(
'124275:3-124275:30',
'124230:26-124230:max',
'124120:1-124120:max',
'124030:1-124030:max',
'124027:24-124027:max',
'124025:3-124025:13',
'124024:2-124024:83',
'124023:38-124023:96',
'124022:66-124022:179',
'124020:12-124020:94',
'124009:1-124009:68',
'124008:1-124008:max',
'124006:1-124006:max',
'123908:2-123908:12',
'123906:18-123906:28',
'123818:2-123818:42',
'123815:8-123815:max',
'123732:56-123732:109',
'123596:2-123596:max',
'123592:2-123592:12'
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR09_R_V5::All')
process.GlobalTag.globaltag = 'GR09_R_34X_V2::All'
#process.GlobalTag.globaltag = cms.string('STARTUP3X_V8I::All')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

# Real data raw to digi
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

# reconstruction sequence for Collisions
process.load('Configuration/StandardSequences/Reconstruction_cff')

# offline DQM
process.load('DQMOffline/Configuration/DQMOffline_cff')

# import of additional standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.localV0Candidates = cms.EDProducer("V0Producer",
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    useSmoothing = cms.bool(True),
    storeSmoothedTracksInRecoVertex = cms.bool(False),
    doPostFitCuts = cms.bool(True),
    doTrackQualityCuts = cms.bool(True),
    # The next parameters are cut values
    # Track quality cuts
    #   Normalized track Chi2:
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track:
    tkNhitsCut = cms.int32(6),

    # Vertex cuts
    vtxChi2Cut = cms.double(7.0),
    collinearityCut = cms.double(0.02),
    #  Setting this one to zero; significance cut is sufficient
    rVtxCut = cms.double(0.0),
#    vtxSignificanceCut = cms.double(15.0),
    vtxSignificanceCut = cms.double(10.0),
    kShortMassCut = cms.double(0.07),
    lambdaMassCut = cms.double(0.05),
    impactParameterSigCut = cms.double(0.5),
    mPiPiCut = cms.double(1.),
    tkDCACut = cms.double(1.),

    # These parameters decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    vertexFitter = cms.InputTag('KalmanVertexFitter')

)

process.mkcands = cms.EDAnalyzer("PiLambda",
  HepMC       = cms.untracked.string('source'),
  VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
  Blist       = cms.untracked.string('Bd')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('PiLambda_BSCNOBEAMHALO-Dec19-bit0-noscrape-bugfix-lumiSel1.root')
)

# this is for filtering on L1 technical trigger bit
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
# bsc minbias in coinidence with bptx and veto on beam halo
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')

#apply the scraping event filter here
process.noScraping= cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

process.ntup = cms.Path( process.hltLevel1GTSeed * process.noScraping * process.localV0Candidates * process.mkcands )

process.schedule = cms.Schedule(process.ntup )
