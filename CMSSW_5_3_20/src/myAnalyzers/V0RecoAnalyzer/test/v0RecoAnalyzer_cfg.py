import FWCore.ParameterSet.Config as cms

process = cms.Process('V0ANALYSIS')

# Standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# Specifically needed for the analyzer
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Geometry_cff')

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
#    vtxSignificanceCut = cms.double(22.0),
    vtxSignificanceCut = cms.double(15.0),
    kShortMassCut = cms.double(0.2),
    lambdaMassCut = cms.double(0.2),
    impactParameterSigCut = cms.double(2.0),
    mPiPiCut = cms.double(0.6),
    tkDCACut = cms.double(1.),

    # These parameters decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    vertexFitter = cms.InputTag('KalmanVertexFitter')

)


process.ksAnalyzer = cms.EDAnalyzer('V0RecoAnalyzer',
    v0Collection = cms.InputTag('generalV0Candidates:Kshort'),
    writeTree = cms.bool(True),
    writeHistos = cms.bool(True),
    ksMassXmin = cms.double(0.44),
    ksMassXmax = cms.double(0.56),
    ksMassNbins = cms.int32(60),
    lamMassXmin = cms.double(1.07),
    lamMassXmax = cms.double(1.16),
    lamMassNbins = cms.int32(60)
)

process.lamAnalyzer = cms.EDAnalyzer('V0RecoAnalyzer',
    v0Collection = cms.InputTag('generalV0Candidates:Lambda'),
    writeTree = cms.bool(True),
    writeHistos = cms.bool(True),
    ksMassXmin = cms.double(0.44),
    ksMassXmax = cms.double(0.56),
    ksMassNbins = cms.int32(60),
    lamMassXmin = cms.double(1.07),
    lamMassXmax = cms.double(1.16),
    lamMassNbins = cms.int32(60)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_3_1_2/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v2/0011/F605AE2E-D990-DE11-BFCB-001D09F2545B.root',
       '/store/relval/CMSSW_3_1_2/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v2/0011/441FC5C8-4C91-DE11-BC78-000423D6BA18.root',
       '/store/relval/CMSSW_3_1_2/RelValQCD_Pt_80_120/GEN-SIM-RECO/MC_31X_V3-v2/0011/1210F6E7-D990-DE11-A327-001D09F2545B.root' 
    )
)

process.GlobalTag.globaltag = 'MC_31X_V3::All'

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('v0analysis.root')
)

# Path definition
process.analysis = cms.Path(process.ksAnalyzer+process.lamAnalyzer)

process.schedule = cms.Schedule(process.analysis)
