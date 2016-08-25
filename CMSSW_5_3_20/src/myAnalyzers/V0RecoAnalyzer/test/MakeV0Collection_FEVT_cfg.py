import FWCore.ParameterSet.Config as cms

process = cms.Process('V0COLL')

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
    tkChi2Cut = cms.double(15.0),
    #   Number of valid hits on track:
    tkNhitsCut = cms.int32(6),

    # Vertex cuts
    vtxChi2Cut = cms.double(7.0),
    collinearityCut = cms.double(0.02),
    #  Setting this one to zero; significance cut is sufficient
    rVtxCut = cms.double(0.0),
#    vtxSignificanceCut = cms.double(22.0),
#    vtxSignificanceCut = cms.double(15.0),
    vtxSignificanceCut = cms.double(5.0),
    kShortMassCut = cms.double(1.),
    lambdaMassCut = cms.double(1.),
    impactParameterSigCut = cms.double(0.5),
    mPiPiCut = cms.double(1.),
    tkDCACut = cms.double(1.),

    # These parameters decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    vertexFitter = cms.InputTag('KalmanVertexFitter')

)

process.v0CollectionMaker = cms.EDProducer("V0CandProducer",
    kShortCollection = cms.InputTag('localV0Candidates:Kshort'),
    lambdaCollection = cms.InputTag('localV0Candidates:Lambda')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    IgnoreCompletely = cms.untracked.vstring('MismatchedInputFiles')
)
     
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/F42A9D46-4FEB-DE11-82D9-001D0967D77E.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/F2D30717-73EA-DE11-A3E9-0024E8768C98.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/E848193C-72EA-DE11-B1CB-0024E8768BFC.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/E0D9773B-72EA-DE11-9011-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/C4781147-72EA-DE11-A119-0024E876A7FA.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/BABAF8C3-71EA-DE11-9D8C-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/BA547B60-72EA-DE11-866D-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/BA133542-72EA-DE11-8AD9-0024E86E8CF1.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/B6B9118D-71EA-DE11-B162-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/AE91A7A8-72EA-DE11-877A-0024E86E8D8D.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/A0666BE7-71EA-DE11-8073-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/9E264E41-72EA-DE11-ACC1-0024E87680CD.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/9A618F92-71EA-DE11-90F1-0024E876841F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/90DABED5-71EA-DE11-B87B-00151796D7F4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/8A7AB41D-72EA-DE11-8410-00151796D7F4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/8A5782D2-71EA-DE11-956F-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/8478D6DB-71EA-DE11-8ED9-0024E876841F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/80E73354-72EA-DE11-ABD3-0024E8768867.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/7E8C1A43-72EA-DE11-B327-001D0967D37D.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/7666748B-71EA-DE11-9BF7-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/70DF824A-72EA-DE11-AF49-0024E8769B1F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/649DF6E2-72EA-DE11-B5F5-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/6081FF95-72EA-DE11-94C4-0024E8768446.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/584EC24A-72EA-DE11-9C7D-0015178C4D14.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/4EC6CE7F-72EA-DE11-BF3D-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/4A7B4571-72EA-DE11-A4A2-0024E8768BFC.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/38520E07-72EA-DE11-9E0A-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/361E798C-71EA-DE11-9FC4-0024E876994B.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/32AE7BB7-71EA-DE11-9F43-0024E876841F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/1A7142E5-72EA-DE11-A19D-0015178C65F4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/1077DAAB-72EA-DE11-8F30-0024E8768867.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/104A4AAA-72EA-DE11-9CA4-00151796D7F4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/087E4821-72EA-DE11-BC72-0015178C65F4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/08346FDF-4EEB-DE11-85A0-001D0967B82E.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/E0811A24-1CEA-DE11-AD64-001D0967DDC3.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/CA8A88FD-27EA-DE11-9305-00151796D760.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/C643DE7F-1DEA-DE11-BE76-001D0967DE13.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/C2B02E96-25EA-DE11-B218-00151796D788.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/B27F9983-1DEA-DE11-A300-0015178C4900.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/4C3C6485-1DEA-DE11-B1B0-001D0967C649.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0101/4AB376B6-1EEA-DE11-A9E9-0015179EDD24.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/F64AD33D-0BEA-DE11-8D9C-0024E8766415.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/F491C10E-0CEA-DE11-85BB-001D0967D49F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/F26FBF26-0BEA-DE11-9AB8-001D0967DB5C.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/EE38F323-0BEA-DE11-A3C4-001D0967D49F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/EE169487-19EA-DE11-88F4-001D0967DA49.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/E81CFDA6-18EA-DE11-898C-001D0967DEF9.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/E0D55622-10EA-DE11-9C3B-001D0967DD0A.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/D2F94DD2-11EA-DE11-84CE-00151796C100.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/D2CCFC84-19EA-DE11-88F1-001D0967DA49.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/CE5ADEDA-0CEA-DE11-B193-00151796D428.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/C81EE936-0BEA-DE11-AED5-0024E8766393.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/C268E61D-10EA-DE11-A0D0-0015178C4B84.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/BE9B3B86-19EA-DE11-9588-001D0967D49F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/B0398B1F-0FEA-DE11-B4AD-00151796D508.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/A6059310-18EA-DE11-90D9-00151796D4B0.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/A4E08436-13EA-DE11-866D-00151796D508.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/98A7CE36-13EA-DE11-ACEC-001D0967DA3A.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/9862E51D-0FEA-DE11-BF05-001D0967DFB7.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/8E97B2F6-0BEA-DE11-91EF-001D0967CF95.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/8E502D21-1CEA-DE11-8B7E-0015178C6704.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/8C89C0F8-0BEA-DE11-9582-00151796D678.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/88681630-17EA-DE11-8259-0024E8768299.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7EEB456F-16EA-DE11-BEDC-001D0967D49F.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7EB6989F-18EA-DE11-A96E-001D0967D37D.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7C40AEFB-14EA-DE11-91FC-0015178C4900.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/72034E8A-19EA-DE11-B723-001D0966E23E.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6EEC7E88-19EA-DE11-A63F-00151796C138.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6AC26C39-17EA-DE11-B954-001D0967CE50.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6A12FB1C-0AEA-DE11-AE97-0024E8768C23.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/66E26D28-10EA-DE11-85C5-00151796D508.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5C8DED2A-0BEA-DE11-8453-001D0967D2DD.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/58987C27-0FEA-DE11-84D9-0024E876808C.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5854992D-17EA-DE11-8E61-00151796D508.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5685AC0D-18EA-DE11-B2F9-001D096760DE.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/4EADFB21-1CEA-DE11-9F0D-00151796C144.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/4AEE4E36-13EA-DE11-9087-00151796D768.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/3AA0B731-17EA-DE11-B535-00151796C1C8.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/30756C31-0BEA-DE11-8BA8-0024E876A82E.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/2E885270-16EA-DE11-BDCD-001D0967CFA9.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/2C0AF991-19EA-DE11-A1F9-0015178C6704.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1E5AD126-0BEA-DE11-9D1D-0015178C48E4.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1CA7FB7D-1BEA-DE11-8EA3-00151796D4B0.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1A55661D-0BEA-DE11-985D-0024E8769B60.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1A14841B-1CEA-DE11-AEEA-0024E8768265.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/10810A7E-16EA-DE11-9FB6-00151796C178.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/0E84AA51-1AEA-DE11-B501-0015179EDF00.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/0AA60F31-17EA-DE11-A6C9-001D0967DA3A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/08AC2C30-0BEA-DE11-87AE-001D0967DB7A.root',
        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/00D4856F-16EA-DE11-A4C8-00151796C0F0.root'
    ),
   inputCommands = cms.untracked.vstring("keep *",
                                         "drop *_gctDigis_*_*",
                                         "drop *_gtDigis_*_*",
					 "drop *_*_*_HLT")
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('123596:1-123596:143')
process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('123592:1-123734:max',
   								    '124009:1-124009:max',
								    '124020:1-124020:max',
								    '124022:1-124024:max',
								    '124027:1-124017:max',
								    '124030:1-124030:max'
								    )

process.GlobalTag.globaltag = 'GR09_R_V4::All'

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *',
					   'keep recoBeamSpot_*_*_*',
					   'keep *_generalTracks_*_*',
					   'keep recoVertexs_*_*_*',
					   'keep recoVertexCompositeCandidates_*_*_*',
                                           'keep patCompositeCandidates_v0CollectionMaker_*_*'),
    fileName = cms.untracked.string('v0Collection-LOOSEVTXSIG.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

process.vee_step = cms.Path(process.localV0Candidates)
process.collect_step = cms.Path(process.v0CollectionMaker)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.vee_step, process.collect_step, process.endjob_step, process.out_step)
