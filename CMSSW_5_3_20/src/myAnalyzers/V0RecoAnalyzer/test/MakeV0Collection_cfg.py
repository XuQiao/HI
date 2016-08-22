import FWCore.ParameterSet.Config as cms

process = cms.Process('V0COLL')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

# Specifically needed for the analyzer
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load('Configuration/StandardSequences/MagneticField_cff')
#process.load('Configuration/StandardSequences/Geometry_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')

#From Eric: Trigger Path
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltSingleTrigger = process.hltHighLevel.clone()
process.hltSingleTrigger.HLTPaths = ["HLT_PAZeroBiasPixel_SingleTrack_v1"]
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')
		
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
#    vtxSignificanceCut = cms.double(15.0),
    vtxSignificanceCut = cms.double(5.0),
    vtxSignificance2DCut = cms.double(15.0),	#by Qiao
    kShortMassCut = cms.double(1.),
    lambdaMassCut = cms.double(1.),
    impactParameterSigCut = cms.double(0.5),
    mPiPiCut = cms.double(1.),
    tkDCACut = cms.double(1.),
trackQualities = cms.vstring('loose'),	#by Qiao

    # These parameters decide whether or not to reconstruct
    #  specific V0 particles
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(True),

    vertexFitter = cms.InputTag('KalmanVertexFitter'),
    
    innerHitPosCut = cms.double(4.)	#by Qiao

)

process.v0CollectionMaker = cms.EDProducer("V0CandProducer",
    kShortCollection = cms.InputTag('localV0Candidates:Kshort'),
    lambdaCollection = cms.InputTag('localV0Candidates:Lambda')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    IgnoreCompletely = cms.untracked.vstring('MismatchedInputFiles')
)
     
process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring(
'/store/hidata/HIRun2013/PAMinBiasUPC/RECO/PromptReco-v1/000/210/818/00000/2C5B71D3-FF68-E211-AAE4-001D09F2441B.root'
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/00D1A776-B477-E211-BE20-5404A63886D6.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/043376E1-BA77-E211-8E53-001D09F28D54.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/04BDF340-B877-E211-8256-003048F0258C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/06B7204E-B477-E211-8D50-002481E94050.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/080794AC-B377-E211-85D8-001D09F24353.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/085B6ECE-BD77-E211-A46A-001D09F23F2A.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/08E090A2-B577-E211-8E54-0019B9F581C9.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0A77959E-B777-E211-B860-003048D2BDD8.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0C1C02A6-B577-E211-858D-002481E0D646.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0C2F54A5-B577-E211-AC4D-003048D2BBF0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0E386ADF-B477-E211-9B3E-001D09F2A690.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0E4A93DE-B677-E211-9015-001D09F24303.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0EA91DA6-B577-E211-BF44-001D09F26509.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/0EC1E550-B777-E211-9B54-001D09F28EA3.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/122EC669-B677-E211-9910-003048F117F6.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/127D0BDF-C077-E211-84FB-BCAEC532970D.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/14586FBD-BA77-E211-BDEA-003048D37580.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/165FAC67-B677-E211-974C-00215AEDFD74.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/16E70AA6-B977-E211-9E01-00237DDBE49C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/1CD92569-B677-E211-9644-002481E0D7C0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/225303E0-B477-E211-AF69-001D09F2305C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/22A0DD18-B977-E211-9577-001D09F28E80.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/26DC5068-B677-E211-BE75-003048F1C58C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/2C4F36A6-B577-E211-AE1B-002481E0D7EC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/2CAB5ADE-C077-E211-B086-BCAEC5329708.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/32C45CE0-B677-E211-ABE6-003048D2BE12.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/34CFB9B0-B377-E211-B7F5-001D09F28F25.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/36091869-B677-E211-B2BD-003048F118C4.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/36E19FA2-B377-E211-8675-BCAEC518FF62.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/384DCF77-B477-E211-9300-5404A63886EE.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/3E9E8F68-B677-E211-98DE-0025B32445E0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/3EAD47A0-B777-E211-8384-003048678098.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/42A11969-B677-E211-84F9-003048D2BA82.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/44DAB8A3-B577-E211-AA9D-003048F118AA.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/4A69F046-B977-E211-BFA5-003048D2BC38.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/4C5B22E2-B677-E211-8B47-001D09F295A1.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/4E6C9E67-B677-E211-AB93-0025B32035BC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/4E9B2D69-B577-E211-86F9-001D09F252DA.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/522ED94C-B777-E211-87BB-001D09F25041.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/52E89CE0-B477-E211-8068-001D09F2960F.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/56D5DDE6-B377-E211-9619-001D09F2525D.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/58575EB0-B377-E211-86DB-001D09F244DE.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5A53ADD6-C277-E211-9697-003048D2BBF0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5ABFB2DF-B677-E211-B5F5-0019B9F581C9.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5AC407E0-B677-E211-ADFE-001D09F28F25.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5AD04F10-C477-E211-A478-BCAEC518FF80.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5C73B7A8-B577-E211-84F1-001D09F24D67.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5CC9F9AB-B577-E211-A8F9-001D09F24682.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/5E015269-B677-E211-BE9D-003048F0258C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/6027DFBB-B377-E211-9CA8-003048D2BB90.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/6203E7E4-B377-E211-8214-BCAEC532971F.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/6289C0E9-B377-E211-9523-001D09F2305C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/648238A9-B377-E211-AABE-BCAEC5329716.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/661C8DE5-B377-E211-ACB8-5404A6388693.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/68CC1DDF-B377-E211-A775-5404A63886D6.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/68ECC377-B477-E211-A8A7-5404A63886BE.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/6C070269-B677-E211-8362-002481E0D7EC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/705737AE-B577-E211-AEFD-001D09F2A49C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/708A345D-B877-E211-84EC-0030486780B4.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/76173F42-B877-E211-99E2-002481E0D646.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/768648DF-C077-E211-A622-0025901D5DB2.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/76F07065-C177-E211-BA80-003048F117EC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/783DBA16-B477-E211-9475-0025901D62A6.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/7AFA86A7-B577-E211-A751-001D09F24763.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/7CC10E76-B477-E211-8966-BCAEC518FF67.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/8023E3E0-B477-E211-89EA-0019B9F72CE5.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/82903476-B477-E211-855F-BCAEC518FF8A.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/86C96BE9-B377-E211-8195-001D09F24D67.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/88BBE8B6-B677-E211-B618-001D09F253C0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/88F3F0E4-C077-E211-B8B7-BCAEC518FF68.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/8A2415A3-B577-E211-91AA-002481E0D73C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/8C5C7BE1-B377-E211-BFC8-0025B320384C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/8CB11BA9-B577-E211-BB68-001D09F2A465.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9061E39E-B777-E211-A73D-003048D3C980.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/92A72E67-B677-E211-9A74-00237DDBE41A.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/962AE477-B477-E211-BC3E-BCAEC532971C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9803E7DF-B477-E211-9AF2-001D09F24DA8.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9838E9E7-B677-E211-BE1A-003048D2BC30.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9848E500-BA77-E211-9175-0030486780B8.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/98C46276-B477-E211-9BFD-BCAEC53296F4.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9C466AA7-B977-E211-8236-00215AEDFD98.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9E0E35DE-C077-E211-AE8C-BCAEC5364C4C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/9EFA66A1-B777-E211-B857-003048D2BC30.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/A2F6EF58-B877-E211-BAEA-003048D374F2.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/A4018E8D-BA77-E211-B72C-001D09F34488.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/A4C20C6B-B677-E211-BE03-00215AEDFCCC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/AA0B123F-BB77-E211-80FB-002481E0D646.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/ACF8F9A5-B977-E211-B9AA-0025B32035BC.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/B01F9EE0-B477-E211-8FD8-001D09F23C73.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/B03052A7-B577-E211-A7A5-003048D2BA82.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/B05211DE-B677-E211-AA2A-001D09F25460.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/B28987B2-B077-E211-8728-5404A640A639.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/B474E8A6-B977-E211-A66B-003048F11942.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/BE402846-B977-E211-B620-001D09F276CF.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/C8101A69-B677-E211-A29A-003048F1BF68.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/C8244B86-B877-E211-9E14-001D09F2983F.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/D0249167-B677-E211-8B5D-0025B320384C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/D275ED5E-BA77-E211-BF47-003048F1BF68.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/D468B609-BB77-E211-BB39-001D09F242EF.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/DA1971A2-B577-E211-AD26-003048D2BEAA.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/E4DFD9AC-B377-E211-B484-001D09F253C0.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/E6127DAF-B077-E211-8927-BCAEC5364C4C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/E6D59C6C-B677-E211-88E0-002481E94050.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/E6F7F174-B477-E211-899F-BCAEC518FF8E.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/EA7ECCA4-B577-E211-9D29-00215AEDFD74.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/EE03F5DF-B477-E211-94A3-001D09F24024.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/EE7C7214-B777-E211-BB5C-5404A63886C4.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/F044D6DE-B477-E211-8334-001D09F26509.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/F4811677-B477-E211-9583-BCAEC518FF74.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/F6554077-B477-E211-9E56-BCAEC5329709.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/FC3A7E6A-B677-E211-B8F3-003048CF94A8.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/FC9F1B68-B677-E211-8297-002481E0E56C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/FE28A21C-B977-E211-9EF0-001D09F2A49C.root',
#'/store/data/Run2013A/PPMinBias/RECO/PromptReco-v1/000/211/797/00000/FE6DB057-BC77-E211-AB43-003048D37538.root',
),
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
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7EEB456F-16EA-DE11-BEDC-001D0967D49F.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7EB6989F-18EA-DE11-A96E-001D0967D37D.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/7C40AEFB-14EA-DE11-91FC-0015178C4900.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/72034E8A-19EA-DE11-B723-001D0966E23E.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6EEC7E88-19EA-DE11-A63F-00151796C138.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6AC26C39-17EA-DE11-B954-001D0967CE50.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/6A12FB1C-0AEA-DE11-AE97-0024E8768C23.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/66E26D28-10EA-DE11-85C5-00151796D508.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5C8DED2A-0BEA-DE11-8453-001D0967D2DD.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/58987C27-0FEA-DE11-84D9-0024E876808C.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5854992D-17EA-DE11-8E61-00151796D508.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/5685AC0D-18EA-DE11-B2F9-001D096760DE.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/4EADFB21-1CEA-DE11-9F0D-00151796C144.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/4AEE4E36-13EA-DE11-9087-00151796D768.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/3AA0B731-17EA-DE11-B535-00151796C1C8.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/30756C31-0BEA-DE11-8BA8-0024E876A82E.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/2E885270-16EA-DE11-BDCD-001D0967CFA9.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/2C0AF991-19EA-DE11-A1F9-0015178C6704.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1E5AD126-0BEA-DE11-9D1D-0015178C48E4.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1CA7FB7D-1BEA-DE11-8EA3-00151796D4B0.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1A55661D-0BEA-DE11-985D-0024E8769B60.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/1A14841B-1CEA-DE11-AEEA-0024E8768265.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/10810A7E-16EA-DE11-9FB6-00151796C178.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/0E84AA51-1AEA-DE11-B501-0015179EDF00.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/0AA60F31-17EA-DE11-A6C9-001D0967DA3A.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/08AC2C30-0BEA-DE11-87AE-001D0967DB7A.root',
#        '/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0100/00D4856F-16EA-DE11-A4C8-00151796C0F0.root'
#    ),
   inputCommands = cms.untracked.vstring("keep *",
                                         "drop *_gctDigis_*_*")
#                                        "drop *_gtDigis_*_*",	#by Qiao
#					 "drop *_*_*_HLT")
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('211797:22-211797:821')

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('123596:1-123596:143')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('123592:1-123734:max',
#   								    '124009:1-124009:max',
#								    '124020:1-124020:max',
#								    '124022:1-124024:max',
#								    '124027:1-124017:max',
#								    '124030:1-124030:max'
#								    )

#process.GlobalTag.globaltag = 'STARTUP3X_V8D::All'
#process.GlobalTag.globaltag = 'GR_R_53_V21::All'
#process.GlobalTag.globaltag = 'START53_V22::All'
process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('drop *',
					   'keep recoBeamSpot_*_*_*',
					   'keep *_generalTracks_*_*',
					   'keep recoVertexs_*_*_*',
					   'keep recoVertexCompositeCandidates_*_*_*',
                                           'keep patCompositeCandidates_v0CollectionMaker_*_*',
					   'keep *_towerMaker_*_*',	#by Qiao
					   'keep *_gtDigis_*_*'),	#by Qiao

	SelectEvents = cms.untracked.PSet(
					SelectEvents = cms.vstring('filter_step')
),

    fileName = cms.untracked.string('/scratch/xuq7/strange_particles/v0Collection.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    )
)

process.reco_step = cms.Path(
                process.hltSingleTrigger *
                process.PAcollisionEventSelection *
                process.pileupVertexFilterCutGplus *
		process.localV0Candidates *
		process.v0CollectionMaker 
)

process.filter_step = cms.Path(
	process.hltSingleTrigger *
	process.PAcollisionEventSelection *
	process.pileupVertexFilterCutGplus
)

#process.vee_step = cms.Path(process.localV0Candidates)
#process.collect_step = cms.Path(process.v0CollectionMaker)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

#Schedule definition
process.schedule = cms.Schedule(
process.filter_step,
process.reco_step,
process.endjob_step,
process.out_step
)

#process.schedule = cms.Schedule(process.vee_step, process.collect_step, process.endjob_step, process.out_step)
