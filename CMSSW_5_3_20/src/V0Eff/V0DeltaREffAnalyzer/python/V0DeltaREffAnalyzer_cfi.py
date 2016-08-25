import FWCore.ParameterSet.Config as cms

matchVeesByDeltaR = cms.EDAnalyzer('V0DeltaREffAnalyzer',
#    ksCollection = cms.InputTag('generalV0Candidates:Kshort'),
#    lamCollection = cms.InputTag('generalV0Candidates:Lambda'),
    hltTag = cms.InputTag('TriggerResults::HLT'),
    v0Collection = cms.InputTag('produceV0PATCands:allVees'),
    pideltaRcut = cms.double(0.1),
    pdeltaRcut = cms.double(0.1),
    deltaLcut = cms.double(10),
    isKshort = cms.bool(True)
)

matchLamsByDeltaR = cms.EDAnalyzer('V0DeltaREffAnalyzer',
    hltTag = cms.InputTag('TriggerResults::HLT'),
    v0Collection = cms.InputTag('produceV0PATCands:allVees'),
    pideltaRcut = cms.double(0.1),
    pdeltaRcut = cms.double(0.1),
    deltaLcut = cms.double(10),
    isKshort = cms.bool(False)
)
