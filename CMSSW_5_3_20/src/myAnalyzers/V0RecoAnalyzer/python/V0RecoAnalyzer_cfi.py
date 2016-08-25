import FWCore.ParameterSet.Config as cms

v0RecoAnalyzer = cms.EDAnalyzer("V0RecoAnalyzer",
    v0Collection = cms.InputTag('generalV0Candidates Kshort')
)
