import FWCore.ParameterSet.Config as cms

process = cms.Process("V0SKIM")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(10000)
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.GlobalTag.globaltag = 'STARTUP3X_V8D::All'

process.veeFilter = cms.EDFilter("V0CandFilter",
#    kShortCollection = cms.InputTag('generalV0Candidates:Kshort'),
    kShortCollection = cms.InputTag('localV0Candidates:Kshort'),
#    lambdaCollection = cms.InputTag('generalV0Candidates:Lambda'),
    lambdaCollection = cms.InputTag('localV0Candidates:Lambda'),
    trackCollection = cms.InputTag('generalTracks')
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:v0Collection123596_gp_goodpixel_FEVT.root'
    )
)

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('123151:3-123151:23')

process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('events_filtered.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('filter_step')
    )
)

process.filter_step = cms.Path(process.veeFilter)
process.endjob_step = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.filter_step, process.endjob_step, process.out_step)

