import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("V0ANALYSIS")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/EndOfProcess_cff')

process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')#by Qiao
process.load('Configuration/Geometry/GeometryIdeal_cff')#by Qiao
process.load('Configuration/StandardSequences/MagneticField_cff')#by Qiao

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')#by Qiao

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(10000)
    input = cms.untracked.int32(-1)
)

#process.options = cms.untracked.PSet(
#    Rethrow = cms.untracked.vstring('ProductNotFound')
#)

#process.GlobalTag.globaltag = 'MC_31X_V3::All'
process.GlobalTag.globaltag = 'START53_V22::All'	#by Qiao
#process.GlobalTag.globaltag = 'GR_R_53_V21::All'	#by Qiao

process.ksAnalyzer = cms.EDAnalyzer("V0RecoAnalyzer",
#    v0Collection = cms.InputTag('generalV0Candidates:Kshort'),
    v0Collection = cms.InputTag('v0CollectionMaker:allVees'),
    pvAvailable = cms.bool(True),
    beamSpotAvailable = cms.bool(True),
    tracksAvailable = cms.bool(True),
    instanceName = cms.string('Kshort'),

hltTag = cms.InputTag('TriggerResults::HLT'),    #by Qiao
pvCollection = cms.InputTag('offlinePrimaryVertices'),	#by Qiao

writeTree = cms.bool(True),
    writeHistos = cms.bool(True),
    ksMassXmin = cms.double(0.305),
    ksMassXmax = cms.double(0.795),
    ksMassNbins = cms.int32(70),
    ksMass_eta_nBins = cms.int32(5),
    ksMass_phi_nBins = cms.int32(6),
#    ksMass_pt_nBins  = cms.int32(5),
    ksMass_pt_nBinSubdiv = cms.int32(1),
    lamMassXmin = cms.double(1.07),
    lamMassXmax = cms.double(1.16),
    lamMassNbins = cms.int32(60),
    lamMass_eta_nBins = cms.int32(5),
    lamMass_phi_nBins = cms.int32(6),
#    lamMass_pt_nBins = cms.int32(5),
    lamMass_pt_nBinSubdiv = cms.int32(1),

    vDau3dIpCut = cms.double(3.),
    v3dIpCut = cms.double(3.),
    otherCandMassSigmaCut = cms.double(2.5)


)




process.lamAnalyzer = cms.EDAnalyzer("V0RecoAnalyzer",
#    v0Collection = cms.InputTag('generalV0Candidates:Lambda'),
    v0Collection = cms.InputTag('v0CollectionMaker:allVees'),
    instanceName = cms.string('Lambda'),

hltTag = cms.InputTag('TriggerResults::HLT'),	#by Qiao
pvCollection = cms.InputTag('offlinePrimaryVertices'),	#by Qiao

    pvAvailable = cms.bool(True),
    beamSpotAvailable = cms.bool(True),
    tracksAvailable = cms.bool(True),
    writeTree = cms.bool(True),
    writeHistos = cms.bool(True),
    ksMassXmin = cms.double(0.44), 
    ksMassXmax = cms.double(0.56),
    ksMassNbins = cms.int32(60),
    ksMass_eta_nBins = cms.int32(5),
    ksMass_phi_nBins = cms.int32(6),
#    ksMass_pt_nBins  = cms.int32(5),
    ksMass_pt_nBinSubdiv = cms.int32(1),
    lamMassXmin = cms.double(1.08),
    lamMassXmax = cms.double(1.16),
    lamMassNbins = cms.int32(40),
    lamMass_eta_nBins = cms.int32(5),
    lamMass_phi_nBins = cms.int32(6),
#    lamMass_pt_nBins = cms.int32(5),
    lamMass_pt_nBinSubdiv = cms.int32(1),

    vDau3dIpCut = cms.double(3.),
    v3dIpCut = cms.double(3.),
    otherCandMassSigmaCut = cms.double(2.5)

)



process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )
mylist = FileUtils.loadListFromFile ('filelist.txt')
for fname in mylist:
 process.source.fileNames.append('file:%s' % (fname))

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#	 'file:/scratch/xuq7/strange_particles/v0Collection.root',
#        'file:v0Collection123596_gp_goodpixel_christophe_FEVT.root'
#        'file:v0CollectionSelfContained.root'
#        'file:v0CollectionSelfContained_2.root',
#        'file:v0CollectionSelfContained_3.root',
#        'file:v0CollectionSelfContained_4.root'
#    )
#)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('/scratch/xuq7/strange_particles/v0analysis.root')
)

process.analysis = cms.Path(
			process.ksAnalyzer*
			process.lamAnalyzer)

process.schedule = cms.Schedule(process.analysis)

