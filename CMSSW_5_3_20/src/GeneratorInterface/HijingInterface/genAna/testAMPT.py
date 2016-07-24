import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.Services_cff")
process.load("GeneratorInterface.HydjetInterface.hydjetDefault_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200)
                                       )
#process.MessageLogger.cerr.FwkReport.reportEvery = 200

#process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
#                                        ignoreTotal=cms.untracked.int32(0),
#                                        oncePerEventMode = cms.untracked.bool(False)
#                                        )

process.source = cms.Source('PoolSource',
                            fileNames = cms.untracked.vstring(),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                            )
mylist = FileUtils.loadListFromFile ('filelist.txt')
for fname in mylist:
 process.source.fileNames.append('file:%s' % (fname))
process.TFileService = cms.Service("TFileService",
                                  fileName = cms.string("pPbAna_jb0.root"))


process.test = cms.EDAnalyzer('AMPTAnalyzer'
                             )

process.p = cms.Path(process.test)




