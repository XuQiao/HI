import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Services_cff import *
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *

from myProducers.V0CandProducer.v0CandProducer_cfi import *
from myAnalyzers.V0RecoAnalyzer.v0RecoAnalyzer_cfi import *

v0analysis = cms.Sequence(localV0Candidates*produceV0PATCands*analyzeKshort*analyzeLambda)

