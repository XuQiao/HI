import FWCore.ParameterSet.Config as cms

from RecoTracker.FinalTrackSelectors.selectHighPurity_cfi import selectHighPurity

highPurityTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('quality("highPurity")')
)

#highPurityTracks = selectHighPurity.clone( copyTrajectories = cms.untracked.bool(True) )

selectTracks = selectHighPurity.clone(
    src = "highPurityTracks",
    keepAllTracks = False,
    # Boolean indicating if adapted primary vertex compatibility cuts are to be applied.
    applyAdaptedPVCuts = cms.bool(True),    
    vertices = "offlinePrimaryVertices",
    vtxNumber = cms.int32(1),  
    vertexCut = "ndof>=2 && abs(z)<=15.0 && abs(position.rho)<=2.0",    
    copyTrajectories = False,             # needed by TrackClusterRemover

    useVertices = True,
    useVtxError = True,
    
    qualityBit = '',                     # set to '' if you don't want to set the bit
    
    min_nhits = cms.uint32(0),
    max_relpterr = 0.1,
    chi2n_par = 9999.9,                  # version with 1D hits modification
    chi2n_no1Dmod_par = 9999.9,          # normalizedChi2 < nLayers * chi2n_par
    d0_par2 = [3.0, 0.0],                # d0E from tk.d0Error
    dz_par2 = [3.0, 0.0],
    
    d0_par1 = [9999., 0.],                 # parameterized nomd0E
    dz_par1 = [9999., 0.],
    res_par = [99999., 99999.],            # residual parameterization (re-check in HI)
    nSigmaZ = 9999.,

    max_z0 = cms.double(1000),
    max_d0 = cms.double(1000),

    minNumberLayers = 0,
    minNumber3DLayers = 0,
    maxNumberLostLayers = 999
    )

nTrkOffTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("selectTracks"),
    cut = cms.string('(pt>0.4) && (abs(eta)<2.4)')
)

nTrkOffFilter120 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(120)
)

nTrkOffFilter150 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(150)
)

nTrkOffFilter185 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(185)
)

nTrkOffFilter220 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(220)
)

nTrkOffFilter260 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(260)
)

nTrkOffFilter300 = cms.EDFilter("TrackCountFilter",
    src = cms.InputTag('nTrkOffTracks'),
    minNumber = cms.uint32(300)
)

ppSingleTrackFilterSequence = cms.Sequence(
                                           highPurityTracks
                                           * selectTracks 
                                           )
