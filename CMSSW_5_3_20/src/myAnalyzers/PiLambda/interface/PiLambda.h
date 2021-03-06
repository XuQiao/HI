// -*- C++ -*-
//
// Package:    PiLambda
// Class:      PiLambda
// 
/**\class PiLambda PiLambda.cc myAnalyzers/PiLambda/src/PiLambda.cc

 Description: <one line class summary>
Make rootTuple for b->s mu mu reconstruction

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer
//         Created:  Mon Apr 21 09:53:19 MDT 2008
// $Id: PiLambda.h,v 1.23 2010/11/16 01:12:06 drell Exp $
//
//

#ifndef _PiLambda_h
#define _PiLambda_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "RecoVertex/V0Producer/interface/V0Producer.h"

#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/Math/interface/angle.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

// For algo triggers
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
// end algo triggers

#include "TFile.h"
#include "TTree.h"

using namespace std;	//by Qiao

//
// class decleration
//

class PiLambda : public edm::EDAnalyzer {
public:
  explicit PiLambda(const edm::ParameterSet&);
  ~PiLambda();
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;
  
  // ----------member data ---------------------------

  InvariantMassFromVertex massCalculator;

  edm::InputTag hltTag;
  std::string hepMC;
  std::string v0Producer;
  std::string v0Type;
  std::string muonType;
  std::string vtxSample;
  std::string v0Collection;
  edm::InputTag trackCollectionName;
  bool doMC;
  bool debug;
  TTree*      tree_;
  int nEvents;
  bool trackInV0(const reco::TrackRef, std::vector<reco::VertexCompositeCandidate>);
  bool aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers);
  bool scrapingTest( const edm::Event&, const edm::EventSetup& );

  unsigned int        nXI, run, event;
  int                 trigHF;
  float               HFpE, HFmE;
  int                 trigTech40_41, trigTech36_39, trigTech34, trigScraping, trigAlgo124, trigHLTminBias;
  int                 nLooseTracks, nHighPurTracks, nLooseTracks_pt_100, nHighPurTracks_pt_100;
  int                 processType;
  float               priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  uint                priNTrk, priNTrk_weightGT_0_5;
  float               bsX, bsY, bsZ, bsXE, bsYE, bsZE;
  vector<float>       *priRefitVtxX, *priRefitVtxY, *priRefitVtxZ, *priRefitVtxXE, *priRefitVtxYE, *priRefitVtxZE, *priRefitVtxCL;
  vector<uint>        *priRefitNTrk; 
  vector<float>       *xiMass, *omMass, *xiVtxCL, *xiPx, *xiPy, *xiPz;
  vector<double>      *xiPxE, *xiPyE, *xiPzE;
  vector<float>       *xictauB, *xictauBE, *xictau3D, *xictau3DE, *xictauMPV, *xictauMPVE;
  vector<float>       *xiDecayVtxX, *xiDecayVtxY, *xiDecayVtxZ;
  vector<double>      *xiDecayVtxXE, *xiDecayVtxYE, *xiDecayVtxZE;
  vector<float>       *VMass, *VCandMass, *VMassError, *VVtxCL, *VPx, *VPy, *VPz, *VPalongXi;
  vector<float>       *VKsMass, *VKsMassError;
  vector<float>       *VDecayVtxX, *VDecayVtxY, *VDecayVtxZ;
  vector<float>       *VDecayVtxXE, *VDecayVtxYE, *VDecayVtxZE;
  vector<float>       *batPiPx, *batPiPy, *batPiPz, *batPiEta, *batPiPhi, *batPiD0, *batPiD0E, *batPiDz, *batPiDzE, *batPiPVweight;
  vector<int>         *batPiQ, *batPiNValidHits, *batPiTrkQual, *nTracks;
  vector<float>       *VTrkPMass, *VTrkPPx, *VTrkPPy, *VTrkPPz, *VTrkPEta, *VTrkPPhi, *VTrkPD0, *VTrkPD0E, *VTrkPDz, *VTrkPDzE, *VTrkPPVweight;
  vector<int>         *VTrkPQ, *VTrkPNValidHits;
  vector<float>       *VTrkPiMass, *VTrkPiPx, *VTrkPiPy, *VTrkPiPz, *VTrkPiEta, *VTrkPiPhi, *VTrkPiD0, *VTrkPiD0E, *VTrkPiDz, *VTrkPiDzE, *VTrkPiPVweight;
  vector<int>         *VTrkPiQ, *VTrkPiNValidHits;
  vector<float>	      *VTransversePCAPrimary, *VTransversePCAPrimaryError, *VLongitudinalPCAPrimary, *VLongitudinalPCAPrimaryError; 
  vector<float>	      *VTransversePCABeamSpot, *VTransversePCABeamSpotError, *VRSig2D, *VRSig3D, *VFLSig2D, *VFLSig3D;
  vector<float>	      *XiTransversePCAPrimary, *XiTransversePCAPrimaryError, *XiLongitudinalPCAPrimary, *XiLongitudinalPCAPrimaryError; 
  vector<float>	      *XiTransversePCABeamSpot, *XiTransversePCABeamSpotError; 
  vector<float>	      *Xi3DIpSig, *XiFLsig3D;
  vector<float>	      *VCosThetaPAndveeVertexToPrimaryVector, *VCosThetaPAndveeVertexToXiVertexVector, *XiCosThetaPAndLineOfFlight;
  vector<double>      *V3dIpWrtPrimary, *V3dIpWrtPrimaryError, *V3dIpWrtPrimarySig;
  vector<double>      *VTrkPi2DIp, *VTrkPi2DIpSig, *VTrkPi3DIp, *VTrkPi3DIpSig;
  vector<double>      *VTrkP2DIp, *VTrkP2DIpSig, *VTrkP3DIp, *VTrkP3DIpSig;
  vector<double>      *batPi2DIp, *batPi2DIpSig, *batPi3DIp, *batPi3DIpSig;
  vector<bool>	      *pionInV0;
 
  unsigned int        nGenXi;
  vector<float>       *genBatPiEta, *genBatPiPhi, *genLamPiEta, *genLamPiPhi, *genLamPEta, *genLamPPhi;
  vector<float>       *genXiP, *genXiEta, *genXiL, *genXiR, *genXiDecayVX, *genXiDecayVY, *genXiDecayVZ;
  vector<float>       *genXiMomL;
  vector<int>         *genXiQ;
  vector<int>         *genXiMotherPDG;
  vector<float>       *genXiProdVtxX, *genXiProdVtxY, *genXiProdVtxZ;
   
};

#endif
