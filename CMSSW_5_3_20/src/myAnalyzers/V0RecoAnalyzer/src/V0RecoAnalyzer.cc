// -*- C++ -*-
//
// Package:    V0RecoAnalyzer
// Class:      V0RecoAnalyzer
// 
/**\class V0RecoAnalyzer V0RecoAnalyzer.cc myAnalyzers/V0RecoAnalyzer/src/V0RecoAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Wed Aug 26 13:34:28 MDT 2009
// $Id: V0RecoAnalyzer.cc,v 1.47 2011/05/23 16:23:05 drell Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"	//by Qiao
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

//#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


#include "DataFormats/Math/interface/angle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

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

#include <typeinfo>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TCanvas.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

// Standard C++ includes
#include <sstream>

// Constants
const double piMass = 0.13957018;
const double piMass2 = piMass * piMass;
const double protonMass = 0.93827203;
const double protonMass2 = protonMass * protonMass;
const double ksMassConst = 0.497648;
const double ksMassConst2 = ksMassConst * ksMassConst;
const double lamMassConst = 1.115683;
const double lamMassConst2 = lamMassConst*lamMassConst;

class V0RecoAnalyzer : public edm::EDAnalyzer {
public:
  explicit V0RecoAnalyzer(const edm::ParameterSet&);
  ~V0RecoAnalyzer();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers);
  bool scrapingTest( const edm::Event&, const edm::EventSetup& );

  edm::InputTag hltTag;

  TTree* theTree;
  TFile* theFile;

  // ***Tree variables***
  // V0 vertex position info
  std::vector<float> *v0VtxX;
  std::vector<float> *v0VtxY;
  std::vector<float> *v0VtxZ;
  std::vector<float> *v0VtxR;
  std::vector<float> *v0VtxChi2;
  std::vector<float> *v0VtxNormChi2;
  /*  std::vector<float> *v0VtxCXX;
  std::vector<float> *v0VtxCYY;
  std::vector<float> *v0VtxCZZ;
  std::vector<float> *v0VtxCYX;
  std::vector<float> *v0VtxCZX;
  std::vector<float> *v0VtxCZY;
  std::vector<float> *v0VtxChi2;
  std::vector<float> *v0VtxNormChi2;*/

  // V0 kinematics
  std::vector<float> *v0pX;
  std::vector<float> *v0pY;
  std::vector<float> *v0pZ;
  std::vector<float> *v0pT;
  std::vector<float> *v0p;
  std::vector<float> *v0Eta;
  std::vector<float> *v0Rapidity;
  std::vector<float> *v0Phi;

  // For the Armenteros plot:
  std::vector<float> *v0pTPlus;
  std::vector<float> *v0pTMinus;
  std::vector<float> *v0pLPlus;
  std::vector<float> *v0pLMinus;
  /*std::vector<float> *v0pXerr;
  std::vector<float> *v0pYerr;
  std::vector<float> *v0pZerr;*/

  // V0 mass and error
  //std::vector<float> *v0Mass;
  //std::vector<float> *v0Masserr;
  //std::vector<float> *v0Massnorm;
  std::vector<float> *v0CandMass;
  std::vector<float> *v0OtherCandMass;
  std::vector<float> *v0CandMassError;
  std::vector<float> *v0OtherCandMassError;
  //std::vector<float> *v0MassFromFit;
  //std::vector<float> *v0MassErrorFromFit;
  //std::vector<float> *v0OtherMassFromFit;
  //std::vector<float> *v0OtherMassErrorFromFit;

  // V0 lifetime
  std::vector<float> *v0Lifetime;
  std::vector<float> *v03DLifetime;
  std::vector<float> *v03DLifetimeError;
  std::vector<float> *v03DLifetimeScaled;
  std::vector<float> *v03DLifetimeScaledError;
  std::vector<float> *v03DLifetimeMPV;
  std::vector<float> *v03DLifetimeMPVError;

  std::vector<int> *v0PDG;
  std::vector<float> *v0VtxSig;
  std::vector<float> *v0VtxSig3D;
  
  // VO point back criteria
  std::vector<float> *v0CosThetaPAndLineOfFlight;
  std::vector<float> *v0PCAPrimary;
  std::vector<float> *v0PCAPrimaryError;
  std::vector<float> *v0PCABeamSpot;
  std::vector<float> *v0PCABeamSpotError;
  std::vector<float> *v3dIpWrtPrimary;
  std::vector<float> *v3dIpWrtPrimaryError;
  std::vector<float> *v3dIpWrtPrimarySig;
  std::vector<float> *thePosDau3DIp;
  std::vector<float> *thePosDau3DIpError;
  std::vector<float> *thePosDau3DIpSig;
  std::vector<float> *theNegDau3DIp;
  std::vector<float> *theNegDau3DIpError;
  std::vector<float> *theNegDau3DIpSig;
  std::vector<int> *thePosDauTkQual;
  std::vector<int> *theNegDauTkQual;
  std::vector<int> *thePosDauTkAlgo;
  std::vector<int> *theNegDauTkAlgo;

  std::vector<int> *posDauNhits;
  std::vector<int> *negDauNhits;

  //std::vector<float> *priVtxChi2;
  //std::vector<float> *priVtxNdof;
  //std::vector<float> *nTracks;
  //std::vector<float> *nV0s;
//   int nLooseTracks;
//   int nHighPurTracks;
//   int nLooseTracks_pt_100;
//   int nHighPurTracks_pt_100;

  int nTracks;
  int nTracks_pt_100;
  int nPromptTracks;
  int nPromptTracks_pt_100;

  int nGenParts_eta_24_pt_05;
  int nGenParts_eta_08_pt_05;
  int nGenParts_eta_24_pt_01;

  int nGenTracks_hp_eta24_pt05_nHitLay6;
  int nGenTracks_hp_eta08_pt05_nHitLay6;
  int nGenTracks_hp_eta24_pt05_nHitLay6_BS;
  int nGenTracks_eta24_pt05;
  int nGenTracks_eta24_pt05_BS;
  int nTracks_eta24_pt01;
  int nTracks_eta24_pt01_BS;

  int nV0s;
  int trigHF;
  float HFpE;
  float HFmE;
  int trigTech40_41;
  int trigTech34;
  int trigTech36_39;
  int trigScraping;
  int trigAlgo124;
  int trigHLTminBias;

  int runNum;
  int evtNum;

  int processType;

  float priVtxChi2;
  float priVtxNdof;
  float priVtxX;
  float priVtxY;
  float priVtxZ;
  int priVtxNtracks;
  int priVtxNtracksRaw;
  int priVtxIsFake;
  int nPriVtxs;
  int nTracksFromOtherPrimary;
  std::vector<int> *priVtxWasRefit;

  std::vector<float> *priRefitVtxX;
  std::vector<float> *priRefitVtxY;
  std::vector<float> *priRefitVtxZ;
  std::vector<float> *priRefitVtxChi2;
  std::vector<float> *priRefitVtxNdof;

  std::vector<float> *v0VtxPriDist;
  std::vector<float> *v0VtxPriDistErr;

  std::vector<float> *v0VtxBSTransDist;
  std::vector<float> *v0VtxBSTransDistErr;

  // V0 daughter track quantities
  std::vector<float> *v0PosDaupT;
  std::vector<float> *v0PosDaup;
  std::vector<float> *v0PosDaupX;
  std::vector<float> *v0PosDaupY;
  std::vector<float> *v0PosDaupZ;
  std::vector<float> *v0PosDauAtVtxpX;
  std::vector<float> *v0PosDauAtVtxpY;
  std::vector<float> *v0PosDauAtVtxpZ;
  std::vector<float> *v0PosDauEta;
  std::vector<float> *v0PosDauPhi;
  //std::vector<float> *v0PosDauNHits;
  std::vector<float> *v0PosDauNValidHits;
  std::vector<float> *v0PosDauNPixelHits;
  std::vector<float> *v0PosDauNStripHits;
  std::vector<float> *v0PosDauChi2;
  std::vector<float> *v0PosDauNormChi2;
  std::vector<float> *v0PosDauD0;

  // V0 daughter track quantities
  std::vector<float> *v0NegDaupT;
  std::vector<float> *v0NegDaup;
  std::vector<float> *v0NegDaupX;
  std::vector<float> *v0NegDaupY;
  std::vector<float> *v0NegDaupZ;
  std::vector<float> *v0NegDauAtVtxpX;
  std::vector<float> *v0NegDauAtVtxpY;
  std::vector<float> *v0NegDauAtVtxpZ;
  std::vector<float> *v0NegDauEta;
  std::vector<float> *v0NegDauPhi;
  //std::vector<float> *v0NegDauNHits;
  std::vector<float> *v0NegDauNValidHits;
  std::vector<float> *v0NegDauNPixelHits;
  std::vector<float> *v0NegDauNStripHits;
  std::vector<float> *v0NegDauChi2;
  std::vector<float> *v0NegDauNormChi2;
  std::vector<float> *v0NegDauD0;

  std::vector<float> *v0DauSumPx;
  std::vector<float> *v0DauSumPy;
  std::vector<float> *v0DauSumPz;
  
  // Quantities associated with daughter track that has lowest pT
  std::vector<float> *v0DaupT_lowestpT;
  std::vector<float> *v0DauEta_lowestpT;
  std::vector<float> *v0DauPhi_lowestpT;

  // Histograms
  TH1F* ksR;
  TH1F* ksRFull;
  TH1F* ksRSignal;
  TH1F* ksRBg;

  TH1F* ksZ;
  TH1F* ksZFull;
  TH1F* ksZSignal;
  TH1F* ksZBg;

  TH2F* ksZvsRFull;
  TH2F* ksZvsRBgSubt;
  TH2F* ksZvsRBg;
  TH2F* ksZvsRSignal;

  TH1F* ksEta;
  TH1F* ksEtaFull;
  TH1F* ksEtaSignal;
  TH1F* ksEtaBg;

  TH1F* ksPhi;
  TH1F* ksPhiFull;
  TH1F* ksPhiSignal;
  TH1F* ksPhiBg;

  TH1F* ksPt;
  TH1F* ksPtFull;
  TH1F* ksPtSignal;
  TH1F* ksPtBg;

  TH1F* ksMass;
  TH1F* ksNormMass;
  TH1F* ksCandMass;

  TH1F* ksDaup;
  TH1F* ksDaupT;
  TH1F* ksDauEta;
  TH1F* ksDauPhi;
  TH1F* ksDauNValidHits;
  TH1F* ksDauNPixelHits;
  TH1F* ksDauNStripHits;
  TH1F* ksDauChi2;
  TH1F* ksDauNormChi2;
  TH1F* ksDauD0;

  TH1F* lamR;
  TH1F* lamRFull;
  TH1F* lamRSignal;
  TH1F* lamRBg;

  TH1F* lamZ;
  TH1F* lamZFull;
  TH1F* lamZSignal;
  TH1F* lamZBg;

  TH2F* lamZvsRFull;
  TH2F* lamZvsRBgSubt;
  TH2F* lamZvsRBg;
  TH2F* lamZvsRSignal;

  TH1F* lamEta;
  TH1F* lamEtaFull;
  TH1F* lamEtaSignal;
  TH1F* lamEtaBg;

  TH1F* lamPhi;
  TH1F* lamPhiFull;
  TH1F* lamPhiSignal;
  TH1F* lamPhiBg;

  TH1F* lamPt;
  TH1F* lamPtFull;
  TH1F* lamPtSignal;
  TH1F* lamPtBg;

  TH1F* lamCandMass;

  TH1F* lamDaup;
  TH1F* lamDaupT;
  TH1F* lamDauEta;
  TH1F* lamDauPhi;
  TH1F* lamDauNValidHits;
  TH1F* lamDauNPixelHits;
  TH1F* lamDauNStripHits;
  TH1F* lamDauChi2;
  TH1F* lamDauNormChi2;
  TH1F* lamDauD0;

  /*TH1F* ksMass_pt_1;
  TH1F* ksMass_pt_2;
  TH1F* ksMass_pt_3;
  TH1F* ksMass_pt_4;
  TH1F* ksMass_pt_5;
  TH1F* ksMass_pt_6;*/
  std::vector<TH1F*> ksMass_pt; 
  std::vector<std::string> ksMass_pt_names_short;
  std::vector<std::string> ksMass_pt_names_long;
  std::vector<double> ksMass_pt_minima;

  /*TH1F* ksMass_phi_1;
  TH1F* ksMass_phi_2;
  TH1F* ksMass_phi_3;
  TH1F* ksMass_phi_4;
  TH1F* ksMass_phi_5;
  TH1F* ksMass_phi_6;*/
  std::vector<TH1F*> ksMass_phi;
  std::vector<std::string> ksMass_phi_names_short;
  std::vector<std::string> ksMass_phi_names_long;
  std::vector<double> ksMass_phi_minima;

  /*TH1F* ksMass_eta_2__2_5;
  TH1F* ksMass_eta_1_5__2;
  TH1F* ksMass_eta_1__1_5;
  TH1F* ksMass_eta_0_5__1;
  TH1F* ksMass_eta_0__0_5;*/
  std::vector<TH1F*> ksMass_eta;
  std::vector<std::string> ksMass_eta_names_short;
  std::vector<std::string> ksMass_eta_names_long;
  std::vector<double> ksMass_eta_minima;

  std::vector<TH1F*> ksMass_eta_phi;
  std::vector<std::string> ksMass_eta_phi_names_short;
  std::vector<std::string> ksMass_eta_phi_names_long;
  std::vector<double> ksMass_eta_phi_etaminima;
  std::vector<double> ksMass_eta_phi_phiminima;

  std::vector<TH1F*> lamMass_pt;
  std::vector<std::string> lamMass_pt_names_short;
  std::vector<std::string> lamMass_pt_names_long;
  std::vector<double> lamMass_pt_minima;

  std::vector<TH1F*> lamMass_phi;
  std::vector<std::string> lamMass_phi_names_short;
  std::vector<std::string> lamMass_phi_names_long;
  std::vector<double> lamMass_phi_minima;

  std::vector<TH1F*> lamMass_eta;
  std::vector<std::string> lamMass_eta_names_short;
  std::vector<std::string> lamMass_eta_names_long;
  std::vector<double> lamMass_eta_minima;

  /*TGraphErrors* ksMassInPtBins;
  TGraphErrors* trueKsMassInPtBins;
  TGraphErrors* ksMassInPhiBins;  for( pat::CompositeCandidateCollection::const_iterator iV0 = patHandle->begin();

  TGraphErrors* trueKsMassInPhiBins;
  TGraphErrors* ksMassInEtaBins;
  TGraphErrors* trueKsMassInEtaBins;*/
  TH1F* ksMassInPtBins;
  TH1F* ksMassInPhiBins;
  TH1F* ksMassInEtaBins;
  TH2F* ksMassBiasEtaVsPhi;

  TH1F* lamMassInPtBins;
  TH1F* lamMassInPhiBins;
  TH1F* lamMassInEtaBins;

  // Strings for algorithm names, etc.
  edm::InputTag v0Collection;

  bool writeTree;
  bool putTreeInTFile;
  bool writeHistos;
  std::string instanceName;

  bool pvAvailable;
  bool beamSpotAvailable;
  bool tracksAvailable;

  edm::InputTag pvColl;

  // Mass plot limits and binning
  double ksMassXmin, ksMassXmax;
  int ksMassNbins;
  double ksMassBinWidth;
  //double ksNormMassXmin, ksNormMassXmax;
  //double ksNormMassBinWidth;

  double lamMassXmin, lamMassXmax;
  int lamMassNbins;
  double lamMassBinWidth;
  //double lamNormMassXmin, lamNormMassXmax;
  //double lamNormMassBinWidth;

  // Parameters for mass bias plots
  int ksMass_eta_nBins;
  int ksMass_phi_nBins;
  //int ksMass_pt_nBins;
  double ksMass_pt_nBinSubdiv;

  int lamMass_eta_nBins;
  int lamMass_phi_nBins;
  //int lamMass_pt_nBins;
  double lamMass_pt_nBinSubdiv;



};


V0RecoAnalyzer::V0RecoAnalyzer(const edm::ParameterSet& iConfig) : 
  hltTag( iConfig.getParameter<edm::InputTag>("hltTag") ),
  v0VtxX(0), v0VtxY(0), v0VtxZ(0), v0VtxR(0), v0VtxChi2(0), v0VtxNormChi2(0),
  /*v0VtxCXX(0), v0VtxCYY(0), v0VtxCZZ(0), v0VtxCYX(0), v0VtxCZX(0), v0VtxCZY(0),
    v0VtxChi2(0), v0VtxNormChi2(0),*/
  v0pX(0), v0pY(0), v0pZ(0), v0pT(0), v0p(0), v0Eta(0), v0Rapidity(0), v0Phi(0), //v0pXerr(0), v0pYerr(0), v0pZerr(0),
  v0pTPlus(0), v0pTMinus(0), v0pLPlus(0), v0pLMinus(0),
  v0CandMass(0), v0OtherCandMass(0), v0CandMassError(0), v0OtherCandMassError(0),
  //v0MassFromFit(0), v0MassErrorFromFit(0),
  //v0OtherMassFromFit(0), v0OtherMassErrorFromFit(0),
  v0Lifetime(0), 
  v03DLifetime(0), v03DLifetimeError(0), v03DLifetimeScaled(0), v03DLifetimeScaledError(0),
  v03DLifetimeMPV(0), v03DLifetimeMPVError(0),
  v0PDG(0), v0VtxSig(0), v0VtxSig3D(0),
  v0CosThetaPAndLineOfFlight(0), v0PCAPrimary(0) ,v0PCAPrimaryError(0),
  v0PCABeamSpot(0),v0PCABeamSpotError(0),
  v3dIpWrtPrimary(0),v3dIpWrtPrimaryError(0),v3dIpWrtPrimarySig(0),
  thePosDau3DIp(0), thePosDau3DIpError(0), thePosDau3DIpSig(0), 
  theNegDau3DIp(0),theNegDau3DIpError(0),theNegDau3DIpSig(0),
  thePosDauTkQual(0), theNegDauTkQual(0),
  thePosDauTkAlgo(0), theNegDauTkAlgo(0),
  posDauNhits(0), negDauNhits(0),
  //nLooseTracks(0), nHighPurTracks(0), nLooseTracks_pt_100(0), nHighPurTracks_pt_100(0),
  nTracks(0), nTracks_pt_100(0), nPromptTracks(0), nPromptTracks_pt_100(0),
  nGenParts_eta_24_pt_05(0), nGenParts_eta_08_pt_05(0), nGenParts_eta_24_pt_01(0),
  nGenTracks_hp_eta24_pt05_nHitLay6(0), nGenTracks_hp_eta08_pt05_nHitLay6(0), 
  nGenTracks_hp_eta24_pt05_nHitLay6_BS(0), nGenTracks_eta24_pt05(0), 
  nGenTracks_eta24_pt05_BS(0), nTracks_eta24_pt01(0), nTracks_eta24_pt01_BS(0), 
  nV0s(0), trigHF(0), HFpE(0), HFmE(0),
  trigTech40_41(0), trigTech34(0), trigTech36_39(0), trigScraping(0), trigAlgo124(0), trigHLTminBias(0),
  runNum(0), evtNum(0), processType(0),
  priVtxChi2(0), priVtxNdof(0), 
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxNtracks(0), priVtxNtracksRaw(0), priVtxIsFake(0), 
  nPriVtxs(0), nTracksFromOtherPrimary(0), priVtxWasRefit(0),
  priRefitVtxX(0), priRefitVtxY(0), priRefitVtxZ(0), priRefitVtxChi2(0), priRefitVtxNdof(0),
  v0VtxPriDist(0), v0VtxPriDistErr(0),
  v0VtxBSTransDist(0), v0VtxBSTransDistErr(0),
  v0PosDaupT(0), v0PosDaup(0), v0PosDaupX(0), v0PosDaupY(0), v0PosDaupZ(0), 
  v0PosDauAtVtxpX(0), v0PosDauAtVtxpY(0), v0PosDauAtVtxpZ(0),
  v0PosDauEta(0), v0PosDauPhi(0), 
  v0PosDauNValidHits(0), v0PosDauNPixelHits(0), v0PosDauNStripHits(0),
  v0PosDauChi2(0), v0PosDauNormChi2(0), v0PosDauD0(0),
  v0NegDaupT(0), v0NegDaup(0), v0NegDaupX(0), v0NegDaupY(0), v0NegDaupZ(0), 
  v0NegDauAtVtxpX(0), v0NegDauAtVtxpY(0), v0NegDauAtVtxpZ(0),
  v0NegDauEta(0), v0NegDauPhi(0), 
  v0NegDauNValidHits(0), v0NegDauNPixelHits(0), v0NegDauNStripHits(0),
  v0NegDauChi2(0), v0NegDauNormChi2(0), v0NegDauD0(0),
  v0DauSumPx(0), v0DauSumPy(0), v0DauSumPz(0),
  v0DaupT_lowestpT(0), v0DauEta_lowestpT(0), v0DauPhi_lowestpT(0),
  v0Collection( iConfig.getParameter<edm::InputTag>("v0Collection") ),
  putTreeInTFile( iConfig.getParameter<bool>("writeTree") ),
  writeHistos( iConfig.getParameter<bool>("writeHistos") ),
  instanceName( iConfig.getParameter<std::string>("instanceName") ),
  pvAvailable( iConfig.getParameter<bool>("pvAvailable") ),
  beamSpotAvailable( iConfig.getParameter<bool>("beamSpotAvailable") ),
  tracksAvailable( iConfig.getParameter<bool>("tracksAvailable") ),
  pvColl( iConfig.getParameter<edm::InputTag>("pvCollection") ),
  ksMassXmin(iConfig.getParameter<double>("ksMassXmin")),
  ksMassXmax(iConfig.getParameter<double>("ksMassXmax")),
  ksMassNbins(iConfig.getParameter<int>("ksMassNbins")),
  ksMassBinWidth(0.),
  lamMassXmin(iConfig.getParameter<double>("lamMassXmin")),
  lamMassXmax(iConfig.getParameter<double>("lamMassXmax")),
  lamMassNbins(iConfig.getParameter<int>("lamMassNbins")),
  lamMassBinWidth(0.),
  ksMass_eta_nBins(iConfig.getParameter<int>("ksMass_eta_nBins")),
  ksMass_phi_nBins(iConfig.getParameter<int>("ksMass_phi_nBins")),
  //ksMass_pt_nBins(iConfig.getParameter<int>("ksMass_pt_nBins")),
  ksMass_pt_nBinSubdiv(iConfig.getParameter<int>("ksMass_pt_nBinSubdiv")),
  lamMass_eta_nBins(iConfig.getParameter<int>("lamMass_eta_nBins")), 
  lamMass_phi_nBins(iConfig.getParameter<int>("lamMass_phi_nBins")),
  //lamMass_pt_nBins(iConfig.getParameter<int>("lamMass_pt_nBins")),
  lamMass_pt_nBinSubdiv(iConfig.getParameter<int>("lamMass_pt_nBinSubdiv")) {
  //std::cout << "Constructor" << std::endl;
  //instanceName = v0Collection.instance();
  writeTree = true;

}

V0RecoAnalyzer::~V0RecoAnalyzer() {

}


void V0RecoAnalyzer::analyze(const edm::Event& iEvent, 
			     const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<HepMCProduct> hepEv;

  GlobalPoint priGenVtx;
  try {
    iEvent.getByType( hepEv );
    processType = hepEv->GetEvent()->signal_process_id();
    //hepEv->GetEvent()->print();
    for(HepMC::GenEvent::vertex_const_iterator iVTX = hepEv->GetEvent()->vertices_begin();
	iVTX != hepEv->GetEvent()->vertices_end();
	iVTX++) {
      if ((*iVTX)->barcode() == -1) {
	priGenVtx = GlobalPoint((*iVTX)->point3d().x(),
				(*iVTX)->point3d().y(),
				(*iVTX)->point3d().z());
      
 //     cout << "GenVertex, bc = " << (*iVTX)->barcode()
 // 	   << ", id = " << (*iVTX)->id()
 //	   << ", position: (" << (*iVTX)->point3d().x()
 //	   << ", " << (*iVTX)->point3d().y()
 //	   << ", " << (*iVTX)->point3d().z() << ")" << endl;
      }
    }
  }
  catch(...) {
  // cout << "Caught" << endl;
    processType = 0;
  }

  Handle<reco::GenParticleCollection> genPH;
  bool isMC = false;//true;
  try {
    iEvent.getByLabel("genParticles", genPH);//PlusSim", genPH);
  }
  catch(...) {
    isMC = false;
  }

  if (isMC) {
    // Loop over genParticles and count them
    int n_eta_24_pt_05 = 0;
    int n_eta_08_pt_05 = 0;
    int n_eta_24_pt_01 = 0;
    for (unsigned int ndx = 0; ndx < genPH->size(); ndx++) {
      const Candidate &theCand = (*genPH)[ndx];
      int pdg = abs(theCand.pdgId());
      if (theCand.numberOfMothers() == 0) {
	//cout << "Particle with ID " << pdg << " has 0 mothers." << endl;
	priGenVtx = GlobalPoint(theCand.daughter(0)->vertex().x(),
				theCand.daughter(0)->vertex().y(),
				theCand.daughter(0)->vertex().z());
	//cout << "GenVtx: " << priGenVtx << endl;
      }

      if (!(pdg == 211
	    || pdg == 321
	    || pdg == 2212
	    || pdg == 11
	    || pdg == 13)) continue;
      if (theCand.status() != 1) continue;

      GlobalPoint genVtx(theCand.vertex().x(),
			 theCand.vertex().y(),
			 theCand.vertex().z());
      GlobalVector distVec(genVtx - priGenVtx);
      //cout << "d0: " << distVec.perp() << ", |dz|: " << fabs(distVec.z()) << endl;
      double dz = fabs(distVec.z());
      double d0 = fabs(distVec.perp());
      double eta = fabs(theCand.eta());
      double pt = theCand.pt();
      if (eta < 2.4 && pt > 0.5 && dz < 0.2 && d0 < 0.2) n_eta_24_pt_05++;
      if (eta < 0.8 && pt > 0.5 && dz < 0.2 && d0 < 0.2) n_eta_08_pt_05++;
      if (eta < 2.4 && pt > 0.1 && dz < 0.2 && d0 < 0.2) n_eta_24_pt_01++;
    }
    nGenParts_eta_24_pt_05 = n_eta_24_pt_05;
    nGenParts_eta_08_pt_05 = n_eta_08_pt_05;
    nGenParts_eta_24_pt_01 = n_eta_24_pt_01;
  }

  Handle<pat::CompositeCandidateCollection> patHandle;
  Handle<reco::BeamSpot> theBeamSpotH;
  Handle<reco::VertexCollection> priVtxs;
  Handle<reco::TrackCollection> theTracks;

  iEvent.getByLabel( v0Collection, patHandle );
  if( putTreeInTFile ) {
    //cout << "Crashy?" << endl;
    //nV0s->push_back( patHandle->size() );
    nV0s = patHandle->size();
    //cout << "Crash?" << endl;
  }
  if( beamSpotAvailable ) {
    iEvent.getByLabel( std::string("offlineBeamSpot"), theBeamSpotH );
  }
  
  ESHandle<TransientTrackBuilder> theTTB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

  reco::Vertex* thePrim = 0;
  reco::Vertex* theOrigPrim = 0;
  if( putTreeInTFile && pvAvailable ) {
    iEvent.getByLabel( pvColl, priVtxs );
    nPriVtxs = priVtxs->size(); 
    if( priVtxs->size() > 0. ) {
      reco::VertexCollection::const_iterator bestVtx = priVtxs->begin();
      if( bestVtx->isValid() ) {// && !bestVtx->isFake() ) {
	thePrim = new reco::Vertex( *bestVtx );
	theOrigPrim = thePrim;
	//priVtxChi2->push_back( bestVtx->normalizedChi2() );
	//priVtxNdof->push_back( bestVtx->ndof() );
	priVtxChi2 = bestVtx->normalizedChi2();
	priVtxNdof = bestVtx->ndof();
	priVtxX = bestVtx->position().x();
	priVtxY = bestVtx->position().y();
	priVtxZ = bestVtx->position().z();
	if( thePrim->isFake() ) priVtxIsFake = 1;
	else priVtxIsFake = 0;
	int nVtxTracks = 0;
	for( Vertex::trackRef_iterator itk = bestVtx->tracks_begin();
	     itk != bestVtx->tracks_end();
	     itk++ ) {
	  if( bestVtx->trackWeight( *itk ) > 0.5 )
	    nVtxTracks++;
	}
	priVtxNtracks = nVtxTracks;
	priVtxNtracksRaw = bestVtx->tracksSize();
      }
    }
  }
  if (priVtxs->size() > 1) {
    for (reco::VertexCollection::const_iterator iVTX = ++priVtxs->begin();
	 iVTX != priVtxs->end();
	 iVTX++) {
      if ((int) iVTX->tracksSize() > nTracksFromOtherPrimary)
	nTracksFromOtherPrimary = iVTX->tracksSize();
    }
  }
  /*if( putTreeInTFile && tracksAvailable ) {
    iEvent.getByLabel( "generalTracks", theTracks );
    //nTracks->push_back( theTracks->size() );
    nTracks = theTracks->size();
    }*/

  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();

  // Check HLT_L1_BscMinBiasOR_BptxPlusORMinus trigger for prescale lumi calculation
  edm::Handle<edm::TriggerResults> hltResults;
  edm::Handle<edm::TriggerNames> trigNames;
  try {
    iEvent.getByLabel(hltTag, hltResults);//edm::InputTag("TriggerResults::REDIGI36"), hltResults);
    //iEvent.getByLabel(edm::InputTag("TriggerNames::HLT"), trigNames);
  }
  catch( ... ) {
    cout << "Couldn't get HLT trigger." << endl;
    return;
  }
  if( !hltResults.isValid() ) {
 //   cout << "No trigger results." << endl;	//by Qiao
  }
  else {
    int ntrigs = hltResults->size();
    //edm::TriggerNames trigNames(*hltResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);
    //trigNames.init(*hltResults);

    for( int ndx = 0; ndx < ntrigs; ndx++ ) {
      TString trigName = trigNames.triggerName(ndx);
      int hltflag = (*hltResults)[ndx].accept();
      //cout << "Trigger " << trigName << " was passed = " << hltflag << endl;
      if( trigName == "HLT_L1_BscMinBiasOR_BptxPlusORMinus" ) {
	trigHLTminBias = hltflag;
      }
    }
  }

  // check HF MinBias trigger
  edm::Handle<CaloTowerCollection> theCaloTowers;
  iEvent.getByLabel("towerMaker",theCaloTowers);
  if (aboveTriggerThreshold( theCaloTowers )) trigHF = 1;
  else trigHF = 0;
   
  //check tech MinBias trigger bits
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel("gtDigis", gtRecord);
  if (gtRecord.isValid()) {
    TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
    if ( ! tw.empty() ) {
      if (tw[40]||tw[41]) 
	trigTech40_41 = 1;
      else 
	trigTech40_41 = 0;
      if( !(tw[36] || tw[37] || tw[38] || tw[39]) )
	trigTech36_39 = 1;
      else
	trigTech36_39 = 0;
      if( tw[34] )
	trigTech34 = 1;
      else
	trigTech34 = 0;
    } 
    else {
      trigTech40_41 = 0;
      trigTech36_39 = 0;
      trigTech34 = 0;
    }
  } 
  else {
    trigTech40_41 = 0; 
    trigTech36_39 = 0;
    trigTech34 = 0;
  }

  // Test noScraping
  bool isNotScraping = scrapingTest( iEvent, iSetup );
  if( isNotScraping ) {
    trigScraping = 1;
  }
  else trigScraping = 0;

  // Test algo bit 124
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd);
  const L1GtTriggerMenu* menu = menuRcd.product();

  string algo("L1_BscMinBiasOR_BptxPlusORMinus");

  const DecisionWord dWord = gtRecord->decisionWord();
  bool bit124passed = false;
  bit124passed |= menu->gtAlgorithmResult(algo, dWord);

  if( bit124passed ) trigAlgo124 = 1;
  else trigAlgo124 = 0;

  GlobalPoint beamSpotPos;
  //GlobalPoint beamSpotPos(theBeamSpotH->position().x(),
  //			  theBeamSpotH->position().y(),
  //			  theBeamSpotH->position().z());
  if( beamSpotAvailable ) {
    beamSpotPos = GlobalPoint(theBeamSpotH->position().x(),
			      theBeamSpotH->position().y(),
			      theBeamSpotH->position().z());
  }
  
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  //cout << "Size: " << patHandle->size() << endl;
  /*for( pat::CompositeCandidateCollection::const_iterator iVEE = patHandle->begin();
       iVEE != patHandle->end();
       iVEE++ ) {
    cout << "PDG: " << iVEE->pdgId() << endl;
    }*/
  
  for( pat::CompositeCandidateCollection::const_iterator iV0 = patHandle->begin();
       iV0 != patHandle->end();
       iV0++ ) {

    if(instanceName == string("Kshort")){
      if(iV0->pdgId() != 310) 
	continue;
    }
    else if(instanceName == string("Lambda")){
      if(abs(iV0->pdgId()) != 3122) 
	continue;
    }
        
    
    reco::Track v0DauLowestpTTrack;

    //set pos daughter to be [0] and negative to [1] 
    std::vector<reco::Track> v0Daughters;
    std::vector<const reco::Candidate*> v0DauCands;
    //v0Daughters.push_back( iV0->getPosTrack() );
    
    reco::Track thePosDauTrk;
    reco::Track theNegDauTrk;
    reco::TrackRef thePosDauTrkRef;
    reco::TrackRef theNegDauTrkRef;
    reco::Candidate* thePosDauCand;
    reco::Candidate* theNegDauCand;
    
    if((dynamic_cast<const reco::RecoCandidate*>(iV0->daughter(0)))->charge() > 0){
      thePosDauTrkRef = ((dynamic_cast<const reco::RecoCandidate*>(iV0->daughter(0)))->track());
      theNegDauTrkRef = ((dynamic_cast<const reco::RecoCandidate*>(iV0->daughter(1)))->track());
      thePosDauTrk = *thePosDauTrkRef;
      theNegDauTrk = *theNegDauTrkRef;
      thePosDauCand = (reco::Candidate*)(iV0->daughter(0));
      theNegDauCand = (reco::Candidate*)(iV0->daughter(1));
    }
    else{
      thePosDauTrkRef = ((dynamic_cast<const reco::RecoCandidate*>(iV0->daughter(1)))->track());
      theNegDauTrkRef = ((dynamic_cast<const reco::RecoCandidate*>(iV0->daughter(0)))->track());
      thePosDauTrk = *thePosDauTrkRef;
      theNegDauTrk = *theNegDauTrkRef;
      thePosDauCand = (reco::Candidate*)(iV0->daughter(1));
      theNegDauCand = (reco::Candidate*)(iV0->daughter(0));
    }

    /*vector<TrackRef> theRefs;
    vector<RecoChargedCandidate> theCCands;
    cout << "Pushing back" << endl;
    theCCands.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
			   (iV0->daughter(0))) );
    cout << "Success 1" << endl;
    theCCands.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
			   (iV0->daughter(1))) );
    cout << "Success 2" << endl;
    for( unsigned int ndx0 = 0; ndx0 < theCCands.size(); ndx0++ ) {
      theRefs.push_back( theCCands[ndx0].track() );
      }*/
    
    v0Daughters.push_back(thePosDauTrk);
    v0Daughters.push_back(theNegDauTrk);
    
    v0DauCands.push_back(thePosDauCand);
    v0DauCands.push_back(theNegDauCand);

    std::vector<reco::TransientTrack> newVtxTracks;
    bool tracksRemoved = false;
    for( std::vector<TrackBaseRef>::const_iterator iTK = theOrigPrim->tracks_begin();
	 iTK != theOrigPrim->tracks_end();
	 iTK++ ) {
      TrackRef currentRef = iTK->castTo<TrackRef>();
      if( !((currentRef == thePosDauTrkRef) || (currentRef == theNegDauTrkRef)) ) {
	newVtxTracks.push_back( theTTB->build( *currentRef ) );
	//cout << "Didn't drop the track.  " << endl;
      }
      else {
	tracksRemoved = true;
	//cout << "Dropped a track in V0s" << endl;
      }
    }

    reco::Vertex refitVtx;
    //reco::Vertex* theOrigPrim = thePrim;
    if( newVtxTracks.size() > 0 ) {
      AdaptiveVertexFitter theAdFitter;
      TransientVertex theNewPrim = theAdFitter.vertex( newVtxTracks );//, *theBeamSpotH );
      if( theNewPrim.isValid() ) {
	refitVtx = theNewPrim;
	if( !tracksRemoved ) 
	  priVtxWasRefit->push_back(0);
	else priVtxWasRefit->push_back(1);
	priRefitVtxChi2->push_back( refitVtx.chi2() );
	priRefitVtxNdof->push_back( refitVtx.ndof() );
	priRefitVtxX->push_back( refitVtx.x() );
	priRefitVtxY->push_back( refitVtx.y() );
	priRefitVtxZ->push_back( refitVtx.z() );
	//priVtxChi2 = refitVtx.normalizedChi2();
	//priVtxNdof = refitVtx.normalizedChi2();
	//priVtxX = refitVtx.x();
	//priVtxY = refitVtx.y();
	//priVtxZ = refitVtx.z();
	//int nTracksInVtx = 0;
	//for( Vertex::trackRef_iterator itk = refitVtx.tracks_begin();
	//     itk != refitVtx.tracks_end();
	//     itk++ ) {
	//  if( refitVtx.trackWeight( *itk ) > 0.5 ) {
	//    nTracksInVtx++;
	//  }
	//}
	//priVtxNtracks = nTracksInVtx;
      }
      else {
	refitVtx = *theOrigPrim;
	if( !tracksRemoved ) priVtxWasRefit->push_back(0);
	else priVtxWasRefit->push_back(-2);
	priRefitVtxChi2->push_back( priVtxChi2 );
	priRefitVtxNdof->push_back( priVtxNdof );
	priRefitVtxX->push_back( priVtxX );
	priRefitVtxY->push_back( priVtxY );
	priRefitVtxZ->push_back( priVtxZ );
      }
    }
    else {
      refitVtx = *theOrigPrim;
      if( !tracksRemoved ) priVtxWasRefit->push_back(0);
      else priVtxWasRefit->push_back(-1);
      priRefitVtxChi2->push_back( priVtxChi2 );
      priRefitVtxNdof->push_back( priVtxNdof );
      priRefitVtxX->push_back( priVtxX );
      priRefitVtxY->push_back( priVtxY );
      priRefitVtxZ->push_back( priVtxZ );
    }
    newVtxTracks.clear();
    thePrim = &refitVtx;

    /*if( priVtxWasRefit && newVtxTracks.size() == thePrim->tracksSize() ) {
      cout << "REFIT: (" << thePrim->x() << ", " << thePrim->y() << ", " 
	   << thePrim->z() << ")" << endl
	   << "ORIG:  (" << theOrigPrim->x() << ", " << theOrigPrim->y() << ", "
	   << theOrigPrim->z() << ")" << endl; 
	   }*/

    if( iV0->daughter(0) ) {
      //cout << "Type of daughter(0): " << typeid(*(iV0->daughter(0))) << endl;
    }
    //cout << "Pushing back cands..." << endl;

    double tkpT = 100000.;
    bool wroteV0TreeHits = false;
    bool wroteKsHistoHits = false;
    //bool wroteLamTreeHits = false;
    bool wroteLamHistoHits = false;

    GlobalVector v0Momentum;
    GlobalVector posDauMomentum;
    GlobalVector negDauMomentum;
    GlobalVector posDauCandMomentum;
    GlobalVector negDauCandMomentum;
    GlobalVector lineOfFlight;
    
    // is there an instance when these are not initialized?
    /*
    posDauMomentum = GlobalVector( 0., 0., 0. );
    negDauMomentum = GlobalVector( 0., 0., 0. );      
    posDauCandMomentum = GlobalVector( 0., 0., 0. );  
    negDauCandMomentum = GlobalVector( 0., 0., 0.);   
    */  
    posDauMomentum = GlobalVector( v0Daughters[0].momentum().x(),      
    				   v0Daughters[0].momentum().y(),
    				   v0Daughters[0].momentum().z() );
    negDauMomentum = GlobalVector( v0Daughters[1].momentum().x(),      
    				   v0Daughters[1].momentum().y(),
    				   v0Daughters[1].momentum().z() );
    posDauCandMomentum = GlobalVector( v0DauCands[0]->momentum().x(),  
    				   v0DauCands[0]->momentum().y(),
    				   v0DauCands[0]->momentum().z() );
    negDauCandMomentum = GlobalVector( v0DauCands[1]->momentum().x(),  
    				   v0DauCands[1]->momentum().y(),
    				   v0DauCands[1]->momentum().z() );

    v0Momentum = GlobalVector( iV0->momentum().x(),
			       iV0->momentum().y(),
			       iV0->momentum().z() );

    double pTPlus = ( posDauCandMomentum.cross( v0Momentum.unit() ) ).mag();
    double pTMinus = ( negDauCandMomentum.cross( v0Momentum.unit() ) ).mag();
    double pLPlus = posDauCandMomentum.dot( v0Momentum.unit() );
    double pLMinus = negDauCandMomentum.dot( v0Momentum.unit() );
    //cout << "pTPlus = " << pTPlus << ", pTMinus = " << pTMinus << endl;

    if(putTreeInTFile){
      v0pTPlus->push_back( pTPlus );
      v0pTMinus->push_back( pTMinus );
      v0pLPlus->push_back( pLPlus );
      v0pLMinus->push_back( pLMinus );
    }
    
    // create a kinematic particle to be used for pca calculation
    
    
    KinematicParticleFactoryFromTransientTrack pFactory;
    vector<RefCountedKinematicParticle> v0Particles;
    vector<double> v0ParticleMasses;
    TransientTrack posDauTT;
    TransientTrack negDauTT;
    //initial chi2 and ndf before kinematic fits.     
    float chi = 0.;				      
    float ndf = 0.;				      
    
    // get kinematic particles          		   
    posDauTT = (*theTTB).build(v0Daughters[0]);   		     
    negDauTT = (*theTTB).build(v0Daughters[1]);   		     
    if( putTreeInTFile){					     
      posDauNhits->push_back( v0Daughters[0].numberOfValidHits() );
      negDauNhits->push_back( v0Daughters[1].numberOfValidHits() );
    }								     

    ClosestApproachInRPhi cApp;
    FreeTrajectoryState posState = posDauTT.impactPointTSCP().theState();
    FreeTrajectoryState negState = negDauTT.impactPointTSCP().theState();

    double sumPx = -1000.;
    double sumPy = -1000.;
    double sumPz = -1000.;

    if (posDauTT.impactPointTSCP().isValid() && negDauTT.impactPointTSCP().isValid()) {
      cApp.calculate(posState, negState);
      if (cApp.status()) {
	GlobalPoint cxPt = cApp.crossingPoint();
	if ((cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) < 120. && abs(cxPt.z()) < 300.) {
	  TrajectoryStateClosestToPoint posTSCP = posDauTT.trajectoryStateClosestToPoint(cxPt);
	  TrajectoryStateClosestToPoint negTSCP = negDauTT.trajectoryStateClosestToPoint(cxPt);
	  if (posTSCP.isValid() && negTSCP.isValid()) {
	    sumPx = posTSCP.momentum().x() + negTSCP.momentum().x();
	    sumPy = posTSCP.momentum().y() + negTSCP.momentum().y();
	    sumPz = posTSCP.momentum().z() + negTSCP.momentum().z();
	  }
	}
      }
    }

    if (putTreeInTFile) {
      v0DauSumPx->push_back(sumPx);
      v0DauSumPy->push_back(sumPy);
      v0DauSumPz->push_back(sumPz);
    }

    float pionSigma = piMass*1.e-6;
    float protonSigma = protonMass*1.e-6;
    
    if(iV0->pdgId() == 310 && instanceName == string("Kshort") && posDauTT.isValid() && negDauTT.isValid() ){
      v0Particles.push_back(pFactory.particle(posDauTT, piMass, chi, ndf, pionSigma));
      v0Particles.push_back(pFactory.particle(negDauTT, piMass, chi, ndf, pionSigma));
      v0ParticleMasses.push_back( piMass );
      v0ParticleMasses.push_back( piMass );
    }
    if(iV0->pdgId() == 3122 && instanceName == string("Lambda") && posDauTT.isValid() && negDauTT.isValid()){
      v0Particles.push_back(pFactory.particle(posDauTT, protonMass, chi, ndf, protonSigma));
      v0Particles.push_back(pFactory.particle(negDauTT, piMass, chi, ndf, pionSigma));
      v0ParticleMasses.push_back( protonMass );
      v0ParticleMasses.push_back( piMass );
    }
    if(iV0->pdgId() == -3122 && instanceName == string("Lambda") && posDauTT.isValid() && negDauTT.isValid()){
      v0Particles.push_back(pFactory.particle(negDauTT, protonMass, chi, ndf, protonSigma));
      v0Particles.push_back(pFactory.particle(posDauTT, piMass, chi, ndf, pionSigma));
      v0ParticleMasses.push_back( protonMass );
      v0ParticleMasses.push_back( piMass );
    }
    
    	   
    KinematicParticleVertexFitter fitter;   							 
    RefCountedKinematicTree v0VertexFitTree;
    bool doKinematicFit = false;
    if(v0Particles.size() == 2){
      doKinematicFit = true;
      v0VertexFitTree = fitter.fit(v0Particles);
    }
    
    if(putTreeInTFile){
      //fill quantities for pos and neg daughters
      v0PosDaupT->push_back(v0Daughters[0].pt());
      v0NegDaupT->push_back(v0Daughters[1].pt());
      
      v0PosDaup->push_back(v0Daughters[0].p());
      v0NegDaup->push_back(v0Daughters[1].p());
      
      v0PosDaupX->push_back(v0Daughters[0].momentum().x());
      v0NegDaupX->push_back(v0Daughters[1].momentum().x());
      
      v0NegDaupY->push_back(v0Daughters[1].momentum().y());
      v0PosDaupY->push_back(v0Daughters[0].momentum().y());
      
      v0PosDaupZ->push_back(v0Daughters[0].momentum().z());
      v0NegDaupZ->push_back(v0Daughters[1].momentum().z());
      
      v0PosDauEta->push_back(v0Daughters[0].eta());
      v0NegDauEta->push_back(v0Daughters[1].eta());
      
      v0PosDauChi2->push_back(v0Daughters[0].chi2());
      v0NegDauChi2->push_back(v0Daughters[1].chi2());
      
      v0PosDauNormChi2->push_back(v0Daughters[0].normalizedChi2());
      v0NegDauNormChi2->push_back(v0Daughters[1].normalizedChi2());
      
      v0PosDauD0->push_back(v0Daughters[0].d0());
      v0NegDauD0->push_back(v0Daughters[1].d0());
      
      if( !wroteV0TreeHits ){
      	v0PosDauNValidHits->push_back( (float) (iV0->userInt("posTkNPixelHits") 
						   + iV0->userInt("posTkNStripHits") ) );
      	v0NegDauNValidHits->push_back( (float) (iV0->userInt("negTkNPixelHits") 
						   + iV0->userInt("negTkNStripHits") ) );
      
      	v0PosDauNPixelHits->push_back( (float) iV0->userInt("posTkNPixelHits") );
      	v0NegDauNPixelHits->push_back( (float) iV0->userInt("negTkNPixelHits") );
	
	v0PosDauNStripHits->push_back( (float) iV0->userInt("posTkNStripHits") );
	v0NegDauNStripHits->push_back( (float) iV0->userInt("negTkNStripHits") );
	wroteV0TreeHits = true;
      }
    }											 
    							
    for( std::vector<reco::Track>::iterator iDAU = v0Daughters.begin();
	 iDAU != v0Daughters.end();
	 iDAU++ ) {
	 
      
      if( iDAU->pt() < tkpT ) {
	v0DauLowestpTTrack = *iDAU;
      }
      tkpT = iDAU->pt();
      // Fill daughter track branches
      if(instanceName == string("Kshort") ) {
	if( (iV0->mass() > ksMassConst - 3*0.0058)
	    && (iV0->mass() < ksMassConst + 3*0.0058) 
	    && writeHistos ) {
	  ksDaupT->Fill( iDAU->pt() );
	  ksDaup->Fill( iDAU->p() );
	  ksDauEta->Fill( iDAU->eta() );
	  ksDauChi2->Fill( iDAU->chi2() );
	  ksDauNormChi2->Fill( iDAU->normalizedChi2() );
	  ksDauD0->Fill( iDAU->d0() );
	  if( !wroteKsHistoHits ) {
	    ksDauNValidHits->Fill( (float) (iV0->userInt("posTkNPixelHits")
				       + iV0->userInt("posTkNStripHits") ) );
	    ksDauNPixelHits->Fill( (float) iV0->userInt("posTkNPixelHits") );
	    ksDauNStripHits->Fill( (float) iV0->userInt("posTkNStripHits") );
	    ksDauNValidHits->Fill( (float) (iV0->userInt("negTkNPixelHits")
				       + iV0->userInt("negTkNStripHits") ) );
	    ksDauNPixelHits->Fill( (float) iV0->userInt("negTkNPixelHits") );
	    ksDauNStripHits->Fill( (float) iV0->userInt("negTkNStripHits") );
	    wroteKsHistoHits = true;
	  }
	}
      }
      // Lambda sigma = 0.002389 ~ 0.0024
      if(instanceName == string("Lambda") ) {
	if( (iV0->mass() > lamMassConst - 3*0.0024)
	    && (iV0->mass() < lamMassConst + 3*0.0024) 
	    && writeHistos ) {
	  lamDaupT->Fill( iDAU->pt() );
	  lamDaup->Fill( iDAU->p() );
	  lamDauEta->Fill( iDAU->eta() );
	  lamDauChi2->Fill( iDAU->chi2() );
	  lamDauNormChi2->Fill( iDAU->normalizedChi2() );
	  lamDauD0->Fill( iDAU->d0() );
	  if( !wroteLamHistoHits ) {
	    lamDauNValidHits->Fill( (float) (iV0->userInt("posTkNPixelHits")
					+ iV0->userInt("posTkNStripHits") ) );
	    lamDauNPixelHits->Fill( (float) iV0->userInt("posTkNPixelHits") );
	    lamDauNStripHits->Fill( (float) iV0->userInt("posTkNStripHits") );
	    lamDauNValidHits->Fill( (float) (iV0->userInt("negTkNPixelHits")
					+ iV0->userInt("negTkNStripHits") ) );
	    lamDauNPixelHits->Fill( (float) iV0->userInt("negTkNPixelHits") );
	    lamDauNStripHits->Fill( (float) iV0->userInt("negTkNStripHits") );
	    wroteLamHistoHits = true;
	  }
	}
      }
    }

    //cout << "They suck." << endl;
    double v0DauPhi_lowestpT_forHistos = 0.;
    // Fill daughter branches at vertex position
    GlobalVector posP( iV0->userFloat("posTkPX"),
		       iV0->userFloat("posTkPY"),
		       iV0->userFloat("posTkPZ") );
    GlobalVector negP( iV0->userFloat("negTkPX"),
		       iV0->userFloat("negTkPY"),
		       iV0->userFloat("negTkPZ") );
    //double protonE = sqrt(posP.mag2() + protonMass2);
    //double antiprotonE = sqrt(negP.mag2() + protonMass2);
    //double posPiE = sqrt(posP.mag2() + piMass2);
    //double negPiE = sqrt(negP.mag2() + piMass2);
    //GlobalVector totalP = posP + negP;
    if(instanceName == string("Kshort") ) {

      //v0VtxChi2->push_back( iV0->UserFloat("vtxChi2") );
      //v0VtxNormChi2->push_back( iV0->UserFloat("vtxNormChi2") );

      //double lamE = protonE + negPiE;
      //double lamBarE = antiprotonE + posPiE;
      //double lamMass = sqrt( lamE*lamE - totalP.mag2() );
      //double lamBarMass = sqrt( lamBarE*lamBarE - totalP.mag2() );
    								        
      if( (iV0->mass() > ksMassConst - 3*0.0058)        	      
          && (iV0->mass() < ksMassConst + 3*0.0058) 
          && writeHistos ) {
        ksDauPhi->Fill( posP.phi() );
        ksDauPhi->Fill( negP.phi() );
      } 					        	      

      /*if( putTreeInTFile ) {
	if( (lamMass - lamMassConst) < (lamBarMass - lamMassConst) ) {
	  //v0OtherCandMass->push_back( lamMass );
	}
	else {
	  //v0OtherCandMass->push_back( lamBarMass );
	}
	}*/
    }
    
    else if(instanceName == string("Lambda") ) {

      //double ksE = posPiE + negPiE;
      //double ksM = sqrt( ksE*ksE - totalP.mag2() );
      
      if( (iV0->mass() > lamMassConst - 3*0.0024)
	  && (iV0->mass() < lamMassConst + 3*0.0024) 
	  && writeHistos) {
	lamDauPhi->Fill( posP.phi() );
	lamDauPhi->Fill( negP.phi() );
      }

      //if( putTreeInTFile )
	//v0OtherCandMass->push_back( ksM );
    }

    if( putTreeInTFile ) {					        
      v0DaupT_lowestpT->push_back( v0DauLowestpTTrack.pt() );
      v0DauEta_lowestpT->push_back( v0DauLowestpTTrack.eta() );
      v0VtxChi2->push_back( iV0->userFloat("vtxChi2") );
      v0VtxNormChi2->push_back( iV0->userFloat("vtxNormChi2") );
    }								        

    if( v0DauLowestpTTrack.charge() > 0. ) {			        
      if( putTreeInTFile ) v0DauPhi_lowestpT->push_back( posP.phi() );
      v0DauPhi_lowestpT_forHistos = posP.phi();
    }								        
    else {							        
      if( putTreeInTFile ) v0DauPhi_lowestpT->push_back( negP.phi() );
      v0DauPhi_lowestpT_forHistos = negP.phi();
    }
    								        
    if( putTreeInTFile ) {					        
      v0PosDauAtVtxpX->push_back( posP.x() );
      v0PosDauAtVtxpY->push_back( posP.y() );
      v0PosDauAtVtxpZ->push_back( posP.z() );
      v0PosDauPhi->push_back( posP.phi() );

      v0NegDauAtVtxpX->push_back( negP.x() );
      v0NegDauAtVtxpY->push_back( negP.y() );
      v0NegDauAtVtxpZ->push_back( negP.z() );
      v0NegDauPhi->push_back( negP.phi() );
    }

    //cout << "They really do." << endl;
    // Fill mass branch
    
    //if( putTreeInTFile ) v0CandMass->push_back( iV0->mass() );
    
    if(instanceName == string("Kshort") && writeHistos) ksCandMass->Fill( iV0->mass() );

    if(instanceName == string("Lambda") && writeHistos )lamCandMass->Fill( iV0->mass() );

    // t = m(L - 15sigma)/p

    // Fill branches vertex quantities
    GlobalVector v0p_( iV0->momentum().x(),
		       iV0->momentum().y(),
		       iV0->momentum().z() );


    

    //set defaults for pca variables...         
    float theV0CosThetaPAndLineOfFlight = -5.;  

    float v3dIpWrtPrimaryValue = -1000.;        
    float v3dIpWrtPrimarySigValue = -1000.;   

    //float v0MassFromFitValue = -1000.;
    //float v0MassErrorFromFitValue = -1000.; 
    float v0MassFromFitValue = iV0->userFloat("massFromFit");
    float v0MassErrorFromFitValue = iV0->userFloat("massErrorFromFit");
    float v0OtherMassFromFitValue = iV0->userFloat("otherMassFromFit");
    float v0OtherMassErrorFromFitValue = iV0->userFloat("otherMassErrorFromFit");

    float v0PCAPrimaryValue = -1000.;	        
    float v0PCAPrimaryErrorValue = -1000.;      
    float v0PCABeamSpotValue = -1000.;	        
    float v0PCABeamSpotErrorValue = -1000.;     
    					        
    float thePosDau3DIpValue = -1000.;	        
    float thePosDau3DIpSigValue = -1000.;       
    float theNegDau3DIpValue = -1000.;	        
    float theNegDau3DIpSigValue = -1000.;
    
    float v0VtxPriDistValue = -1000.;     
    float v0VtxPriDistErrValue = -1000.;     
    float ltime3D = -1000.;
    float ltime3DErr = -1000.;
    float ltime3D_scaled = -1000.;
    float ltime3D_scaledErr = -1000.;
    float ltime3D_MPV = -1000.;
    float ltime3D_MPVErr = -1000.;
    
    // Lifetime stuff
    if( beamSpotAvailable ) {    
      double cov00 = iV0->userFloat("vtxCov00");
      double cov01 = iV0->userFloat("vtxCov01");
      double cov11 = iV0->userFloat("vtxCov11");
      double cov02 = iV0->userFloat("vtxCov02");
      double cov12 = iV0->userFloat("vtxCov12");
      double cov22 = iV0->userFloat("vtxCov22");
      std::vector<double> covv;
      covv.push_back( cov00 );
      covv.push_back( cov01 );
      covv.push_back( cov11 );
      covv.push_back( cov02 );
      covv.push_back( cov12 );
      covv.push_back( cov22 );

      SMatrixSym3D vtxCov( covv.begin(), covv.end() );

      // Adding stuff for 3D lifetime.  This is the Vertex.
      //  Chi2 is BOGUS.  It shouldn't be needed for VertexDistance,
      //  so don't expect it to be right for any other purpose.
      reco::Vertex v0Vtx(iV0->vertex(), vtxCov, 1., 1., 2);
      //v0VtxPriDist(Err)

      /*bool isKs = false;
      bool isLam = false;
      if( putTreeInTFile && iV0->pdgId() == 310 
	  && instanceName == string("Kshort") ) {
	isKs = true;
      }
      if( putTreeInTFile && abs(iV0->pdgId()) == 3122 
	  && instanceName == string("Lambda") ) {
	isLam = true;
	}*/
	
      if( pvAvailable && thePrim ) {
	//cout << "In new PV stuff" << endl;
	VertexDistance3D theDistCalculator;

	double ksBetaGamma  = v0p_.mag() / ksMassConst;
	double lamBetaGamma = v0p_.mag() / lamMassConst;

	GlobalPoint priVtxPos( thePrim->position().x(),
			       thePrim->position().y(),
			       thePrim->position().z() );
	GlobalPoint vVtxPos( v0Vtx.position().x(),
			     v0Vtx.position().y(),
			     v0Vtx.position().z() );

	GlobalVector sep3D = vVtxPos - priVtxPos;

	float scale = fabs( sep3D.dot(v0p_) / sep3D.mag() / v0p_.mag() ); 

	Measurement1D theDist1D =
	  theDistCalculator.signedDistance(*thePrim, v0Vtx, v0p_);
	  
	v0VtxPriDistValue = theDist1D.value();
	v0VtxPriDistErrValue = theDist1D.error();

	// Now calculate by matrices
	AlgebraicMatrix31 pB;
	pB(0,0) = v0p_.x();
	pB(1,0) = v0p_.y();
	pB(2,0) = v0p_.z();

	AlgebraicMatrix13 pBT;
	pBT(0,0) = v0p_.x();
	pBT(0,1) = v0p_.y();
	pBT(0,2) = v0p_.z();

	AlgebraicMatrix31 PV;
	PV(0,0) = thePrim->position().x();
	PV(1,0) = thePrim->position().y();
	PV(2,0) = thePrim->position().z();

	AlgebraicMatrix31 BV;
	BV(0,0) = v0Vtx.position().z();
	BV(1,0) = v0Vtx.position().y();
	BV(2,0) = v0Vtx.position().z();
	AlgebraicMatrix31 lxyz = BV - PV;
	AlgebraicMatrix33 PVError( thePrim->error() );
	AlgebraicMatrix33 BVError( v0Vtx.error() );
	AlgebraicMatrix33 lxyzError = PVError + BVError;
	lxyzError.Invert();
	AlgebraicMatrix11 a = pBT * lxyzError * pB;
	AlgebraicMatrix11 b = pBT * lxyzError * lxyz;
	double num = b(0,0);
	double deno = a(0,0);
	
	//cout << "Just did distance calculation" << endl;
	//ltime3D = fabs( theDist1D.value() ) * iV0->mass() / v0p_.mag();
	if( iV0->pdgId() == 310 ) {
	  ltime3D = fabs( theDist1D.value() ) / ksBetaGamma;
	  ltime3DErr = fabs( theDist1D.error() ) / ksBetaGamma;
	  ltime3D_scaled = sep3D.dot(v0p_) / ksBetaGamma / v0p_.mag();
	  ltime3D_scaledErr = fabs( theDist1D.error() ) * scale / ksBetaGamma;
	  ltime3D_MPV = num * ksMassConst / deno;
	  ltime3D_MPVErr = ksMassConst / sqrt(deno);
	}
	else if( abs( iV0->pdgId() ) == 3122 ) {
	  ltime3D = fabs( theDist1D.value() ) / lamBetaGamma;
	  ltime3DErr = fabs( theDist1D.error() ) / lamBetaGamma;
	  ltime3D_scaled = sep3D.dot(v0p_) / lamBetaGamma / v0p_.mag();
	  ltime3D_scaledErr = fabs( theDist1D.error() ) * scale / lamBetaGamma;
	  ltime3D_MPV = num * lamMassConst / deno;
	  ltime3D_MPVErr = lamMassConst / sqrt(deno);
	}

	//cout << "Found 3D lifetime: " << ltime3D << endl;
	
	// calculate line of flight between primary and secondary
	lineOfFlight = GlobalVector (iV0->vertex().x() - priVtxX,
				     iV0->vertex().y() - priVtxY,
				     iV0->vertex().z() - priVtxZ);
	
	// get angle between momentum vector and line of flight			     
	float theAngle = angle(lineOfFlight.x(), 
			      lineOfFlight.y(), 
	  		      lineOfFlight.z(),
			      v0Momentum.x(),
			      v0Momentum.y(),
			      v0Momentum.z());
	
	theV0CosThetaPAndLineOfFlight = cos(theAngle);
	
	// calulate pca for both daughter wrt to primary
	pair<bool,Measurement1D> thePosDau3DIpPair = IPTools::absoluteImpactParameter3D(posDauTT, (*thePrim));
	if(thePosDau3DIpPair.first){
	  thePosDau3DIpValue = thePosDau3DIpPair.second.value();
	  thePosDau3DIpSigValue = thePosDau3DIpPair.second.significance();
	}

	pair<bool,Measurement1D> theNegDau3DIpPair = IPTools::absoluteImpactParameter3D(negDauTT, (*thePrim));
	if(theNegDau3DIpPair.first){
	  theNegDau3DIpValue = theNegDau3DIpPair.second.value();
	  theNegDau3DIpSigValue = theNegDau3DIpPair.second.significance();
	}
	
	// calculate the pca of the V0 wrt to the primary (beamspot)
	
	const GlobalPoint thePrimaryVertexPoint(priVtxX, priVtxY, priVtxZ);
	  
	if( doKinematicFit && v0VertexFitTree->isValid()){
	  
      	  v0VertexFitTree->movePointerToTheTop();						       

          RefCountedKinematicParticle v0FitKinematicParticle = v0VertexFitTree->currentParticle();     
      	  RefCountedKinematicVertex v0FitKinematicVertex = v0VertexFitTree->currentDecayVertex();
	  /*InvariantMassFromVertex massWithErrCalculator;
	  Measurement1D massWithError 
	    = massWithErrCalculator.invariantMass(v0FitKinematicVertex, v0ParticleMasses);

	  v0MassFromFitValue = massWithError.value();
	  v0MassErrorFromFitValue = massWithError.error();*/
          
	  KinematicState theCurrentKinematicState = v0FitKinematicParticle->currentState();

	  /*KinematicParameters theKinParams = theCurrentKinematicState.kinematicParameters();
	  KinematicParametersError theKinParamErrors 
	    = theCurrentKinematicState.kinematicParametersError();
	  AlgebraicSymMatrix77 mtrx = theKinParamErrors.matrix();
	  double kinParams66 = sqrt( mtrx(6,6) );

	  if( theKinParamErrors.isValid() && theKinParams.isValid() ) {
	    v0MassFromFitValue = theKinParams.mass();
	    v0MassErrorFromFitValue = kinParams66;
	    }*/

          //Get trajectory state from kinematic state
          FreeTrajectoryState theV0FTS = theCurrentKinematicState.freeTrajectoryState();
          //Get TransientTrack from FTS
          TransientTrack v0TT = (*theTTB).build(theV0FTS);
	 
	  pair<bool,Measurement1D> theV03DIpPair = IPTools::signedImpactParameter3D(v0TT, v0Momentum, (*thePrim));
	  if(theV03DIpPair.first ){  // set to false if the TT is invalid...
	     v3dIpWrtPrimaryValue = theV03DIpPair.second.value();
	     v3dIpWrtPrimarySigValue = theV03DIpPair.second.significance();
	  }

	  if(v0TT.isValid()){
      	    TrajectoryStateClosestToPoint  theV0TrjPCAPrimary = v0TT.trajectoryStateClosestToPoint(thePrimaryVertexPoint);
	    if(theV0TrjPCAPrimary.isValid()){
      	      v0PCAPrimaryValue = theV0TrjPCAPrimary.perigeeParameters().transverseImpactParameter();
      	      if(theV0TrjPCAPrimary.hasError())v0PCAPrimaryErrorValue = theV0TrjPCAPrimary.perigeeError().transverseImpactParameterError();
      	    }  

	    TrajectoryStateClosestToPoint theV0TrjPCABeamSpot = v0TT.trajectoryStateClosestToPoint(beamSpotPos);
	    if(theV0TrjPCABeamSpot.isValid()){
      	      v0PCABeamSpotValue = theV0TrjPCABeamSpot.perigeeParameters().transverseImpactParameter();
      	      if(theV0TrjPCABeamSpot.hasError()) v0PCABeamSpotErrorValue = theV0TrjPCABeamSpot.perigeeError().transverseImpactParameterError();
      	    }
      	  }
	}
	//cout << "Deleting primary we created before" << endl;
	//delete thePrim;
      }

      GlobalPoint primPos(thePrim->position().x(),
			  thePrim->position().y(),
			  thePrim->position().z());
      //if (!priVtxIsFake) beamSpotPos = primPos;

      SMatrixSym3D totalCov = theBeamSpotH->covariance3D() + vtxCov;
      SMatrixSym3D totalCov3D = vtxCov + thePrim->covariance();
      SVector3 distanceVector( iV0->vertex().x() - beamSpotPos.x(),
			       iV0->vertex().y() - beamSpotPos.y(),
			       0. );
      //SVector3 distanceVector( iV0->vertex().x() - thePrim->position().x(),
      //		       iV0->vertex().y() - thePrim->position().y(),
      //		       0. );
      double rVtxMag = ROOT::Math::Mag(distanceVector);
      double sigmaRvtxMag = 
	sqrt(ROOT::Math::Similarity(totalCov, distanceVector)) / rVtxMag;

      SVector3 distanceVector3D( iV0->vertex().x() - thePrim->position().x(),
				 iV0->vertex().y() - thePrim->position().y(),
				 iV0->vertex().z() - thePrim->position().z() );
      double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
      double sigmaRvtxMag3D =
	sqrt(ROOT::Math::Similarity(totalCov3D, distanceVector3D)) / rVtxMag3D;



      //double sigmaR = (cov00*(iV0->vertex().x()*iV0->vertex().x())
      //	       + cov11*(iV0->vertex().y()*iV0->vertex().y())
      //	       + cov01*(iV0->vertex().x()*iV0->vertex().y())*2 );
      //double lifetime = v0p_.mag() * (iV0->vertex().Rho())/iV0->mass();
      //double ltime = (iV0->vertex().Rho() - 15.*sigmaR) 
      double ltime = ( rVtxMag - 15.*sigmaRvtxMag )
	//double ltime = (iV0->vertex().Rho())
	* iV0->mass()/v0p_.transverse();
      if( iV0->pdgId() == 310 && instanceName == string("Kshort") ) {
// 	cout << "V0RecoAnalyzer: 2d = " << rVtxMag/sigmaRvtxMag
// 	     << ", 3d = " << rVtxMag3D/sigmaRvtxMag3D << endl;
	if( putTreeInTFile ) {
	  v0Lifetime->push_back( ltime );
	  v0VtxSig->push_back( rVtxMag/sigmaRvtxMag );
	  v0VtxSig3D->push_back( rVtxMag3D/sigmaRvtxMag3D );
	  v0VtxBSTransDist->push_back( rVtxMag );
	  v0VtxBSTransDistErr->push_back( sigmaRvtxMag );
	  thePosDauTkQual->push_back( iV0->userInt("posQual") );
	  theNegDauTkQual->push_back( iV0->userInt("negQual") );
	  thePosDauTkAlgo->push_back( iV0->userInt("posAlgo") );
	  theNegDauTkAlgo->push_back( iV0->userInt("negAlgo") );
	}
	if( writeHistos ) {
	  //  Nothing yet.
	}
      }
      if( abs(iV0->pdgId()) == 3122 && instanceName == string("Lambda") ) {
	if( putTreeInTFile ) {
	  v0Lifetime->push_back( ltime );
	  v0VtxSig->push_back( rVtxMag/sigmaRvtxMag );
	  v0VtxSig3D->push_back( rVtxMag3D/sigmaRvtxMag3D );
	  v0VtxBSTransDist->push_back( rVtxMag );
	  v0VtxBSTransDistErr->push_back( sigmaRvtxMag );
	  thePosDauTkQual->push_back( iV0->userInt("posQual") );
	  theNegDauTkQual->push_back( iV0->userInt("negQual") );
	  thePosDauTkAlgo->push_back( iV0->userInt("posAlgo") );
	  theNegDauTkAlgo->push_back( iV0->userInt("negAlgo") );
	}
      }
    }//if(beamSpotAvailable)
    
    if(putTreeInTFile){ //fill pca variables						       
      v0CosThetaPAndLineOfFlight->push_back(theV0CosThetaPAndLineOfFlight);
      
      v3dIpWrtPrimary->push_back(v3dIpWrtPrimaryValue);
      v3dIpWrtPrimarySig->push_back(v3dIpWrtPrimarySigValue);
      
      v0PCAPrimary->push_back(v0PCAPrimaryValue);
      v0PCAPrimaryError->push_back(v0PCAPrimaryErrorValue);
      v0PCABeamSpot->push_back(v0PCABeamSpotValue);
      v0PCABeamSpotError->push_back(v0PCABeamSpotErrorValue);
      
      thePosDau3DIp->push_back(thePosDau3DIpValue);
      thePosDau3DIpSig->push_back(thePosDau3DIpSigValue);
      theNegDau3DIp->push_back(theNegDau3DIpValue);
      theNegDau3DIpSig->push_back(theNegDau3DIpSigValue);
      
      v0VtxPriDist->push_back(v0VtxPriDistValue);
      v0VtxPriDistErr->push_back(v0VtxPriDistErrValue);
      v03DLifetime->push_back( ltime3D );
      v03DLifetimeError->push_back( ltime3DErr );
      v03DLifetimeScaled->push_back( ltime3D_scaled );
      v03DLifetimeScaledError->push_back( ltime3D_scaledErr );
      v03DLifetimeMPV->push_back( ltime3D_MPV );
      v03DLifetimeMPVError->push_back( ltime3D_MPVErr );
      
      v0VtxX->push_back( iV0->vertex().x() );
      v0VtxY->push_back( iV0->vertex().y() );
      v0VtxZ->push_back( iV0->vertex().z() );
      v0VtxR->push_back( iV0->vertex().Rho() );
      v0pX->push_back( v0p_.x() );
      v0pY->push_back( v0p_.y() );
      v0pZ->push_back( v0p_.z() );
      v0pT->push_back( v0p_.transverse() );
      v0p->push_back( v0p_.mag() );
      v0Eta->push_back( v0p_.eta() );

      double v0E_ = 0.;
      if( iV0->pdgId() == 310 ) {
	v0E_ = sqrt( v0p_.mag()*v0p_.mag() + ksMassConst2 );
      }
      else {
	v0E_ = sqrt( v0p_.mag()*v0p_.mag() + lamMassConst2 );
      }
      double v0y = 0.5*log((v0E_ + v0p_.z()) / (v0E_ - v0p_.z()));
      v0Rapidity->push_back( v0y );

      v0Phi->push_back( v0p_.phi() );
      v0PDG->push_back( iV0->pdgId() );

      //v0MassFromFit->push_back( v0MassFromFitValue );
      v0CandMass->push_back( v0MassFromFitValue );
      //v0MassErrorFromFit->push_back( v0MassErrorFromFitValue );
      v0CandMassError->push_back( v0MassErrorFromFitValue );
      //v0OtherMassFromFit->push_back( v0OtherMassFromFitValue );
      //v0OtherMassErrorFromFit->push_back( v0OtherMassErrorFromFitValue );
      v0OtherCandMass->push_back( v0OtherMassFromFitValue );
      v0OtherCandMassError->push_back( v0OtherMassErrorFromFitValue );
    }


    // Fill branches with vertex quantities
    if(instanceName == string("Kshort") && writeHistos ) {
      ksRFull->Fill( iV0->vertex().Rho() );
      ksZFull->Fill( iV0->vertex().z() );
      ksZvsRFull->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
      ksEtaFull->Fill( v0p_.eta() );
      ksPhiFull->Fill( v0p_.phi() );
      ksPtFull->Fill( v0p_.transverse() );
      if( iV0->mass() > (ksMassConst - 0.01) 
          && iV0->mass() < (ksMassConst + 0.01) ) {
        ksRSignal->Fill( iV0->vertex().Rho() );
        ksZSignal->Fill( iV0->vertex().z() );
        ksZvsRSignal->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
        ksEtaSignal->Fill( v0p_.eta() );
        ksPhiSignal->Fill( v0p_.phi() );
        ksPtSignal->Fill( v0p_.transverse() );
      }
      if( (iV0->mass() > 0.445 && iV0->mass() < 0.465) ||
          (iV0->mass() > 0.55  && iV0->mass() < 0.57 ) ) {
        ksRBg->Fill( iV0->vertex().Rho() );
        ksZBg->Fill( iV0->vertex().z() );
        //ksZvsRBg
        ksEtaBg->Fill( v0p_.eta() );
        ksPhiBg->Fill( v0p_.phi() );
        ksPtBg->Fill( v0p_.transverse() );
      }
      if( (iV0->mass() < (ksMassConst - 0.01))
          || (iV0->mass() > (ksMassConst + 0.01)) ) {
        ksZvsRBg->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
      }
    }
    if(instanceName == string("Lambda") && writeHistos ) {
      lamRFull->Fill( iV0->vertex().Rho() );
      lamZFull->Fill( iV0->vertex().z() );
      lamZvsRFull->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
      lamEtaFull->Fill( v0p_.eta() );
      lamPhiFull->Fill( v0p_.phi() );
      lamPtFull->Fill( v0p_.transverse() );
      if( iV0->mass() > (lamMassConst - 0.01) 
          && iV0->mass() < (lamMassConst + 0.01) ) {
        lamRSignal->Fill( iV0->vertex().Rho() );
        lamZSignal->Fill( iV0->vertex().z() );
        lamZvsRSignal->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
        lamEtaSignal->Fill( v0p_.eta() );
        lamPhiSignal->Fill( v0p_.phi() );
        lamPtSignal->Fill( v0p_.transverse() );
      }
      if( (iV0->mass() > 0.445 && iV0->mass() < 0.465) ||
          (iV0->mass() > 0.55  && iV0->mass() < 0.57 ) ) {
        lamRBg->Fill( iV0->vertex().Rho() );
        lamZBg->Fill( iV0->vertex().z() );
        //lamZvsRBg
        lamEtaBg->Fill( v0p_.eta() );
        lamPhiBg->Fill( v0p_.phi() );
        lamPtBg->Fill( v0p_.transverse() );
      }
      if( (iV0->mass() < (lamMassConst - 0.01))
          || (iV0->mass() > (lamMassConst + 0.01)) ) {
        lamZvsRBg->Fill( iV0->vertex().z(), iV0->vertex().Rho(), 1. );
      }
    }

    if( iV0->pdgId() == 310 
	&& instanceName == string("Kshort") 
	&& writeHistos) {
      double ksMass_eta_binWidth = 5. / (double) (ksMass_eta_nBins*2);
      double ksMass_eta_min = -2.5;
      for( vector<TH1F*>::iterator iETA = ksMass_eta.begin();
	   iETA != ksMass_eta.end();
	   iETA++ ) {
	if( v0DauLowestpTTrack.eta() > ksMass_eta_min 
	    && v0DauLowestpTTrack.eta() 
	    < (ksMass_eta_min + ksMass_eta_binWidth) ) {
	  (*iETA)->Fill( iV0->mass() );
	}
	ksMass_eta_min += ksMass_eta_binWidth;
      }

      double ksMass_phi_binWidth = 2*M_PI / (double) ksMass_phi_nBins;
      double ksMass_phi_min = -M_PI;
      for( vector<TH1F*>::iterator iPHI = ksMass_phi.begin();
	   iPHI != ksMass_phi.end();
	   iPHI++ ) {
	if( v0DauPhi_lowestpT_forHistos > ksMass_phi_min
	    && v0DauPhi_lowestpT_forHistos 
	    < (ksMass_phi_min + ksMass_phi_binWidth) ) {
	  (*iPHI)->Fill( iV0->mass() );
	}
	ksMass_phi_min += ksMass_phi_binWidth;
      }

      ksMass_phi_binWidth = 2*ksMass_phi_binWidth;
      ksMass_phi_min = -M_PI;
      ksMass_eta_min = -2.5;
      int phiCount = 0;
      int nPhi = ksMass_phi_nBins / 2;
      for( vector<TH1F*>::iterator iEP = ksMass_eta_phi.begin();
	   iEP != ksMass_eta_phi.end();
	   iEP++ ) {
	if( (v0DauPhi_lowestpT_forHistos > ksMass_phi_min
	     && v0DauPhi_lowestpT_forHistos < (ksMass_phi_min + ksMass_phi_binWidth))
	    && (v0DauLowestpTTrack.eta() > ksMass_eta_min
		&& v0DauLowestpTTrack.eta() < (ksMass_eta_min + ksMass_eta_binWidth) ) ) {
	  (*iEP)->Fill( iV0->mass() );
	}
	phiCount++;
	if( phiCount < nPhi ) {
	  ksMass_phi_min += ksMass_phi_binWidth;
	}
	else {
	  ksMass_phi_min = -M_PI;
	  ksMass_eta_min += ksMass_eta_binWidth;
	  phiCount = 0;
	}
      }

      // Need to do pT plots here
      int ptCounter = 0;
      for( vector<TH1F*>::iterator iPT = ksMass_pt.begin();
	   iPT != ksMass_pt.end();
	   iPT++ ) {
	if( ksMass_pt_minima[ptCounter] < 3. ) {
	  if( v0DauLowestpTTrack.pt() > ksMass_pt_minima[ptCounter]
	      && v0DauLowestpTTrack.pt() <= ksMass_pt_minima[ptCounter+1] ) {
	    (*iPT)->Fill( iV0->mass() );
	  }
	}
	else {
	  if( v0DauLowestpTTrack.pt() > 3. ) {
	    (*iPT)->Fill( iV0->mass() );
	  }
	}
	ptCounter++;
      }
    }

    if( abs(iV0->pdgId()) == 3122 
	&& instanceName == string("Lambda") 
	&& writeHistos ) {
      double lamMass_eta_binWidth = 5. / (double) (lamMass_eta_nBins*2);
      double lamMass_eta_min = -2.5;
      for( vector<TH1F*>::iterator iETA = lamMass_eta.begin();
	   iETA != lamMass_eta.end();
	   iETA++ ) {
	if( v0DauLowestpTTrack.eta() > lamMass_eta_min 
	    && v0DauLowestpTTrack.eta() 
	    < (lamMass_eta_min + lamMass_eta_binWidth) ) {
	  (*iETA)->Fill( iV0->mass() );
	}
	lamMass_eta_min += lamMass_eta_binWidth;
      }

      double lamMass_phi_binWidth = 2*M_PI / (double) lamMass_phi_nBins;
      double lamMass_phi_min = -M_PI;
      for( vector<TH1F*>::iterator iPHI = lamMass_phi.begin();
	   iPHI != lamMass_phi.end();
	   iPHI++ ) {
	if( v0DauPhi_lowestpT_forHistos > lamMass_phi_min
	    && v0DauPhi_lowestpT_forHistos 
	    < (lamMass_phi_min + lamMass_phi_binWidth) ) {
	  (*iPHI)->Fill( iV0->mass() );
	}
	lamMass_phi_min += lamMass_phi_binWidth;
      }

      // Need to do pT plots here
      int ptCounter2 = 0;
      for( vector<TH1F*>::iterator iPT = lamMass_pt.begin();
	   iPT != lamMass_pt.end();
	   iPT++ ) {
	if( lamMass_pt_minima[ptCounter2] < 3. ) {
	  if( v0DauLowestpTTrack.pt() > lamMass_pt_minima[ptCounter2]
	      && v0DauLowestpTTrack.pt() <= lamMass_pt_minima[ptCounter2+1] ) {
	    (*iPT)->Fill( iV0->mass() );
	  }
	}
	else {
	  if( v0DauLowestpTTrack.pt() > 3. ) {
	    (*iPT)->Fill( iV0->mass() );
	  }
	}
	ptCounter2++;
      }
    }
    //thePrim = theOrigPrim;
  }

  if( theOrigPrim ) delete theOrigPrim;

  // Fill the tree
  if( putTreeInTFile ) {
    theTree->Fill();
  }

  if( putTreeInTFile ) {
    v0VtxX->clear(); v0VtxY->clear(); v0VtxZ->clear(); v0VtxR->clear(); v0VtxChi2->clear(); v0VtxNormChi2->clear();
    //v0VtxCXX->clear(); v0VtxCYY->clear(); v0VtxCZZ->clear(); v0VtxCYX->clear(); v0VtxCZX->clear(); v0VtxCZY->clear();
    //v0VtxChi2->clear(); v0VtxNormChi2->clear();
    
    v0pX->clear(); v0pY->clear(); v0pZ->clear(); v0pT->clear(); v0p->clear(); v0Eta->clear(); v0Phi->clear();
    v0pTPlus->clear(); v0pTMinus->clear(); v0pLPlus->clear(); v0pLMinus->clear();
    //v0pXerr->clear(); v0pYerr->clear(); v0pZerr->clear();
    
    v0CandMass->clear(); v0OtherCandMass->clear();
    v0CandMassError->clear(); v0OtherCandMassError->clear();
    //v0MassFromFit->clear(); v0MassErrorFromFit->clear();
    //v0OtherMassFromFit->clear(); v0OtherMassErrorFromFit->clear();
    v0Lifetime->clear(); 
    v03DLifetime->clear(); v03DLifetimeError->clear(); v03DLifetimeScaled->clear(); 
    v03DLifetimeScaledError->clear(); v03DLifetimeMPV->clear(); v03DLifetimeMPVError->clear();
    v0PDG->clear(); v0VtxSig->clear(); v0VtxSig3D->clear();
    v0CosThetaPAndLineOfFlight->clear(); v0PCAPrimary->clear(); v0PCAPrimaryError->clear(); 
    v0PCABeamSpot->clear(); v0PCABeamSpotError->clear();
    v3dIpWrtPrimary->clear(); v3dIpWrtPrimaryError->clear(); v3dIpWrtPrimarySig->clear();
    thePosDau3DIp->clear();thePosDau3DIpError->clear();thePosDau3DIpSig->clear();
    theNegDau3DIp->clear();theNegDau3DIpError->clear();theNegDau3DIpSig->clear();
    thePosDauTkQual->clear(); theNegDauTkQual->clear();
    thePosDauTkAlgo->clear(); theNegDauTkAlgo->clear();

    posDauNhits->clear(); negDauNhits->clear();
    //priVtxChi2->clear(); 
    //priVtxNdof->clear(); nTracks->clear(); nV0s->clear();

    v0VtxPriDist->clear(); v0VtxPriDistErr->clear();
    v0VtxBSTransDist->clear(); v0VtxBSTransDistErr->clear();
    
    v0PosDaupT->clear(); v0PosDaup->clear(); v0PosDaupX->clear(); v0PosDaupY->clear(); v0PosDaupZ->clear();
    v0PosDauAtVtxpX->clear(); v0PosDauAtVtxpY->clear(); v0PosDauAtVtxpZ->clear();
    v0PosDauEta->clear(); v0PosDauPhi->clear();
    v0PosDauNValidHits->clear(); v0PosDauNStripHits->clear(); v0PosDauNPixelHits->clear();
    v0PosDauChi2->clear(); v0PosDauNormChi2->clear(); v0PosDauD0->clear();
    
    v0NegDaupT->clear(); v0NegDaup->clear(); v0NegDaupX->clear(); v0NegDaupY->clear(); v0NegDaupZ->clear();
    v0NegDauAtVtxpX->clear(); v0NegDauAtVtxpY->clear(); v0NegDauAtVtxpZ->clear();
    v0NegDauEta->clear(); v0NegDauPhi->clear();
    v0NegDauNValidHits->clear(); v0NegDauNStripHits->clear(); v0NegDauNPixelHits->clear();
    v0NegDauChi2->clear(); v0NegDauNormChi2->clear(); v0NegDauD0->clear();
    
    v0DaupT_lowestpT->clear(); v0DauEta_lowestpT->clear(); v0DauPhi_lowestpT->clear();

    v0DauSumPx->clear(); v0DauSumPy->clear(); v0DauSumPz->clear();

    //nLooseTracks = nHighPurTracks = nLooseTracks_pt_100 = nHighPurTracks_pt_100 = 0;
    nTracks = nTracks_pt_100 = nPromptTracks = nPromptTracks_pt_100 = 0;
    nV0s = 0;
    trigHF = trigTech40_41 = trigTech36_39 = trigTech34 = trigScraping = trigAlgo124 = trigHLTminBias = 0;
    HFpE = HFmE = 0.;
    priVtxChi2 = priVtxNdof = priVtxX = priVtxY = priVtxZ = 0.;
    priVtxNtracks = priVtxIsFake = priVtxNtracksRaw = 0;
    priVtxWasRefit->clear();
    nPriVtxs = 0; nTracksFromOtherPrimary = 0;
    runNum = evtNum = 0;

    nGenParts_eta_24_pt_05 = nGenParts_eta_08_pt_05 = nGenParts_eta_24_pt_01 = 0;

    nGenTracks_hp_eta24_pt05_nHitLay6 = nGenTracks_hp_eta08_pt05_nHitLay6 = 
      nGenTracks_hp_eta24_pt05_nHitLay6_BS = nGenTracks_eta24_pt05 = 
      nGenTracks_eta24_pt05_BS = nTracks_eta24_pt01 = nTracks_eta24_pt01_BS = 0;

    priRefitVtxChi2->clear(); priRefitVtxNdof->clear(); priRefitVtxX->clear(); priRefitVtxY->clear(); priRefitVtxZ->clear();
    v0Rapidity->clear();
  }
}


void V0RecoAnalyzer::beginJob() {
  using namespace std;

  edm::Service<TFileService> fs;

  if( putTreeInTFile ) {
    theTree = fs->make<TTree>("ntuple", "vee ntuple");
    theTree->Branch("v0VtxX", &v0VtxX);
    theTree->Branch("v0VtxY", &v0VtxY);
    theTree->Branch("v0VtxZ", &v0VtxZ);
    theTree->Branch("v0VtxR", &v0VtxR);
    //theTree->Branch("v0VtxChi2", &v0VtxChi2);
    v0VtxChi2 = new std::vector<float>;
    theTree->Branch("v0VtxNormChi2", &v0VtxNormChi2);
  /*theTree->Branch("v0VtxCXX", &v0VtxCXX);
    theTree->Branch("v0VtxCYY", &v0VtxCYY);
    theTree->Branch("v0VtxCZZ", &v0VtxCZZ);
    theTree->Branch("v0VtxCYX", &v0VtxCYX);
    theTree->Branch("v0VtxCZX", &v0VtxCZX);
    theTree->Branch("v0VtxCZY", &v0VtxCZY);
    theTree->Branch("v0VtxChi2", &v0VtxChi2);
    theTree->Branch("v0VtxNormChi2", &v0VtxNormChi2);*/
    
    theTree->Branch("v0pX", &v0pX);
    theTree->Branch("v0pY", &v0pY);
    theTree->Branch("v0pZ", &v0pZ);
    theTree->Branch("v0pT", &v0pT);
    //theTree->Branch("v0p", &v0p);
    v0p = new std::vector<float>;
    theTree->Branch("v0Eta", &v0Eta);
    theTree->Branch("v0Rapidity", &v0Rapidity);
    theTree->Branch("v0Phi", &v0Phi);
  /*theTree->Branch("v0pXerr", &v0pXerr);
    theTree->Branch("v0pYerr", &v0pYerr);
    theTree->Branch("v0pZerr", &v0pZerr);*/

    /*theTree->Branch("v0pTPlus", &v0pTPlus);
    theTree->Branch("v0pTMinus", &v0pTMinus);
    theTree->Branch("v0pLPlus", &v0pLPlus);
    theTree->Branch("v0pLMinus", &v0pLMinus);*/
    v0pTPlus = new std::vector<float>;
    v0pLPlus = new std::vector<float>;
    v0pTMinus = new std::vector<float>;
    v0pLMinus = new std::vector<float>;
  
    theTree->Branch("v0CandMass", &v0CandMass);
    theTree->Branch("v0OtherCandMass", &v0OtherCandMass);
    theTree->Branch("v0CandMassError", &v0CandMassError);
    theTree->Branch("v0OtherCandMassError", &v0OtherCandMassError);

    /*theTree->Branch("v0MassFromFit", &v0MassFromFit);
    theTree->Branch("v0MassErrorFromFit", &v0MassErrorFromFit);
    theTree->Branch("v0OtherMassFromFit", &v0OtherMassFromFit);
    theTree->Branch("v0OtherMassErrorFromFit", &v0OtherMassErrorFromFit);*/

    theTree->Branch("v0Lifetime", &v0Lifetime);

    theTree->Branch("v03DLifetime", &v03DLifetime);
    theTree->Branch("v03DLifetimeError", &v03DLifetimeError);
    theTree->Branch("v03DLifetimeScaled", &v03DLifetimeScaled);
    theTree->Branch("v03DLifetimeScaledError", &v03DLifetimeScaledError);
    theTree->Branch("v03DLifetimeMPV", &v03DLifetimeMPV);
    theTree->Branch("v03DLifetimeMPVError", &v03DLifetimeMPVError);

    theTree->Branch("v0PDG", &v0PDG);

    theTree->Branch("v0VtxSig", &v0VtxSig);
    theTree->Branch("v0VtxSig3D", &v0VtxSig3D);
    
    theTree->Branch("v0CosThetaPAndLineOfFlight", &v0CosThetaPAndLineOfFlight);
    theTree->Branch("v0PCAPrimary", &v0PCAPrimary);
    theTree->Branch("v0PCAPrimaryError", &v0PCAPrimaryError);
    theTree->Branch("v0PCABeamSpot", &v0PCABeamSpot);
    theTree->Branch("v0PCABeamSpotError", &v0PCABeamSpotError);
    theTree->Branch("v3dIpWrtPrimary", &v3dIpWrtPrimary);
    //theTree->Branch("v3dIpWrtPrimaryError", &v3dIpWrtPrimaryError);
    v3dIpWrtPrimaryError = new std::vector<float>;
    theTree->Branch("v3dIpWrtPrimarySig", &v3dIpWrtPrimarySig);
    theTree->Branch("thePosDau3DIp", &thePosDau3DIp);
    //theTree->Branch("thePosDau3DIpError", &thePosDau3DIpError);
    thePosDau3DIpError = new std::vector<float>;
    theTree->Branch("thePosDau3DIpSig", &thePosDau3DIpSig);
    theTree->Branch("theNegDau3DIp", &theNegDau3DIp);
    //theTree->Branch("theNegDau3DIpError", &theNegDau3DIpError);
    theNegDau3DIpError = new std::vector<float>;
    theTree->Branch("theNegDau3DIpSig", &theNegDau3DIpSig);
    theTree->Branch("thePosDauTkQual", &thePosDauTkQual);
    theTree->Branch("theNegDauTkQual", &theNegDauTkQual);
    theTree->Branch("thePosDauTkAlgo", &thePosDauTkAlgo);
    theTree->Branch("theNegDauTkAlgo", &theNegDauTkAlgo);

    theTree->Branch("posDauNhits", &posDauNhits);
    theTree->Branch("negDauNhits", &negDauNhits);

    //theTree->Branch("priVtxChi2", &priVtxChi2);
    //theTree->Branch("priVtxNdof", &priVtxNdof);
    //theTree->Branch("nTracks", &nTracks);
    //theTree->Branch("nV0s", &nV0s);
    theTree->Branch("priVtxChi2", &priVtxChi2, "priVtxChi2/F");
    theTree->Branch("privtxNdof", &priVtxNdof, "priVtxNdof/F");
    theTree->Branch("priVtxX", &priVtxX, "priVtxX/F");
    theTree->Branch("priVtxY", &priVtxY, "priVtxY/F");
    theTree->Branch("priVtxZ", &priVtxZ, "priVtxZ/F");
    theTree->Branch("priVtxNtracks", &priVtxNtracks, "priVtxNtracks/I");
    theTree->Branch("priVtxNtracksRaw", &priVtxNtracksRaw, "priVtxNtracksRaw/I");
    //theTree->Branch("priVtxIsFake", &priVtxIsFake, "priVtxIsFake/I");
    theTree->Branch("priVtxWasRefit", &priVtxWasRefit);
    theTree->Branch("nPriVtxs", &nPriVtxs, "nPriVtxs/I");
    theTree->Branch("nTracksFromOtherPrimary", &nTracksFromOtherPrimary, "nTracksFromOtherPrimary/I");

//     theTree->Branch("nLooseTracks", &nLooseTracks, "nLooseTracks/I");
//     theTree->Branch("nHighPurTracks", &nHighPurTracks, "nHighPurTracks/I");
//     theTree->Branch("nLooseTracks_pt_100", &nLooseTracks_pt_100, "nLooseTracks_pt_100/I");
//     theTree->Branch("nHighPurTracks_pt_100", &nHighPurTracks_pt_100, "nHighPurTracks_pt_100/I");

    theTree->Branch("nTracks", &nTracks, "nTracks/I");
    theTree->Branch("nTracks_pt_100", &nTracks_pt_100, "nTracks_pt_100/I");
    theTree->Branch("nPromptTracks", &nPromptTracks, "nPromptTracks/I");
    theTree->Branch("nPromptTracks_pt_100", &nPromptTracks_pt_100, "nPromptTracks_pt_100/I");

    theTree->Branch("nGenParts_eta_24_pt_05", &nGenParts_eta_24_pt_05, "nGenParts_eta_24_pt_05/I");
    theTree->Branch("nGenParts_eta_08_pt_05", &nGenParts_eta_08_pt_05, "nGenParts_eta_08_pt_05/I");
    theTree->Branch("nGenParts_eta_24_pt_01", &nGenParts_eta_24_pt_01, "nGenParts_eta_24_pt_01/I");

    theTree->Branch("nGenTracks_hp_eta24_pt05_nHitLay6", 
		    &nGenTracks_hp_eta24_pt05_nHitLay6, 
		    "nGenTracks_hp_eta24_pt05_nHitLay6/I");
    theTree->Branch("nGenTracks_hp_eta08_pt05_nHitLay6",
		    &nGenTracks_hp_eta08_pt05_nHitLay6,
		    "nGenTracks_hp_eta08_pt05_nHitLay6/I");
    theTree->Branch("nGenTracks_hp_eta24_pt05_nHitLay6_BS",
		    &nGenTracks_hp_eta24_pt05_nHitLay6_BS,
		    "nGenTracks_hp_eta24_pt05_nHitLay6_BS/I");
    theTree->Branch("nGenTracks_eta24_pt05",
		    &nGenTracks_eta24_pt05,
		    "nGenTracks_eta24_pt05/I");
    theTree->Branch("nGenTracks_eta24_pt05_BS",
		    &nGenTracks_eta24_pt05_BS,
		    "nGenTracks_eta24_pt05_BS/I");
    theTree->Branch("nTracks_eta24_pt01",
		    &nTracks_eta24_pt01,
		    "nTracks_eta24_pt01/I");
    theTree->Branch("nTracks_eta24_pt01_BS",
		    &nTracks_eta24_pt01_BS,
		    "nTracks_eta24_pt01_BS/I");

    theTree->Branch("nV0s", &nV0s, "nV0s/I");
    theTree->Branch("trigHF", &trigHF, "trigHF/I");
    theTree->Branch("HFpE", &HFpE, "HFpE/F");
    theTree->Branch("HFmE", &HFmE, "HFmE/F");
    theTree->Branch("trigTech40_41", &trigTech40_41, "trigTech40_41/I");
    theTree->Branch("trigTech34", &trigTech34, "trigTech34/I");
    theTree->Branch("trigTech36_39", &trigTech36_39, "trigTech36_39/I");
    theTree->Branch("trigScraping", &trigScraping, "trigScraping/I");
    theTree->Branch("trigAlgo124", &trigAlgo124, "trigAlgo124/I");
    theTree->Branch("trigHLTminBias", &trigHLTminBias, "trigHLTminBias/I");

    theTree->Branch("runNum", &runNum, "runNum/I");
    theTree->Branch("evtNum", &evtNum, "evtNum/I");

    theTree->Branch("processType", &processType, "processType/I");

    theTree->Branch("priRefitVtxChi2", &priRefitVtxChi2);
    theTree->Branch("priRefitVtxNdof", &priRefitVtxNdof);
    theTree->Branch("priRefitVtxX", &priRefitVtxX);
    theTree->Branch("priRefitVtxY", &priRefitVtxY);
    theTree->Branch("priRefitVtxZ", &priRefitVtxZ);

    theTree->Branch("v0VtxPriDist", &v0VtxPriDist);
    theTree->Branch("v0VtxPriDistErr", &v0VtxPriDistErr);
    theTree->Branch("v0VtxBSTransDist", &v0VtxBSTransDist);
    theTree->Branch("v0VtxBSTransDistErr", &v0VtxBSTransDistErr);
  
    //theTree->Branch("v0PosDaupT", &v0PosDaupT);
    //theTree->Branch("v0PosDaup", &v0PosDaup);
    v0PosDaupT = new std::vector<float>;
    v0PosDaup = new std::vector<float>;
    theTree->Branch("v0PosDaupX", &v0PosDaupX);
    theTree->Branch("v0PosDaupY", &v0PosDaupY);
    theTree->Branch("v0PosDaupZ", &v0PosDaupZ);
    theTree->Branch("v0PosDauAtVtxpX", &v0PosDauAtVtxpX);
    theTree->Branch("v0PosDauAtVtxpY", &v0PosDauAtVtxpY);
    theTree->Branch("v0PosDauAtVtxpZ", &v0PosDauAtVtxpZ);
    theTree->Branch("v0PosDauEta", &v0PosDauEta);
    //theTree->Branch("v0PosDauPhi", &v0PosDauPhi);
    v0PosDauPhi = new std::vector<float>;
    //theTree->Branch("v0PosDauNValidHits", &v0PosDauNValidHits);
    v0PosDauNValidHits = new std::vector<float>;
    theTree->Branch("v0PosDauNPixelHits", &v0PosDauNPixelHits);
    theTree->Branch("v0PosDauNStripHits", &v0PosDauNStripHits);
    theTree->Branch("v0PosDauChi2", &v0PosDauChi2);
    theTree->Branch("v0PosDauNormChi2", &v0PosDauNormChi2);
    theTree->Branch("v0PosDauD0", &v0PosDauD0);

    theTree->Branch("v0DaupT_lowestpT", &v0DaupT_lowestpT);
    
    //theTree->Branch("v0NegDaupT", &v0NegDaupT);
    //theTree->Branch("v0NegDaup", &v0NegDaup);
    v0NegDaupT = new std::vector<float>;
    v0NegDaup = new std::vector<float>;
    theTree->Branch("v0NegDaupX", &v0NegDaupX);
    theTree->Branch("v0NegDaupY", &v0NegDaupY);
    theTree->Branch("v0NegDaupZ", &v0NegDaupZ);
    theTree->Branch("v0NegDauAtVtxpX", &v0NegDauAtVtxpX);
    theTree->Branch("v0NegDauAtVtxpY", &v0NegDauAtVtxpY);
    theTree->Branch("v0NegDauAtVtxpZ", &v0NegDauAtVtxpZ);
    theTree->Branch("v0NegDauEta", &v0NegDauEta);
    //theTree->Branch("v0NegDauPhi", &v0NegDauPhi);
    v0NegDauPhi = new std::vector<float>;
    //theTree->Branch("v0NegDauNValidHits", &v0NegDauNValidHits);
    v0NegDauNValidHits = new std::vector<float>;
    theTree->Branch("v0NegDauNPixelHits", &v0NegDauNPixelHits);
    theTree->Branch("v0NegDauNStripHits", &v0NegDauNStripHits);
    theTree->Branch("v0NegDauChi2", &v0NegDauChi2);
    theTree->Branch("v0NegDauNormChi2", &v0NegDauNormChi2);
    theTree->Branch("v0NegDauD0", &v0NegDauD0);

    //theTree->Branch("v0DauSumPx", &v0DauSumPx);
    //theTree->Branch("v0DauSumPy", &v0DauSumPy);
    //theTree->Branch("v0DauSumPz", &v0DauSumPz);
    v0DauSumPx = new std::vector<float>;
    v0DauSumPy = new std::vector<float>;
    v0DauSumPz = new std::vector<float>;

    theTree->Branch("v0DauEta_lowestpT", &v0DauEta_lowestpT);
    theTree->Branch("v0DauPhi_lowestpT", &v0DauPhi_lowestpT);
  }

    ksMassBinWidth = (ksMassXmax - ksMassXmin) / (double) ksMassNbins;
    lamMassBinWidth = (lamMassXmax - lamMassXmin) / (double) lamMassNbins;

  if( writeHistos ) {
    if( instanceName == string("Kshort") ) {
      ksR = fs->make<TH1F>("ksR", "K^{0}_{S} radial distance from beam line",
			   50, 0., 40.);
      ksRFull = fs->make<TH1F>("ksRFull", "K^{0}_{S} radial distance from beam line",
			       50, 0., 40.);
      ksRSignal = fs->make<TH1F>("ksRSignal", "K^{0}_{S} radial distance from beam line",
				 50, 0., 40.);
      ksRBg = fs->make<TH1F>("ksRBg", "K^{0}_{S} radial distance from beam line",
			     50, 0., 40.);
      
      ksZ = fs->make<TH1F>("ksZ", "K^{0}_{S} z position",
			   50, -150., 150.);
      ksZFull = fs->make<TH1F>("ksZFull", "K^{0}_{S} z position",
			       50, -150., 150.);
      ksZSignal = fs->make<TH1F>("ksZSignal", "K^{0}_{S} z position",
				 50, -150., 150.);
      ksZBg = fs->make<TH1F>("ksZBg", "K^{0}_{S} z position",
			     50, -150., 150.);

      ksZvsRFull = fs->make<TH2F>("ksZvsRFull", "K^{0}_{S} z position vs. decay radius",
				  120, -120., 120., 40, 0., 40.);
      ksZvsRBgSubt = fs->make<TH2F>("ksZvsRBgSubt", "K^{0}_{S} z position vs. decay radius",
				    120, -120., 120., 40, 0., 40.);
      ksZvsRBg = fs->make<TH2F>("ksZvsRBg", "K^{0}_{S} z position vs. decay radius",
				120, -120., 120., 40, 0., 40.);
      ksZvsRSignal = fs->make<TH2F>("ksZvsRSignal", "K^{0}_{S} z position vs. decay radius",
				    120, -120., 120., 40, 0., 40.);
    
      ksEta = fs->make<TH1F>("ksEta", "K^{0}_{S} momentum #eta",
			     40, -2.5, 2.5);
      ksEtaFull = fs->make<TH1F>("ksEtaFull", "K^{0}_{S} momentum #eta",
				 40, -2.5, 2.5);
      ksEtaSignal = fs->make<TH1F>("ksEtaSignal", "K^{0}_{S} momentum #eta",
				   40, -2.5, 2.5);
      ksEtaBg = fs->make<TH1F>("ksEtaBg", "K^{0}_{S} momentum #eta",
			       40, -2.5, 2.5);
    
      ksPhi = fs->make<TH1F>("ksPhi", "K^{0}_{S} momentum #phi",
			     60, -M_PI, M_PI);
      ksPhiFull = fs->make<TH1F>("ksPhiFull", "K^{0}_{S} momentum #phi",
				 60, -M_PI, M_PI);
      ksPhiSignal = fs->make<TH1F>("ksPhiSignal", "K^{0}_{S} momentum #phi",
				   60, -M_PI, M_PI);
      ksPhiBg = fs->make<TH1F>("ksPhiBg", "K^{0}_{S} momentum #phi",
			       60, -M_PI, M_PI);
      
      ksPt = fs->make<TH1F>("ksPt", "K^{0}_{S} p_{T}",
			    70, 0., 20.);
      ksPtFull = fs->make<TH1F>("ksPtFull", "K^{0}_{S} p_{T}",
				70, 0., 20.);
      ksPtSignal = fs->make<TH1F>("ksPtSignal", "K^{0}_{S} p_{T}",
				  70, 0., 20.);
      ksPtBg = fs->make<TH1F>("ksPtBg", "K^{0}_{S} p_{T}",
			      70, 0., 20.);

      ksCandMass = fs->make<TH1F>("ksCandMass", "K^{0}_{S} mass from Candidate",
				  ksMassNbins, ksMassXmin, ksMassXmax);
	
      ksDaup = fs->make<TH1F>("ksDaup", "p of K^{0}_{S} daughter tracks",
			      100, 0., 30.);
      ksDaupT = fs->make<TH1F>("ksDaupT", "p_{T} of K^{0}_{S} daughter tracks",
			       100, 0., 20.);
      ksDauEta = fs->make<TH1F>("ksDauEta", "#eta of K^{0}_{S} daughter tracks",
				100, -2.5, 2.5);
      ksDauPhi = fs->make<TH1F>("ksDauPhi", "#phi of K^{0}_{S} daughter track at vertex position",
				70, -M_PI, M_PI);
      ksDauNValidHits = fs->make<TH1F>("ksDauNValHits", 
				       "Number of valid hits on K^{0}_{S} daughter tracks", 
				       80, 0., 40.);
      ksDauNPixelHits = fs->make<TH1F>("ksDauNPixelHits", 
				       "Number of valid pixel hits on K^{0}_{S} daughter tracks",
				       30 , 0., 15.);
      ksDauNStripHits = fs->make<TH1F>("ksDauNStripHits", 
				       "Number of valid strip hits on K^{0}_{S} daughter tracks",
				       80, 0., 40.);
      ksDauChi2 = fs->make<TH1F>("ksDauChi2", "#chi^{2} of K^{0}_{S} daughter tracks",
				 100, 0., 100.);
      ksDauNormChi2 = fs->make<TH1F>("ksDauNormChi2", "#chi^{2}/ndof of K^{0}_{S} daughter tracks",
				     100, 0., 5.);
      ksDauD0 = fs->make<TH1F>("ksDauD0", "abs(d_{0}) of K^{0}_{S} daughter tracks",
			       100, 0., 15.);

      double ksMass_eta_binWidth = 5. / (double) (ksMass_eta_nBins*2);
      for( double ksMass_eta_min = -2.5;
	   ksMass_eta_min < 2.5;
	   ksMass_eta_min += ksMass_eta_binWidth ) {
	ostringstream etaNameShort;
	ostringstream etaNameLong;
	etaNameShort << "ksMass_eta_" << ksMass_eta_min;
	etaNameLong << "K^{0}_{S} invariant mass " << ksMass_eta_min
		    << " < #eta < " << ksMass_eta_min + ksMass_eta_binWidth;
	/*ksMass_eta.push_back( fs->make<TH1F>((string("ksMass_eta_") 
					      + string(ksMass_eta_min)).c_str(),
					     (string("K^{0}_{S} invariant mass eta min ")
					     + string(ksMass_eta_min)).c_str(),*/
	ksMass_eta.push_back( fs->make<TH1F>(etaNameShort.str().c_str(),
					     etaNameLong.str().c_str(),
					     ksMassNbins, ksMassXmin, ksMassXmax) );
	ksMass_eta_names_short.push_back( etaNameShort.str() );
	ksMass_eta_names_long.push_back( etaNameLong.str() );
	ksMass_eta_minima.push_back( ksMass_eta_min );
      }
      ksMassInEtaBins = fs->make<TH1F>("KsMassBiasVsEta", 
				       "K^{0}_{S} mass bias vs. #eta of lowest-p_{T} track",
				       ksMass_eta_nBins*2, -2.5, 2.5); 

      double ksMass_phi_binWidth = 2*M_PI / (double) ksMass_phi_nBins;
      for( double ksMass_phi_min = -M_PI;
	   ksMass_phi_min < M_PI - 0.0001;
	   ksMass_phi_min += ksMass_phi_binWidth ) {
	ostringstream phiNameShort;
	ostringstream phiNameLong;
	phiNameShort << "ksMass_phi_" << ksMass_phi_min;
	phiNameLong << "K^{0}_{S} invariant mass " << ksMass_phi_min
		    << " < #phi < " << ksMass_phi_min + ksMass_phi_binWidth;
	/*ksMass_phi.push_back( fs->make<TH1F>((string("ksMass_phi_")
					      + string(ksMass_phi_min)).c_str(),
					     (string("K^{0}_{S} invariant mass phi min ")
					     + string(ksMass_phi_min)).c_str(),*/
	ksMass_phi.push_back( fs->make<TH1F>(phiNameShort.str().c_str(),
					     phiNameLong.str().c_str(),
					     ksMassNbins, ksMassXmin, ksMassXmax) );
	ksMass_phi_names_short.push_back( phiNameShort.str() );
	ksMass_phi_names_long.push_back( phiNameLong.str() );
	ksMass_phi_minima.push_back( ksMass_phi_min );
      }
      ksMassInPhiBins = fs->make<TH1F>("KsMassBiasVsPhi", 
				       "K^{0}_{S} mass bias vs. #phi of lowest-p_{T} track",
				       ksMass_phi_nBins, -M_PI, M_PI); 

      for(double ksMass_eta_min = -2.5;
	  ksMass_eta_min < 2.5;
	  ksMass_eta_min += ksMass_eta_binWidth ) {
	for(double ksMass_phi_min = -M_PI;
	    ksMass_phi_min < M_PI - 0.0001;
	    ksMass_phi_min += ksMass_phi_binWidth*2 ) {
	  ostringstream shortNames2D;
	  ostringstream longNames2D;
	  shortNames2D << "ksMass_2D_" << ksMass_eta_min << "_" << ksMass_phi_min;
	  longNames2D << "Ks invariant mass, " << ksMass_eta_min << " < #eta < "
		      << ksMass_eta_min + ksMass_eta_binWidth << ", "
		      << ksMass_phi_min << " < #phi < " << ksMass_phi_min + ksMass_phi_binWidth;
	  ksMass_eta_phi.push_back( fs->make<TH1F>(shortNames2D.str().c_str(),
						   longNames2D.str().c_str(),
						   ksMassNbins/2, ksMassXmin, ksMassXmax) );
	  ksMass_eta_phi_names_short.push_back( shortNames2D.str() );
	  ksMass_eta_phi_names_long.push_back( longNames2D.str() );
	  ksMass_eta_phi_etaminima.push_back( ksMass_eta_min );
	  ksMass_eta_phi_phiminima.push_back( ksMass_phi_min );
	}
      }
      ksMassBiasEtaVsPhi = fs->make<TH2F>("KsMassBiasEtaVsPhi",
					  "K^{0}_{S} mass bias #eta vs. #phi",
					  ksMass_phi_nBins/2,
					  -M_PI, M_PI,
					  ksMass_eta_nBins*2,
					  -2.5, 2.5);

      double ksMass_pt_binWidth1 = 0.25 / (double) ksMass_pt_nBinSubdiv;//0. < pT < 0.25
      double ksMass_pt_binWidth2 = ksMass_pt_binWidth1;//0.25 < pT < 0.5
      double ksMass_pt_binWidth3 = 0.5 / (double) ksMass_pt_nBinSubdiv;//0.5 < pT < 1.
      double ksMass_pt_binWidth4 = ksMass_pt_binWidth3;//1. < pT < 1.5
      double ksMass_pt_binWidth5 = ksMass_pt_binWidth3;//1.5 < pT < 2.
      double ksMass_pt_binWidth6 = 1.0 / (double) ksMass_pt_nBinSubdiv;//2. < pT < 3.
      double ksMass_pt_min = 0.;
      //cout << "Making pT histos" << endl;
      while( ksMass_pt_min <= 3. ) {
	//cout << "ksMass_pt_min: " << ksMass_pt_min << endl;
	ostringstream ptNameShort;
	ostringstream ptNameLong;
	ptNameShort << "ksMass_pt_" << ksMass_pt_min;
	ptNameLong << "K^{0}_{S} invariant mass, p_{T} bin minimum = "
		   << ksMass_pt_min;
	ksMass_pt.push_back( fs->make<TH1F>(ptNameShort.str().c_str(),
					    ptNameLong.str().c_str(),
					    ksMassNbins, ksMassXmin, ksMassXmax) );
	ksMass_pt_names_short.push_back( ptNameShort.str() );
	ksMass_pt_names_long.push_back( ptNameLong.str() );
	ksMass_pt_minima.push_back( ksMass_pt_min );
	//cout << "Minima size is now " << ksMas
	if( ksMass_pt_min < 0.25
	    && (ksMass_pt_min + ksMass_pt_binWidth1) <= 0.25)  {
	  ksMass_pt_min += ksMass_pt_binWidth1;
	}
	else if( ksMass_pt_min < 0.25
		 && (ksMass_pt_min + ksMass_pt_binWidth1) > 0.25 ) {
	  ksMass_pt_min += ksMass_pt_binWidth2;
	}
	else if( ksMass_pt_min < 0.5
		 && (ksMass_pt_min + ksMass_pt_binWidth2) <= 0.5 ) {
	  ksMass_pt_min += ksMass_pt_binWidth2;
	}
	else if( ksMass_pt_min < 0.5
		 && (ksMass_pt_min + ksMass_pt_binWidth2) > 0.5 ) {
	  ksMass_pt_min += ksMass_pt_binWidth3;
	}
	else if( ksMass_pt_min < 1.
		 && (ksMass_pt_min + ksMass_pt_binWidth3) <= 1. ) {
	  ksMass_pt_min += ksMass_pt_binWidth3;
	}
	else if( ksMass_pt_min < 1.
		 && (ksMass_pt_min + ksMass_pt_binWidth3) > 1. ) {
	  ksMass_pt_min += ksMass_pt_binWidth4;
	}
	else if( ksMass_pt_min < 1.5
		 && (ksMass_pt_min + ksMass_pt_binWidth4) <= 1.5 ) {
	  ksMass_pt_min += ksMass_pt_binWidth4;
	}
	else if( ksMass_pt_min < 1.5
		 && (ksMass_pt_min + ksMass_pt_binWidth4) > 1.5 ) {
	  ksMass_pt_min += ksMass_pt_binWidth5;
	}
	else if( ksMass_pt_min < 2.
		 && (ksMass_pt_min + ksMass_pt_binWidth5) <= 2. ) {
	  ksMass_pt_min += ksMass_pt_binWidth5;
	}
	else if( ksMass_pt_min < 2.
		 && (ksMass_pt_min + ksMass_pt_binWidth5) > 2. ) {
	  ksMass_pt_min += ksMass_pt_binWidth6;
	}
	else if( ksMass_pt_min < 3.
		 && (ksMass_pt_min + ksMass_pt_binWidth6) <= 3. ) {
	  ksMass_pt_min += ksMass_pt_binWidth6;
	}
	else if( ksMass_pt_min < 3.) {
	  ksMass_pt_min = 3.;
	}
	else {
	  ksMass_pt_min = 4.;
	}

      }

      vector<float> minimaForTH1;
      unsigned int ndx = 0;
      for( vector<double>::iterator iMIN = ksMass_pt_minima.begin();
	   iMIN != ksMass_pt_minima.end();
	   iMIN++ ) {
	minimaForTH1.push_back(*iMIN);
	if( ndx == ksMass_pt_minima.size() - 1 ) {
	  //cout << "Pushed back a 15. " << endl;
	  minimaForTH1.push_back( 6. );
	}
	ndx++;
      }
      ksMassInPtBins = fs->make<TH1F>("KsMassBiasVsPt",
				      "K^{0}_{S} mass bias vs. p_{T} of lowest-p_{T} track",
				      ksMass_pt_minima.size(), &minimaForTH1[0]);

    }
    if( instanceName == string("Lambda") ) {
      lamR = fs->make<TH1F>("lamR", "#Lambda^{0} radial distance from beam line",
			    50, 0., 40.);
      lamRFull = fs->make<TH1F>("lamRFull", "#Lambda^{0} radial distance from beam line",
				50, 0., 40.);
      lamRSignal = fs->make<TH1F>("lamRSignal", "#Lambda^{0} radial distance from beam line",
				  50, 0., 40.);
      lamRBg = fs->make<TH1F>("lamRBg", "#Lambda^{0} radial distance from beam line",
			      50, 0., 40.);
      
      lamZ = fs->make<TH1F>("lamZ", "#Lambda^{0} z position",
			    50, -150., 150.);
      lamZFull = fs->make<TH1F>("lamZFull", "#Lambda^{0} z position",
				50, -150., 150.);
      lamZSignal = fs->make<TH1F>("lamZSignal", "#Lambda^{0} z position",
				  50, -150., 150.);
      lamZBg = fs->make<TH1F>("lamZBg", "#Lambda^{0} z position",
			      50, -150., 150.);

      lamZvsRFull = fs->make<TH2F>("lamZvsRFull", "#Lambda^{0} z position vs. decay radius",
				   120, -120., 120., 40, 0., 40.);
      lamZvsRBgSubt = fs->make<TH2F>("lamZvsRBgSubt", "#Lambda^{0} z position vs. decay radius",
				     120, -120., 120., 40, 0., 40.);
      lamZvsRBg = fs->make<TH2F>("lamZvsRBg", "#Lambda^{0} z position vs. decay radius",
				 120, -120., 120., 40, 0., 40.);
      lamZvsRSignal = fs->make<TH2F>("lamZvsRSignal", "#Lambda^{0} z position vs. decay radius",
				     120, -120., 120., 40, 0., 40.);
    
      lamEta = fs->make<TH1F>("lamEta", "#Lambda^{0} momentum #eta",
			      40, -2.5, 2.5);
      lamEtaFull = fs->make<TH1F>("lamEtaFull", "#Lambda^{0} momentum #eta",
				  40, -2.5, 2.5);
      lamEtaSignal = fs->make<TH1F>("lamEtaSignal", "#Lambda^{0} momentum #eta",
				    40, -2.5, 2.5);
      lamEtaBg = fs->make<TH1F>("lamEtaBg", "#Lambda^{0} momentum #eta",
				40, -2.5, 2.5);
      
      lamPhi = fs->make<TH1F>("lamPhi", "#Lambda^{0} momentum #phi",
			      60, -M_PI, M_PI);
      lamPhiFull = fs->make<TH1F>("lamPhiFull", "#Lambda^{0} momentum #phi",
				  60, -M_PI, M_PI);
      lamPhiSignal = fs->make<TH1F>("lamPhiSignal", "#Lambda^{0} momentum #phi",
				    60, -M_PI, M_PI);
      lamPhiBg = fs->make<TH1F>("lamPhiBg", "#Lambda^{0} momentum #phi",
				60, -M_PI, M_PI);
    
      lamPt = fs->make<TH1F>("lamPt", "#Lambda^{0} p_{T}",
			     70, 0., 20.);
      lamPtFull = fs->make<TH1F>("lamPtFull", "#Lambda^{0} p_{T}",
				 70, 0., 20.);
      lamPtSignal = fs->make<TH1F>("lamPtSignal", "#Lambda^{0} p_{T}",
				   70, 0., 20.);
      lamPtBg = fs->make<TH1F>("lamPtBg", "#Lambda^{0} p_{T}",
			       70, 0., 20.);

      lamCandMass = fs->make<TH1F>("lamCandMass", "#Lambda^{0} mass from Candidate",
				   lamMassNbins, lamMassXmin, lamMassXmax);
      
      lamDaup = fs->make<TH1F>("lamDaup", "p of #Lambda^{0} daughter tracks",
			       100, 0., 30.);
      lamDaupT = fs->make<TH1F>("lamDaupT", "p_{T} of #Lambda^{0} daughter tracks",
				100, 0., 20.);
      lamDauEta = fs->make<TH1F>("lamDauEta", "#eta of #Lambda^{0} daughter tracks",
				 100, -2.5, 2.5);
      lamDauPhi = fs->make<TH1F>("lamDauPhi", "#phi of #Lambda^{0} daughter track at vertex position",
				 70, -M_PI, M_PI);
      lamDauNValidHits = fs->make<TH1F>("lamDauNValHits", 
					"Number of valid hits on #Lambda^{0} daughter tracks", 
					80, 0., 40.);
      lamDauNPixelHits = fs->make<TH1F>("lamDauNPixelHits", 
					"Number of valid pixel hits on #Lambda^{0} daughter tracks",
					30 , 0., 15.);
      lamDauNStripHits = fs->make<TH1F>("lamDauNStripHits", 
					"Number of valid strip hits on #Lambda^{0} daughter tracks",
					80, 0., 40.);
      lamDauChi2 = fs->make<TH1F>("lamDauChi2", "#chi^{2} of #Lambda^{0} daughter tracks",
				  100, 0., 100.);
      lamDauNormChi2 = fs->make<TH1F>("lamDauNormChi2", "#chi^{2}/ndof of #Lambda^{0} daughter tracks",
				      100, 0., 5.);
      lamDauD0 = fs->make<TH1F>("lamDauD0", "abs(d_{0}) of #Lambda^{0} daughter tracks",
				100, 0., 15.);

      double lamMass_eta_binWidth = 5. / (double) (lamMass_eta_nBins*2);
      for( double lamMass_eta_min = -2.5;
	   lamMass_eta_min < 2.5;
	   lamMass_eta_min += lamMass_eta_binWidth ) {
	ostringstream etaNameShort;
	ostringstream etaNameLong;
	etaNameShort << "lamMass_eta_" << lamMass_eta_min;
	etaNameLong << "#Lambda^{0} invariant mass " << lamMass_eta_min
		    << " < #eta < " << lamMass_eta_min + lamMass_eta_binWidth;
	/*lamMass_eta.push_back( fs->make<TH1F>((string("lamMass_eta_") 
					       + string(lamMass_eta_min)).c_str(),
					      (string("#Lambda^{0} invariant mass eta min ")
					      + string(lamMass_eta_min)).c_str(),*/
	lamMass_eta.push_back( fs->make<TH1F>(etaNameShort.str().c_str(),
					      etaNameLong.str().c_str(),
					      lamMassNbins, lamMassXmin, lamMassXmax) );
	lamMass_eta_names_short.push_back( etaNameShort.str() );
	lamMass_eta_names_long.push_back( etaNameLong.str() );
	lamMass_eta_minima.push_back( lamMass_eta_min );
      }
      lamMassInEtaBins = fs->make<TH1F>("LamMassBiasVsEta", 
					"#Lambda^{0} mass bias vs. #eta of lowest-p_{T} track",
					lamMass_eta_nBins*2, -2.5, 2.5); 

      double lamMass_phi_binWidth = 2*M_PI / (double) lamMass_phi_nBins;
      for( double lamMass_phi_min = -M_PI;
	   lamMass_phi_min < M_PI - 0.001;
	   lamMass_phi_min += lamMass_phi_binWidth ) {
	ostringstream phiNameShort;
	ostringstream phiNameLong;
	phiNameShort << "lamMass_phi_" << lamMass_phi_min;
	phiNameLong << "#Lambda^{0} invariant mass " << lamMass_phi_min
		    << " < #phi < " << lamMass_phi_min + lamMass_phi_binWidth;
	/*lamMass_phi.push_back( fs->make<TH1F>((string("lamMass_phi_")
					       + string(lamMass_phi_min)).c_str(),
					      (string("#Lambda^{0} invariant mass phi min ")
					      + string(lamMass_phi_min)).c_str(),*/
	lamMass_phi.push_back( fs->make<TH1F>(phiNameShort.str().c_str(),
					      phiNameLong.str().c_str(),
					      lamMassNbins, lamMassXmin, lamMassXmax) );
	lamMass_phi_names_short.push_back( phiNameShort.str() );
	lamMass_phi_names_long.push_back( phiNameLong.str() );
	lamMass_phi_minima.push_back( lamMass_phi_min );
      }
      lamMassInPhiBins = fs->make<TH1F>("LamMassBiasVsPhi",
					"#Lambda^{0} mass bias vs. #phi of lowest-p_{T} track",
					lamMass_phi_nBins, -M_PI, M_PI);

      double lamMass_pt_binWidth1 = 0.25 / (double) lamMass_pt_nBinSubdiv;//0. < pT < 0.25
      double lamMass_pt_binWidth2 = lamMass_pt_binWidth1;//0.25 < pT < 0.5
      double lamMass_pt_binWidth3 = 0.5 / (double) lamMass_pt_nBinSubdiv;//0.5 < pT < 1.
      double lamMass_pt_binWidth4 = lamMass_pt_binWidth3;//1. < pT < 1.5
      double lamMass_pt_binWidth5 = lamMass_pt_binWidth3;//1.5 < pT < 2.
      double lamMass_pt_binWidth6 = 1.0 / (double) lamMass_pt_nBinSubdiv;//2. < pT < 3.
      double lamMass_pt_min = 0.;
      //cout << "Making pT histos" << endl;
      while( lamMass_pt_min <= 3. ) {
	//cout << "lamMass_pt_min: " << lamMass_pt_min << endl;
	ostringstream ptNameShort;
	ostringstream ptNameLong;
	ptNameShort << "lamMass_pt_" << lamMass_pt_min;
	ptNameLong << "#Lambda^{0} invariant mass, p_{T} bin minimum = "
		   << lamMass_pt_min;
	lamMass_pt.push_back( fs->make<TH1F>(ptNameShort.str().c_str(),
					    ptNameLong.str().c_str(),
					    lamMassNbins, lamMassXmin, lamMassXmax) );
	lamMass_pt_names_short.push_back( ptNameShort.str() );
	lamMass_pt_names_long.push_back( ptNameLong.str() );
	lamMass_pt_minima.push_back( lamMass_pt_min );
	//cout << "Minima size is now " << lamMas
	if( lamMass_pt_min < 0.25
	    && (lamMass_pt_min + lamMass_pt_binWidth1) <= 0.25)  {
	  lamMass_pt_min += lamMass_pt_binWidth1;
	}
	else if( lamMass_pt_min < 0.25
		 && (lamMass_pt_min + lamMass_pt_binWidth1) > 0.25 ) {
	  lamMass_pt_min += lamMass_pt_binWidth2;
	}
	else if( lamMass_pt_min < 0.5
		 && (lamMass_pt_min + lamMass_pt_binWidth2) <= 0.5 ) {
	  lamMass_pt_min += lamMass_pt_binWidth2;
	}
	else if( lamMass_pt_min < 0.5
		 && (lamMass_pt_min + lamMass_pt_binWidth2) > 0.5 ) {
	  lamMass_pt_min += lamMass_pt_binWidth3;
	}
	else if( lamMass_pt_min < 1.
		 && (lamMass_pt_min + lamMass_pt_binWidth3) <= 1. ) {
	  lamMass_pt_min += lamMass_pt_binWidth3;
	}
	else if( lamMass_pt_min < 1.
		 && (lamMass_pt_min + lamMass_pt_binWidth3) > 1. ) {
	  lamMass_pt_min += lamMass_pt_binWidth4;
	}
	else if( lamMass_pt_min < 1.5
		 && (lamMass_pt_min + lamMass_pt_binWidth4) <= 1.5 ) {
	  lamMass_pt_min += lamMass_pt_binWidth4;
	}
	else if( lamMass_pt_min < 1.5
		 && (lamMass_pt_min + lamMass_pt_binWidth4) > 1.5 ) {
	  lamMass_pt_min += lamMass_pt_binWidth5;
	}
	else if( lamMass_pt_min < 2.
		 && (lamMass_pt_min + lamMass_pt_binWidth5) <= 2. ) {
	  lamMass_pt_min += lamMass_pt_binWidth5;
	}
	else if( lamMass_pt_min < 2.
		 && (lamMass_pt_min + lamMass_pt_binWidth5) > 2. ) {
	  lamMass_pt_min += lamMass_pt_binWidth6;
	}
	else if( lamMass_pt_min < 3.
		 && (lamMass_pt_min + lamMass_pt_binWidth6) <= 3. ) {
	  lamMass_pt_min += lamMass_pt_binWidth6;
	}
	else if( lamMass_pt_min < 3.) {
	  lamMass_pt_min = 3.;
	}
	else {
	  lamMass_pt_min = 4.;
	}

      }

      vector<float> minimaForTH1;
      unsigned int ndx = 0;
      for( vector<double>::iterator iMIN = lamMass_pt_minima.begin();
	   iMIN != lamMass_pt_minima.end();
	   iMIN++ ) {
	minimaForTH1.push_back(*iMIN);
	if( ndx == lamMass_pt_minima.size() - 1 ) {
	  //cout << "Pushed back a 15. " << endl;
	  minimaForTH1.push_back( 6. );
	}
	ndx++;
      }
      lamMassInPtBins = fs->make<TH1F>("LamMassBiasVsPt",
				      "#Lambda^{0} mass bias vs. p_{T} of lowest-p_{T} track",
				      lamMass_pt_minima.size(), &minimaForTH1[0]);
      
    }
  }
}

void V0RecoAnalyzer::endJob() {
  using namespace std;

  if( writeHistos && instanceName == string("Kshort") ) {
    ostringstream ksoss;
    ostringstream ksoss_doubgaus;
    ksoss << ksMassBinWidth << "*gausn(0) + [3] + [4]*(x-0.49767)"
	  << " + [5]*(x-0.49767)^2";
    // p0*exp(-0.5*((x-p1)/p2)^2)
    /*ksoss_doubgaus << ksMassBinWidth 
		   << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2) + "
		   << ksMassBinWidth 
		   << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2) + "
		   << "[5] + [6]*(x-0.49767)";*/
    double ksMassBinWidth_DG = ksMassBinWidth / sqrt(2*3.141592654);
    ksoss_doubgaus << ksMassBinWidth_DG 	
		   << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		   << ksMassBinWidth_DG 
		   << "*[0]*(1 - [3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		   << "[5] + [6]*(x-[1]) + [7]*(x-[1])^2";

    TF1 *ksFit_doubGaus = new TF1("ksFit_doubGaus", 
				  ksoss_doubgaus.str().c_str(),
				  ksMassXmin, ksMassXmax);
    ksFit_doubGaus->SetParName( 0, "yield" );
    ksFit_doubGaus->SetParName( 1, "mean" );
    ksFit_doubGaus->SetParName( 2, "sigma1" );
    ksFit_doubGaus->SetParName( 3, "fraction" );
    ksFit_doubGaus->SetParLimits( 3, 0., 1. );
    ksFit_doubGaus->SetParName( 4, "sigma2" );
    ksFit_doubGaus->SetParName( 5, "const" );
    ksFit_doubGaus->SetParName( 6, "slope" );
    ksFit_doubGaus->SetParName( 7, "quadconst" );

    ksFit_doubGaus->SetParameter( 0, 4.4444*ksCandMass->GetMaximum() );
    ksFit_doubGaus->SetParameter( 1, 0.49767 );
    ksFit_doubGaus->SetParameter( 2, 0.0045 );
    ksFit_doubGaus->SetParameter( 3, 0.75 );
    ksFit_doubGaus->SetParameter( 4, 0.011 );
    ksFit_doubGaus->SetParameter( 5, 0.1*ksCandMass->GetMaximum() );
    ksFit_doubGaus->SetParameter( 6, 1. );
    ksFit_doubGaus->SetParameter( 7, -1. );

    TF1 *ksFit = new TF1("ksFit", ksoss.str().c_str(), ksMassXmin, ksMassXmax);
    ksFit->SetParName( 0, "yield" );
    ksFit->SetParName( 1, "mean" );
    ksFit->SetParName( 2, "sigma" );
    ksFit->SetParName( 3, "const" );
    ksFit->SetParName( 4, "slope" );
    ksFit->SetParName( 5, "quadconst" );

    ksFit->SetParameter( 0, 4.4444*ksCandMass->GetMaximum() );
    ksFit->SetParameter( 1, 0.49767 );
    ksFit->SetParameter( 2, 0.005 );
    ksFit->SetParLimits( 2, 0., 0.010 );
    ksFit->SetParameter( 3, 0.1*ksCandMass->GetMaximum() );
    ksFit->SetParameter( 4, 1. );
    ksFit->SetParameter( 5, -1. );

    /*
    TF1 *ksFit_doubGaus = new TF1("ksFitDG", ksoss_doubgaus.str().c_str(),
				  ksMassXmin, ksMassXmax);
    ksFit_doubGaus->SetParName( 0, "yield" );
    ksFit_doubGaus->SetParName( 1, "mean" );
    ksFit_doubGaus->SetParName( 2, "sigma1" );
    ksFit_doubGaus->SetParName( 3, "fraction" );
    ksFit_doubGaus->SetParName( 4, "sigma2" );
    ksFit_doubGaus->SetParName( 5, "const" );
    ksFit_doubGaus->SetParName( 6, "slope" );

    ksFit_doubGaus->SetParameter( 0, 4.4444*ksCandMass->GetMaximum() );
    ksFit_doubGaus->SetParameter( 1, 0.49767 );
    ksFit_doubGaus->SetParameter( 2, 0.005 );
    ksFit_doubGaus->SetParameter( 3, 0.9 );
    ksFit_doubGaus->SetParameter( 4, 0.007 );
    ksFit_doubGaus->SetParameter( 5, 0.1*ksCandMass->GetMaximum() );
    ksFit_doubGaus->SetParameter( 6, 0. );
    */

    ksR->Add( ksRSignal, ksRBg, 1., -0.5 );
    ksZ->Add( ksZSignal, ksZBg, 1., -0.5 );
    ksZvsRBgSubt->Add( ksZvsRFull, ksZvsRBg, 1., -0.5);
    ksEta->Add( ksEtaSignal, ksEtaBg, 1., -0.5 );
    ksPhi->Add( ksPhiSignal, ksPhiBg, 1., -0.5 );
    ksPt->Add( ksPtSignal, ksPtBg, 1., -0.5 );

    ksCandMass->Fit("ksFit_doubGaus", "RLE");

    // Now do fits for mass bias plots
    int ksEtaCount = 0;
    for( vector<TH1F*>::iterator iETA = ksMass_eta.begin();
	 iETA != ksMass_eta.end();
	 iETA++ ) {
      //cout << "---> Fitting " << ksMass_eta_names_short[ksEtaCount] << endl;
      ksFit->SetParameter( 0, 4.4444*(*iETA)->GetMaximum() );
      ksFit->SetParameter( 1, 0.49767 );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*(*iETA)->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = (*iETA)->Fit("ksFit", "RLE");
      if( !fitStatus ) {
	ksMassInEtaBins->SetBinContent( ksEtaCount + 1, ksFit->GetParameter(1) );
	ksMassInEtaBins->SetBinError( ksEtaCount + 1, ksFit->GetParError(1) );
      }
      else {
	ksMassInEtaBins->SetBinContent( ksEtaCount + 1, ksMassConst );
	ksMassInEtaBins->SetBinError( ksEtaCount + 1, 0. );
      }
      ksEtaCount++;
    }

    //cout << "Doing Ks phi, count = " << ksMass_phi.size() << endl;
    int ksPhiCount = 0;
    for( vector<TH1F*>::iterator iPHI = ksMass_phi.begin();
	 iPHI != ksMass_phi.end();
	 iPHI++ ) {
      cout << "---> Fitting " << ksMass_phi_names_short[ksPhiCount] << endl;
      ksFit->SetParameter( 0, 4.4444*(*iPHI)->GetMaximum() );
      ksFit->SetParameter( 1, 0.49767 );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*(*iPHI)->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = (*iPHI)->Fit("ksFit", "RLE");
      if( !fitStatus ) {
	ksMassInPhiBins->SetBinContent( ksPhiCount + 1, ksFit->GetParameter(1) );
	ksMassInPhiBins->SetBinError( ksPhiCount + 1, ksFit->GetParError(1) );
      }
      else {
	ksMassInPhiBins->SetBinContent( ksPhiCount + 1, ksMassConst );
	ksMassInPhiBins->SetBinError( ksPhiCount + 1, 0. );
      }
      ksPhiCount++;
    }

    ksEtaCount = ksPhiCount = 0;
    //int numEta = ksMass_eta_nBins * 2;
    int numPhi = ksMass_phi_nBins / 2;
    for( vector<TH1F*>::iterator iEP = ksMass_eta_phi.begin();
	 iEP != ksMass_eta_phi.end();
	 iEP++ ) {
      ksFit->SetParameter( 0, 4.4444*(*iEP)->GetMaximum() );
      ksFit->SetParameter( 1, 0.49767 );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*(*iEP)->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = (*iEP)->Fit("ksFit", "RLE");
      if( !fitStatus ) {
	ksMassBiasEtaVsPhi->SetBinContent( ksPhiCount + 1, 
					   ksEtaCount + 1, ksFit->GetParameter(1) );
	ksMassBiasEtaVsPhi->SetBinError( ksPhiCount + 1,
					 ksEtaCount + 1, ksFit->GetParError(1) );
      }
      else {
	ksMassBiasEtaVsPhi->SetBinContent( ksPhiCount + 1,
					ksEtaCount + 1, 0.497648 );
	ksMassBiasEtaVsPhi->SetBinError( ksPhiCount + 1,
					 ksEtaCount + 1, 0. );
      }
      ksPhiCount++;
      if( ksPhiCount == numPhi ) {
	ksEtaCount++;
	ksPhiCount = 0;
      }
    }


    int ksPtCount = 0;
    for( vector<TH1F*>::iterator iPT = ksMass_pt.begin();
	 iPT != ksMass_pt.end();
	 iPT++ ) {
      //cout << "---> Fitting " << ksMass_pt_names_short[ksPtCount] << endl;
      ksFit->SetParameter( 0, 4.4444*(*iPT)->GetMaximum() );
      ksFit->SetParameter( 1, 0.49767 );
      ksFit->SetParameter( 2, 0.005 );
      ksFit->SetParameter( 3, 0.1*(*iPT)->GetMaximum() );
      ksFit->SetParameter( 4, 1. );
      ksFit->SetParameter( 5, -1. );

      int fitStatus = (*iPT)->Fit("ksFit", "RLE");
      if( !fitStatus ) {
	ksMassInPtBins->SetBinContent( ksPtCount + 1, ksFit->GetParameter(1) );
	ksMassInPtBins->SetBinError( ksPtCount + 1, ksFit->GetParError(1) );
      }
      else {
	ksMassInPtBins->SetBinContent( ksPtCount + 1, ksMassConst );
	ksMassInPtBins->SetBinError( ksPtCount + 1, 0. );
      }
      ksPtCount++;
    }
  }

  if( writeHistos && instanceName == string("Lambda") ) {
    ostringstream lamoss;
    lamoss << lamMassBinWidth << "*gausn(0) + [3]*(x - "
	   << piMass + protonMass << ")^(1/2) + [4]*(x - "
	   << piMass + protonMass << ")^(3/2)";
    ostringstream lamfitoss_doubGaus;
    double lamMassBinWidth_DG = lamMassBinWidth/sqrt(2*3.141592654);
    lamfitoss_doubGaus << lamMassBinWidth_DG
		       << "*[0]*[3]*exp(-0.5*((x-[1])/[2])^2)/[2] + "
		       << lamMassBinWidth_DG
		       << "*[0]*(1-[3])*exp(-0.5*((x-[1])/[4])^2)/[4] + "
		       << "[5]*(x - " << piMass + protonMass << ")^(1/2) + [6]*(x-"
		       << piMass + protonMass << ")^(3/2)";
    TF1 *lamFit_doubGaus = new TF1("lamFit_doubGaus", lamfitoss_doubGaus.str().c_str(),
				   piMass+protonMass, lamMassXmax);
    lamFit_doubGaus->SetParName( 0, "yield" );
    lamFit_doubGaus->SetParName( 1, "mean" );
    lamFit_doubGaus->SetParName( 2, "sigma1" );
    lamFit_doubGaus->SetParName( 3, "fraction" );
    lamFit_doubGaus->SetParLimits( 3, 0., 1. );
    lamFit_doubGaus->SetParName( 4, "sigma2" );
    lamFit_doubGaus->SetParName( 5, "sqrt" );
    lamFit_doubGaus->SetParName( 6, "sqrtcube" );
    lamFit_doubGaus->SetParameter( 0, 8.*lamCandMass->GetMaximum() );
    lamFit_doubGaus->SetParameter( 1, 1.1159 );
    lamFit_doubGaus->SetParameter( 2, 0.0016 );
    lamFit_doubGaus->SetParameter( 3, 0.44 );
    lamFit_doubGaus->SetParameter( 4, 0.0037);
    lamFit_doubGaus->SetParameter( 5, 0.001*lamCandMass->GetMaximum()/lamMassBinWidth );
    lamFit_doubGaus->SetParameter( 6, 0.001*lamCandMass->GetMaximum()/lamMassBinWidth );

    TF1 *lamFit = new TF1("lamFit", lamoss.str().c_str(), piMass+protonMass, lamMassXmax);
    lamFit->SetParName( 0, "yield" );
    lamFit->SetParName( 1, "GausMean" );
    lamFit->SetParName( 2, "GausSigma" );
    lamFit->SetParName( 3, "c1" );
    lamFit->SetParName( 4, "c2" );

    lamFit->SetParameter( 0, 2.6*lamCandMass->GetMaximum() );
    lamFit->SetParameter( 1, 1.1156 );
    lamFit->SetParameter( 2, 0.002 );
    lamFit->SetParameter( 3, 0.75*lamCandMass->GetMaximum() );
    lamFit->SetParameter( 4, 1.85*lamCandMass->GetMaximum() );

    lamR->Add( lamRSignal, lamRBg, 1., -0.5 );
    lamZ->Add( lamZSignal, lamZBg, 1., -0.5 );
    lamZvsRBgSubt->Add( lamZvsRFull, lamZvsRBg, 1., -0.5);
    lamEta->Add( lamEtaSignal, lamEtaBg, 1., -0.5 );
    lamPhi->Add( lamPhiSignal, lamPhiBg, 1., -0.5 );
    lamPt->Add( lamPtSignal, lamPtBg, 1., -0.5 );

    lamCandMass->Fit("lamFit_doubGaus", "RLE");

    int lamEtaCount = 0;
    for( vector<TH1F*>::iterator iETA = lamMass_eta.begin();
	 iETA != lamMass_eta.end();
	 iETA++ ) {
      //cout << "---> Fitting " << lamMass_eta_names_short[lamEtaCount] << endl;
      lamFit->SetParameter( 0, 2.6*(*iETA)->GetMaximum() );
      lamFit->SetParameter( 1, 1.1156 );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*(*iETA)->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*(*iETA)->GetMaximum() );

      int fitStatus = (*iETA)->Fit("lamFit", "RLE");
      if( !fitStatus ) {
	lamMassInEtaBins->SetBinContent( lamEtaCount + 1, lamFit->GetParameter(1) );
	lamMassInEtaBins->SetBinError( lamEtaCount + 1, lamFit->GetParError(1) );
      }
      else {
	lamMassInEtaBins->SetBinContent( lamEtaCount + 1, lamMassConst );
	lamMassInEtaBins->SetBinError( lamEtaCount + 1, 0. );
      }
      lamEtaCount++;
    }

    int lamPhiCount = 0;
    for( vector<TH1F*>::iterator iPHI = lamMass_phi.begin();
	 iPHI != lamMass_phi.end();
	 iPHI++ ) {
      //cout << "---> Fitting " << lamMass_phi_names_short[lamPhiCount] << endl;
      lamFit->SetParameter( 0, 2.6*(*iPHI)->GetMaximum() );
      lamFit->SetParameter( 1, 1.1156 );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*(*iPHI)->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*(*iPHI)->GetMaximum() );

      int fitStatus = (*iPHI)->Fit("lamFit", "RLE");
      if( !fitStatus ) {
	lamMassInPhiBins->SetBinContent( lamPhiCount + 1, lamFit->GetParameter(1) );
	lamMassInPhiBins->SetBinError( lamPhiCount + 1, lamFit->GetParError(1) );
      }
      else {
	lamMassInPhiBins->SetBinContent( lamPhiCount + 1, lamMassConst );
	lamMassInPhiBins->SetBinError( lamPhiCount + 1, 0. );
      }
      lamPhiCount++;
    }

    lamEtaCount = lamPhiCount = 0;

    int lamPtCount = 0;
    for( vector<TH1F*>::iterator iPT = lamMass_pt.begin();
	 iPT != lamMass_pt.end();
	 iPT++ ) {
      //cout << "---> Fitting " << lamMass_pt_names_short[lamPtCount] << endl;
      lamFit->SetParameter( 0, 2.6*(*iPT)->GetMaximum() );
      lamFit->SetParameter( 1, 1.1156 );
      lamFit->SetParameter( 2, 0.002 );
      lamFit->SetParameter( 3, 0.75*(*iPT)->GetMaximum() );
      lamFit->SetParameter( 4, 1.85*(*iPT)->GetMaximum() );

      int fitStatus = (*iPT)->Fit("lamFit", "RLE");
      //cout << "Fit returned " << isGoodFit << endl;
      if( !fitStatus ) {
	lamMassInPtBins->SetBinContent( lamPtCount + 1, lamFit->GetParameter(1) );
	lamMassInPtBins->SetBinError( lamPtCount + 1, lamFit->GetParError(1) );
      }
      else {
	lamMassInPtBins->SetBinContent( lamPtCount + 1, lamMassConst );
	lamMassInPtBins->SetBinError( lamPtCount + 1, 0. );
      }
      lamPtCount++;
    }
  }
}

bool V0RecoAnalyzer::aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers)
{
 int negTowers = 0;
 int posTowers = 0;
 for(CaloTowerCollection::const_iterator cal = towers->begin(); cal != towers->end(); ++cal) {
    for(unsigned int i = 0; i < cal->constituentsSize(); i++) {
       const DetId id = cal->constituent(i);
       if(id.det() == DetId::Hcal) {
         HcalSubdetector subdet=(HcalSubdetector(id.subdetId()));
         if(subdet == HcalForward) {
// 	   if (fabs(cal->eta()) < 3.)
// 	     cout << "cal->eta() = " << cal->eta() << endl;
           if(cal->eta()<-3.) {
	     if (cal->energy() > 3.)
	       negTowers++;
	     if (cal->energy() > HFmE)
	       HFmE = cal->energy();
	   }
           if(cal->eta()>3.) {
	     if (cal->energy() > 3.)
	       posTowers++;
	     if (cal->energy() > HFpE)
	       HFpE = cal->energy();
	   }
         }
      }
    }
 }
 if(negTowers>0 && posTowers>0)
   return true;

 return false;
}

bool V0RecoAnalyzer::scrapingTest( const edm::Event& iEvent, 
					const edm::EventSetup& iSetup ) {
  bool accepted;
  float fraction;
  unsigned int numtrack = 10;
  float thresh = 0.2;
  reco::TrackBase::TrackQuality _trackQuality;
  // get GeneralTracks collection

  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks",tkRef);    
  const reco::TrackCollection* tkColl = tkRef.product();

  //std::cout << "Total Number of Tracks " << tkColl->size() << endl;
  
  int numhighpurity=0;
  int numhighpurity_pt_100 = 0;
  int numloose_pt_100 = 0;
  _trackQuality = reco::TrackBase::qualityByName("highPurity");

  edm::Handle<reco::BeamSpot> bsH;
  iEvent.getByLabel("offlineBeamSpot", bsH);

  edm::Handle<reco::VertexCollection> thePrims;
 // iEvent.getByLabel("offlinePrimaryVerticesWithMBTracks", thePrims);
  iEvent.getByLabel("offlinePrimaryVertices",thePrims);	//by Qiao
  reco::TrackBase::Point priVtxPos(thePrims->begin()->x(),
				   thePrims->begin()->y(),
				   thePrims->begin()->z());
  reco::TrackBase::Point bsPos(bsH->position().x(), bsH->position().y(), bsH->position().z());

  reco::TrackCollection::const_iterator itk = tkColl->begin();
  reco::TrackCollection::const_iterator itk_e = tkColl->end();
  for(;itk!=itk_e;++itk){
    // std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
    if(itk->quality(_trackQuality)) numhighpurity++;
    if(itk->quality(_trackQuality) && itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numhighpurity_pt_100++;
    if(itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numloose_pt_100++;
    if(itk->quality(_trackQuality) 
       && fabs(itk->eta()) < 2.4 
       && itk->pt() > 0.5 
       && itk->hitPattern().trackerLayersWithMeasurement() >= 6
       && itk->dxy(priVtxPos) < 0.2
       && itk->dz(priVtxPos) < 0.2) 
      nGenTracks_hp_eta24_pt05_nHitLay6++;
    if(itk->quality(_trackQuality)
       && fabs(itk->eta()) < 0.8
       && itk->pt() > 0.5
       && itk->hitPattern().trackerLayersWithMeasurement() >= 6
       && itk->dxy(priVtxPos) < 0.2
       && itk->dz(priVtxPos) < 0.2) 
      nGenTracks_hp_eta08_pt05_nHitLay6++;
    if(itk->quality(_trackQuality)
       && fabs(itk->eta()) < 2.4
       && itk->pt() > 0.5
       && itk->hitPattern().trackerLayersWithMeasurement() >= 6
       && itk->dxy(bsPos) < 0.2)
      nGenTracks_hp_eta24_pt05_nHitLay6_BS++;
    if(fabs(itk->eta()) < 2.4
       && itk->pt() > 0.1
       && itk->dxy(priVtxPos) < 0.2
       && itk->dz(priVtxPos) < 0.2)
      nGenTracks_eta24_pt05++;
    if(fabs(itk->eta()) < 2.4
       && itk->pt() > 0.1
       && itk->dxy(bsPos) < 0.2)
      nGenTracks_eta24_pt05_BS++;
  }
  fraction = (float)numhighpurity/(float)tkColl->size();

  accepted = (tkColl->size()<=numtrack)||(fraction>thresh);

  edm::Handle<reco::TrackCollection> tkH;
  edm::ESHandle<TransientTrackBuilder> theTTB;
  edm::ESHandle<MagneticField> bFieldH;
  edm::ESHandle<GlobalTrackingGeometry> globTkGeomH;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTB);
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomH);
  try{
 //   iEvent.getByLabel("generalPlusMinBiasTracks", tkH);
    iEvent.getByLabel("generalTracks", tkH);
}
  catch(...) {
//    std::cout << "Problem finding track collection." << endl;
    std::cout << "Problem finding track collection." << std::endl;	//by Qiao
    return accepted;
    return accepted;
  }

  using namespace reco;
  TSCBLBuilderNoMaterial blsBuilder;
  for (unsigned int ndx = 0; ndx < tkH->size(); ndx++) {
    TrackRef tmpRef(tkH, ndx);
    TransientTrack tmpTk( *tmpRef, &(*bFieldH), globTkGeomH);
//    TrajectoryStateTransform theTransform;
//    FreeTrajectoryState initialFTS = theTransform.initialFreeState(*tmpRef, bFieldH.product());
    FreeTrajectoryState initialFTS = trajectoryStateTransform::initialFreeState(*tmpRef, bFieldH.product());
    TrajectoryStateClosestToBeamLine tscb(blsBuilder(initialFTS, *bsH));
//    TrajectoryStateClosestToBeamLine tscb(blsBuilder(initialFTS, *bsH));
//    TrajectoryStateClosestToBeamLine tscb(blsBuilder(initialFTS, *bsH));
    if (tscb.isValid()) {
      double d0 = tscb.transverseImpactParameter().value();
      if (fabs(d0) < 0.4) {
	nPromptTracks++;
	if (tmpRef->pt() > 0.1 && fabs(tmpRef->eta()) < 2.5) nPromptTracks_pt_100++;
      }
    }
    nTracks++;
    if (tmpRef->pt() > 0.1 && fabs(tmpRef->eta()) < 2.5) nTracks_pt_100++;

    if (fabs(tmpRef->eta()) < 2.4
	&& tmpRef->pt() > 0.1
	&& tmpRef->dxy(priVtxPos) < 0.2
	&& tmpRef->dz(priVtxPos) < 0.2)
      nTracks_eta24_pt01++;
    if (fabs(tmpRef->eta()) < 2.4
	&& tmpRef->pt() > 0.1
	&& tmpRef->dxy(bsPos) < 0.2)
      nTracks_eta24_pt01_BS++;
      
      
  }

//   nLooseTracks = tkColl->size();
//   nHighPurTracks = numhighpurity;
//   nLooseTracks_pt_100 = numloose_pt_100;
//   nHighPurTracks_pt_100 = numhighpurity_pt_100;




  /*if (debugOn) {
    int ievt = iEvent.id().event();
    int irun = iEvent.id().run();
    int ils = iEvent.luminosityBlock();
    int bx = iEvent.bunchCrossing();
    
    cout << "FilterOutScraping_debug: Run " << irun << " Event " << ievt << " Lumi Block " << ils << " Bunch Crossing " << bx << " Fraction " << fraction << " NTracks " << tkColl->size() << " Accepted " << accepted << endl;
    }*/
 
  //if (applyfilter)
    return accepted;
    //else
    //return true;
}


/*
bool V0RecoAnalyzer::scrapingTest( const edm::Event& iEvent, 
				   const edm::EventSetup& iSetup ) {
  bool accepted = false;
  float fraction = 0;  

  // These are the parameters that are passed to the noScraping filter.
  unsigned int numtrack = 10;
  double thresh = 0.2;

  reco::TrackBase::TrackQuality _trackQuality;
  // get GeneralTracks collection

  edm::Handle<reco::TrackCollection> tkRef;
  iEvent.getByLabel("generalTracks",tkRef);    
  const reco::TrackCollection* tkColl = tkRef.product();

  //std::cout << "Total Number of Tracks " << tkColl->size() << endl;
  
  int numhighpurity=0;
  int numhighpurity_pt_100 = 0;
  int numloose_pt_100 = 0;
  _trackQuality = reco::TrackBase::qualityByName("highPurity");

  if(tkColl->size()>numtrack){ 
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      // std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
      if(itk->quality(_trackQuality)) numhighpurity++;
      if(itk->quality(_trackQuality) && itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numhighpurity_pt_100++;
      if(itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numloose_pt_100++;
    }
    fraction = (float)numhighpurity/(float)tkColl->size();
    if(fraction>thresh) accepted=true;

    nLooseTracks = tkColl->size();
    nHighPurTracks = numhighpurity;
    nLooseTracks_pt_100 = numloose_pt_100;
    nHighPurTracks_pt_100 = numhighpurity_pt_100;
  }
  else{
    //if less than 10 Tracks accept the event anyway    
    accepted= true;
    reco::TrackCollection::const_iterator itk = tkColl->begin();
    reco::TrackCollection::const_iterator itk_e = tkColl->end();
    for(;itk!=itk_e;++itk){
      // std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
      if(itk->quality(_trackQuality)) numhighpurity++;
      if(itk->quality(_trackQuality) && itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numhighpurity_pt_100++;
      if(itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numloose_pt_100++;
    }
    nLooseTracks = tkColl->size();
    nHighPurTracks = numhighpurity;
    nLooseTracks_pt_100 = numloose_pt_100;
    nHighPurTracks_pt_100 = numhighpurity_pt_100;
  }
  
  
//   if (debugOn) {
//     int ievt = iEvent.id().event();
//     int irun = iEvent.id().run();
//     int ils = iEvent.luminosityBlock();
//     int bx = iEvent.bunchCrossing();
    
//     cout << "FilterOutScraping_debug: Run " << irun << " Event " << ievt << " Lumi Block " << ils << " Bunch Crossing " << bx << " Fraction " << fraction << " NTracks " << tkColl->size() << " Accepted " << accepted << endl;
//     }
 
  //if (applyfilter)
    return accepted;
    //else
    //return true;
}
*/

DEFINE_FWK_MODULE(V0RecoAnalyzer);
