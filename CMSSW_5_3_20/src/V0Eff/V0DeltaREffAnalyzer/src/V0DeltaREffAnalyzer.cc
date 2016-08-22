// -*- C++ -*-
//
// Package:    V0DeltaREffAnalyzer
// Class:      V0DeltaREffAnalyzer
// 
/**\class V0DeltaREffAnalyzer V0DeltaREffAnalyzer.cc V0Eff/V0DeltaREffAnalyzer/src/V0DeltaREffAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Brian Drell
//         Created:  Thu Apr 22 15:12:27 MDT 2010
// $Id: V0DeltaREffAnalyzer.cc,v 1.27 2010/10/01 19:20:10 drell Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

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

#include "TMath.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TSystem.h"

//
// class declaration
//

class V0DeltaREffAnalyzer : public edm::EDAnalyzer {
public:
  explicit V0DeltaREffAnalyzer(const edm::ParameterSet&);
  ~V0DeltaREffAnalyzer();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers);
  bool scrapingTest(const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------

  edm::InputTag hltTag;

  InvariantMassFromVertex massCalculator;

  //edm::InputTag ksCollection;
  //edm::InputTag lamCollection;
  edm::InputTag v0Collection;
  edm::InputTag pvCollection;
  double pideltaRcut;
  double pdeltaRcut;
  double deltaLcut;

  bool kshort;

  // Tree
  TTree* eff;
  //TTree* lameff;

  // Branches, single variable
  int trigHF;
  int trigTech40_41;
  int trigScraping;
  int trigTech36_39;
  int trigTech34;
  int trigAlgo124;

  int trigHLTminBias;
  float HFpE;
  float HFmE;

//   int nLooseTracks;
//   int nHighPurTracks;
//   int nLooseTracks_pt_100;
//   int nHighPurTracks_pt_100;

  int nTracks;
  int nTracks_pt_100;
  int nPromptTracks;
  int nPromptTracks_pt_100;

  int nGenParts;
  int nGenParts_pt_100;

  int runNum;
  int evtNum;
  int processType;

  float priVtxX;
  float priVtxY;
  float priVtxZ;
  float priVtxChi2;
  float priVtxNdof;
  int priVtxNtracks;
  int priVtxIsFake;
  int nPriVtxs;

  int priVtxWasRefit;
  float priRefitVtxX;
  float priRefitVtxY;
  float priRefitVtxZ;
  float priRefitVtxChi2;
  float priRefitVtxNdof;

  // Branch vectors
  float genEta;
  float genRapidity;
  float genpT;
  float genp;
  float genPhi;
  float genRho;
  float genR;
  float genMass;
  float genCtau;

  float gendauMaxEta;

  float recoEta;
  float recoRapidity;
  float recopT;
  float recop;
  float recoPhi;
  float recoRho;
  float recoR;
  float recoMass;
  float recoCtau;

  //float recoMassFromFit;
  //float recoMassErrorFromFit;
  float recoMassError;

  // NEW
  float recoVtxX;
  float recoVtxY;
  float recoVtxZ;
  float recoVtxSig;
  float recoVtxSig3D;
  float recoVtxChi2;
  float recoVtxNormChi2;
  float recopX;
  float recopY;
  float recopZ;
  float recoOtherCandMass;

  //float recoOtherMassFromFit;
  //float recoOtherMassErrorFromFit;
  float recoOtherCandMassError;

  float recoPosDauChi2;
  float recoPosDauNormChi2;
  float recoPosDauD0;
  float recoNegDauChi2;
  float recoNegDauNormChi2;
  float recoNegDauD0;

  float recoPosDauPx;
  float recoPosDauPy;
  float recoPosDauPz;
  float recoNegDauPx;
  float recoNegDauPy;
  float recoNegDauPz;

  float v3dIpWrtPrimary;
  float v3dIpWrtPrimaryError;
  float v3dIpWrtPrimarySig;

  float posDau3DIp;
  float posDau3DIpError;
  float posDau3DIpSig;
  float negDau3DIp;
  float negDau3DIpError;
  float negDau3DIpSig;

  int posDauNhits;
  int negDauNhits;

  // OLD
  float posDeltaR;
  float posDeltaEta;
  float posDeltaPhi;
  float pospT;
  float posp;
  float negDeltaR;
  float negDeltaEta;
  float negDeltaPhi;
  float negpT;
  float negp;
  float deltaL;

  float genProdVtxX;
  float genProdVtxY;
  float genProdVtxZ;

  float genDecVtxX;
  float genDecVtxY;
  float genDecVtxZ;

  int v0MatchStatus;

  int genPDG;
  int genPosDauPDG;
  int genNegDauPDG;

  float genPosDauPx;
  float genPosDauPy;
  float genPosDauPz;
  float genNegDauPx;
  float genNegDauPy;
  float genNegDauPz;


  bool debug;


};

//
// constants, enums and typedefs
//

const double ksMassConst = 0.497614;
const double ksMassConst2 = ksMassConst * ksMassConst;
const double lamMassConst = 1.115683;
const double lamMassConst2 = lamMassConst * lamMassConst;
const double piMassConst = 0.13957018;
const double protonMassConst = 0.93827203;
float pionSigma = piMassConst*1.e-6;
float protonSigma = protonMassConst*1.e-6;

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

enum VeeTypes{none=0, ks=1, lam=2};

//
// static data member definitions
//

//
// constructors and destructor
//
V0DeltaREffAnalyzer::V0DeltaREffAnalyzer(const edm::ParameterSet& iConfig):
  //ksCollection(iConfig.getParameter<edm::InputTag>("ksCollection")),
  //lamCollection(iConfig.getParameter<edm::InputTag>("lamCollection")),
  hltTag(iConfig.getParameter<edm::InputTag>("hltTag")),
  v0Collection(iConfig.getParameter<edm::InputTag>("v0Collection")),
  pvCollection(iConfig.getUntrackedParameter<edm::InputTag>("pvCollection", edm::InputTag("offlinePrimaryVertices"))),
  pideltaRcut(iConfig.getParameter<double>("pideltaRcut")),
  pdeltaRcut(iConfig.getParameter<double>("pdeltaRcut")),
  deltaLcut(iConfig.getParameter<double>("deltaLcut")),
  kshort(iConfig.getParameter<bool>("isKshort")),
  trigHF(0), trigTech40_41(0), trigScraping(0), trigTech36_39(0), trigTech34(0), trigAlgo124(0),
  trigHLTminBias(0), HFpE(0), HFmE(0),
  //nLooseTracks(0), nHighPurTracks(0), nLooseTracks_pt_100(0), nHighPurTracks_pt_100(0),
  nTracks(0), nTracks_pt_100(0), nPromptTracks(0), nPromptTracks_pt_100(0),
  nGenParts(0), nGenParts_pt_100(0),
  runNum(0), evtNum(0), processType(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxChi2(0), priVtxNdof(0), priVtxNtracks(0), priVtxIsFake(0),
  nPriVtxs(0),
  genEta(0), genRapidity(0), 
  genpT(0), genp(0), genPhi(0), genRho(0), genR(0), 
  genMass(0), genCtau(0),
  gendauMaxEta(0),
  recoEta(0), recoRapidity(0),
  recopT(0), recop(0), recoPhi(0), recoRho(0), recoR(0),
  recoMass(0), recoCtau(0), recoMassError(0),
  //recoMassFromFit(0), recoMassErrorFromFit(0),
  recoVtxX(0), recoVtxY(0), recoVtxZ(0), recoVtxSig(0), recoVtxChi2(0), recoVtxNormChi2(0),
  recopX(0), recopY(0), recopZ(0), recoOtherCandMass(0), recoOtherCandMassError(0),
  //recoOtherMassFromFit(0), recoOtherMassErrorFromFit(0),
  recoPosDauChi2(0), recoPosDauNormChi2(0), recoPosDauD0(0),
  recoNegDauChi2(0), recoNegDauNormChi2(0), recoNegDauD0(0),
  recoPosDauPx(0), recoPosDauPy(0), recoPosDauPz(0),
  recoNegDauPx(0), recoNegDauPy(0), recoNegDauPz(0),
  v3dIpWrtPrimary(0), v3dIpWrtPrimaryError(0), v3dIpWrtPrimarySig(0),
  posDau3DIp(0), posDau3DIpError(0), posDau3DIpSig(0),
  negDau3DIp(0), negDau3DIpError(0), negDau3DIpSig(0),
  posDauNhits(0), negDauNhits(0),
  posDeltaR(0), posDeltaEta(0), posDeltaPhi(0),
  pospT(0), posp(0), 
  negDeltaR(0), negDeltaEta(0), negDeltaPhi(0),
  negpT(0), negp(0), deltaL(0),
  genProdVtxX(0), genProdVtxY(0), genProdVtxZ(0),
  genDecVtxX(0), genDecVtxY(0), genDecVtxZ(0),
  v0MatchStatus(0),
  genPDG(0), genPosDauPDG(0), genNegDauPDG(0),
  genPosDauPx(0), genPosDauPy(0), genPosDauPz(0),
  genNegDauPx(0), genNegDauPy(0), genNegDauPz(0),

  debug(false) {

}


V0DeltaREffAnalyzer::~V0DeltaREffAnalyzer() {

}


//
// member functions
//



// ------------ method called to for each event  ------------


void V0DeltaREffAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  runNum = iEvent.id().run();
  runNum = iEvent.id().event();

  Handle<edm::HepMCProduct> hepEv;
  iEvent.getByType( hepEv );
  processType = hepEv->GetEvent()->signal_process_id();

  Handle<pat::CompositeCandidateCollection> patH;
  Handle<reco::VertexCollection> thePrimH;
  Handle<reco::GenParticleCollection> genPH;
  Handle<reco::BeamSpot> theBeamSpotH;

  iEvent.getByLabel("genParticlesPlusSim", genPH);
  iEvent.getByLabel(v0Collection, patH);
  iEvent.getByLabel(pvCollection, thePrimH);
  iEvent.getByLabel(string("offlineBeamSpot"), theBeamSpotH);

  int n = 0;
  int n_pt_100 = 0;
  for (unsigned int ndx = 0; ndx < genPH->size(); ndx++) {
    const Candidate &theCand = (*genPH)[ndx];
    int pdg = abs(theCand.pdgId());
    if (!(pdg == 211
	  || pdg == 321
	  || pdg == 2212
	  || pdg == 11
	  || pdg == 13)) continue;
    double x = theCand.vertex().x();
    double y = theCand.vertex().y();
    double z = theCand.vertex().z();
    double r = sqrt(x*x + y*y);
    if (r > 3 || fabs(z) > 30) continue;
    n++;
    if (theCand.pt() > 0.1 && theCand.eta() < 2.5) n_pt_100++;
  }
  nGenParts = n;
  nGenParts_pt_100 = n_pt_100;

  ESHandle<TransientTrackBuilder> theTTB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTB);

  // Check HLT_L1_BscMinBiasOR_BptxPlusORMinus trigger for prescale lumi calculation
  edm::Handle<edm::TriggerResults> hltResults;
  try {
    iEvent.getByLabel(hltTag, hltResults);
  }
  catch( ... ) {
    cout << "Couldn't get HLT trigger." << endl;
    return;
  }
  if( !hltResults.isValid() ) {
    cout << "No trigger results." << endl;
  }
  else {
    int ntrigs = hltResults->size();
    //edm::TriggerNames trigNames(*hltResults);
    //trigNames.init(*hltResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);

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
   
  //check Tech MinBias trigger bit (40||41) (trigTech40_41) 
  //  and also NOT (36||37||38||39) (Tech2)
  edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
  iEvent.getByLabel("gtDigis", gtRecord);
  if( gtRecord.isValid() ) {
    TechnicalTriggerWord tw = gtRecord->technicalTriggerWord();
    if( !tw.empty() ) {
      if( tw[40] || tw[41] ) 
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

  bool isNotScraping = scrapingTest(iEvent, iSetup);
  if (isNotScraping) trigScraping = 1;
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

  GlobalPoint beamSpotPos( theBeamSpotH->position().x(),
			   theBeamSpotH->position().y(),
			   theBeamSpotH->position().z() );

  nPriVtxs = thePrimH->size();

  reco::VertexCollection::const_iterator bestVtx = thePrimH->begin();
  if( !bestVtx->isValid() )//|| bestVtx->isFake() ) 
    return;

  if( bestVtx->isFake() ) priVtxIsFake = 1;
  else priVtxIsFake = 0;
  priVtxX = bestVtx->position().x();
  priVtxY = bestVtx->position().y();
  priVtxZ = bestVtx->position().z();
  priVtxChi2 = bestVtx->chi2();
  priVtxNdof = bestVtx->ndof();
  int nVtxTracks = 0;
  for( Vertex::trackRef_iterator itk = bestVtx->tracks_begin();
       itk != bestVtx->tracks_end();
       itk++ ) {
    if( bestVtx->trackWeight( *itk ) > 0.5 )
      nVtxTracks++;
  }
  priVtxNtracks = nVtxTracks;

  reco::Vertex* thePrim = 0;
  reco::Vertex* theOrigPrim = 0;

  thePrim = new reco::Vertex(*bestVtx);
  theOrigPrim = thePrim;


  bool writeTree = false;

  // Search GenParticles for K0S, lambda and match daughters
  const Candidate* genPosDau = 0;
  const Candidate* genNegDau = 0;
  for (unsigned int ndx = 0; ndx < genPH->size(); ndx++) {
    const Candidate & theVeeCand = (*genPH)[ndx];
    VeeTypes veeType = none;
    bool isAntiparticle = false;
    if (kshort && theVeeCand.pdgId() == 310) veeType = ks;
    if (!kshort && abs(theVeeCand.pdgId()) == 3122) {
      veeType = lam;
      if (theVeeCand.pdgId() < 0) isAntiparticle = true;
    }
    if (veeType == none) continue;

    bool skipMatch = false;
    bool isGoodDecay = true;
    if (theVeeCand.numberOfDaughters() != 2) skipMatch = true;

    GlobalPoint genVtx(0., 0., 0.);
    if (!skipMatch) {
      unsigned int posGenDauNdx = 0;
      unsigned int negGenDauNdx = 1;
      if (theVeeCand.daughter(0)->charge() < 0) {
	posGenDauNdx = 1;
	negGenDauNdx = 0;
      }
      const Candidate* genPosDau_ = theVeeCand.daughter(posGenDauNdx);
      const Candidate* genNegDau_ = theVeeCand.daughter(negGenDauNdx);

      genPosDau = genPosDau_;
      genNegDau = genNegDau_;

      if (veeType == ks
	  && ( abs(genPosDau->pdgId()) != 211 || abs(genNegDau->pdgId()) != 211
	       || genPosDau->vertex() != genNegDau->vertex() )) {
	isGoodDecay = false;
      }
      if (veeType == lam) {
	if (!isAntiparticle) {
	  if( abs(genPosDau->pdgId()) != 2212 || abs(genNegDau->pdgId()) != 211
	      || genPosDau->vertex() != genNegDau->vertex() ) {
	    isGoodDecay = false;
	  }
	}
	else {
	  if( abs(genPosDau->pdgId()) != 211 || abs(genNegDau->pdgId()) != 2212
	      || genPosDau->vertex() != genNegDau->vertex() ) {
	    isGoodDecay = false;
	  }
	}
      }
      genEta = theVeeCand.momentum().eta();
      double genE = -1.;
      if (veeType == ks) {
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + ksMassConst2);
      }
      else {
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + lamMassConst2);
      }
      double geny = 0.5*log((genE + theVeeCand.momentum().z())
		     / (genE - theVeeCand.momentum().z()));

      genRapidity = geny;
      genpT = theVeeCand.pt();
      genp = theVeeCand.p();
      genPhi = theVeeCand.momentum().phi();
      genMass = theVeeCand.mass();
      genPDG = theVeeCand.pdgId();
      genPosDauPDG = genPosDau->pdgId();
      genNegDauPDG = genNegDau->pdgId();
      genPosDauPx = genPosDau->px();
      genPosDauPy = genPosDau->py();
      genPosDauPz = genPosDau->pz();
      genNegDauPx = genNegDau->px();
      genNegDauPy = genNegDau->py();
      genNegDauPz = genNegDau->pz();

      double maxEta = genPosDau->momentum().eta();
      if(abs(genNegDau->momentum().eta()) > abs(maxEta)) {
	maxEta = genNegDau->momentum().eta();
      }
      gendauMaxEta = maxEta;
      genVtx = GlobalPoint( genPosDau->vertex().x(),
			    genPosDau->vertex().y(),
			    genPosDau->vertex().z() );
      GlobalPoint genProdVtx( theVeeCand.vertex().x(),
			      theVeeCand.vertex().y(),
			      theVeeCand.vertex().z() );

      genProdVtxX = theVeeCand.vertex().x();
      genProdVtxY = theVeeCand.vertex().y();
      genProdVtxZ = theVeeCand.vertex().z();

      genDecVtxX = genPosDau->vertex().x();
      genDecVtxY = genPosDau->vertex().y();
      genDecVtxZ = genPosDau->vertex().z();

      genRho = (genVtx - genProdVtx).perp();
      genR = (genVtx - genProdVtx).mag();
      genCtau = genR * genMass / genp;
      writeTree = true;
    }
    else if (theVeeCand.numberOfDaughters() ) {
      const Candidate* genDau = theVeeCand.daughter(0);
      genEta = theVeeCand.momentum().eta();

      double genE = 0.;
      if (veeType == ks) {
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + ksMassConst2);
      } else {
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + lamMassConst2);
      }
      double geny = 0.5*log((genE + theVeeCand.momentum().z())
			    / (genE - theVeeCand.momentum().z()));
      genRapidity = geny;
      genpT = theVeeCand.pt();
      genp = theVeeCand.p();
      genPhi = theVeeCand.momentum().phi();
      genMass = theVeeCand.mass();
      genPDG = theVeeCand.pdgId();
      genPosDauPDG = -1000;
      genNegDauPDG = -1000;
      gendauMaxEta = -1000.;

      genVtx = GlobalPoint(genDau->vertex().x(),
			   genDau->vertex().y(),
			   genDau->vertex().z());
      GlobalPoint genProdVtx(theVeeCand.vertex().x(),
			     theVeeCand.vertex().y(),
			     theVeeCand.vertex().z());
      genProdVtxX = genProdVtx.x();
      genProdVtxY = genProdVtx.y();
      genProdVtxZ = genProdVtx.z();

      genDecVtxX = genVtx.x();
      genDecVtxY = genVtx.y();
      genDecVtxZ = genVtx.z();

      genRho = (genVtx - genProdVtx).perp();
      genR = (genVtx - genProdVtx).mag();
      genCtau = genR * genMass / genp;
      writeTree = true;
    } else {
      genEta = theVeeCand.momentum().eta();

      double genE = 0.;
      if (veeType == ks) 
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + ksMassConst2);
      else
	genE = sqrt(theVeeCand.p()*theVeeCand.p() + lamMassConst2);
      double geny = 0.5*log((genE + theVeeCand.momentum().z())
		     / (genE - theVeeCand.momentum().z()));

      genRapidity = geny;
      genpT = theVeeCand.pt();
      genp = theVeeCand.p();
      genPhi = theVeeCand.momentum().phi();
      genMass = theVeeCand.mass();
      genPDG = theVeeCand.pdgId();
      genPosDauPDG = -1000;
      genNegDauPDG = -1000;
      gendauMaxEta = -1000.;

      genProdVtxX = theVeeCand.vertex().x();
      genProdVtxY = theVeeCand.vertex().y();
      genProdVtxZ = theVeeCand.vertex().z();

      genDecVtxX = -1000.;
      genDecVtxY = -1000.;
      genDecVtxZ = -1000.;

      genRho = -1000.;
      genR = -1000.;
      genCtau = -1000.;
      writeTree = true;
    }

    // Match reco now
      int matchStatus = 0;
      int lastMatchStatus = 0;
      double posDeltaRVal = -1.;
      double negDeltaRVal = -1.;

      double posDeltaEtaVal = -1.;
      double posDeltaPhiVal = -1.;
      double negDeltaEtaVal = -1.;
      double negDeltaPhiVal = -1.;

      double lastposDeltaEtaVal = -1.;
      double lastposDeltaPhiVal = -1.;
      double lastnegDeltaEtaVal = -1.;
      double lastnegDeltaPhiVal = -1.;

      double lastposDeltaRVal = -1.;
      double lastnegDeltaRVal = -1.;
      double deltaLVal = -1.;
      double lastdeltaLVal = -1.;

      for (pat::CompositeCandidateCollection::const_iterator iKS = patH->begin();
	   iKS != patH->end();
	   iKS++) {
	if (skipMatch || !isGoodDecay) break;
	if (veeType == ks && iKS->pdgId() != 310) continue;
	if (veeType == lam && abs(iKS->pdgId()) != 3122) continue;

	double recoE = 0.;
	if (veeType == ks)
	  recoE = sqrt(iKS->p()*iKS->p() + ksMassConst2);
	else
	  recoE = sqrt(iKS->p()*iKS->p() + lamMassConst2);
	double recoy = 0.5*log((recoE + iKS->momentum().z())
			/ (recoE - iKS->momentum().z()));
	
	bool newMatchThisIteration = false;

	GlobalPoint recoVtx(iKS->vertex().x(),
			    iKS->vertex().y(),
			    iKS->vertex().z());
	GlobalPoint priPos(bestVtx->position().x(),
			   bestVtx->position().y(),
			   bestVtx->position().z());
	GlobalVector lineOfFlight(recoVtx.x() - priPos.x(),
				  recoVtx.y() - priPos.y(),
				  recoVtx.z() - priPos.z());

	double recoDecayRadius = (recoVtx - priPos).perp();
	double recoDecayDistance = (recoVtx - priPos).mag();

	deltaLVal = (genVtx - recoVtx).mag();

	unsigned int posRecoDauNdx = 0;
	unsigned int negRecoDauNdx = 1;
	if (iKS->daughter(0)->charge() < 0) {
	  posRecoDauNdx = 1;
	  negRecoDauNdx = 0;
	}

	const Candidate* recoPosDau = iKS->daughter(posRecoDauNdx);
	const Candidate* recoNegDau = iKS->daughter(negRecoDauNdx);
	reco::TrackRef posDauTrackRef 
	  = ((dynamic_cast<const reco::RecoCandidate*>(iKS->daughter(posRecoDauNdx)))->track());
	reco::TrackRef negDauTrackRef
	  = ((dynamic_cast<const reco::RecoCandidate*>(iKS->daughter(negRecoDauNdx)))->track());
	reco::Track posDauTrack = *posDauTrackRef;
	reco::Track negDauTrack = *negDauTrackRef;

	// Kinematic fit for IP calculation
	KinematicParticleFactoryFromTransientTrack pFactory;
	vector<RefCountedKinematicParticle> v0Particles;
	vector<double> v0ParticleMasses;
	TransientTrack posDauTT = theTTB->build(posDauTrack);
	TransientTrack negDauTT = theTTB->build(negDauTrack);
	float chi = 0.;
	float ndf = 0.;

	std::vector<reco::TransientTrack> newVtxTracks;
	bool tracksRemoved = false;
	for (std::vector<TrackBaseRef>::const_iterator iTK = theOrigPrim->tracks_begin();
	     iTK != theOrigPrim->tracks_end();
	     iTK++) {
	  TrackRef currentRef = iTK->castTo<TrackRef>();
	  if (!((currentRef == posDauTrackRef) || (currentRef == negDauTrackRef))) {
	    newVtxTracks.push_back(theTTB->build(*currentRef));
	  } else {
	    tracksRemoved = true;
	  }
	}

	reco::Vertex refitVtx;
	if (newVtxTracks.size() > 0) {
	  AdaptiveVertexFitter theAdFitter;
	  TransientVertex theNewPrim = theAdFitter.vertex(newVtxTracks);
	  if (theNewPrim.isValid()) {
	    refitVtx = theNewPrim;
	    if (!tracksRemoved) priVtxWasRefit = 0;
	    else priVtxWasRefit = 1;
	    priRefitVtxChi2 = refitVtx.chi2();
	    priRefitVtxNdof = refitVtx.ndof();
	    priRefitVtxX = refitVtx.x();
	    priRefitVtxY = refitVtx.y();
	    priRefitVtxZ = refitVtx.z();
	  } else {
	    refitVtx = *theOrigPrim;
	    if (!tracksRemoved) priVtxWasRefit = 0;
	    else priVtxWasRefit = -2;
	    priRefitVtxChi2 = refitVtx.chi2();
	    priRefitVtxNdof = refitVtx.ndof();
	    priRefitVtxX = refitVtx.x();
	    priRefitVtxY = refitVtx.y();
	    priRefitVtxZ = refitVtx.z();
	  }
	} else {
	  refitVtx = *theOrigPrim;
	  if (!tracksRemoved) priVtxWasRefit = 0;
	  else priVtxWasRefit = -1;
	  priRefitVtxChi2 = refitVtx.chi2();
	  priRefitVtxNdof = refitVtx.ndof();
	  priRefitVtxX = refitVtx.x();
	  priRefitVtxY = refitVtx.y();
	  priRefitVtxZ = refitVtx.z();
	}
	newVtxTracks.clear();
	thePrim = &refitVtx;



	if (veeType == ks) {
	  v0Particles.push_back(pFactory.particle(posDauTT, piMassConst, chi, ndf, pionSigma));
	  v0ParticleMasses.push_back( piMassConst );
	  v0Particles.push_back(pFactory.particle(negDauTT, piMassConst, chi, ndf, pionSigma));
	  v0ParticleMasses.push_back( piMassConst );
	} else {
	  if (iKS->pdgId() > 0) {//is a lambda
	    v0Particles.push_back(pFactory.particle(posDauTT, protonMassConst, chi, ndf, protonSigma));
	    v0Particles.push_back(pFactory.particle(negDauTT, piMassConst, chi, ndf, pionSigma));
	  } else {
	    v0Particles.push_back(pFactory.particle(posDauTT, piMassConst, chi, ndf, pionSigma));
	    v0Particles.push_back(pFactory.particle(negDauTT, protonMassConst, chi, ndf, protonSigma));
	  }
	}

	// Continue from line 815
	KinematicParticleVertexFitter theFitter;
	RefCountedKinematicTree v0VtxFitTree;
	bool doKinFit = false;
	if( v0Particles.size() == 2) {
	  doKinFit = true;
	  v0VtxFitTree = theFitter.fit(v0Particles);
	}

	pair<bool,Measurement1D> thePosDau3DIpPair = 
	  IPTools::absoluteImpactParameter3D(posDauTT, (*bestVtx));
	pair<bool,Measurement1D> theNegDau3DIpPair = 
	  IPTools::absoluteImpactParameter3D(negDauTT, (*bestVtx));

	double posDau3DIpVal = -1000., posDau3DIpErrorVal = -1000., posDau3DIpSigVal = -1000.;
	double negDau3DIpVal = -1000., negDau3DIpErrorVal = -1000., negDau3DIpSigVal = -1000.;

	if (thePosDau3DIpPair.first) {
	  posDau3DIpVal = thePosDau3DIpPair.second.value();
	  posDau3DIpErrorVal = thePosDau3DIpPair.second.error();
	  posDau3DIpSigVal = thePosDau3DIpPair.second.significance();
	}
	if (theNegDau3DIpPair.first) {
	  negDau3DIpVal = theNegDau3DIpPair.second.value();
	  negDau3DIpErrorVal = theNegDau3DIpPair.second.error();
	  negDau3DIpSigVal = theNegDau3DIpPair.second.significance();
	} 

	double v3dIpWrtPrimaryVal = -1000.;
	double v3dIpWrtPrimaryErrorVal = -1000.;
	double v3dIpWrtPrimarySigVal = -1000.;
	if( doKinFit && v0VtxFitTree->isValid() ) {

	  v0VtxFitTree->movePointerToTheTop();

	  RefCountedKinematicParticle v0FitKinPart = v0VtxFitTree->currentParticle();
	  RefCountedKinematicVertex v0FitKinVtx = v0VtxFitTree->currentDecayVertex();

	  KinematicState theCurrentKinematicState = v0FitKinPart->currentState();

	  FreeTrajectoryState theV0FTS = theCurrentKinematicState.freeTrajectoryState();

	  TransientTrack v0TT = theTTB->build(theV0FTS);

	  pair<bool,Measurement1D> theV03dIpPair = IPTools::signedImpactParameter3D(v0TT,
										    lineOfFlight,
										    *bestVtx);
	  if (theV03dIpPair.first) {
	    v3dIpWrtPrimaryVal = theV03dIpPair.second.value();
	    v3dIpWrtPrimaryErrorVal = theV03dIpPair.second.error();
	    v3dIpWrtPrimarySigVal = theV03dIpPair.second.significance();
	  }
	}
	GlobalVector posP(iKS->userFloat("posTkPX"),
			  iKS->userFloat("posTkPY"),
			  iKS->userFloat("posTkPZ") );
	GlobalVector negP(iKS->userFloat("negTkPX"),
			  iKS->userFloat("negTkPY"),
			  iKS->userFloat("negTkPZ") );
	double protonE = sqrt( posP.mag2() + protonMassConst*protonMassConst );
	double antiProtonE = sqrt( negP.mag2() + protonMassConst*protonMassConst );
	double posPiE = sqrt( posP.mag2() + piMassConst*piMassConst );
	double negPiE = sqrt( negP.mag2() + piMassConst*piMassConst );
	GlobalVector totalP = posP + negP;
	double ksE = posPiE + negPiE;
	double lamE = protonE + negPiE;
	double lamBarE = antiProtonE + posPiE;
	double ksMass = sqrt( ksE*ksE - totalP.mag2() );
	double lamMass = sqrt( lamE*lamE - totalP.mag2() );
	double lamBarMass = sqrt( lamBarE*lamBarE - totalP.mag2() );
	double otherCandMassVal = -1000.;
	if (veeType == ks) {
	  if (protonE > antiProtonE) {
	    otherCandMassVal = lamMass;
	  } else {
	    otherCandMassVal = lamBarMass;
	  }
	} else {
	  otherCandMassVal = ksMass;
	}

	std::vector<double> covv;
	covv.push_back( iKS->userFloat("vtxCov00") );
	covv.push_back( iKS->userFloat("vtxCov01") );
	covv.push_back( iKS->userFloat("vtxCov11") );
	covv.push_back( iKS->userFloat("vtxCov02") );
	covv.push_back( iKS->userFloat("vtxCov12") );
	covv.push_back( iKS->userFloat("vtxCov22") );

	SMatrixSym3D vtxCov( covv.begin(), covv.end() );
	covv.clear();

	GlobalPoint primPos(thePrim->position().x(),
			    thePrim->position().y(),
			    thePrim->position().z());
	//if (!priVtxIsFake) beamSpotPos = primPos;

	SMatrixSym3D totalCov2D = theBeamSpotH->covariance3D() + vtxCov;
	SMatrixSym3D totalCov = thePrim->covariance() + vtxCov;
	SVector3 distanceVector( iKS->vertex().x() - beamSpotPos.x(),
				 iKS->vertex().y() - beamSpotPos.y(),
				 0. );
	//SVector3 distanceVector( iKS->vertex().x() - bestVtx->position().x(),
	//			 iKS->vertex().y() - bestVtx->position().y(),
	//			 0. );
	double rVtxMag = ROOT::Math::Mag(distanceVector);
	double sigmaRvtxMag = 
	  sqrt(ROOT::Math::Similarity( totalCov2D, distanceVector )) / rVtxMag;	
	SVector3 distanceVector3D( iKS->vertex().x() - thePrim->position().x(),
				   iKS->vertex().y() - thePrim->position().y(),
				   iKS->vertex().z() - thePrim->position().z() );
	double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
	double sigmaRvtxMag3D =
	  sqrt(ROOT::Math::Similarity( totalCov, distanceVector3D )) / rVtxMag3D;

	reco::Track posDauTk( *((dynamic_cast<const reco::RecoCandidate*>(recoPosDau))->track()));
	reco::Track negDauTk( *((dynamic_cast<const reco::RecoCandidate*>(recoNegDau))->track()));

	double posDauChi2 = -1000., posDauNormChi2 = -1000., posDauD0 = -1000.;
	double negDauChi2 = -1000., negDauNormChi2 = -1000.;//, negDauD0 = -1000.;
	int posDauNhitsVal = -1000, negDauNhitsVal = -1000;
	posDauChi2 = posDauTk.chi2();
	posDauNormChi2 = posDauTk.normalizedChi2();
	posDauD0 = posDauTk.d0();
	negDauChi2 = negDauTk.chi2();
        negDauNormChi2 = negDauTk.normalizedChi2();
	//negDauD0 = negDauTk.d0();
	posDauNhitsVal = posDauTk.numberOfValidHits();
	negDauNhitsVal = negDauTk.numberOfValidHits();

	// Cut on DeltaR
	double genPosDauPhi = genPosDau->momentum().phi();
	double recoPosDauPhi = recoPosDau->momentum().phi();

	double genPosDauEta = genPosDau->momentum().eta();
	double recoPosDauEta = recoPosDau->momentum().eta();

	double genNegDauPhi = genNegDau->momentum().phi();
	double recoNegDauPhi = recoNegDau->momentum().phi();

	double genNegDauEta = genNegDau->momentum().eta();
	double recoNegDauEta = recoNegDau->momentum().eta();

	double posDeltaPhi = reco::deltaPhi( genPosDauPhi, recoPosDauPhi);
	posDeltaPhiVal = posDeltaPhi;
	double negDeltaPhi = reco::deltaPhi( genNegDauPhi, recoNegDauPhi);
	negDeltaPhiVal = negDeltaPhi;
	double posDeltaEta = genPosDauEta - recoPosDauEta;
	posDeltaEtaVal = posDeltaEta;
	double negDeltaEta = genNegDauEta - recoNegDauEta;
	negDeltaEtaVal = negDeltaEta;

	while( posDeltaPhi > M_PI ) posDeltaPhi -= 2*M_PI;
	while( posDeltaPhi < -M_PI ) posDeltaPhi += 2*M_PI;
	while( negDeltaPhi > M_PI ) negDeltaPhi -= 2*M_PI;
	while( negDeltaPhi < -M_PI ) negDeltaPhi += 2*M_PI;

	posDeltaRVal = sqrt( posDeltaPhi*posDeltaPhi + posDeltaEta*posDeltaEta );
	negDeltaRVal = sqrt( negDeltaPhi*negDeltaPhi + negDeltaEta*negDeltaEta );


	if( posDeltaRVal < pideltaRcut && negDeltaRVal < pideltaRcut && deltaLVal < deltaLcut) {
	  // Good deltaR match for both tracks, good deltaL match for vertices
	  matchStatus = 1;
	  newMatchThisIteration = true;
	  if(debug) cout << "Found good match.  genEta: " << theVeeCand.momentum().eta()
			 << ", recoEta: " << iKS->momentum().eta() << endl;
	}
	else if ( posDeltaRVal < pideltaRcut && negDeltaRVal < pideltaRcut ) {
	  matchStatus = 2;
	  newMatchThisIteration = true;
	}
	else if ( (posDeltaRVal < pideltaRcut || negDeltaRVal < pideltaRcut) 
		  && deltaLVal < deltaLcut) {
	  // Good deltaR match for only 1 track, good vertex match
	  matchStatus = 3;
	  newMatchThisIteration = true;
	}
	else if ( (posDeltaRVal < pideltaRcut || negDeltaRVal < pideltaRcut) ) {
	  // Good match for 1 track, no vertex match
	  matchStatus = 4;
	  newMatchThisIteration = true;
	}


	 
	if( lastMatchStatus > 0 && matchStatus >= lastMatchStatus ) {
	  matchStatus = lastMatchStatus;
	  posDeltaRVal = lastposDeltaRVal;
	  negDeltaRVal = lastnegDeltaRVal;
	  deltaLVal = lastdeltaLVal;
	  posDeltaEtaVal = lastposDeltaEtaVal;
	  posDeltaPhiVal = lastposDeltaPhiVal;
	  negDeltaEtaVal = lastnegDeltaEtaVal;
	  negDeltaPhiVal = lastnegDeltaPhiVal;
	} else if (matchStatus && newMatchThisIteration) {
	  recoEta = iKS->momentum().eta();
	  recoRapidity = recoy;
	  recopT = iKS->pt();
	  recop = iKS->p();
	  recoPhi = iKS->momentum().phi();
	  recoRho = recoDecayRadius;
	  recoR = recoDecayDistance;
	  recoCtau = recoDecayDistance * iKS->mass() / iKS->p();
	  recoVtxX = iKS->vertex().x();
	  recoVtxY = iKS->vertex().y();
	  recoVtxZ = iKS->vertex().z();
	  recoVtxSig = rVtxMag / sigmaRvtxMag;
	  recoVtxSig3D = rVtxMag3D / sigmaRvtxMag3D;
	  recoVtxChi2 = iKS->userFloat("vtxChi2");
	  recoVtxNormChi2 = iKS->userFloat("vtxNormChi2");
	  recopX = iKS->momentum().x();
	  recopY = iKS->momentum().y();
	  recopZ = iKS->momentum().z();
	  recoOtherCandMass = otherCandMassVal;
	  recoOtherCandMassError = iKS->userFloat("otherMassErrorFromFit");
	  recoMass = iKS->mass();
	  recoMassError = iKS->userFloat("massErrorFromFit");
	  recoPosDauChi2 = posDauChi2;
	  recoPosDauNormChi2 = posDauNormChi2;
	  recoPosDauD0 = posDauD0;
	  recoNegDauChi2 = negDauChi2;
	  recoNegDauNormChi2 = negDauNormChi2;
	  v3dIpWrtPrimary = v3dIpWrtPrimaryVal;
	  v3dIpWrtPrimaryError = v3dIpWrtPrimaryErrorVal;
	  v3dIpWrtPrimarySig = v3dIpWrtPrimarySigVal;
	  posDauNhits = posDauNhitsVal;
	  negDauNhits = negDauNhitsVal;
	  posDau3DIp = posDau3DIpVal;
	  posDau3DIpError = posDau3DIpErrorVal;
	  posDau3DIpSig = posDau3DIpSigVal;
	  negDau3DIp = negDau3DIpVal;
	  negDau3DIpError = negDau3DIpErrorVal;
	  negDau3DIpSig = negDau3DIpSigVal;
	  recoPosDauPx = recoPosDau->px();
	  recoPosDauPy = recoPosDau->py();
	  recoPosDauPz = recoPosDau->pz();
	  recoNegDauPx = recoNegDau->px();
	  recoNegDauPy = recoNegDau->py();
	  recoNegDauPz = recoNegDau->pz();

	}
	lastMatchStatus = matchStatus;
	lastposDeltaRVal = posDeltaRVal;
	lastnegDeltaRVal = negDeltaRVal;
	lastposDeltaEtaVal = posDeltaEtaVal;
	lastposDeltaPhiVal = posDeltaPhiVal;
	lastnegDeltaEtaVal = negDeltaEtaVal;
	lastnegDeltaPhiVal = negDeltaPhiVal;
	lastdeltaLVal = deltaLVal;
      }
      v0MatchStatus = matchStatus;
      posDeltaR = posDeltaRVal;
      posDeltaEta = posDeltaEtaVal;
      posDeltaPhi = posDeltaPhiVal;
      if (!skipMatch && isGoodDecay) {
	pospT = genPosDau->pt();
	posp = genPosDau->p();
      } else {
	pospT = -1000.;
	posp = -1000.;
      }
      negDeltaR = negDeltaRVal;
      negDeltaEta = negDeltaEtaVal;
      negDeltaPhi = negDeltaPhiVal;
      if (!skipMatch && isGoodDecay) {
	negpT = genNegDau->pt();
	negp = genNegDau->p();
      } else {
	negpT = -1000.;
	negp = -1000.;
      }
      deltaL = deltaLVal;
      if (writeTree) {
	eff->Fill();
      }

      genEta = genRapidity = genpT = genp = genPhi = genRho = genR = genMass = genCtau = 0.;
      gendauMaxEta = 0.;
      recoEta = recoRapidity = recopT = recop = recoPhi = recoRho = recoR = recoMass = recoMassError = recoCtau = 0.;
      recoVtxX = recoVtxY = recoVtxZ = recoVtxSig = recoVtxSig3D = recoVtxChi2 = recoVtxNormChi2 
	= recopX = recopY = recopZ = recoOtherCandMass = recoOtherCandMassError = 0.;
      recoPosDauChi2 = recoPosDauNormChi2 = recoPosDauD0
	= recoNegDauChi2 = recoNegDauNormChi2 = recoNegDauD0 = 0.;
      recoPosDauPx = recoPosDauPy = recoPosDauPz = 0;
      recoNegDauPx = recoNegDauPy = recoNegDauPz = 0;
      v3dIpWrtPrimary = v3dIpWrtPrimaryError = v3dIpWrtPrimarySig = 0.;
      posDau3DIp = posDau3DIpError = posDau3DIpSig
	= negDau3DIp = negDau3DIpError = negDau3DIpSig = 0.;
      posDauNhits = negDauNhits = 0;
      posDeltaR = posDeltaEta = posDeltaPhi = pospT = posp
	= negDeltaR = negDeltaEta = negDeltaPhi = negpT = negp = 0.;
      genProdVtxX = genProdVtxY = genProdVtxZ
	= genDecVtxX = genDecVtxY = genDecVtxZ = 0.;
      v0MatchStatus = 0;
      genPDG = genPosDauPDG = genNegDauPDG = 0;
      genPosDauPx = genPosDauPy = genPosDauPz = 0;
      genNegDauPx = genNegDauPy = genNegDauPz = 0;

  }
  // Clear event-wide variables
  trigHF = trigTech40_41 = trigScraping = trigTech36_39 = trigTech34 = trigAlgo124 = trigHLTminBias = 0;
  HFpE = HFmE = 0.;
  //nLooseTracks = nHighPurTracks = nLooseTracks_pt_100 = nHighPurTracks_pt_100 = 0;
  nTracks = nTracks_pt_100 = nPromptTracks = nPromptTracks_pt_100 = 0;
  nGenParts = nGenParts_pt_100 = 0;
  runNum = evtNum = processType = 0;
  priVtxX = priVtxY = priVtxZ = priVtxChi2 = priVtxNdof = priVtxNtracks = priVtxIsFake = 0;
  nPriVtxs = 0;


}

void V0DeltaREffAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  if (kshort) {
    eff = fs->make<TTree>("kseff", "ks efficiency");
    eff->Branch("kstrigHF", &trigHF, "trigHF/I");
    eff->Branch("HFpE", &HFpE, "HFpE/F");
    eff->Branch("HFmE", &HFmE, "HFmE/F");
    eff->Branch("kstrigTech40_41", &trigTech40_41, "trigTech40_41/I");
    eff->Branch("kstrigScraping", &trigScraping, "trigScraping/I");
    eff->Branch("kstrigTech36_39", &trigTech36_39, "trigTech36_39/I");
    eff->Branch("kstrigTech34", &trigTech34, "trigTech34/I");
    eff->Branch("trigHLTminBias", &trigHLTminBias, "trigHLTminBias/I");
    eff->Branch("kstrigAlgo124", &trigAlgo124, "trigAlgo124/I");

 //    eff->Branch("nLooseTracks", &nLooseTracks, "nLooseTracks/I");
//     eff->Branch("nHighPurTracks", &nHighPurTracks, "nHighPurTracks/I");
//     eff->Branch("nLooseTracks_pt_100", &nLooseTracks_pt_100, "nLooseTracks_pt_100/I");
//     eff->Branch("nHighPurTracks_pt_100", &nHighPurTracks_pt_100, "nHighPurTracks_pt_100/I");

    eff->Branch("nTracks", &nTracks, "nTracks/I");
    eff->Branch("nTracks_pt_100", &nTracks_pt_100, "nTracks_pt_100/I");
    eff->Branch("nPromptTracks", &nPromptTracks, "nPromptTracks/I");
    eff->Branch("nPromptTracks_pt_100", &nPromptTracks_pt_100, "nPromptTracks_pt_100/I");

    eff->Branch("nGenParts", &nGenParts, "nGenParts/I");
    eff->Branch("nGenParts_pt_100", &nGenParts_pt_100, "nGenParts_pt_100/I");

    eff->Branch("runNum", &runNum, "runNum/I");
    eff->Branch("evtNum", &evtNum, "evtNum/I");
    eff->Branch("processType", &processType, "processType/I");

    eff->Branch("priVtxX", &priVtxX, "priVtxX/F");
    eff->Branch("priVtxY", &priVtxY, "priVtxY/F");
    eff->Branch("priVtxZ", &priVtxZ, "priVtxZ/F");
    eff->Branch("priVtxChi2", &priVtxChi2, "priVtxChi2/F");
    eff->Branch("priVtxNdof", &priVtxNdof, "priVtxNdof/F");
    eff->Branch("priVtxNtracks", &priVtxNtracks, "priVtxNtracks/I");
    //eff->Branch("priVtxIsFake", &priVtxIsFake, "priVtxIsFake/I");
    eff->Branch("nPriVtxs", &nPriVtxs, "nPriVtxs/I");

    eff->Branch("priVtxWasRefit", &priVtxWasRefit, "priVtxWasRefit/I");
    eff->Branch("priRefitVtxChi2", &priRefitVtxChi2, "priRefitVtxChi2/F");
    eff->Branch("priRefitVtxNdof", &priRefitVtxNdof, "priRefitVtxNdof/F");
    eff->Branch("priRefitVtxX", &priRefitVtxX, "priRefitVtxX/F");
    eff->Branch("priRefitVtxY", &priRefitVtxY, "priRefitVtxY/F");
    eff->Branch("priRefitVtxZ", &priRefitVtxZ, "priRefitVtxZ/F");

    eff->Branch("ksgenEta", &genEta, "genEta/F");
    eff->Branch("ksgenRapidity", &genRapidity, "genRapidity/F");
    eff->Branch("ksgenpT", &genpT, "genpT/F");
    eff->Branch("ksgenp", &genp, "genp/F");
    eff->Branch("ksgenPhi", &genPhi, "genPhi/F");
    eff->Branch("ksgenRho", &genRho, "genRho/F");
    eff->Branch("ksgenR", &genR, "genR/F");
    eff->Branch("ksgenMass", &genMass, "genMass/F");
    eff->Branch("ksgenCtau", &genCtau, "genCtau/F");
    eff->Branch("ksgendauMaxEta", &gendauMaxEta, "gendauMaxEta/F");

    eff->Branch("ksrecoEta", &recoEta, "recoEta/F");
    eff->Branch("ksrecoRapidity", &recoRapidity, "recoRapidity/F");
    eff->Branch("ksrecopT", &recopT, "recopT/F");
    eff->Branch("ksrecop", &recop, "recop/F");
    eff->Branch("ksrecoPhi", &recoPhi, "recoPhi/F");
    eff->Branch("ksrecoRho", &recoRho, "recoRho/F");
    eff->Branch("ksrecoR", &recoR, "recoR/F");
    eff->Branch("ksrecoMass", &recoMass, "recoMass/F");
    eff->Branch("ksrecoMassError", &recoMassError, "recoMassError/F");
    eff->Branch("ksrecoCtau", &recoCtau, "recoCtau/F");

    eff->Branch("ksrecoVtxX", &recoVtxX, "recoVtxX/F");
    eff->Branch("ksrecoVtxY", &recoVtxY, "recoVtxY/F");
    eff->Branch("ksrecoVtxZ", &recoVtxZ, "recoVtxZ/F");
    eff->Branch("ksrecoVtxSig", &recoVtxSig, "recoVtxSig/F");
    eff->Branch("ksrecoVtxSig3D", &recoVtxSig3D, "recoVtxSig3D/F");
    eff->Branch("ksrecoVtxChi2", &recoVtxChi2, "recoVtxChi2/F");
    //eff->Branch("ksrecoVtxNormChi2", &recoVtxNormChi2, "recoVtxNormChi2/F");
    eff->Branch("ksrecopX", &recopX, "recopX/F");
    eff->Branch("ksrecopY", &recopY, "recopY/F");
    eff->Branch("ksrecopZ", &recopZ, "recopZ/F");
    eff->Branch("ksrecoOtherCandMass", &recoOtherCandMass, "recoOtherCandMass/F");
    eff->Branch("ksrecoOtherCandMassError", &recoOtherCandMassError, "recoOtherCandMassError/F");

    eff->Branch("ksrecoPosDauChi2", &recoPosDauChi2, "recoPosDauChi2/F");
    eff->Branch("ksrecoPosDauNormChi2", &recoPosDauNormChi2, "recoPosDauNormChi2/F");
    eff->Branch("ksrecoPosDauD0", &recoPosDauD0, "recoPosDauD0/F");
    eff->Branch("ksrecoNegDauChi2", &recoNegDauChi2, "recoNegDauChi2/F");
    eff->Branch("ksrecoNegDauNormChi2", &recoNegDauNormChi2, "recoNegDauNormChi2/F");
    eff->Branch("ksrecoNegDauD0", &recoNegDauD0, "recoNegDauD0/F");
    eff->Branch("ksrecoPosDauPx", &recoPosDauPx, "recoPosDauPx/F");
    eff->Branch("ksrecoPosDauPy", &recoPosDauPy, "recoPosDauPy/F");
    eff->Branch("ksrecoPosDauPz", &recoPosDauPz, "recoPosDauPz/F");
    eff->Branch("ksrecoNegDauPx", &recoNegDauPx, "recoNegDauPx/F");
    eff->Branch("ksrecoNegDauPy", &recoNegDauPy, "recoNegDauPy/F");
    eff->Branch("ksrecoNegDauPz", &recoNegDauPz, "recoNegDauPz/F");

    eff->Branch("ks3dIpWrtPrimary", &v3dIpWrtPrimary, "v3dIpWrtPrimary/F");
    //eff->Branch("ks3dIpWrtPrimaryError", &v3dIpWrtPrimaryError, "v3dIpWrtPrimaryError/F");
    eff->Branch("ks3dIpWrtPrimarySig", &v3dIpWrtPrimarySig, "v3dIpWrtPrimarySig/F");

    eff->Branch("ksposDau3DIp", &posDau3DIp, "posDau3DIp/F");
    //eff->Branch("ksposDau3DIpError", &posDau3DIpError, "posDau3DIpError/F");
    eff->Branch("ksposDau3DIpSig", &posDau3DIpSig, "posDau3DIpSig/F");
    eff->Branch("ksnegDau3DIp", &negDau3DIp, "negDau3DIp/F");
    //eff->Branch("ksnegDau3DIpError", &negDau3DIpError, "negDau3DIpError/F");
    eff->Branch("ksnegDau3DIpSig", &negDau3DIpSig, "negDau3DIpSig/F");

    eff->Branch("ksposDauNhits", &posDauNhits, "posDauNhits/I");
    eff->Branch("ksnegDauNhits", &negDauNhits, "negDauNhits/I");

    eff->Branch("ksposDeltaR", &posDeltaR, "posDeltaR/F");
    eff->Branch("ksposDeltaEta", &posDeltaEta, "posDeltaEta/F");
    eff->Branch("ksposDeltaPhi", &posDeltaPhi, "posDeltaPhi/F");
    eff->Branch("kspospT", &pospT, "pospT/F");
    eff->Branch("ksposp", &posp, "posp/F");
    eff->Branch("ksnegDeltaR", &negDeltaR, "negDeltaR/F");
    eff->Branch("ksnegDeltaEta", &negDeltaEta, "negDeltaEta/F");
    eff->Branch("ksnegDeltaPhi", &negDeltaPhi, "negDeltaPhi/F");
    eff->Branch("ksnegpT", &negpT, "negpT/F");
    eff->Branch("ksnegp", &negp, "negp/F");
    eff->Branch("ksdeltaL", &deltaL, "deltaL/F");

    eff->Branch("ksgenProdVtxX", &genProdVtxX, "genProdVtxX/F");
    eff->Branch("ksgenProdVtxY", &genProdVtxY, "genProdVtxY/F");
    eff->Branch("ksgenProdVtxZ", &genProdVtxZ, "genProdVtxZ/F");

    eff->Branch("ksgenDecVtxX", &genDecVtxX, "genDecVtxX/F");
    eff->Branch("ksgenDecVtxY", &genDecVtxY, "genDecVtxY/F");
    eff->Branch("ksgenDecVtxZ", &genDecVtxZ, "genDecVtxZ/F");

    eff->Branch("ksv0MatchStatus", &v0MatchStatus, "v0MatchStatus/I");

    eff->Branch("ksgenPDG", &genPDG, "genPDG/I");
    eff->Branch("ksgenPosDauPDG", &genPosDauPDG, "genPosDauPDG/I");
    eff->Branch("ksgenNegDauPDG", &genNegDauPDG, "genNegDauPDG/I");
    eff->Branch("ksgenPosDauPx", &genPosDauPx, "genPosDauPx/F");
    eff->Branch("ksgenPosDauPy", &genPosDauPy, "genPosDauPy/F");
    eff->Branch("ksgenPosDauPz", &genPosDauPz, "genPosDauPz/F");
    eff->Branch("ksgenNegDauPx", &genNegDauPx, "genNegDauPx/F");
    eff->Branch("ksgenNegDauPy", &genNegDauPy, "genNegDauPy/F");
    eff->Branch("ksgenNegDauPz", &genNegDauPz, "genNegDauPz/F");
  } else {
    eff = fs->make<TTree>("lameff", "lam efficiency");
    eff->Branch("lamtrigHF", &trigHF, "trigHF/I");
    eff->Branch("HFpE", &HFpE, "HFpE/F");
    eff->Branch("HFmE", &HFmE, "HFmE/F");
    eff->Branch("lamtrigTech40_41", &trigTech40_41, "trigTech40_41/I");
    eff->Branch("lamtrigScraping", &trigScraping, "trigScraping/I");
    eff->Branch("lamtrigTech36_39", &trigTech36_39, "trigTech36_39/I");
    eff->Branch("lamtrigTech34", &trigTech34, "trigTech34/I");
    eff->Branch("trigHLTminBias", &trigHLTminBias, "trigHLTminBias/I");
    eff->Branch("lamtrigAlgo124", &trigAlgo124, "trigAlgo124/I");

 //    eff->Branch("nLooseTracks", &nLooseTracks, "nLooseTracks/I");
//     eff->Branch("nHighPurTracks", &nHighPurTracks, "nHighPurTracks/I");
//     eff->Branch("nLooseTracks_pt_100", &nLooseTracks_pt_100, "nLooseTracks_pt_100/I");
//     eff->Branch("nHighPurTracks_pt_100", &nHighPurTracks_pt_100, "nHighPurTracks_pt_100/I");

    eff->Branch("nTracks", &nTracks, "nTracks/I");
    eff->Branch("nTracks_pt_100", &nTracks_pt_100, "nTracks_pt_100/I");
    eff->Branch("nPromptTracks", &nPromptTracks, "nPromptTracks/I");
    eff->Branch("nPromptTracks_pt_100", &nPromptTracks_pt_100, "nPromptTracks_pt_100/I");

    eff->Branch("nGenParts", &nGenParts, "nGenParts/I");
    eff->Branch("nGenParts_pt_100", &nGenParts_pt_100, "nGenParts_pt_100/I");

    eff->Branch("runNum", &runNum, "runNum/I");
    eff->Branch("evtNum", &evtNum, "evtNum/I");
    eff->Branch("processType", &processType, "processType/I");

    eff->Branch("priVtxX", &priVtxX, "priVtxX/F");
    eff->Branch("priVtxY", &priVtxY, "priVtxY/F");
    eff->Branch("priVtxZ", &priVtxZ, "priVtxZ/F");
    eff->Branch("priVtxChi2", &priVtxChi2, "priVtxChi2/F");
    eff->Branch("priVtxNdof", &priVtxNdof, "priVtxNdof/F");
    eff->Branch("priVtxNtracks", &priVtxNtracks, "priVtxNtracks/I");
    //eff->Branch("priVtxIsFake", &priVtxIsFake, "priVtxIsFake/I");
    eff->Branch("nPriVtxs", &nPriVtxs, "nPriVtxs/I");

    eff->Branch("priVtxWasRefit", &priVtxWasRefit, "priVtxWasRefit/I");
    eff->Branch("priRefitVtxChi2", &priRefitVtxChi2, "priRefitVtxChi2/F");
    eff->Branch("priRefitVtxNdof", &priRefitVtxNdof, "priRefitVtxNdof/F");
    eff->Branch("priRefitVtxX", &priRefitVtxX, "priRefitVtxX/F");
    eff->Branch("priRefitVtxY", &priRefitVtxY, "priRefitVtxY/F");
    eff->Branch("priRefitVtxZ", &priRefitVtxZ, "priRefitVtxZ/F");

    eff->Branch("lamgenEta", &genEta, "genEta/F");
    eff->Branch("lamgenRapidity", &genRapidity, "genRapidity/F");
    eff->Branch("lamgenpT", &genpT, "genpT/F");
    eff->Branch("lamgenp", &genp, "genp/F");
    eff->Branch("lamgenPhi", &genPhi, "genPhi/F");
    eff->Branch("lamgenRho", &genRho, "genRho/F");
    eff->Branch("lamgenR", &genR, "genR/F");
    eff->Branch("lamgenMass", &genMass, "genMass/F");
    eff->Branch("lamgenCtau", &genCtau, "genCtau/F");
    eff->Branch("lamgendauMaxEta", &gendauMaxEta, "gendauMaxEta/F");

    eff->Branch("lamrecoEta", &recoEta, "recoEta/F");
    eff->Branch("lamrecoRapidity", &recoRapidity, "recoRapidity/F");
    eff->Branch("lamrecopT", &recopT, "recopT/F");
    eff->Branch("lamrecop", &recop, "recop/F");
    eff->Branch("lamrecoPhi", &recoPhi, "recoPhi/F");
    eff->Branch("lamrecoRho", &recoRho, "recoRho/F");
    eff->Branch("lamrecoR", &recoR, "recoR/F");
    eff->Branch("lamrecoMass", &recoMass, "recoMass/F");
    eff->Branch("lamrecoMassError", &recoMassError, "recoMassError/F");
    eff->Branch("lamrecoCtau", &recoCtau, "recoCtau/F");

    eff->Branch("lamrecoVtxX", &recoVtxX, "recoVtxX/F");
    eff->Branch("lamrecoVtxY", &recoVtxY, "recoVtxY/F");
    eff->Branch("lamrecoVtxZ", &recoVtxZ, "recoVtxZ/F");
    eff->Branch("lamrecoVtxSig", &recoVtxSig, "recoVtxSig/F");
    eff->Branch("lamrecoVtxSig3D", &recoVtxSig3D, "recoVtxSig3D/F");
    eff->Branch("lamrecoVtxChi2", &recoVtxChi2, "recoVtxChi2/F");
    //eff->Branch("lamrecoVtxNormChi2", &recoVtxNormChi2, "recoVtxNormChi2/F");
    eff->Branch("lamrecopX", &recopX, "recopX/F");
    eff->Branch("lamrecopY", &recopY, "recopY/F");
    eff->Branch("lamrecopZ", &recopZ, "recopZ/F");
    eff->Branch("lamrecoOtherCandMass", &recoOtherCandMass, "recoOtherCandMass/F");
    eff->Branch("lamrecoOtherCandMassError", &recoOtherCandMassError, "recoOtherCandMassError/F");

    eff->Branch("lamrecoPosDauChi2", &recoPosDauChi2, "recoPosDauChi2/F");
    eff->Branch("lamrecoPosDauNormChi2", &recoPosDauNormChi2, "recoPosDauNormChi2/F");
    eff->Branch("lamrecoPosDauD0", &recoPosDauD0, "recoPosDauD0/F");
    eff->Branch("lamrecoNegDauChi2", &recoNegDauChi2, "recoNegDauChi2/F");
    eff->Branch("lamrecoNegDauNormChi2", &recoNegDauNormChi2, "recoNegDauNormChi2/F");
    eff->Branch("lamrecoNegDauD0", &recoNegDauD0, "recoNegDauD0/F");
    eff->Branch("lamrecoPosDauPx", &recoPosDauPx, "recoPosDauPx/F");
    eff->Branch("lamrecoPosDauPy", &recoPosDauPy, "recoPosDauPy/F");
    eff->Branch("lamrecoPosDauPz", &recoPosDauPz, "recoPosDauPz/F");
    eff->Branch("lamrecoNegDauPx", &recoNegDauPx, "recoNegDauPx/F");
    eff->Branch("lamrecoNegDauPy", &recoNegDauPy, "recoNegDauPy/F");
    eff->Branch("lamrecoNegDauPz", &recoNegDauPz, "recoNegDauPz/F");

    eff->Branch("lam3dIpWrtPrimary", &v3dIpWrtPrimary, "v3dIpWrtPrimary/F");
    //eff->Branch("lam3dIpWrtPrimaryError", &v3dIpWrtPrimaryError, "v3dIpWrtPrimaryError/F");
    eff->Branch("lam3dIpWrtPrimarySig", &v3dIpWrtPrimarySig, "v3dIpWrtPrimarySig/F");

    eff->Branch("lamposDau3DIp", &posDau3DIp, "posDau3DIp/F");
    //eff->Branch("lamposDau3DIpError", &posDau3DIpError, "posDau3DIpError/F");
    eff->Branch("lamposDau3DIpSig", &posDau3DIpSig, "posDau3DIpSig/F");
    eff->Branch("lamnegDau3DIp", &negDau3DIp, "negDau3DIp/F");
    //eff->Branch("lamnegDau3DIpError", &negDau3DIpError, "negDau3DIpError/F");
    eff->Branch("lamnegDau3DIpSig", &negDau3DIpSig, "negDau3DIpSig/F");

    eff->Branch("lamposDauNhits", &posDauNhits, "posDauNhits/I");
    eff->Branch("lamnegDauNhits", &negDauNhits, "negDauNhits/I");

    eff->Branch("lamposDeltaR", &posDeltaR, "posDeltaR/F");
    eff->Branch("lamposDeltaEta", &posDeltaEta, "posDeltaEta/F");
    eff->Branch("lamposDeltaPhi", &posDeltaPhi, "posDeltaPhi/F");
    eff->Branch("lampospT", &pospT, "pospT/F");
    eff->Branch("lamposp", &posp, "posp/F");
    eff->Branch("lamnegDeltaR", &negDeltaR, "negDeltaR/F");
    eff->Branch("lamnegDeltaEta", &negDeltaEta, "negDeltaEta/F");
    eff->Branch("lamnegDeltaPhi", &negDeltaPhi, "negDeltaPhi/F");
    eff->Branch("lamnegpT", &negpT, "negpT/F");
    eff->Branch("lamnegp", &negp, "negp/F");
    eff->Branch("lamdeltaL", &deltaL, "deltaL/F");

    eff->Branch("lamgenProdVtxX", &genProdVtxX, "genProdVtxX/F");
    eff->Branch("lamgenProdVtxY", &genProdVtxY, "genProdVtxY/F");
    eff->Branch("lamgenProdVtxZ", &genProdVtxZ, "genProdVtxZ/F");

    eff->Branch("lamgenDecVtxX", &genDecVtxX, "genDecVtxX/F");
    eff->Branch("lamgenDecVtxY", &genDecVtxY, "genDecVtxY/F");
    eff->Branch("lamgenDecVtxZ", &genDecVtxZ, "genDecVtxZ/F");

    eff->Branch("lamv0MatchStatus", &v0MatchStatus, "v0MatchStatus/I");

    eff->Branch("lamgenPDG", &genPDG, "genPDG/I");
    eff->Branch("lamgenPosDauPDG", &genPosDauPDG, "genPosDauPDG/I");
    eff->Branch("lamgenNegDauPDG", &genNegDauPDG, "genNegDauPDG/I");
    eff->Branch("lamgenPosDauPx", &genPosDauPx, "genPosDauPx/F");
    eff->Branch("lamgenPosDauPy", &genPosDauPy, "genPosDauPy/F");
    eff->Branch("lamgenPosDauPz", &genPosDauPz, "genPosDauPz/F");
    eff->Branch("lamgenNegDauPx", &genNegDauPx, "genNegDauPx/F");
    eff->Branch("lamgenNegDauPy", &genNegDauPy, "genNegDauPy/F");
    eff->Branch("lamgenNegDauPz", &genNegDauPz, "genNegDauPz/F");
  }
}


// ------------ method called once each job just after ending the event loop  ------------
void V0DeltaREffAnalyzer::endJob() {
}

bool V0DeltaREffAnalyzer::aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers)
{
 int negTowers = 0;
 int posTowers = 0;
 for(CaloTowerCollection::const_iterator cal = towers->begin(); cal != towers->end(); ++cal) {
    for(unsigned int i = 0; i < cal->constituentsSize(); i++) {
       const DetId id = cal->constituent(i);
       if(id.det() == DetId::Hcal) {
         HcalSubdetector subdet=(HcalSubdetector(id.subdetId()));
         if(subdet == HcalForward) {
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

bool V0DeltaREffAnalyzer::scrapingTest( const edm::Event& iEvent, 
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

  reco::TrackCollection::const_iterator itk = tkColl->begin();
  reco::TrackCollection::const_iterator itk_e = tkColl->end();
  for(;itk!=itk_e;++itk){
    // std::cout << "HighPurity?  " << itk->quality(_trackQuality) << std::endl;
    if(itk->quality(_trackQuality)) numhighpurity++;
    if(itk->quality(_trackQuality) && itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numhighpurity_pt_100++;
    if(itk->pt() > 0.1 && fabs(itk->eta()) < 2.5) numloose_pt_100++;
  }
  fraction = (float)numhighpurity/(float)tkColl->size();

  accepted = (tkColl->size()<=numtrack)||(fraction>thresh);

  edm::Handle<reco::TrackCollection> tkH;
  edm::Handle<reco::BeamSpot> bsH;
  iEvent.getByLabel("offlineBeamSpot", bsH);
  edm::ESHandle<TransientTrackBuilder> theTTB;
  edm::ESHandle<MagneticField> bFieldH;
  edm::ESHandle<GlobalTrackingGeometry> globTkGeomH;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTB);
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomH);
  try{
    iEvent.getByLabel("generalPlusMinBiasTracks", tkH);
  }
  catch(...) {

   // std::cout << "Problem finding track collection." << endl;
    std::cout << "Problem finding track collection." << std::endl;	//by Qiao
    return accepted;
  }

  using namespace reco;
  TSCBLBuilderNoMaterial blsBuilder;
  for (unsigned int ndx = 0; ndx < tkH->size(); ndx++) {
    TrackRef tmpRef(tkH, ndx);
    TransientTrack tmpTk( *tmpRef, &(*bFieldH), globTkGeomH);
 //   TrajectoryStateTransform theTransform;
 //   FreeTrajectoryState initialFTS = theTransform.initialFreeState(*tmpRef, bFieldH.product());
    FreeTrajectoryState initialFTS =trajectoryStateTransform::initialFreeState(*tmpRef, bFieldH.product());	//by Qiao
    TrajectoryStateClosestToBeamLine tscb(blsBuilder(initialFTS, *bsH));
    if (tscb.isValid()) {
      double d0 = tscb.transverseImpactParameter().value();
      if (d0 < 0.4) {
	nPromptTracks++;
	if (tmpRef->pt() > 0.1 && fabs(tmpRef->eta()) < 2.5) nPromptTracks_pt_100++;
      }
    }
    nTracks++;
    if (tmpRef->pt() > 0.1 && fabs(tmpRef->eta()) < 2.5) nTracks_pt_100++;
      
      
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

//define this as a plug-in
DEFINE_FWK_MODULE(V0DeltaREffAnalyzer);
