// -*- C++ -*-
//
// Package:    PiLambda
// Class:      PiLambda
// 
/**\class PiLambda PiLambda.cc myAnalyzers/PiLambda/src/PiLambda.cc

 Description: <one line class summary>
Make rootTuple for Lambda+track (Xi or Omega) analyses

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Keith Ulmer
//         Created:  Wed May  7 13:15:04 MDT 2008
// $Id: PiLambda.cc,v 1.32 2010/11/16 01:12:07 drell Exp $
//
//


// system include files
#include <memory>

// user include files
//#include "../interface/PiLambda.h" //by Qiao
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "TFile.h"
#include "TTree.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include <vector>
#include <utility>

#include "myAnalyzers/PiLambda/interface/PiLambda.h"
//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;
const double piMassConst_ = 0.13957018;
const double protonMassConst_ = 0.93827203;
const double ksMassConst_ = 0.497614;

//
// static data member definitions
//

//
// constructors and destructor
//
PiLambda::PiLambda(const edm::ParameterSet& iConfig)
  :
  hltTag(iConfig.getParameter<edm::InputTag>("hltTag")),
  hepMC(iConfig.getUntrackedParameter<std::string>("HepMC",std::string("source")) ),
  muonType(iConfig.getUntrackedParameter("muonType",std::string("muons"))),
  vtxSample( iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices")) ), 
  v0Collection( iConfig.getUntrackedParameter<std::string>("V0Collection",std::string("generalV0Candidates")) ),
  trackCollectionName( iConfig.getUntrackedParameter<edm::InputTag>("trackCollection", edm::InputTag("generalTracks")) ),
  doMC(iConfig.getUntrackedParameter<bool>("DoMC",false) ),
  debug(iConfig.getUntrackedParameter<bool>("Debug",false) ),

  tree_(0), nEvents(0), nXI(0), run(0), event(0), 
  trigHF(0), HFpE(0), HFmE(0), trigTech40_41(0), trigTech36_39(0), trigTech34(0), trigScraping(0), trigAlgo124(0), trigHLTminBias(0),
  nLooseTracks(0), nHighPurTracks(0), nLooseTracks_pt_100(0), nHighPurTracks_pt_100(0),
  processType(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0), priNTrk(0),
  priNTrk_weightGT_0_5(0),
  bsX(0), bsY(0), bsZ(0), bsXE(0), bsYE(0), bsZE(0),
  priRefitVtxX(0), priRefitVtxY(0), priRefitVtxZ(0), priRefitVtxXE(0), priRefitVtxYE(0), priRefitVtxZE(0), priRefitVtxCL(0), priRefitNTrk(0),
  xiMass(0), omMass(0), xiVtxCL(0), xiPx(0), xiPy(0), xiPz(0), xiPxE(0), xiPyE(0), xiPzE(0), 
  xictauB(0), xictauBE(0), xictau3D(0), xictau3DE(0), xictauMPV(0), xictauMPVE(0),
  xiDecayVtxX(0), xiDecayVtxY(0), xiDecayVtxZ(0), xiDecayVtxXE(0), xiDecayVtxYE(0), xiDecayVtxZE(0),
  VMass(0), VCandMass(0), VMassError(0), VVtxCL(0), VPx(0), VPy(0), VPz(0), VPalongXi(0),
  VKsMass(0), VKsMassError(0),
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0),
  VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  batPiPx(0), batPiPy(0), batPiPz(0), batPiEta(0), batPiPhi(0), batPiD0(0), batPiD0E(0), batPiDz(0), batPiDzE(0), batPiPVweight(0), batPiQ(0), 
  batPiNValidHits(0), batPiTrkQual(0), nTracks(0), 
  VTrkPMass(0), VTrkPPx(0), VTrkPPy(0), VTrkPPz(0), VTrkPEta(0), VTrkPPhi(0),
  VTrkPD0(0), VTrkPD0E(0), VTrkPDz(0), VTrkPDzE(0), VTrkPPVweight(0), VTrkPQ(0), VTrkPNValidHits(0),
  VTrkPiMass(0), VTrkPiPx(0), VTrkPiPy(0), VTrkPiPz(0), VTrkPiEta(0), VTrkPiPhi(0),
  VTrkPiD0(0), VTrkPiD0E(0), VTrkPiDz(0), VTrkPiDzE(0), VTrkPiPVweight(0), VTrkPiQ(0), VTrkPiNValidHits(0),
  VTransversePCAPrimary(0), VTransversePCAPrimaryError(0), VLongitudinalPCAPrimary(0), VLongitudinalPCAPrimaryError(0), 
  VTransversePCABeamSpot(0), VTransversePCABeamSpotError(0), VRSig2D(0), VRSig3D(0), VFLSig2D(0), VFLSig3D(0),
  XiTransversePCAPrimary(0), XiTransversePCAPrimaryError(0), XiLongitudinalPCAPrimary(0), XiLongitudinalPCAPrimaryError(0),
  XiTransversePCABeamSpot(0), XiTransversePCABeamSpotError(0), 
  Xi3DIpSig(0), XiFLsig3D(0),
  VCosThetaPAndveeVertexToPrimaryVector(0), VCosThetaPAndveeVertexToXiVertexVector(0), XiCosThetaPAndLineOfFlight(0),
  V3dIpWrtPrimary(0),V3dIpWrtPrimarySig(0), 
  VTrkPi2DIp(0),VTrkPi2DIpSig(0),VTrkPi3DIp(0),VTrkPi3DIpSig(0),
  VTrkP2DIp(0),VTrkP2DIpSig(0),VTrkP3DIp(0),VTrkP3DIpSig(0),
  batPi2DIp(0),batPi2DIpSig(0),batPi3DIp(0),batPi3DIpSig(0),
  pionInV0(0),
  
  nGenXi(0), genBatPiEta(0), genBatPiPhi(0), genLamPiEta(0), genLamPiPhi(0), genLamPEta(0), genLamPPhi(0),
  genXiP(0), genXiEta(0), genXiL(0), genXiR(0), genXiDecayVX(0), genXiDecayVY(0),
  genXiDecayVZ(0), genXiMomL(0), genXiQ(0),
  genXiMotherPDG(0), genXiProdVtxX(0), genXiProdVtxY(0), genXiProdVtxZ(0)

{
   //now do what ever initialization is needed
}


PiLambda::~PiLambda()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void


PiLambda::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  //The mass of a muon and the insignificant mass sigma        
  //to avoid singularities in the covariance matrix.	       
  ParticleMass pion_mass = 0.13957018;			       
  ParticleMass proton_mass = 0.938272013;		       
  ParticleMass kaon_mass = 0.493677;			       
  float pion_sigma = pion_mass*1.e-6;			       
  float proton_sigma = proton_mass*1.e-6;		       
  float kaon_sigma = kaon_mass*1.e-6;			       
  ParticleMass lambda_mass = 1.115683;			       
  float lambda_sigma = 0.000006;			       

  if( doMC ) {
    Handle<HepMCProduct> hepEv;
    iEvent.getByType( hepEv );
    processType = hepEv->GetEvent()->signal_process_id();
  }
  
  // get event content information
  
  nEvents++;
  
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  
  ESHandle<TransientTrackBuilder> theTTB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB);
  
   // get primary vertex
   Handle<reco::VertexCollection> recVtxs;
   iEvent.getByLabel(vtxSample, recVtxs);
   unsigned int nVtxTrks = 0;
   reco::VertexCollection::const_iterator bestVtx = recVtxs->begin();
   //if(debug) cout <<"total number of primary vertices = " << recVtxs->size() << endl;
   for (reco::VertexCollection::const_iterator vtx = recVtxs->begin();
	vtx != recVtxs->end(); ++vtx){
     if (nVtxTrks < vtx->tracksSize() ) {
       nVtxTrks = vtx->tracksSize();
       bestVtx = vtx;
     }
   }

   int nVtxTk_weight = 0;
   for( Vertex::trackRef_iterator itk = bestVtx->tracks_begin();
	itk != bestVtx->tracks_end();
	itk++ ) {
     if( bestVtx->trackWeight( *itk ) > 0.5 )
       nVtxTk_weight++;
   }
   priNTrk_weightGT_0_5 = nVtxTk_weight;

   //get primary with beamspot constraint
   Handle<reco::VertexCollection> recVtxsBS;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", recVtxsBS);
   
   nVtxTrks = 0;
   reco::VertexCollection::const_iterator bestVtxBS = recVtxsBS->begin();
   for (reco::VertexCollection::const_iterator vtxBS = recVtxsBS->begin();
	vtxBS != recVtxsBS->end(); ++vtxBS){
     if (nVtxTrks < vtxBS->tracksSize() ) {
       nVtxTrks = vtxBS->tracksSize();
       bestVtxBS = vtxBS;
     }
   }

   //get beamspot
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
   if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
   else  cout <<"No beam spot available from EventSetup" << endl;

   const GlobalPoint theBeamSpot = GlobalPoint(beamSpotHandle->position().x(),
 	  				       beamSpotHandle->position().y(),
	  				       beamSpotHandle->position().z());

   run = iEvent.id().run();
   event = iEvent.id().event();

   if(debug) cout <<"processing event = " << event << ", run = " << run << endl;

   //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   // try reconstruction without fitting modules
   //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Handle<vector<VertexCompositeCandidate> > theV0Handle;
   iEvent.getByLabel(v0Collection, "Lambda", theV0Handle);

   //CombinatoriaTrack
   edm::Handle<reco::TrackCollection> trackCollection;
   iEvent.getByLabel(trackCollectionName,trackCollection);

   // first check HF MinBias trigger
   // get calo towers
   edm::Handle<CaloTowerCollection> theCaloTowers;
   iEvent.getByLabel("towerMaker",theCaloTowers);
   if (aboveTriggerThreshold( theCaloTowers )) trigHF = 1;
   else trigHF = 0;

  // Check HLT_L1_BscMinBiasOR_BptxPlusORMinus trigger for prescale lumi calculation
  edm::Handle<edm::TriggerResults> hltResults;
  try {
    iEvent.getByLabel(hltTag, hltResults);//edm::InputTag("TriggerResults::REDIGI36"), hltResults);
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
    
    const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltResults);

    for( int ndx = 0; ndx < ntrigs; ndx++ ) {
      TString trigName = trigNames.triggerName(ndx);
      int hltflag = (*hltResults)[ndx].accept();
      //cout << "Trigger " << trigName << " was passed = " << hltflag << endl;
      if( trigName == "HLT_L1_BscMinBiasOR_BptxPlusORMinus" ) {
	trigHLTminBias = hltflag;
      }
    }
  }
   
   // next check BSC MinBias trigger bit (40||41)

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

   bool isNotScraping = scrapingTest( iEvent, iSetup );
   if( isNotScraping ) {
     trigScraping = 1;
   }
   else
     trigScraping = 0;

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

   if(debug) cout <<"size of Lambda collection = " << theV0Handle->size() << " and size of track collection = " << trackCollection->size() << endl;

   if ( theV0Handle->size()>0 && trackCollection->size()>=1 ) {
     for ( vector<VertexCompositeCandidate>::const_iterator iVee = theV0Handle->begin();
	   iVee != theV0Handle->end(); ++iVee ) {
	   
       // apply some early cuts (looser than final selection) to reduce ntuple size
	
       // check lambda mass to be within +- 20MeV...
       float massWindow = 0.020;
       if(iVee->mass() > lambda_mass + massWindow || iVee->mass() < lambda_mass - massWindow) continue;
       // after finding values also cut on xi vertex probability > 0.0001 and batPion Dxyz/sigma > 0.5

       int lamPDG = iVee->pdgId();
       bool lamIsParticle = (lamPDG > 0);
       
       //get Ks tracks from V0 candidate
       vector<RecoChargedCandidate> v0daughters;
       vector<TrackRef> theDaughterTracks;

       //check momentum and add pion first, proton second
       if (iVee->daughter(0)->momentum().mag2() < iVee->daughter(1)->momentum().mag2()){
	 v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				  (iVee->daughter(0))) );
	 v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				  (iVee->daughter(1))) );
       } else {
	 v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				  (iVee->daughter(1))) );
	 v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
				  (iVee->daughter(0))) );
       }
       for(unsigned int j = 0; j < v0daughters.size(); ++j) {
	 theDaughterTracks.push_back(v0daughters[j].track());
       }

       vector<TransientTrack> tracksForKalmanFit;
       for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
	   tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
       }

       KalmanVertexFitter theKsFitter(true);
       CachingVertex<5> ksVtxForInvMass = theKsFitter.vertex(tracksForKalmanFit);
       vector<double> ksDauMasses;
       for (int ndx = 0; ndx < 2; ndx++) {
	   ksDauMasses.push_back(piMassConst_);
       }

       vector<double> lamDauMasses;
       if (lamIsParticle && tracksForKalmanFit[0].charge() > 0) {
	   lamDauMasses.push_back(protonMassConst_);
	   lamDauMasses.push_back(piMassConst_);
       } else if (!lamIsParticle && tracksForKalmanFit[0].charge() < 0) {
	   lamDauMasses.push_back(protonMassConst_);
	   lamDauMasses.push_back(piMassConst_);
       } else {
	   lamDauMasses.push_back(piMassConst_);
	   lamDauMasses.push_back(protonMassConst_);
       }

       Measurement1D ksMassWErr = massCalculator.invariantMass(ksVtxForInvMass, ksDauMasses);
       double lamOtherMassVal = ksMassWErr.value();
       double lamOtherMassErrorVal = ksMassWErr.error();

       Measurement1D lamMassWErr = massCalculator.invariantMass(ksVtxForInvMass, lamDauMasses);
       double lamMassVal = lamMassWErr.value();
       double lamMassErrorVal = lamMassWErr.error();

       //if(debug) cout <<"using V0 with mass = " << iVee->mass() << endl;

       for(unsigned int indx = 0; indx < trackCollection->size(); indx++) {
	 TrackRef glbTrack( trackCollection, indx );
         // get at least 3 hits for tracks
	 
	 if ( glbTrack.isNull() ) {
	   if(debug) cout <<"continue due to no track ref" << endl;
	   continue;
	 }
	
	 if(glbTrack->numberOfValidHits() < 3) continue;
	 
	 double weightVP = 0;
	 double weightVPi = 0;
	 double weightXiPi = 0;
	 int trkQuality  = -1;

	 // check Track Quality:
	 reco::TrackBase::TrackQuality trackQualityUndef         =  reco::TrackBase::qualityByName("undefQuality");
	 reco::TrackBase::TrackQuality trackQualityLoose         =  reco::TrackBase::qualityByName("loose");
	 reco::TrackBase::TrackQuality trackQualityTight         =  reco::TrackBase::qualityByName("tight");
	 reco::TrackBase::TrackQuality trackQualityhighPur       =  reco::TrackBase::qualityByName("highPurity");

	 if(glbTrack->quality(trackQualityUndef))           trkQuality = 5;
	 if(glbTrack->quality(trackQualityLoose))          trkQuality = 0; 
	 if(glbTrack->quality(trackQualityTight))          trkQuality = 1; 
	 if(glbTrack->quality(trackQualityhighPur))        trkQuality = 2; 
 
	 bool match = false;
	 // check if the pion is used in any v0 candidate in the collection and flag it
	 for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
	   if (glbTrack->charge() == theDaughterTracks[j]->charge() && 
	       glbTrack->momentum() == theDaughterTracks[j]->momentum() ) {
	     //if(debug) cout <<"Match found between batPi and V0 track early with p = " << glbTrack->momentum() << " and " << theDaughterTracks[j]->momentum() << endl;
	     match = true;
	   }
	   if (match) break;
	 }
	 if (match) { 
	   if(debug) cout <<"Match found between batPi and V0 track" << endl;
	   continue; 
	 } // Track is already used in making the V0
	
	 pionInV0->push_back(trackInV0(glbTrack, (*theV0Handle))); // check if pion is in *any* good V0
	 //if(debug) cout <<"have 3 good tracks" << endl;

	 TransientTrack pionTT(theDaughterTracks[0], &(*bFieldHandle) );
	 TransientTrack protonTT(theDaughterTracks[1], &(*bFieldHandle) );
	 TransientTrack batPionTT(glbTrack, &(*bFieldHandle) );

	 if (!batPionTT.isValid()) {
	   if(debug) cout <<"continue from invalid batPionTT" << endl;
	   continue;
	 }

       	 //Creating a KinematicParticleFactory
	 KinematicParticleFactoryFromTransientTrack pFactory;
	 
	 //initial chi2 and ndf before kinematic fits.
	 float chi = 0.;
	 float ndf = 0.;
	 vector<RefCountedKinematicParticle> VeeParticles;
	 VeeParticles.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));
	 VeeParticles.push_back(pFactory.particle(protonTT,proton_mass,chi,ndf,proton_sigma));
	 
	 KinematicParticleVertexFitter fitter;   
	 RefCountedKinematicTree veeVertexFitTree;
	 veeVertexFitTree = fitter.fit(VeeParticles); 
	 if (!veeVertexFitTree->isValid()) {
	   if(debug) cout <<"invalid vertex from the ks vertex fit" << std::endl;
	   continue; 
	 }
	 veeVertexFitTree->movePointerToTheTop();

	 RefCountedKinematicParticle vee_vFit_noMC = veeVertexFitTree->currentParticle();
	 RefCountedKinematicVertex vee_vFit_vertex_noMC = veeVertexFitTree->currentDecayVertex();
	 
	 veeVertexFitTree->movePointerToTheFirstChild();
	 RefCountedKinematicParticle veePi1 = veeVertexFitTree->currentParticle();
	 veeVertexFitTree->movePointerToTheNextChild();
	 RefCountedKinematicParticle veePi2 = veeVertexFitTree->currentParticle();

	 KinematicParameters veePiKP = veePi1->currentState().kinematicParameters();
	 KinematicParameters veePKP = veePi2->currentState().kinematicParameters();

	 // now apply Lambda mass constraint
	 KinematicParticleFitter csFitterVee;
	 KinematicConstraint * lam_c = new MassKinematicConstraint(lambda_mass,lambda_sigma);
	 // add mass constraint to the lambda fit to do a constrained fit:  

	 veeVertexFitTree->movePointerToTheTop();
	 veeVertexFitTree = csFitterVee.fit(lam_c,veeVertexFitTree);
	 if (!veeVertexFitTree->isValid()){
	   if(debug) cout <<"caught an exception in the lambda mass constraint fit" << std::endl;
	   continue; 
	 }
	 veeVertexFitTree->movePointerToTheTop();
	 RefCountedKinematicParticle vee_vFit_withMC = veeVertexFitTree->currentParticle();
	 
	 if (!vee_vFit_withMC->currentState().isValid()) {
	   if(debug) cout <<"continue from invalid state from vee fit" << endl;
	   continue;
	 }

	 vector<RefCountedKinematicParticle> xiFitParticles;
	 vector<RefCountedKinematicParticle> omFitParticles;
	
	 xiFitParticles.push_back(pFactory.particle(batPionTT,pion_mass,chi,ndf,pion_sigma)); 
	 omFitParticles.push_back(pFactory.particle(batPionTT,kaon_mass,chi,ndf,kaon_sigma));
	 
	 xiFitParticles.push_back(vee_vFit_withMC);
	 omFitParticles.push_back(vee_vFit_withMC);
	 
	 RefCountedKinematicTree xiFitTree = fitter.fit(xiFitParticles);  
	 if (!xiFitTree->isValid()) {
	   if(debug) cout <<"caught an exception in the Xi vertex fit" << std::endl;
	   continue;
	 }
	 
	 
	 RefCountedKinematicTree omFitTree = fitter.fit(omFitParticles);  
	 if (!omFitTree->isValid()) {
	   if(debug) cout <<"caught an exception in the Xi vertex fit" << std::endl;
	   continue;
	 }

	 xiFitTree->movePointerToTheTop();
	 RefCountedKinematicParticle xiCand = xiFitTree->currentParticle();
	 
	 omFitTree->movePointerToTheTop();
	 RefCountedKinematicParticle omCand = omFitTree->currentParticle();
	 
	 RefCountedKinematicVertex xiDecayVertex = xiFitTree->currentDecayVertex();
	 if (!xiDecayVertex->vertexIsValid()){
	   if(debug) cout <<"Xi MC fit vertex is not valid" << endl;
	   continue;
	 }
	 
	 if ( xiDecayVertex->chiSquared()<0 || xiDecayVertex->chiSquared()>1000 ) {
	   if(debug) cout <<" failed chi2 cut in fit with value = " << xiDecayVertex->chiSquared() << endl;
	   continue;
	 }
	 
         float xiC2Prob =
	    ChiSquaredProbability((double)(xiDecayVertex->chiSquared()),(double)(xiDecayVertex->degreesOfFreedom()));
         if (xiC2Prob < 0.0001) continue;
	 
	 if ( xiCand->currentState().mass() > 3 ) {
	   if(debug) cout <<"continue from mass > 3 in Xi fit with mass = " << xiCand->currentState().mass() << endl;
	   continue;
	 }
	 
	 // get children from final Xi fit
	 xiFitTree->movePointerToTheFirstChild();
	 RefCountedKinematicParticle batPionCand = xiFitTree->currentParticle();
	 xiFitTree->movePointerToTheNextChild();
	 RefCountedKinematicParticle ksCandMC = xiFitTree->currentParticle();
	 
	 // get batchlor pion and V0 parameters from Xi fit
	 KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
	 
	 KinematicParameters VCandKP = ksCandMC->currentState().kinematicParameters();

	 // check for Xi tracks in the primary vertex 
	 vector<reco::TransientTrack> newVertexTracks;	 
	 
	 for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtx->tracks_begin();
	       iTrack != bestVtx->tracks_end(); ++iTrack) {
	   // compare primary tracks to check for matches with Xi cand
	   TrackRef trackRef = iTrack->castTo<TrackRef>();
	   
	   // the 3 tracks in the Xi cand are theDaughterTracks[0] theDaughterTracks[1] glbTrack
	   if ( theDaughterTracks[0]==trackRef ) {
	     weightVP = bestVtx->trackWeight ( *iTrack );
	   }
	   if ( theDaughterTracks[1]==trackRef ) {
	     weightVPi = bestVtx->trackWeight ( *iTrack );
	   }
	   if (	glbTrack==trackRef ) {
	     // let's get the weight for the track.
	     weightXiPi = bestVtx->trackWeight ( *iTrack );
	   }
	   
	   // try refitting the primary without the tracks in the Xi reco candidate	   
	   if ( !((theDaughterTracks[0]==trackRef) || (theDaughterTracks[1]==trackRef) ||
		  (glbTrack==trackRef)) ) {
	     TransientTrack tt(trackRef, &(*bFieldHandle) );
	     newVertexTracks.push_back(tt);
	   }
	   

	 }

	 //if(debug) cout <<"number of tracks in primary = " << bestVtx->tracksSize() << endl;
	 //if(debug) cout <<"number of tracks from vertex tt vector = " << newVertexTracks.size() << endl;

	 Vertex refitVtx;

	 // if no tracks in primary or no reco track included in primary then don't do anything
	 GlobalPoint vert;
	 if (  newVertexTracks.size()>0 && (bestVtx->tracksSize()!=newVertexTracks.size()) ) {
	   
	   AdaptiveVertexFitter theFitter;
	   TransientVertex v = theFitter.vertex(newVertexTracks);//, beamSpot);
	   if ( v.isValid() ) {
	     // calculate the impact parameters in 2D and 3D for the Xi candidate tracks wrt the new primary vertex
	     vert = v.position();
	     refitVtx = v;
	   } else {
	     // use original PV called bestVtx, which is a reco::Vertex
	     vert = GlobalPoint(bestVtx->x(),bestVtx->y(),bestVtx->z());
	     refitVtx = *bestVtx;
	   }
	 } else {
	   // use original PV called bestVtx, which is a reco::Vertex
	   vert = GlobalPoint(bestVtx->x(),bestVtx->y(),bestVtx->z());
	   refitVtx = *bestVtx;
	 }

         // get vertex displacement here 
         std::vector<double> vVtxEVec;
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().cxx() );
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().cyx() );
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().cyy() );
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().czx() );
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().czy() );
         vVtxEVec.push_back( vee_vFit_vertex_noMC->error().czz() );
         SMatrixSym3D vVtxCov(vVtxEVec.begin(), vVtxEVec.end());

         std::vector<double> xiVtxEVec;
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().cxx() );
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().cyx() );
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().cyy() );
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().czx() );
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().czy() );
         xiVtxEVec.push_back( vee_vFit_vertex_noMC->error().czz() );
         SMatrixSym3D xiVtxCov(xiVtxEVec.begin(), xiVtxEVec.end());

         SMatrixSym3D totalCov = vVtxCov + refitVtx.covariance();
	 SMatrixSym3D lamCov = vVtxCov + xiVtxCov;
	 SMatrixSym3D xiCov = xiVtxCov + refitVtx.covariance();

         SVector3 distanceVector2D(  vee_vFit_vertex_noMC->position().x() - refitVtx.position().x(),
	  		             vee_vFit_vertex_noMC->position().y() - refitVtx.position().y(),
			             0. );
         double rVtxMag2D = ROOT::Math::Mag(distanceVector2D);
         double sigmaRvtxMag2D = 
	   sqrt(ROOT::Math::Similarity(totalCov, distanceVector2D)) / rVtxMag2D;

         SVector3 distanceVector3D( vee_vFit_vertex_noMC->position().x() - refitVtx.position().x(),
	  			    vee_vFit_vertex_noMC->position().y() - refitVtx.position().y(),
				    vee_vFit_vertex_noMC->position().z() - refitVtx.position().z() );
         double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
         double sigmaRvtxMag3D =
	   sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / rVtxMag3D;

	 TrajectoryStateClosestToPoint  trajPion = pionTT.trajectoryStateClosestToPoint( vert );
	 double pionD0 = trajPion.perigeeParameters().transverseImpactParameter();
	 double pionD0E = trajPion.perigeeError().transverseImpactParameterError();
	 double pionDz = trajPion.perigeeParameters().longitudinalImpactParameter();
	 double pionDzE = trajPion.perigeeError().longitudinalImpactParameterError();
	 pair<bool,Measurement1D> thePion3DIpPair = IPTools::absoluteImpactParameter3D(pionTT, refitVtx);
	 float VTrkPi2DIpValue = -1000.;
	 float VTrkPi2DIpSigValue = -1000.;
	 float VTrkPi3DIpValue = -1000.;
	 float VTrkPi3DIpSigValue = -1000.;
	 if(thePion3DIpPair.first){
	   VTrkPi3DIpValue = thePion3DIpPair.second.value();
	   VTrkPi3DIpSigValue = thePion3DIpPair.second.significance();
	 }

	 TrajectoryStateClosestToPoint trajPionBS = pionTT.trajectoryStateClosestToPoint(theBeamSpot);
	 if(trajPionBS.isValid()){
	   VTrkPi2DIpValue = trajPionBS.perigeeParameters().transverseImpactParameter();
	   if(trajPionBS.hasError()) VTrkPi2DIpSigValue = VTrkPi2DIpValue/trajPionBS.perigeeError().transverseImpactParameterError();
	 }

	 TrajectoryStateClosestToPoint  trajProton = protonTT.trajectoryStateClosestToPoint( vert );
	 double protonD0 = trajProton.perigeeParameters().transverseImpactParameter();
	 double protonD0E = trajProton.perigeeError().transverseImpactParameterError();
	 double protonDz = trajProton.perigeeParameters().longitudinalImpactParameter();
	 double protonDzE = trajProton.perigeeError().longitudinalImpactParameterError();
	 pair<bool,Measurement1D> theProton3DIpPair = IPTools::absoluteImpactParameter3D(protonTT, refitVtx);
	 float VTrkP2DIpValue = -1000.;
	 float VTrkP2DIpSigValue = -1000.;
	 float VTrkP3DIpValue = -1000.;
	 float VTrkP3DIpSigValue = -1000.;
	 if(theProton3DIpPair.first){
	   VTrkP3DIpValue = theProton3DIpPair.second.value();
	   VTrkP3DIpSigValue = theProton3DIpPair.second.significance();
	 } 

	 TrajectoryStateClosestToPoint trajProtonBS = protonTT.trajectoryStateClosestToPoint(theBeamSpot);
	 if(trajProtonBS.isValid()){
	   VTrkP2DIpValue = trajProtonBS.perigeeParameters().transverseImpactParameter();
	   if(trajProtonBS.hasError()) VTrkP2DIpSigValue = VTrkP2DIpValue/trajProtonBS.perigeeError().transverseImpactParameterError();
	 }
	 
	 TrajectoryStateClosestToPoint  trajBatPion = batPionTT.trajectoryStateClosestToPoint( vert );
	 double batPionD0 = trajBatPion.perigeeParameters().transverseImpactParameter();
	 double batPionD0E = trajBatPion.perigeeError().transverseImpactParameterError();
	 double batPionDz = trajBatPion.perigeeParameters().longitudinalImpactParameter();
	 double batPionDzE = trajBatPion.perigeeError().longitudinalImpactParameterError();
	 pair<bool,Measurement1D> theBatPi3DIpPair = IPTools::absoluteImpactParameter3D(batPionTT, refitVtx);
	 float batPi2DIpValue = -1000.;
	 float batPi2DIpSigValue = -1000.;
	 float batPi3DIpValue = -1000.;
	 float batPi3DIpSigValue = -1000.;
	 if(theBatPi3DIpPair.first){
	   batPi3DIpValue = theBatPi3DIpPair.second.value();
	   batPi3DIpSigValue = theBatPi3DIpPair.second.significance();
	 }

         if (batPi3DIpSigValue < 1.0) continue;

	 TrajectoryStateClosestToPoint trajBatPionBS = batPionTT.trajectoryStateClosestToPoint(theBeamSpot);
	 if(trajBatPionBS.isValid()){
	   batPi2DIpValue = trajBatPionBS.perigeeParameters().transverseImpactParameter();
	   if(trajBatPionBS.hasError()) batPi2DIpSigValue = batPi2DIpValue/trajBatPionBS.perigeeError().transverseImpactParameterError();
	 }

	 //Get kinematic state from kinematic particle for V
	 KinematicState veeKinematicState = vee_vFit_noMC->currentState();
	 //Get trajectory state from kinematic state
	 FreeTrajectoryState theVFTS = veeKinematicState.freeTrajectoryState();
	 //Get TransientTrack from FTS
	 TransientTrack VTT = (*theTTB).build(theVFTS);
	 // get 3d ip wrt to vertex, the sign is wrt to line of flight vector
	 GlobalVector veeLineOfFlight = GlobalVector((*xiDecayVertex).position().x()  - iVee->vertex().x(),
						     (*xiDecayVertex).position().y()  - iVee->vertex().y(),
						     (*xiDecayVertex).position().z()  - iVee->vertex().z());
	 pair<bool,Measurement1D> theVee3DIpPair = IPTools::signedImpactParameter3D(VTT, veeLineOfFlight, refitVtx);
	 float V3dIpWrtPrimaryValue = -1000.;
	 float V3dIpWrtPrimarySigValue = -1000.;
	 if(theVee3DIpPair.first){
	    V3dIpWrtPrimaryValue = theVee3DIpPair.second.value();
	    V3dIpWrtPrimarySigValue = theVee3DIpPair.second.significance();
	 }

	 
	 //Get TrajectoryStateClosestToPoint from TT
	 const GlobalPoint refitVertexGPoint = GlobalPoint(refitVtx.x(), refitVtx.y(), refitVtx.z());
	 
	 float VTransversePCAPrimaryValue = -1000.;
	 float VTransversePCAPrimaryErrorValue = -1000.;
	 float VLongitudinalPCAPrimaryValue = -1000.;
	 float VLongitudinalPCAPrimaryErrorValue = -1000.;
	 
	 float VTransversePCABeamSpotValue = -1000.;
	 float VTransversePCABeamSpotErrorValue = -1000.;
	  
	 if(VTT.isValid()){
	   TrajectoryStateClosestToPoint  theVTrjPCAPrimary = VTT.trajectoryStateClosestToPoint(refitVertexGPoint);
	   if(theVTrjPCAPrimary.isValid()){
	     VTransversePCAPrimaryValue = theVTrjPCAPrimary.perigeeParameters().transverseImpactParameter();
	     VLongitudinalPCAPrimaryValue = theVTrjPCAPrimary.perigeeParameters().longitudinalImpactParameter();
	     if(theVTrjPCAPrimary.hasError()){
	       VTransversePCAPrimaryErrorValue = theVTrjPCAPrimary.perigeeError().transverseImpactParameterError();
	       VLongitudinalPCAPrimaryErrorValue = theVTrjPCAPrimary.perigeeError().longitudinalImpactParameterError();
	     }
	   }
	   TrajectoryStateClosestToPoint theVTrjPCABeamSpot = VTT.trajectoryStateClosestToPoint(theBeamSpot);
	   if(theVTrjPCABeamSpot.isValid()){
	     VTransversePCABeamSpotValue = theVTrjPCABeamSpot.perigeeParameters().transverseImpactParameter();
	     if(theVTrjPCABeamSpot.hasError())VTransversePCABeamSpotErrorValue = theVTrjPCABeamSpot.perigeeError().transverseImpactParameterError();
	   }
	 }
	 
         // Calculate some variables to determine if the Lam vertex is downstream of the Xi vertex
	 
	 // calculate line of flight between primary and secondary
	 GlobalVector veeVertexToPrimaryVector = GlobalVector ( vee_vFit_vertex_noMC->position().x() - refitVtx.x(),
	  			      		                vee_vFit_vertex_noMC->position().y() - refitVtx.y(),
	 			      		                vee_vFit_vertex_noMC->position().z() - refitVtx.z());
	
	 // get angle between momentum vector and line of flight		      
	 float theAngle = angle(veeVertexToPrimaryVector.x(),veeVertexToPrimaryVector.y(),veeVertexToPrimaryVector.z(),
 	 		        float(VCandKP.momentum().x()),float(VCandKP.momentum().y()),float(VCandKP.momentum().z()));	
	  
	 SVector3 lamFlightVector( vee_vFit_vertex_noMC->position().x() - xiDecayVertex->position().x(),
	  		           vee_vFit_vertex_noMC->position().y() - xiDecayVertex->position().y(),
			           vee_vFit_vertex_noMC->position().z() - xiDecayVertex->position().z() );
         double lamFlightMag = ROOT::Math::Mag(lamFlightVector);
         double lamFlightSigma = sqrt(ROOT::Math::Similarity(lamCov, lamFlightVector)) / lamFlightMag;

         float lamDotXi = xiCand->currentState().globalMomentum().dot( VCandKP.momentum() );
         float lamPalongXi = lamDotXi/xiCand->currentState().globalMomentum().mag();
         
	 float VFLSig;
         if ( lamDotXi > 0 ) VFLSig = lamFlightMag/lamFlightSigma ;
         else VFLSig = -1*lamFlightMag/lamFlightSigma ;

	 // get angle between momentum vector and line of flight		      
	 float theAngle2 = angle(float(lamFlightVector[0]),float(lamFlightVector[1]),float(lamFlightVector[2]),
 	 		        float(VCandKP.momentum().x()),float(VCandKP.momentum().y()),float(VCandKP.momentum().z()));	 

         // get Xi flight length significance					       
	 SVector3 xiFlightVector( xiDecayVertex->position().x() - refitVtx.x(),
	  		          xiDecayVertex->position().y() - refitVtx.y(),
			          xiDecayVertex->position().z() - refitVtx.z() );
         double xiFlightMag = ROOT::Math::Mag(xiFlightVector);
         double xiFlightSigma = sqrt(ROOT::Math::Similarity(xiCov, xiFlightVector)) / xiFlightMag;

         ////////////////////////////////////////////////////////////////
	 // fill candidate variables now
	 xiMass->push_back(xiCand->currentState().mass());
	 xiPx->push_back(xiCand->currentState().globalMomentum().x());
	 xiPy->push_back(xiCand->currentState().globalMomentum().y());
	 xiPz->push_back(xiCand->currentState().globalMomentum().z());
	 
	 omMass->push_back(omCand->currentState().mass());

         KinematicState theCurrentXiCandKinematicState = xiCand->currentState();
	 FreeTrajectoryState theXiFTS = theCurrentXiCandKinematicState.freeTrajectoryState();
	 TransientTrack xiTT = (*theTTB).build(theXiFTS);

	 float XiTransversePCAPrimaryValue = -1000.;
	 float XiTransversePCAPrimaryErrorValue = -1000.;
	 float XiLongitudinalPCAPrimaryValue = -1000.;
	 float XiLongitudinalPCAPrimaryErrorValue = -1000.;
	 float XiTransversePCABeamSpotValue = -1000.;
	 float XiTransversePCABeamSpotErrorValue = -1000.;
         float Xi3DIpSigValue = -1000;

	 if(xiTT.isValid()){
	   TrajectoryStateClosestToPoint theXiTrjPCAPrimary = xiTT.trajectoryStateClosestToPoint(refitVertexGPoint);
	   if(theXiTrjPCAPrimary.isValid()){
	     XiTransversePCAPrimaryValue = theXiTrjPCAPrimary.perigeeParameters().transverseImpactParameter();
	     XiLongitudinalPCAPrimaryValue = theXiTrjPCAPrimary.perigeeParameters().longitudinalImpactParameter();
	     if(theXiTrjPCAPrimary.hasError()){
	       XiTransversePCAPrimaryErrorValue = theXiTrjPCAPrimary.perigeeError().transverseImpactParameterError();
	       XiLongitudinalPCAPrimaryErrorValue = theXiTrjPCAPrimary.perigeeError().longitudinalImpactParameterError();
	     }
	   }
 	   pair<bool,Measurement1D> xi3DIpPair = IPTools::absoluteImpactParameter3D(xiTT, refitVtx);
	   if(xi3DIpPair.first){
	     //Xi3DIpValue = xi3DIpPair.second.value();
	     Xi3DIpSigValue = xi3DIpPair.second.significance();
	   }
	   TrajectoryStateClosestToPoint theXiTrjPCABeamSpot = xiTT.trajectoryStateClosestToPoint(theBeamSpot);
	   if(theXiTrjPCABeamSpot.isValid()){
	     XiTransversePCABeamSpotValue  = theXiTrjPCABeamSpot.perigeeParameters().transverseImpactParameter();
	     if(theXiTrjPCABeamSpot.hasError())XiTransversePCABeamSpotErrorValue = theXiTrjPCABeamSpot.perigeeError().transverseImpactParameterError();
	   }
	 }
	 
	 GlobalVector XiLineOfFlight = GlobalVector ((*xiDecayVertex).position().x()  - refitVtx.x(),
						     (*xiDecayVertex).position().y()  - refitVtx.y(),
						     (*xiDecayVertex).position().z()  - refitVtx.z());
	 float XiAngle = (XiLineOfFlight.x(),XiLineOfFlight.y(),XiLineOfFlight.z(),
	                 xiCand->currentState().globalMomentum().x(),
			 xiCand->currentState().globalMomentum().y(),
			 xiCand->currentState().globalMomentum().z());
	 
	 xiPxE->push_back( sqrt( xiCand->currentState().kinematicParametersError().matrix()(3,3) ) );
	 xiPyE->push_back( sqrt( xiCand->currentState().kinematicParametersError().matrix()(4,4) ) );
	 xiPzE->push_back( sqrt( xiCand->currentState().kinematicParametersError().matrix()(5,5) ) );
	 
	     
	 xiVtxCL->push_back( xiC2Prob );
	 xiDecayVtxX->push_back((*xiDecayVertex).position().x());
	 xiDecayVtxY->push_back((*xiDecayVertex).position().y());
	 xiDecayVtxZ->push_back((*xiDecayVertex).position().z());
	 xiDecayVtxXE->push_back( sqrt((*xiDecayVertex).error().cxx()) );
	 xiDecayVtxYE->push_back( sqrt((*xiDecayVertex).error().cyy()) );
	 xiDecayVtxZE->push_back( sqrt((*xiDecayVertex).error().czz()) );
	 
	 VDecayVtxX->push_back( vee_vFit_vertex_noMC->position().x() );
	 VDecayVtxY->push_back( vee_vFit_vertex_noMC->position().y() );
	 VDecayVtxZ->push_back( vee_vFit_vertex_noMC->position().z() );
	 
	 VDecayVtxXE->push_back( sqrt(vee_vFit_vertex_noMC->error().cxx()) );
	 VDecayVtxYE->push_back( sqrt(vee_vFit_vertex_noMC->error().cyy()) );
	 VDecayVtxZE->push_back( sqrt(vee_vFit_vertex_noMC->error().czz()) );
	 VVtxCL->push_back( ChiSquaredProbability((double)(vee_vFit_vertex_noMC->chiSquared()),(double)(vee_vFit_vertex_noMC->degreesOfFreedom())) );
	 
	 VPx->push_back( VCandKP.momentum().x() );
	 VPy->push_back( VCandKP.momentum().y() );
	 VPz->push_back( VCandKP.momentum().z() );
	 VPalongXi->push_back( lamPalongXi );
	 
	 VMass->push_back( veeKinematicState.mass() ); 

	 VCandMass->push_back(lamMassVal);
	 VMassError->push_back(lamMassErrorVal);

	 VKsMass->push_back(lamOtherMassVal);
	 VKsMassError->push_back(lamOtherMassErrorVal);
	 
	 VTrkPPx->push_back(veePKP.momentum().x());
	 VTrkPPy->push_back(veePKP.momentum().y());
	 VTrkPPz->push_back(veePKP.momentum().z());
	 VTrkPEta->push_back(veePKP.momentum().eta());
	 VTrkPPhi->push_back(veePKP.momentum().phi());
	 VTrkPMass->push_back(veePKP.mass());
	 VTrkPD0->push_back( protonD0 );
	 VTrkPD0E->push_back( protonD0E ); 
	 VTrkPDz->push_back( protonDz );
	 VTrkPDzE->push_back( protonDzE ); 
	 VTrkPQ->push_back( theDaughterTracks[1]->charge() );
	 VTrkPPVweight->push_back( weightVP );
	 VTrkPNValidHits->push_back( theDaughterTracks[1]->numberOfValidHits() );

	 VTrkPiPx->push_back(veePiKP.momentum().x());
	 VTrkPiPy->push_back(veePiKP.momentum().y());
	 VTrkPiPz->push_back(veePiKP.momentum().z());
	 VTrkPiEta->push_back(veePiKP.momentum().eta());
	 VTrkPiPhi->push_back(veePiKP.momentum().phi());
	 VTrkPiMass->push_back(veePiKP.mass());
	 VTrkPiD0->push_back( pionD0 );
	 VTrkPiD0E->push_back( pionD0E ); 
	 VTrkPiDz->push_back( pionDz );
	 VTrkPiDzE->push_back( pionDzE ); 
	 VTrkPiQ->push_back( theDaughterTracks[0]->charge() ); 
	 VTrkPiPVweight->push_back( weightVPi );
	 VTrkPiNValidHits->push_back( theDaughterTracks[0]->numberOfValidHits() );
	 
	 batPiPx->push_back(batPionKP.momentum().x());
	 batPiPy->push_back(batPionKP.momentum().y());
	 batPiPz->push_back(batPionKP.momentum().z());
         batPiEta->push_back(batPionKP.momentum().eta());
	 batPiPhi->push_back(batPionKP.momentum().phi());
	 batPiD0->push_back( batPionD0 );
	 batPiD0E->push_back( batPionD0E );
	 batPiDz->push_back( batPionDz );
	 batPiDzE->push_back( batPionDzE );
	 batPiPVweight->push_back( weightXiPi );
	 batPiQ->push_back(glbTrack->charge() );
	 batPiNValidHits->push_back(glbTrack->numberOfValidHits() );
	 batPiTrkQual->push_back( trkQuality );

	 nTracks->push_back( trackCollection->size() );

	 VCosThetaPAndveeVertexToPrimaryVector->push_back(cos(theAngle));
	 VCosThetaPAndveeVertexToXiVertexVector->push_back(cos(theAngle2));
	 XiCosThetaPAndLineOfFlight->push_back(cos(XiAngle));
	 
	 VTransversePCAPrimary->push_back(VTransversePCAPrimaryValue);
	 VTransversePCAPrimaryError->push_back(VTransversePCAPrimaryErrorValue);
	 VLongitudinalPCAPrimary->push_back(VLongitudinalPCAPrimaryValue);
	 VLongitudinalPCAPrimaryError->push_back(VLongitudinalPCAPrimaryErrorValue);
	 VTransversePCABeamSpot->push_back(VTransversePCABeamSpotValue);
	 VTransversePCABeamSpotError->push_back(VTransversePCABeamSpotErrorValue);
	 VRSig2D->push_back( rVtxMag2D/sigmaRvtxMag2D );
	 VRSig3D->push_back( rVtxMag3D/sigmaRvtxMag3D );
	 VFLSig2D->push_back( 0.0 );
	 VFLSig3D->push_back( VFLSig );
	 XiTransversePCAPrimary->push_back(XiTransversePCAPrimaryValue);
	 XiTransversePCAPrimaryError->push_back(XiTransversePCAPrimaryErrorValue);
	 XiLongitudinalPCAPrimary->push_back(XiLongitudinalPCAPrimaryValue);
	 XiLongitudinalPCAPrimaryError->push_back(XiLongitudinalPCAPrimaryErrorValue);
	 XiTransversePCABeamSpot->push_back(XiTransversePCABeamSpotValue);
	 XiTransversePCABeamSpotError->push_back(XiTransversePCABeamSpotErrorValue);
         Xi3DIpSig->push_back( Xi3DIpSigValue );
	 XiFLsig3D->push_back( xiFlightMag/xiFlightSigma );

         V3dIpWrtPrimary->push_back(V3dIpWrtPrimaryValue);
         V3dIpWrtPrimarySig->push_back(V3dIpWrtPrimarySigValue);
	 VTrkPi2DIp->push_back(VTrkPi2DIpValue);
	 VTrkPi2DIpSig->push_back(VTrkPi2DIpSigValue);
	 VTrkPi3DIp->push_back(VTrkPi3DIpValue);
	 VTrkPi3DIpSig->push_back(VTrkPi3DIpSigValue);
	 VTrkP2DIp->push_back(VTrkP2DIpValue);
	 VTrkP2DIpSig->push_back(VTrkP2DIpSigValue);	 
	 VTrkP3DIp->push_back(VTrkP3DIpValue);
	 VTrkP3DIpSig->push_back(VTrkP3DIpSigValue);
	 batPi2DIp->push_back(batPi2DIpValue);
	 batPi2DIpSig->push_back(batPi2DIpSigValue);
	 batPi3DIp->push_back(batPi3DIpValue);
	 batPi3DIpSig->push_back(batPi3DIpSigValue);


         // calculate Xi lifetime in various ways
         // first calculate ctau with ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)
	 float ximass = 1.32171;
	 float betagamma = (xiCand->currentState().globalMomentum().mag()/ximass);
	 GlobalPoint XiVP = GlobalPoint( xiDecayVertex->position() );
	 GlobalPoint PVP = GlobalPoint( refitVtx.position().x(), refitVtx.position().y(), refitVtx.position().z() );
	 GlobalVector sep3D = XiVP-PVP;
	 	       
	 GlobalVector pXiV = xiCand->currentState().globalMomentum();	      
	 float xictau_temp = (ximass*(sep3D.dot(pXiV)))/(pXiV.dot(pXiV));   
	 xictau3D->push_back( xictau_temp );

	 // calculate ctau error. Momentum error is negligible compared to the vertex errors, so don't worry about it
	 GlobalError XiVE = xiDecayVertex->error();
	 GlobalError PVE = GlobalError( refitVtx.error() );
	 VertexDistance3D theVertexDistance3D; 
         Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(XiVP, XiVE), VertexState(PVP, PVE) );
	 double myError = TheMeasurement.error();	 
	 
	 //  ctau is defined by the portion of the flight distance along the compoenent of the B momementum, so only
	 // consider the error of that component, too, which is accomplished by scaling by ((VB-VP)(dot)PB)/|VB-VP|*|PB|	 
         float scale = abs( (sep3D.dot(pXiV))/(sep3D.mag()*pXiV.mag()) ); 		 
         float xictauE_temp =  (myError*scale)/betagamma;
	 xictau3DE->push_back( xictauE_temp );

	 // next do ctau in the lazy way of XiVtx-priVtx
         xictauB->push_back( TheMeasurement.value()/betagamma );
	 xictauBE->push_back( TheMeasurement.error()/betagamma );

	 //calculate most probable value for ctau in 3D
	 
	 AlgebraicMatrix31 pB;
	 pB(0,0) = xiCand->currentState().globalMomentum().x();
	 pB(1,0) = xiCand->currentState().globalMomentum().y();
	 pB(2,0) = xiCand->currentState().globalMomentum().z();
	 
	 
	 AlgebraicMatrix13 pBT;
	 pBT(0,0) = xiCand->currentState().globalMomentum().x();
	 pBT(0,1) = xiCand->currentState().globalMomentum().y();
	 pBT(0,2) = xiCand->currentState().globalMomentum().z();
	 
	 AlgebraicMatrix31 PV;
	 PV(0,0) = refitVtx.position().x();
	 PV(1,0) = refitVtx.position().y();
	 PV(2,0) = refitVtx.position().z();
	 AlgebraicMatrix31 BV;
	 BV(0,0) = xiDecayVertex->position().x();
	 BV(1,0) = xiDecayVertex->position().y();
	 BV(2,0) = xiDecayVertex->position().z();
	 AlgebraicMatrix31 lxyz = BV-PV;
	 AlgebraicMatrix33 PVError(refitVtx.error());
	 AlgebraicMatrix33 BVError(xiDecayVertex->error().matrix_new());
	 AlgebraicMatrix33 lxyzError = PVError + BVError;
	 lxyzError.Invert();
	 
	 float xictauMPV_temp;
	 AlgebraicMatrix11 a = pBT * lxyzError * pB ;
	 AlgebraicMatrix11 b = pBT * lxyzError * lxyz;
	 double num(b(0,0));
	 double deno(a(0,0));
	 xictauMPV_temp = (num*ximass)/(deno);
	 xictauMPV->push_back(xictauMPV_temp);     

	 float xictauMPVE_temp = ximass*(sqrt(1/(deno)));
	 xictauMPVE->push_back( xictauMPVE_temp );

	 priVtxX = bestVtx->x();
	 priVtxY = bestVtx->y();
	 priVtxZ = bestVtx->z();
	 priVtxXE = bestVtx->xError();
	 priVtxYE = bestVtx->yError();
	 priVtxZE = bestVtx->zError();
	 priVtxCL = ChiSquaredProbability((double)(bestVtx->chi2()),(double)(bestVtx->ndof())); 
	 priNTrk = bestVtx->tracksSize();

	 priRefitVtxX->push_back(refitVtx.x());
	 priRefitVtxY->push_back(refitVtx.y());
	 priRefitVtxZ->push_back(refitVtx.z());
	 priRefitVtxXE->push_back(refitVtx.xError());
	 priRefitVtxYE->push_back(refitVtx.yError());
	 priRefitVtxZE->push_back(refitVtx.zError());
	 priRefitVtxCL->push_back(ChiSquaredProbability((double)(refitVtx.chi2()),(double)(refitVtx.ndof()))); 
	 priRefitNTrk->push_back(refitVtx.tracksSize());

	 bsX = beamSpot.x0();
	 bsY = beamSpot.y0();
	 bsZ = beamSpot.z0();
	 bsXE = beamSpot.x0Error();
	 bsYE = beamSpot.y0Error();
	 bsZE = beamSpot.z0Error();

	 nXI++;
	 
	 VeeParticles.clear();
	 xiFitParticles.clear();
	 omFitParticles.clear();
	 
       } 
     } 
   } // if V0Handle > 0 and muHandle > 1
   
     // if MC, check for generated Xi and Om
   if (doMC) {
     Handle<GenParticleCollection> genParticles;
     iEvent.getByLabel("genParticlesPlusSim", genParticles);

     for( size_t k = 0; k < genParticles->size(); k++ ) {
       const Candidate &Cand = (*genParticles)[ k ];
       float batPiGenEta = -99.0;
       float batPiGenPhi = -99.0;
       float lamPiGenEta = -99.0;
       float lamPiGenPhi = -99.0;
       float lamPGenEta = -99.0;
       float lamPGenPhi = -99.0;
       float xiGenP = -99.0;
       float xiGenEta = -99.0;
       float xiGenL = -99.0;
       float xiGenR = -99.0;
       float xiGenProdVX = -99.0;
       float xiGenProdVY = -99.0;
       float xiGenProdVZ = -99.0;
       float xiGenDecayVX = -99.0;
       float xiGenDecayVY = -99.0;
       float xiGenDecayVZ = -99.0;
       float xiGenMomL = -1.0;
       int xiGenMommyPDG = -99;
       if ( abs(Cand.pdgId())==3312 ) {
         if(debug) cout << "Found a Xi" << endl;
	 xiGenMommyPDG = Cand.mother(0)->pdgId();
         xiGenP = Cand.p();
         xiGenEta = Cand.eta();
         xiGenProdVX = Cand.vertex().x();
         xiGenProdVY = Cand.vertex().y();
         xiGenProdVZ = Cand.vertex().z();
	 if (Cand.numberOfMothers()>0) {
	   xiGenMomL = sqrt( (Cand.mother(0)->vertex().x()-xiGenProdVX) * (Cand.mother(0)->vertex().x()-xiGenProdVX) +
        	             (Cand.mother(0)->vertex().y()-xiGenProdVY) * (Cand.mother(0)->vertex().y()-xiGenProdVY) +
         	             (Cand.mother(0)->vertex().z()-xiGenProdVZ) * (Cand.mother(0)->vertex().z()-xiGenProdVZ) );
         }
	 for ( uint i = 0; i < Cand.numberOfDaughters(); i++){
           const Candidate * xiDau = Cand.daughter(i);
           if (debug) cout << "Xi daughter " << i << " has id = " << xiDau->pdgId() << endl;
           if ( abs(xiDau->pdgId())==211 ) {
	     if (debug) cout << "found pion from Xi" << endl;
	     batPiGenEta = xiDau->eta();
	     batPiGenPhi = xiDau->phi();
             xiGenDecayVX = xiDau->vertex().x();
             xiGenDecayVY = xiDau->vertex().y();
             xiGenDecayVZ = xiDau->vertex().z();
	     xiGenR = sqrt( (xiGenDecayVX-xiGenProdVX)*(xiGenDecayVX-xiGenProdVX) + 
	                    (xiGenDecayVY-xiGenProdVY)*(xiGenDecayVY-xiGenProdVY) ); 
	     xiGenL = sqrt( xiGenR*xiGenR + (xiGenDecayVZ-xiGenProdVZ)*(xiGenDecayVZ-xiGenProdVZ)	);
	   }
	   if ( abs(xiDau->pdgId())==3122 ) {
  	     if (debug) cout << "found lambda from Xi" << endl;
	     xiGenDecayVX = xiDau->vertex().x();
             xiGenDecayVY = xiDau->vertex().y();
             xiGenDecayVZ = xiDau->vertex().z();
	     xiGenR = sqrt( (xiGenDecayVX-xiGenProdVX)*(xiGenDecayVX-xiGenProdVX) + 
	                    (xiGenDecayVY-xiGenProdVY)*(xiGenDecayVY-xiGenProdVY) ); 
	     xiGenL = sqrt( xiGenR*xiGenR + (xiGenDecayVZ-xiGenProdVZ)*(xiGenDecayVZ-xiGenProdVZ)	);
	     for ( uint j = 0; j < xiDau->numberOfDaughters(); j++){
               const Candidate * xiGDau = xiDau->daughter(j);
               if (debug) cout << "Xi grand daughter " << j << " has id = " << xiGDau->pdgId() << endl;
	       if ( abs(xiGDau->pdgId())==211 ) {
	         if (debug) cout << "found pion from lambda from Xi" << endl;
	         lamPiGenEta = xiGDau->eta();
	         lamPiGenPhi = xiGDau->phi();
	       }
	       if ( abs(xiGDau->pdgId())==2212 ) {
	         if (debug) cout << "found proton from lambda from Xi" << endl;
	         lamPGenEta = xiGDau->eta();
	         lamPGenPhi = xiGDau->phi();
	       }
             }
	   }
         }
         genBatPiEta->push_back(batPiGenEta); 
         genBatPiPhi->push_back(batPiGenPhi); 
         genLamPiEta->push_back(lamPiGenEta); 
         genLamPiPhi->push_back(lamPiGenPhi); 
         genLamPEta->push_back(lamPGenEta); 
         genLamPPhi->push_back(lamPGenPhi);   
         genXiP->push_back(xiGenP);
         genXiEta->push_back(xiGenEta);
         genXiL->push_back(xiGenL);
         genXiR->push_back(xiGenR);
       
         genXiDecayVX->push_back(xiGenDecayVX);    
         genXiDecayVY->push_back(xiGenDecayVY);   
         genXiDecayVZ->push_back(xiGenDecayVZ);   
	 genXiProdVtxX->push_back(xiGenProdVX);
	 genXiProdVtxY->push_back(xiGenProdVY);
	 genXiProdVtxZ->push_back(xiGenProdVZ);
	 genXiMomL->push_back(xiGenMomL);
	 genXiQ->push_back(Cand.charge());
	 genXiMotherPDG->push_back(xiGenMommyPDG);
         nGenXi++; 
       }
       if ( abs(Cand.pdgId())==3334 ) {
         if (debug) cout << "Found an Om" << endl;
       }
     }
   }
   //fill the tree and clear the vectors
   if ( nXI > 0 || nGenXi>0 ) {
     //if(debug) cout <<"filling tree" << endl;
     tree_->Fill();
   }
   nXI = 0; run = 0; event = 0; trigHF = 0; HFpE = 0; HFmE = 0; trigTech40_41 = 0; trigTech36_39 = 0; trigTech34 = 0; trigScraping = 0; trigAlgo124 = 0; trigHLTminBias = 0; processType = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0; priNTrk = 0;
   bsX = 0;; bsY = 0;; bsZ = 0;; bsXE = 0;; bsYE = 0;; bsZE = 0;;
   priRefitVtxX->clear(); priRefitVtxY->clear(); priRefitVtxZ->clear(); priRefitVtxXE->clear(); priRefitVtxYE->clear(); priRefitVtxZE->clear(); priRefitVtxCL->clear(); priRefitNTrk->clear();
   xiMass->clear(); omMass->clear(); xiVtxCL->clear(); xiPx->clear(); xiPy->clear(); xiPz->clear(); 
   xiPxE->clear(); xiPyE->clear(); xiPzE->clear();
   xictauB->clear(); xictauBE->clear(); xictau3D->clear(); xictau3DE->clear(); xictauMPV->clear(); xictauMPVE->clear();
   xiDecayVtxX->clear(); xiDecayVtxY->clear(); xiDecayVtxZ->clear(); 
   xiDecayVtxXE->clear(); xiDecayVtxYE->clear(); xiDecayVtxZE->clear(); 
   VMass->clear(); VCandMass->clear(); VMassError->clear(); VKsMass->clear(); VKsMassError->clear();
   VVtxCL->clear(); VPx->clear(); VPy->clear(); VPz->clear(); VPalongXi->clear();
   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   batPiPx->clear(); batPiPy->clear(); batPiPz->clear(); 
   batPiEta->clear(); batPiPhi->clear(); batPiD0->clear(); batPiD0E->clear(); 
   batPiDz->clear(); batPiDzE->clear(); batPiPVweight->clear();
   batPiQ->clear(); batPiNValidHits->clear(); batPiTrkQual->clear(); nTracks->clear();
   VTrkPMass->clear(); VTrkPPx->clear(); VTrkPPy->clear(); VTrkPPz->clear(); VTrkPEta->clear(); VTrkPPhi->clear();
   VTrkPD0->clear(); VTrkPD0E->clear(); VTrkPDz->clear(); VTrkPDzE->clear(); VTrkPQ->clear(); VTrkPPVweight->clear();
   VTrkPNValidHits->clear();
   VTrkPiMass->clear(); VTrkPiPx->clear(); VTrkPiPy->clear(); VTrkPiPz->clear(); 
   VTrkPiEta->clear(); VTrkPiPhi->clear();
   VTrkPiNValidHits->clear();
   VTrkPiD0->clear(); VTrkPiD0E->clear(); VTrkPiDz->clear(); VTrkPiDzE->clear(); VTrkPiQ->clear(); VTrkPiPVweight->clear();
   VTransversePCAPrimary->clear(); VTransversePCAPrimaryError->clear(); VLongitudinalPCAPrimary->clear(); VLongitudinalPCAPrimaryError->clear();
   VTransversePCABeamSpot->clear(); VTransversePCABeamSpotError->clear();VCosThetaPAndveeVertexToPrimaryVector->clear();
   VCosThetaPAndveeVertexToXiVertexVector->clear();
   VRSig2D->clear(); VRSig3D->clear(); VFLSig2D->clear(); VFLSig3D->clear();
   XiTransversePCAPrimary->clear(); XiTransversePCAPrimaryError->clear(); XiLongitudinalPCAPrimary->clear(); XiLongitudinalPCAPrimaryError->clear();
   XiTransversePCABeamSpot->clear(); XiTransversePCABeamSpotError->clear();XiCosThetaPAndLineOfFlight->clear();
   Xi3DIpSig->clear(); XiFLsig3D->clear(); 
   V3dIpWrtPrimary->clear(); V3dIpWrtPrimarySig->clear();
   VTrkPi2DIp->clear();VTrkPi2DIpSig->clear();VTrkPi3DIp->clear();VTrkPi3DIpSig->clear();
   VTrkP2DIp->clear();VTrkP2DIpSig->clear();VTrkP3DIp->clear();VTrkP3DIpSig->clear();
   batPi2DIp->clear();batPi2DIpSig->clear();batPi3DIp->clear();batPi3DIpSig->clear();
   pionInV0->clear();

   nGenXi = 0; genBatPiEta->clear(); genBatPiPhi->clear(); genLamPiEta->clear(); genLamPiPhi->clear();
   genLamPEta->clear(); genLamPPhi->clear(); genXiP->clear(); genXiEta->clear();
   genXiL->clear(); genXiR->clear(); genXiDecayVX->clear(); genXiDecayVY->clear(); genXiDecayVZ->clear(); 
   genXiMomL->clear(); genXiQ->clear();
   genXiMotherPDG->clear();
   genXiProdVtxX->clear(); genXiProdVtxY->clear(); genXiProdVtxZ->clear();
 
}

// ------------ method called once each job just before starting event loop  ------------

void 
PiLambda::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","btosmumu ntuple");

  tree_->Branch("nXI",&nXI,"nXI/i");
  tree_->Branch("run",&run,"run/i");
  tree_->Branch("event",&event,"event/i");
  tree_->Branch("trigHF",&trigHF,"trigHF/I");
  tree_->Branch("HFpE",&HFpE,"HFpE/f");
  tree_->Branch("HFmE",&HFmE,"HFmE/f");
  tree_->Branch("trigTech40_41",&trigTech40_41,"trigTech40_41/I");
  tree_->Branch("trigTech36_39",&trigTech36_39,"trigTech36_39/I");
  tree_->Branch("trigTech34",&trigTech34,"trigTech34/I");
  tree_->Branch("trigScraping", &trigScraping, "trigScraping/I");
  tree_->Branch("trigAlgo124", &trigAlgo124, "trigAlgo124/I");
  tree_->Branch("trigHLTminBias", &trigHLTminBias, "trigHLTminBias/I");
  tree_->Branch("nLooseTracks", &nLooseTracks, "nLooseTracks/I");
  tree_->Branch("nHighPurTracks", &nHighPurTracks, "nHighPurTracks/I");
  tree_->Branch("nLooseTracks_pt_100", &nLooseTracks_pt_100, "nLooseTracks_pt_100/I");
  tree_->Branch("nHighPurTracks_pt_100", &nHighPurTracks_pt_100, "nHighPurTracks_pt_100/I");
  tree_->Branch("processType", &processType, "processType/I");
  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");
  tree_->Branch("priNTrk",&priNTrk, "priNTrk/i");
  tree_->Branch("priNTrk_weightGT_0_5", &priNTrk_weightGT_0_5, "priNTrk_weightGT_0_5/I");
  tree_->Branch("bsX",&bsX, "bsX/f");
  tree_->Branch("bsY",&bsY, "bsY/f");
  tree_->Branch("bsZ",&bsZ, "bsZ/f");
  tree_->Branch("bsXE",&bsXE, "bsXE/f");
  tree_->Branch("bsYE",&bsYE, "bsYE/f");
  tree_->Branch("bsZE",&bsZE, "bsZE/f");
  tree_->Branch("priRefitVtxX",&priRefitVtxX);
  tree_->Branch("priRefitVtxY",&priRefitVtxY);
  tree_->Branch("priRefitVtxZ",&priRefitVtxZ);
  tree_->Branch("priRefitVtxXE",&priRefitVtxXE);
  tree_->Branch("priRefitVtxYE",&priRefitVtxYE);
  tree_->Branch("priRefitVtxZE",&priRefitVtxZE);
  tree_->Branch("priRefitVtxCL",&priRefitVtxCL);
  tree_->Branch("priRefitNTrk",&priRefitNTrk);
  tree_->Branch("xiMass",&xiMass);
  tree_->Branch("omMass",&omMass);
  tree_->Branch("xiVtxCL",&xiVtxCL);
  tree_->Branch("xiPx",&xiPx);
  tree_->Branch("xiPy",&xiPy);
  tree_->Branch("xiPz",&xiPz);
  tree_->Branch("xiPxE",&xiPxE);
  tree_->Branch("xiPyE",&xiPyE);
  tree_->Branch("xiPzE",&xiPzE);
  tree_->Branch("xiDecayVtxX",&xiDecayVtxX);
  tree_->Branch("xiDecayVtxY",&xiDecayVtxY);
  tree_->Branch("xiDecayVtxZ",&xiDecayVtxZ);
  tree_->Branch("xiDecayVtxXE",&xiDecayVtxXE);
  tree_->Branch("xiDecayVtxYE",&xiDecayVtxYE);
  tree_->Branch("xiDecayVtxZE",&xiDecayVtxZE);
  tree_->Branch("xictauB",&xictauB);
  tree_->Branch("xictauBE",&xictauBE); 
  tree_->Branch("xictau3D",&xictau3D);
  tree_->Branch("xictau3DE",&xictau3DE);  
  tree_->Branch("xictauMPV",&xictauMPV);
  tree_->Branch("xictauMPVE",&xictauMPVE);
  tree_->Branch("VMass",&VMass);
  tree_->Branch("VCandMass", &VCandMass);
  tree_->Branch("VMassError", &VMassError);
  tree_->Branch("VKsMass", &VKsMass);
  tree_->Branch("VKsMassError", &VKsMassError);
  tree_->Branch("VVtxCL",&VVtxCL);
  tree_->Branch("VPx",&VPx);
  tree_->Branch("VPy",&VPy);
  tree_->Branch("VPz",&VPz);
  tree_->Branch("VPalongXi",&VPalongXi);
  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("batPiPx",&batPiPx);
  tree_->Branch("batPiPy",&batPiPy);
  tree_->Branch("batPiPz",&batPiPz);
  tree_->Branch("batPiEta",&batPiEta);
  tree_->Branch("batPiPhi",&batPiPhi);
  tree_->Branch("batPiD0",&batPiD0);
  tree_->Branch("batPiD0E",&batPiD0E);
  tree_->Branch("batPiDz",&batPiDz);
  tree_->Branch("batPiDzE",&batPiDzE);
  tree_->Branch("batPiPVweight",&batPiPVweight);
  tree_->Branch("batPiQ",&batPiQ);
  tree_->Branch("batPiNValidHits",&batPiNValidHits);
  tree_->Branch("batPiTrkQual",&batPiTrkQual);
  tree_->Branch("nTracks",&nTracks);
  tree_->Branch("VTrkPTrk1Mass",&VTrkPMass);
  tree_->Branch("VTrkPPx",&VTrkPPx);
  tree_->Branch("VTrkPPy",&VTrkPPy);
  tree_->Branch("VTrkPPz",&VTrkPPz);
  tree_->Branch("VTrkPEta",&VTrkPEta);
  tree_->Branch("VTrkPPhi",&VTrkPPhi);
  tree_->Branch("VTrkPD0",&VTrkPD0);
  tree_->Branch("VTrkPD0E",&VTrkPD0E);
  tree_->Branch("VTrkPDz",&VTrkPDz);
  tree_->Branch("VTrkPDzE",&VTrkPDzE);
  tree_->Branch("VTrkPQ",&VTrkPQ);
  tree_->Branch("VTrkPPVweight",&VTrkPPVweight);
  tree_->Branch("VTrkPNValidHits",&VTrkPNValidHits);
  tree_->Branch("VTrkPiMass",&VTrkPiMass);
  tree_->Branch("VTrkPiPx",&VTrkPiPx);
  tree_->Branch("VTrkPiPy",&VTrkPiPy);
  tree_->Branch("VTrkPiPz",&VTrkPiPz);
  tree_->Branch("VTrkPiEta",&VTrkPiEta);
  tree_->Branch("VTrkPiPhi",&VTrkPiPhi);
  tree_->Branch("VTrkPiD0",&VTrkPiD0);
  tree_->Branch("VTrkPiD0E",&VTrkPiD0E);
  tree_->Branch("VTrkPiDz",&VTrkPiDz);
  tree_->Branch("VTrkPiDzE",&VTrkPiDzE);
  tree_->Branch("VTrkPiQ",&VTrkPiQ);
  tree_->Branch("VTrkPiPVweight",&VTrkPiPVweight);
  tree_->Branch("VTrkPiNValidHits",&VTrkPiNValidHits);
  tree_->Branch("VTransversePCAPrimary",&VTransversePCAPrimary);
  tree_->Branch("VTransversePCAPrimaryError",&VTransversePCAPrimaryError);
  tree_->Branch("VLongitudinalPCAPrimary",&VLongitudinalPCAPrimary);
  tree_->Branch("VLongitudinalPCAPrimaryError",&VLongitudinalPCAPrimaryError);
  tree_->Branch("VTransversePCABeamSpot",&VTransversePCABeamSpot);
  tree_->Branch("VTransversePCABeamSpotError",&VTransversePCABeamSpotError);
  tree_->Branch("VRSig2D", &VRSig2D);
  tree_->Branch("VRSig3D", &VRSig3D);  
  tree_->Branch("VFLSig2D", &VFLSig2D);
  tree_->Branch("VFLSig3D", &VFLSig3D);  
  tree_->Branch("VCosThetaPAndveeVertexToPrimaryVector",&VCosThetaPAndveeVertexToPrimaryVector);
  tree_->Branch("VCosThetaPAndveeVertexToXiVertexVector",&VCosThetaPAndveeVertexToXiVertexVector); 
  tree_->Branch("XiTransversePCAPrimary",&XiTransversePCAPrimary);
  tree_->Branch("XiTransversePCAPrimaryError",&XiTransversePCAPrimaryError);
  tree_->Branch("XiLongitudinalPCAPrimary",&XiLongitudinalPCAPrimary);
  tree_->Branch("XiLongitudinalPCAPrimaryError",&XiLongitudinalPCAPrimaryError);
  tree_->Branch("XiTransversePCABeamSpot",&XiTransversePCABeamSpot);
  tree_->Branch("XiTransversePCABeamSpotError",&XiTransversePCABeamSpotError);
  tree_->Branch("XiCosThetaPAndLineOfFlight",&XiCosThetaPAndLineOfFlight);
  tree_->Branch("Xi3DIpSig",&Xi3DIpSig);
  tree_->Branch("XiFLsig3D",&XiFLsig3D);
  tree_->Branch("V3dIpWrtPrimary",&V3dIpWrtPrimary);
  tree_->Branch("V3dIpWrtPrimarySig",&V3dIpWrtPrimarySig);
  tree_->Branch("VTrkPi2DIp",&VTrkPi2DIp);
  tree_->Branch("VTrkPi2DIpSig",&VTrkPi2DIpSig);
  tree_->Branch("VTrkPi3DIp",&VTrkPi3DIp);
  tree_->Branch("VTrkPi3DIpSig",&VTrkPi3DIpSig);
  tree_->Branch("VTrkP2DIp",&VTrkP2DIp);
  tree_->Branch("VTrkP2DIpSig",&VTrkP2DIpSig);
  tree_->Branch("VTrkP3DIp",&VTrkP3DIp);
  tree_->Branch("VTrkP3DIpSig",&VTrkP3DIpSig);
  tree_->Branch("batPi2DIp",&batPi2DIp);
  tree_->Branch("batPi2DIpSig",&batPi2DIpSig);
  tree_->Branch("batPi3DIp",&batPi3DIp);
  tree_->Branch("batPi3DIpSig",&batPi3DIpSig);
  tree_->Branch("pionInV0",&pionInV0);
  
  tree_->Branch("nGenXi",&nGenXi,"nGenXi/i");
  tree_->Branch("genBatPiEta",&genBatPiEta);
  tree_->Branch("genBatPiPhi",&genBatPiPhi);
  tree_->Branch("genLamPiEta",&genLamPiEta);
  tree_->Branch("genLamPiPhi",&genLamPiPhi);
  tree_->Branch("genLamPEta",&genLamPEta);
  tree_->Branch("genLamPPhi",&genLamPPhi);

  tree_->Branch("genXiP",&genXiP);
  tree_->Branch("genXiEta",&genXiEta);
  tree_->Branch("genXiL",&genXiL);
  tree_->Branch("genXiR",&genXiR);
  tree_->Branch("genXiDecayVX",&genXiDecayVX);
  tree_->Branch("genXiDecayVY",&genXiDecayVY);
  tree_->Branch("genXiDecayVZ",&genXiDecayVZ);
  tree_->Branch("genXiMomL",&genXiMomL);
  tree_->Branch("genXiQ",&genXiQ);

  tree_->Branch("genXiMotherPDG", &genXiMotherPDG);
  tree_->Branch("genXiProdVtxX", &genXiProdVtxX);
  tree_->Branch("genXiProdVtxY", &genXiProdVtxY);
  tree_->Branch("genXiProdVtxZ", &genXiProdVtxZ);

}


// ------------ method called once each job just after ending the event loop  ------------
void 
PiLambda::endJob() {
  cout << "Ending job with nEvents = " << nEvents << endl;
  //tree_->GetDirectory()->cd();
  //tree_->Write();
}

// checks if global track is part of a good V0
bool PiLambda::trackInV0(const reco::TrackRef theTrk, std::vector<reco::VertexCompositeCandidate> theV0Collection){
  bool trkIsMatched = false;
  ParticleMass lambdaMass = 1.115683;
  //loop over v0 collection
  for(std::vector<reco::VertexCompositeCandidate>::const_iterator theV0It = theV0Collection.begin(); 
      theV0It != theV0Collection.end(); ++theV0It){

    if( (theV0It->mass() > (lambdaMass + 0.008)) || (theV0It->mass() < (lambdaMass - 0.008)) )continue;  
    
    // get daughter tracks... doing same as above
    std::vector<reco::RecoChargedCandidate> v0daughters;						 

    //check momentum and add pion first, proton second  				 
    if (theV0It->daughter(0)->momentum().mag2() < theV0It->daughter(1)->momentum().mag2()){	 
      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
    			       (theV0It->daughter(0))) );
      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
    			       (theV0It->daughter(1))) );
    } else {										 
      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
    			       (theV0It->daughter(1))) );
      v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
    			       (theV0It->daughter(0))) );
    }											 
    for(unsigned int j = 0; j < v0daughters.size(); ++j) {
      const reco::TrackRef theDaughterTrk = v0daughters[j].track();				 
      if(theDaughterTrk->charge() == theTrk->charge() && theDaughterTrk->momentum() == theTrk->momentum()){
        //cout << "found a track match in track in v0 test" << endl;
	trkIsMatched = true;
	break;
      }
    }
    if(trkIsMatched)break;											 
  }
  //check if track is compatible with any of the doughter tracks
  return trkIsMatched;
}

bool PiLambda::aboveTriggerThreshold(edm::Handle<CaloTowerCollection> towers)
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

bool PiLambda::scrapingTest( const edm::Event& iEvent, 
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
DEFINE_FWK_MODULE(PiLambda);

