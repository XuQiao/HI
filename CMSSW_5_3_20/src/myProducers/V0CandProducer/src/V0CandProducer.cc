// -*- C++ -*-
//
// Package:    V0CandProducer
// Class:      V0CandProducer
// 
/**\class V0CandProducer V0CandProducer.cc myProducers/V0CandProducer/src/V0CandProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Mon Sep 28 16:32:58 MDT 2009
// $Id: V0CandProducer.cc,v 1.11 2011/05/23 16:24:27 drell Exp $
//
//


// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "DataFormats/V0Candidate/interface/V0Candidate.h"

//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"   //by Qiao
//#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" 	//by Qiao

#include "myProducers/V0CandProducer/interface/V0Candidate.h"   //by Qiao

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"


class V0CandProducer : public edm::EDProducer {
public:
  explicit V0CandProducer(const edm::ParameterSet&);
  ~V0CandProducer();
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  InvariantMassFromVertex massCalculator;
      
  edm::InputTag k0sCollection;
  edm::InputTag lamCollection;

  bool embedTracks;
};

const double piMassConst = 0.13957018;
const double protonMassConst = 0.93827203;
const double lamMassConst = 1.115683;

V0CandProducer::V0CandProducer(const edm::ParameterSet& iConfig):
  k0sCollection( iConfig.getParameter<edm::InputTag>("kShortCollection")),
  lamCollection( iConfig.getParameter<edm::InputTag>("lambdaCollection")),
  embedTracks( iConfig.getUntrackedParameter<bool>("embedTracks", true))
{

  produces< std::vector<reco::V0Candidate> >("allVees");
  produces< pat::CompositeCandidateCollection >("allVees");
  
}


V0CandProducer::~V0CandProducer() {
}

void V0CandProducer::produce(edm::Event& iEvent, 
			     const edm::EventSetup& iSetup) {
   using namespace edm;
   using namespace std;

   ESHandle<TransientTrackBuilder> ttBuilder;

   Handle<reco::VertexCompositeCandidateCollection> theKshorts;
   Handle<reco::VertexCompositeCandidateCollection> theLambdas;

   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttBuilder);

   iEvent.getByLabel(k0sCollection, theKshorts);
   iEvent.getByLabel(lamCollection, theLambdas);

   std::vector<reco::V0Candidate> theVeeCands;
   std::auto_ptr<pat::CompositeCandidateCollection> 
     patOut(new pat::CompositeCandidateCollection);

 /*  reco::TrackBase::TrackQuality loose_ = reco::TrackBase::qualityByName("loose");
   reco::TrackBase::TrackQuality tight_ = reco::TrackBase::qualityByName("tight");
   reco::TrackBase::TrackQuality highPurity_ = reco::TrackBase::qualityByName("highPurity");
   reco::TrackBase::TrackQuality confirmed_ = reco::TrackBase::qualityByName("confirmed");
   reco::TrackBase::TrackQuality goodIterative_ = reco::TrackBase::qualityByName("goodIterative");
   reco::TrackBase::TrackQuality qualSize = reco::TrackBase::qualityByName("qualitySize");
*/
//    cout << "loose:" << loose_
// 	<< " tight:" << tight_
// 	<< " highPurity:" << highPurity_
// 	<< " confirmed:" << confirmed_
// 	<< " goodIterative:" << goodIterative_
// 	<< " qualitySize:" << qualSize
// 	<< endl;

   // Loop over the K0s collection, find which tracks have which charge, and push back
   //  V0Candidate objects
   for( reco::VertexCompositeCandidateCollection::const_iterator iKS = theKshorts->begin();
	iKS != theKshorts->end();
	iKS++ ) {
     reco::VertexCompositeCandidate theKS(*iKS);
     reco::VertexCompositeCandidate* iKS_ = &theKS;
     reco::CompositeCandidate* iKSptr 
       = dynamic_cast<reco::CompositeCandidate *> (iKS_);
     pat::CompositeCandidate ksPatCand(*iKSptr);
     ksPatCand.clearDaughters();
 //     if (embedTracks) {
//        ksPatCand.daughter(0)->embedTrack();
//        ksPatCand.daughter(1)->embedTrack();
//      }

     int numOfPosDaughter = -1;
     int numOfNegDaughter = -1;
     if( iKS->daughter(0)->charge() > 0. && iKS->daughter(1)->charge() < 0. ) {
       numOfPosDaughter = 0;
       numOfNegDaughter = 1;
     }
     else if( iKS->daughter(1)->charge() > 0. && iKS->daughter(0)->charge() < 0. ){
       numOfPosDaughter = 1;
       numOfNegDaughter = 0;
     }
     if( numOfPosDaughter >= 0 && numOfNegDaughter >= 0 ) {
       // PAT
       pat::GenericParticle patPosDau( *(iKS->daughter(numOfPosDaughter)) );
       pat::GenericParticle patNegDau( *(iKS->daughter(numOfNegDaughter)) );
       //const reco::RecoChargedCandidate* posDauCand = dynamic_cast<const reco::RecoChargedCandidate*>(iKS->daughter(numOfPosDaughter));
       //const reco::RecoChargedCandidate* negDauCand = dynamic_cast<const reco::RecoChargedCandidate*>(iKS->daughter(numOfNegDaughter));
       //reco::RecoChargedCandidate posDauCand_ = *posDauCand;
       //reco::CandidatePtr test(posDauCand_);

       //cout << "patPosDau_eta: " << patPosDau.eta() << endl;
       //cout << "patPosDau_track_chi2: " << patPosDau.track()->chi2() << endl;
//        if (embedTracks) {
// 	 patPosDau.embedTrack();
// 	 patNegDau.embedTrack();
//        }
       //ksPatCand.addDaughter( patPosDau, "posTk" );
       //ksPatCand.addDaughter( patNegDau, "negTk" );

       //pat::GenericParticle patPosDau();
       //pat::GenericParticle patNegDau();

       // Non-PAT
       reco::RecoChargedCandidate posCand1 
	 = *(dynamic_cast<const reco::RecoChargedCandidate*>(iKS->daughter(numOfPosDaughter)));
       reco::RecoChargedCandidate negCand1
	 = *(dynamic_cast<const reco::RecoChargedCandidate*>(iKS->daughter(numOfNegDaughter)));
       patPosDau.setTrack( posCand1.track(), embedTracks );
       patNegDau.setTrack( negCand1.track(), embedTracks );

       //int posQual = 0;
       //int negQual = 0;
       int posAlgo = -1;
       int negAlgo = -1;
       //int posLoose = 0, posTight = 0, posHighPur = 0, posConfirmed = 0,
       //	 posGoodIter = 0, posQualSize = 0;
       posAlgo = posCand1.track()->algo();
       negAlgo = negCand1.track()->algo();
//        if (posCand1.track()->quality(loose_)) posLoose = 1;
//        if (negCand1.track()->quality(loose_)) negQual = 1;
//        if (posCand1.track()->quality(tight_)) posTight = 1;
//        if (posCand1.track()->quality(highPurity_)) posHighPur = 1;
//        if (posCand1.track()->quality(confirmed_)) posConfirmed = 1;
//        if (negCand1.track()->quality(confirmed_)) negQual = 2;
//        if (posCand1.track()->quality(goodIterative_)) posGoodIter = 1;
//        if (negCand1.track()->quality(goodIterative_)) negQual = 3;
//        if (posCand1.track()->quality(qualSize)) posQualSize = 1;
//        if (negCand1.track()->quality(qualSize)) negQual = 4;

//        posLoose = posCand1.track()->quality(loose_);
//        posTight = posCand1.track()->quality(tight_);
//        posHighPur = posCand1.track()->quality(highPurity_);
//        posConfirmed = posCand1.track()->quality(confirmed_);
//        posGoodIter = posCand1.track()->quality(goodIterative_);
//        posQualSize = posCand1.track()->quality(qualSize);

//        cout << k0sCollection.label() << ":"<< k0sCollection.instance() << endl;
//        cout << "mask:" << posCand1.track()->qualityMask()
// 	    << " loose:" << posLoose
// 	    << " tight:" << posTight
// 	    << " highPurity:" << posHighPur
// 	    << " confirmed:" << posConfirmed
// 	    << " goodIterative:" << posGoodIter
// 	    << " qualitySize(" << qualSize << "):" << posQualSize << endl;
//        cout << "Mask should have been "
// 	    << posLoose*1 + posTight*2 + posHighPur*4 + posConfirmed*8 + posGoodIter*16 + posQualSize*32
// 	    << endl;

       ksPatCand.addUserInt("posQual", posCand1.track()->qualityMask());//posQual);
       ksPatCand.addUserInt("negQual", negCand1.track()->qualityMask());//negQual);
       ksPatCand.addUserInt("posAlgo", posAlgo);
       ksPatCand.addUserInt("negAlgo", negAlgo);

       ksPatCand.addDaughter( patPosDau );
       ksPatCand.addDaughter( patNegDau );

       //cout << "Tk Chi2: " 
       //   << (dynamic_cast<const reco::RecoCandidate*>(ksPatCand.daughter(0)))->track()->chi2() << endl;
       reco::V0Candidate temp( *iKS, 
			       *(posCand1.track()), 
			       *(negCand1.track()) );

       const reco::HitPattern& posHitPatt = posCand1.track()->hitPattern();
       const reco::HitPattern& negHitPatt = negCand1.track()->hitPattern();
       reco::TransientTrack posDaughterTT 
	 = ttBuilder->build( temp.getPosTrack() );
       reco::TransientTrack negDaughterTT 
	 = ttBuilder->build( temp.getNegTrack() );
       GlobalPoint vtxpos( temp.getVtxCC().vx(), 
			   temp.getVtxCC().vy(), 
			   temp.getVtxCC().vz() );

       KalmanVertexFitter theKsFitter(true);
       vector<reco::TransientTrack> ksTracks;
       ksTracks.push_back( posDaughterTT );
       ksTracks.push_back( negDaughterTT );
       CachingVertex<5> ksVtxForInvMass = theKsFitter.vertex( ksTracks );
       vector<double> ksDauMasses;
       ksDauMasses.push_back( piMassConst );
       ksDauMasses.push_back( piMassConst );

       vector<double> ksOtherDauMasses;
       ksOtherDauMasses.push_back( protonMassConst );
       ksOtherDauMasses.push_back( piMassConst );

       vector<double> ksOtherDauMasses2;
       ksOtherDauMasses2.push_back( piMassConst );
       ksOtherDauMasses2.push_back( protonMassConst );

       Measurement1D ksMassWErr 
	 = massCalculator.invariantMass( ksVtxForInvMass, ksDauMasses );
       Measurement1D ksOtherMassWErr1
	 = massCalculator.invariantMass( ksVtxForInvMass, ksOtherDauMasses );
       Measurement1D ksOtherMassWErr2
	 = massCalculator.invariantMass( ksVtxForInvMass, ksOtherDauMasses2 );

       /*cout << "k0: " << "REAL: " << iKS->mass() << "; "
	    << ksMassWErr.value() << ", " << ksMassWErr.error() << endl;
       cout << ksOtherMassWErr1.value() << ", " << ksOtherMassWErr1.error()
	    << endl;
       cout << ksOtherMassWErr2.value() << ", " << ksOtherMassWErr2.error()
	    << endl;
	    cout << "-----------" << endl;*/


       ksPatCand.addUserFloat( "massFromFit",
			       ksMassWErr.value() );
       ksPatCand.addUserFloat( "massErrorFromFit",
			       ksMassWErr.error() );
       if( (ksOtherMassWErr1.value() - lamMassConst)
	   < (ksOtherMassWErr2.value() - lamMassConst) ) {
	 ksPatCand.addUserFloat( "otherMassFromFit",
				 ksOtherMassWErr1.value() );
	 ksPatCand.addUserFloat( "otherMassErrorFromFit",
				 ksOtherMassWErr1.error() );
       }
       else {
	 ksPatCand.addUserFloat( "otherMassFromFit",
				 ksOtherMassWErr2.value() );
	 ksPatCand.addUserFloat( "otherMassErrorFromFit",
				 ksOtherMassWErr2.error() );
       }


       temp.setPosNPixelHits( posHitPatt.numberOfValidPixelHits() );
       temp.setNegNPixelHits( negHitPatt.numberOfValidPixelHits() );
       temp.setPosNStripHits( posHitPatt.numberOfValidStripHits() );
       temp.setNegNStripHits( negHitPatt.numberOfValidStripHits() );

       ksPatCand.addUserInt( "posTkNPixelHits",
			     posHitPatt.numberOfValidPixelHits() );
       ksPatCand.addUserInt( "negTkNPixelHits",
			     negHitPatt.numberOfValidPixelHits() );
       ksPatCand.addUserInt( "posTkNStripHits", 
			     posHitPatt.numberOfValidStripHits() );
       ksPatCand.addUserInt( "negTkNStripHits",
			     negHitPatt.numberOfValidStripHits() );

       TrajectoryStateClosestToPoint posDaughterTSCP 
	 = posDaughterTT.trajectoryStateClosestToPoint(vtxpos);
       TrajectoryStateClosestToPoint negDaughterTSCP
	 = negDaughterTT.trajectoryStateClosestToPoint(vtxpos);
       temp.setPosMomentum( posDaughterTSCP.momentum() );
       temp.setNegMomentum( negDaughterTSCP.momentum() );

       ksPatCand.addUserFloat( "posTkPX",
			       posDaughterTSCP.momentum().x() );
       ksPatCand.addUserFloat( "posTkPY",
			       posDaughterTSCP.momentum().y() );
       ksPatCand.addUserFloat( "posTkPZ",
			       posDaughterTSCP.momentum().z() );
       ksPatCand.addUserFloat( "negTkPX",
			       negDaughterTSCP.momentum().x() );
       ksPatCand.addUserFloat( "negTkPY",
			       negDaughterTSCP.momentum().y() );
       ksPatCand.addUserFloat( "negTkPZ",
			       negDaughterTSCP.momentum().z() );

       ksPatCand.addUserFloat( "vtxCov00",
			       iKS->vertexCovariance(0,0) );
       ksPatCand.addUserFloat( "vtxCov01",
			       iKS->vertexCovariance(0,1) );
       ksPatCand.addUserFloat( "vtxCov11",
			       iKS->vertexCovariance(1,1) );
       ksPatCand.addUserFloat( "vtxCov02",
			       iKS->vertexCovariance(0,2) );
       ksPatCand.addUserFloat( "vtxCov12",
			       iKS->vertexCovariance(1,2) );
       ksPatCand.addUserFloat( "vtxCov22",
			       iKS->vertexCovariance(2,2) );
       ksPatCand.addUserFloat( "vtxChi2",
			       iKS->vertexChi2() );
       ksPatCand.addUserFloat( "vtxNormChi2",
			       iKS->vertexNormalizedChi2() );

       theVeeCands.push_back( temp );
       patOut->push_back( ksPatCand );
     }
   }

   // Rinse and repeat for Lambda collection
   for( reco::VertexCompositeCandidateCollection::const_iterator iLM = theLambdas->begin();
	iLM != theLambdas->end();
	iLM++ ) {
     int lamPDG = iLM->pdgId();
     // PAT stuff
     reco::VertexCompositeCandidate theLM(*iLM);
     reco::VertexCompositeCandidate* iLM_ = &theLM;
     reco::CompositeCandidate* iLMptr 
       = dynamic_cast<reco::CompositeCandidate *> (iLM_);
     pat::CompositeCandidate lamPatCand(*iLMptr);
     lamPatCand.clearDaughters();
//      if (embedTracks) {
//        lamPatCand.daughter(0)->embedTrack();
//        lamPatCand.daughter(1)->embedTrack();
//      }

     // Non-PAT
     int numOfPosDaughter = -1;
     int numOfNegDaughter = -1;
     if( iLM->daughter(0)->charge() > 0. && iLM->daughter(1)->charge() < 0. ) {
       numOfPosDaughter = 0;
       numOfNegDaughter = 1;
     }
     else if( iLM->daughter(1)->charge() > 0. && iLM->daughter(0)->charge() < 0. ) {
       numOfPosDaughter = 1;
       numOfNegDaughter = 0;
     }
     if( numOfPosDaughter >= 0 && numOfNegDaughter >= 0 ) {
       // PAT
       pat::GenericParticle patPosDau( *(iLM->daughter(numOfPosDaughter)) );
       pat::GenericParticle patNegDau( *(iLM->daughter(numOfNegDaughter)) );
//        if (embedTracks) {
// 	 patPosDau.embedTrack();
// 	 patNegDau.embedTrack();
//        }
       //lamPatCand.addDaughter( patPosDau, "posTk" );
       //lamPatCand.addDaughter( patNegDau, "negTk" );
       //pat::GenericParticle patPosDau();
       //pat::GenericParticle patNegDau();

       // Non-PAT
       reco::RecoChargedCandidate posCand2 
	 = *(dynamic_cast<const reco::RecoChargedCandidate*>(iLM->daughter(numOfPosDaughter)));
       reco::RecoChargedCandidate negCand2
	 = *(dynamic_cast<const reco::RecoChargedCandidate*>(iLM->daughter(numOfNegDaughter)));

       patPosDau.setTrack( posCand2.track(), embedTracks );
       patNegDau.setTrack( negCand2.track(), embedTracks );

       //int posQual = 0;
       //int negQual = 0;
       int posAlgo = posCand2.track()->algo();
       int negAlgo = negCand2.track()->algo();
//        if (posCand2.track()->quality(loose_)) posQual = 1;
//        if (negCand2.track()->quality(loose_)) negQual = 1;
//        if (posCand2.track()->quality(confirmed_)) posQual = 2;
//        if (negCand2.track()->quality(confirmed_)) negQual = 2;
//        if (posCand2.track()->quality(goodIterative_)) posQual = 3;
//        if (negCand2.track()->quality(goodIterative_)) negQual = 3;
//        if (posCand2.track()->quality(qualSize)) posQual = 4;
//        if (negCand2.track()->quality(qualSize)) negQual = 4;
       lamPatCand.addUserInt("posQual", posCand2.track()->qualityMask());//posQual);
       lamPatCand.addUserInt("negQual", negCand2.track()->qualityMask());//negQual);
       lamPatCand.addUserInt("posAlgo", posAlgo);
       lamPatCand.addUserInt("negAlgo", negAlgo);

       lamPatCand.addDaughter( patPosDau );
       lamPatCand.addDaughter( patNegDau );

       reco::V0Candidate temp( *iLM,
			       *(posCand2.track()),
			       *(negCand2.track()) );
       const reco::HitPattern& posHitPatt = posCand2.track()->hitPattern();
       const reco::HitPattern& negHitPatt = negCand2.track()->hitPattern();
       reco::TransientTrack posDaughterTT 
	 = ttBuilder->build( temp.getPosTrack() );
       reco::TransientTrack negDaughterTT 
	 = ttBuilder->build( temp.getNegTrack() );
       GlobalPoint vtxpos( temp.getVtxCC().vx(), 
			   temp.getVtxCC().vy(), 
			   temp.getVtxCC().vz() );

       KalmanVertexFitter theLamFitter(true);
       vector<reco::TransientTrack> lamTracks;
       lamTracks.push_back( posDaughterTT );
       lamTracks.push_back( negDaughterTT );
       CachingVertex<5> lamVtxForInvMass = theLamFitter.vertex( lamTracks );
       vector<double> lamDauMasses;
       if( lamPDG > 0 ) {
	 lamDauMasses.push_back( protonMassConst );
	 lamDauMasses.push_back( piMassConst );
       }
       else {
	 lamDauMasses.push_back( piMassConst );
	 lamDauMasses.push_back( protonMassConst );
       }
       vector<double> lamOtherDauMasses;
       for( int ii = 0; ii < 2; ii++ ) {
	 lamOtherDauMasses.push_back( piMassConst );
       }

       Measurement1D lamMassWErr
	 = massCalculator.invariantMass( lamVtxForInvMass, lamDauMasses );
       Measurement1D lamOtherMassWErr
	 = massCalculator.invariantMass( lamVtxForInvMass, lamOtherDauMasses );

       //cout << "Cand mass: " << iLM->mass() << ", Fit mass: "
       //   << lamMassWErr.value() << endl;

       lamPatCand.addUserFloat( "massFromFit",
				lamMassWErr.value() );
       lamPatCand.addUserFloat( "massErrorFromFit",
				lamMassWErr.error() );
       lamPatCand.addUserFloat( "otherMassFromFit",
				lamOtherMassWErr.value() );
       lamPatCand.addUserFloat( "otherMassErrorFromFit",
				lamOtherMassWErr.error() );


       temp.setPosNPixelHits( posHitPatt.numberOfValidPixelHits() );
       temp.setNegNPixelHits( negHitPatt.numberOfValidPixelHits() );
       temp.setPosNStripHits( posHitPatt.numberOfValidStripHits() );
       temp.setNegNStripHits( negHitPatt.numberOfValidStripHits() );

       lamPatCand.addUserInt( "posTkNPixelHits",
			      posHitPatt.numberOfValidPixelHits() );
       lamPatCand.addUserInt( "negTkNPixelHits",
			      negHitPatt.numberOfValidPixelHits() );
       lamPatCand.addUserInt( "posTkNStripHits", 
			      posHitPatt.numberOfValidStripHits() );
       lamPatCand.addUserInt( "negTkNStripHits",
			      negHitPatt.numberOfValidStripHits() );

       TrajectoryStateClosestToPoint posDaughterTSCP 
	 = posDaughterTT.trajectoryStateClosestToPoint(vtxpos);
       TrajectoryStateClosestToPoint negDaughterTSCP
	 = negDaughterTT.trajectoryStateClosestToPoint(vtxpos);
       temp.setPosMomentum( posDaughterTSCP.momentum() );
       temp.setNegMomentum( negDaughterTSCP.momentum() );

       lamPatCand.addUserFloat( "posTkPX",
				posDaughterTSCP.momentum().x() );
       lamPatCand.addUserFloat( "posTkPY",
				posDaughterTSCP.momentum().y() );
       lamPatCand.addUserFloat( "posTkPZ",
				posDaughterTSCP.momentum().z() );
       lamPatCand.addUserFloat( "negTkPX",
				negDaughterTSCP.momentum().x() );
       lamPatCand.addUserFloat( "negTkPY",
				negDaughterTSCP.momentum().y() );
       lamPatCand.addUserFloat( "negTkPZ",
				negDaughterTSCP.momentum().z() );

       lamPatCand.addUserFloat( "vtxCov00",
			       iLM->vertexCovariance(0,0) );
       lamPatCand.addUserFloat( "vtxCov01",
			       iLM->vertexCovariance(0,1) );
       lamPatCand.addUserFloat( "vtxCov11",
			       iLM->vertexCovariance(1,1) );
       lamPatCand.addUserFloat( "vtxCov02",
			       iLM->vertexCovariance(0,2) );
       lamPatCand.addUserFloat( "vtxCov12",
			       iLM->vertexCovariance(1,2) );
       lamPatCand.addUserFloat( "vtxCov22",
			       iLM->vertexCovariance(2,2) );
       lamPatCand.addUserFloat( "vtxChi2",
				iLM->vertexChi2() );
       lamPatCand.addUserFloat( "vtxNormChi2",
				iLM->vertexNormalizedChi2() );

       theVeeCands.push_back( temp );
       patOut->push_back( lamPatCand );
     }
   }

   // Allocate vector and put it into the event
   std::auto_ptr< std::vector<reco::V0Candidate> > 
     veesOut( new std::vector<reco::V0Candidate> );
   veesOut->reserve( theVeeCands.size() );
   std::copy( theVeeCands.begin(),
	      theVeeCands.end(),
	      std::back_inserter( *veesOut ) );
   iEvent.put( veesOut, std::string("allVees") );

   iEvent.put( patOut, std::string("allVees") );
 
}

void V0CandProducer::beginJob() {
  // if (embedTracks) std::cout << "Embedding tracks in cand objects." << endl;
 if (embedTracks) std::cout << "Embedding tracks in cand objects." << std:: endl;	//by Qiao
}

void V0CandProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0CandProducer);
