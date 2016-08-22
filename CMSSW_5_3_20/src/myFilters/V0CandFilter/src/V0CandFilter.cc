// -*- C++ -*-
//
// Package:    V0CandFilter
// Class:      V0CandFilter
// 
/**\class V0CandFilter V0CandFilter.cc myFilters/V0CandFilter/src/V0CandFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Sat Dec  5 09:48:43 MST 2009
// $Id: V0CandFilter.cc,v 1.3 2009/12/07 02:20:35 drell Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class V0CandFilter : public edm::EDFilter {
   public:
      explicit V0CandFilter(const edm::ParameterSet&);
      ~V0CandFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

  edm::InputTag k0sCollection;
  edm::InputTag lamCollection;
  edm::InputTag trackCollection;
};

//
// constants, enums and typedefs
//

const double ksMass = 0.49765;
const double lamMass = 1.11568;

//
// static data member definitions
//

//
// constructors and destructor
//
V0CandFilter::V0CandFilter(const edm::ParameterSet& iConfig):
  k0sCollection( iConfig.getParameter<edm::InputTag>("kShortCollection") ),
  lamCollection( iConfig.getParameter<edm::InputTag>("lambdaCollection") ),
  trackCollection( iConfig.getParameter<edm::InputTag>("trackCollection") ) {

}


V0CandFilter::~V0CandFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
V0CandFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<reco::VertexCompositeCandidateCollection> theKshorts;
   Handle<reco::VertexCompositeCandidateCollection> theLambdas;
   Handle<reco::TrackCollection> theTracks;

   iEvent.getByLabel( k0sCollection, theKshorts );
   iEvent.getByLabel( lamCollection, theLambdas );
   iEvent.getByLabel( trackCollection, theTracks );

   int tkSize = theTracks->size();

   bool filtersel = false;
   // Loop over K0S:
   for( reco::VertexCompositeCandidateCollection::const_iterator iKS = theKshorts->begin();
	iKS != theKshorts->end();
	iKS++ ) {
     bool ksmsel = false;
     double kschi2 = iKS->vertexNormalizedChi2();
     double ksm = iKS->mass();
     double ksx = iKS->vertex().x();
     double ksy = iKS->vertex().y();
     double kseta = iKS->momentum().eta();
     double ksr = sqrt(ksx*ksx + ksy*ksy);
     if( ksm > (ksMass - 0.005) && ksm < (ksMass + 0.005) ) {
       ksmsel = true;
     }
     if( ksmsel && kschi2 < 3. && fabs(kseta) < 2. 
	 && (ksr > 9. && ksr < 20.) 
	 && tkSize < 30)
       filtersel = true;
   }

   for( reco::VertexCompositeCandidateCollection::const_iterator iLAM = theLambdas->begin();
	iLAM != theLambdas->end();
	iLAM++ ) {
     bool lammsel = false;
     double lamchi2 = iLAM->vertexNormalizedChi2();
     double lamm = iLAM->mass();
     double lamx = iLAM->vertex().x();
     double lamy = iLAM->vertex().y();
     double lamr = sqrt(lamx*lamx + lamy*lamy);
     double lameta = iLAM->momentum().eta();
     if( lamm > (lamMass - 0.005) && lamm < (lamMass + 0.005) ) {
       lammsel = true;
     }
     if( lammsel && lamchi2 < 3. && fabs(lameta) < 2. 
	 && (lamr > 9. && lamr < 20.) 
	 && tkSize < 30)
       filtersel = true;
   }
   //if( filtersel ) std::cout << "Filter found a good K0S or Lambda" << std::endl;
   //else if (!filtersel) std::cout << "Filter didn't find any good ones." << std::endl; 

   return filtersel;
}

// ------------ method called once each job just before starting event loop  ------------
void 
V0CandFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
V0CandFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0CandFilter);
