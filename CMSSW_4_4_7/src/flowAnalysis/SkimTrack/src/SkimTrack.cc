// -*- C++ -*-
//
// Package:    SkimTrack
// Class:      SkimTrack
// 
/**\class SkimTrack SkimTrack.cc Analysis/SkimTrack/src/SkimTrack.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qiao Xu
//         Created:  Tue Aug  5 15:08:18 CDT 2014
// $Id$
//
//


// system include files
#include <memory>
#include "TComplex.h"
#include "TString.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

//
// class declaration
//

class SkimTrack : public edm::EDAnalyzer {
   public:
      explicit SkimTrack(const edm::ParameterSet&);
      ~SkimTrack();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
	
     	edm::Service<TFileService> fs;
	edm::InputTag srcTrack;
	double vzcutmin;
	double vzcutmax;
	double etaMin;
	double etaMax;
	double ptMin;
	double ptMax;
	int quality_;
	double ptEcut_, dzcut_, dxycut_;
	int mult, nTracks;
	std::vector<int> trkbin;
        double phi[10000], eta[10000], pt[10000];
	TTree *tree_;
};

//
// constructors and destructor
//
SkimTrack::SkimTrack(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
	srcTrack=iConfig.getUntrackedParameter<edm::InputTag>("srcTrack");
	quality_=iConfig.getUntrackedParameter<int>("quality");
	dzcut_=iConfig.getUntrackedParameter<double>("dzcut");
	dxycut_=iConfig.getUntrackedParameter<double>("dxycut");
	ptEcut_=iConfig.getUntrackedParameter<double>("ptEcut");
  	vzcutmin=iConfig.getUntrackedParameter<double>("vzcutmin");
  	vzcutmax=iConfig.getUntrackedParameter<double>("vzcutmax");
  	etaMin=iConfig.getUntrackedParameter<double>("etaMin");
  	etaMax=iConfig.getUntrackedParameter<double>("etaMax");
  	ptMin=iConfig.getUntrackedParameter<double>("ptMin");
  	ptMax=iConfig.getUntrackedParameter<double>("ptMax");
	trkbin=iConfig.getUntrackedParameter<std::vector<int> >("trkbin");
}


SkimTrack::~SkimTrack()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SkimTrack::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
    
  mult=0;	nTracks=0;

  double vx=-999.;
  double vy=-999.;
  double vz=-999.;
  double vxError=-999.;
  double vyError=-999.;
  double vzError=-999.;
  math::XYZVector vtxPos(0,0,0);
    // only looks at the first vertex collection
    const reco::VertexCollection * recoVertices;
    edm::Handle<reco::VertexCollection> vertexCollection;
    iEvent.getByLabel("offlinePrimaryVertices",vertexCollection);
    recoVertices = vertexCollection.product();
    unsigned int daughter = 0;
    int greatestvtx = 0;
    
    for (unsigned int j = 0 ; j< recoVertices->size(); ++j){
      daughter = (*recoVertices)[j].tracksSize();
      if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = j;
    }

    if(recoVertices->size()>0){
      vx = (*recoVertices)[greatestvtx].position().x();
      vy = (*recoVertices)[greatestvtx].position().y();
      vz = (*recoVertices)[greatestvtx].position().z();
      vxError = (*recoVertices)[greatestvtx].xError();
      vyError = (*recoVertices)[greatestvtx].yError();
      vzError = (*recoVertices)[greatestvtx].zError();
    }

       math::XYZPoint v1(vx,vy, vz);
    
    if(fabs(vz)<vzcutmax){
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(srcTrack,tracks);

    for(unsigned it=0; it<tracks->size(); ++it){
      const reco::Track & track = (*tracks)[it];
       
       double dz= track.dz(v1);
       double dzsigma = sqrt(track.dzError()*track.dzError()+vzError*vzError);
       double dxy= track.dxy(v1);
       double dxysigma = sqrt(track.dxyError()*track.dxyError()+vxError*vyError);

       if(track.quality(reco::TrackBase::highPurity) && track.pt()>0.4 && fabs(track.eta())<2.4 && track.ptError()/track.pt()<0.1 && fabs(dz/dzsigma)<3.0 && fabs(dxy/dxysigma)<3.0)
		nTracks++;
	if(track.pt()<ptMin || track.pt()>ptMax)	continue;
	if(track.eta()<etaMin || track.eta()>etaMax)	continue;
	if(quality_==0){
		if(!track.quality(reco::TrackBase::highPurity))	continue;
	}	
	else if(quality_==1){
		if(!track.quality(reco::TrackBase::loose))	continue;
	}
	else{
		if(!track.quality(reco::TrackBase::tight))	continue;
	}
       if(track.pt()>0.1  && track.ptError()/track.pt()< ptEcut_ && fabs(dz/dzsigma)<dzcut_ && fabs(dxy/dxysigma)<dxycut_) {
	pt[mult]=track.pt();	eta[mult]=track.eta();	phi[mult]=track.phi();
	mult++;
	}

   	} // tracks;
	if(nTracks<trkbin[0] && nTracks>=trkbin[trkbin.size()-1])
		tree_->Fill();
	}//vz
}


// ------------ method called once each job just before starting event loop  ------------
void 
SkimTrack::beginJob()
{
  tree_=fs->make<TTree>("TrackTree","track skim tree");
  tree_->Branch("mult",&(this->mult),"mult/I");
  tree_->Branch("ntrk",&(this->nTracks),"ntrk/I");
  tree_->Branch("pt",this->pt,"pt[mult]/D");
  tree_->Branch("eta",this->eta,"eta[mult]/D");
  tree_->Branch("phi",this->phi,"phi[mult]/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SkimTrack::endJob() 
{
        
}

// ------------ method called when starting to processes a run  ------------
void 
SkimTrack::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SkimTrack::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SkimTrack::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SkimTrack::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SkimTrack::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SkimTrack);
