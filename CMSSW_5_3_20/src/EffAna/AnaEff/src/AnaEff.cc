// -*- C++ -*-
//
// Package:    AnaEff
// Class:      AnaEff
// 
/**\class AnaEff AnaEff.cc EffAna/AnaEff/src/AnaEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qiao Xu
//         Created:  Mon Feb 24 16:46:35 CST 2014
// $Id$
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

#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CastorReco/interface/CastorTower.h"
#include "DataFormats/CastorReco/interface/CastorCluster.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "CalibCalorimetry/CastorCalib/interface/CastorDbASCIIIO.h"
#include "CondFormats/CastorObjects/interface/CastorChannelQuality.h"
#include "CondFormats/DataRecord/interface/CastorChannelQualityRcd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include "TRandom3.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
//

class AnaEff : public edm::EDAnalyzer {
   public:
      explicit AnaEff(const edm::ParameterSet&);
      ~AnaEff();

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
   TH1D* hvz;
   TH1D* hMultDS;
   TH1D* hHFp4DS;
   TH1D* hHFpDS;
   TH1D* hHFm4DS;
   TH1D* hHFmDS;
   TH1D* hHF4DS;
   TH1D* hpi0;
   TH1D* hpip;
   TF1* fVz;
   edm::ESHandle < ParticleDataTable > pdt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AnaEff::AnaEff(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


AnaEff::~AnaEff()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AnaEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
   using namespace HepMC;

   iSetup.getData(pdt);

   double nPi0=0;
   double nPip=0;
   bool posDS = false; bool negDS = false;
	double weight=1.;
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
	
      Handle<HepMCProduct> mc;
      iEvent.getByLabel("generator",mc);
      const GenEvent* evt;
      evt = mc->GetEvent();

   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){

         //if((*it)->status() != 1) cout<<(*it)->status()<<endl;
         int pdg_id = (*it)->pdg_id();
         HepPDT::ParticleID particleID(pdg_id);
         double eta = (*it)->momentum().eta();
         //double phi = (*it)->momentum().phi();
         double pt = (*it)->momentum().perp();
         double energy = (*it)->momentum().e();
         const ParticleData * part = pdt->particle(pdg_id );
         double tau =  pdt->particle(pdg_id)->lifetime();
         //if(pt>0.3&&pdg_id==111) cout<<"tau="<<tau<<", eta="<<eta<<", status="<<(*it)->status()<<endl; 
         //if(pdg_id==111) { nPi0 +=1; hpi0->Fill(pt); cout<<"energy= "<<energy<<", sqrt(p^2+m^2)="<<sqrt(pt*pt+(*it)->momentum().pz()*(*it)->momentum().pz()+ (*it)->momentum().m2())<<endl; }
         if(pdg_id==111) { nPi0 +=1; hpi0->Fill(pt); }
         if(pdg_id==211) { nPip +=1; hpip->Fill(pt); }

         if(particleID.isValid() && part){
           if ( tau  > 1e-18 || tau == 0.0 )
           {
             if( energy > 3.0 && eta > 3.0 && eta < 5.0 ) posDS = true;
             if( energy > 3.0 && eta < -3.0 && eta > -5.0 ) negDS = true;
           }
 
         }

   }
   //cout<<"xxxxxxxxx per event, pi0="<<nPi0<<", pip="<<nPip<<endl;
   //if(posDS&&negDS) hMultDS->Fill(mult);

  // vertex and generalTracks
     double vx[1000];
     double vy[1000];
     double vz[1000];
     double vxError[1000];
     double vyError[1000];
     double vzError[1000];

       const reco::VertexCollection * recoVertices;
       edm::Handle<reco::VertexCollection> vertexCollection;
       //cout <<"srcVertex_[iv]= " <<srcVertex_[iv]<<endl;
       iEvent.getByLabel("offlinePrimaryVertices",vertexCollection);
       recoVertices = vertexCollection.product();
       unsigned int daughter = 0;
       //int nVertex = 0;
       int greatestvtx = 0;

       //nVertex = recoVertices->size();
       for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
         daughter = (*recoVertices)[i].tracksSize();
         if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = i;
           //cout <<"Vertex: "<< (*recoVertices)[i].position().z()<<" , tracksSize= "<<daughter<<",  vertex size="<<recoVertices->size()<<endl;
       }

       if(recoVertices->size()>0){
         vx[1] = (*recoVertices)[greatestvtx].position().x();
         vy[1] = (*recoVertices)[greatestvtx].position().y();
         vz[1] = (*recoVertices)[greatestvtx].position().z();
         vxError[1] = (*recoVertices)[greatestvtx].xError();
         vyError[1] = (*recoVertices)[greatestvtx].yError();
         vzError[1] = (*recoVertices)[greatestvtx].zError();
       }else{
         vx[1] =  -99;
         vy[1] =  -99;
         vz[1] =  -99;
         vxError[1] =  -99;
         vyError[1] =  -99;
         vzError[1] =  -99;
       }

    math::XYZPoint v1(vx[1],vy[1], vz[1]);

    int nTracks = 0;
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel("generalTracks",tracks);
    for(unsigned it=0; it<tracks->size(); ++it){
      const reco::Track & track = (*tracks)[it];
      //double eta = track.eta();
      //double phi = track.phi();
      //double pt  = track.pt();
       double dz = track.dz(v1);
       double dzsigma = sqrt(track.dzError()*track.dzError()+vzError[1]*vzError[1]);
       double dxy= track.dxy(v1);
       double dxysigma = sqrt(track.dxyError()*track.dxyError()+vxError[1]*vyError[1]);
       if(!( track.quality(reco::TrackBase::highPurity)!=1 || track.ptError()/track.pt()>0.1 || fabs(dz/dzsigma)>3.0 || fabs(dxy/dxysigma)>3.0)) {
         nTracks++;
        }
    } // tracks;

   //edm::Handle<reco::Centrality> centrality;
   //iEvent.getByLabel("pACentrality",centrality);

     double etHFtowerSumPlus=0; 
     double etHFtowerSumPlusEta4=0;
     double etHFtowerSumMinus=0; 
     double etHFtowerSumMinusEta4=0;
     //double etHFtowerSumEta4=0;
     Handle<CaloTowerCollection> towers;
     iEvent.getByLabel("towerMaker",towers);
        for( size_t i = 0; i<towers->size(); ++ i){
           const CaloTower & tower = (*towers)[ i ];
           double eta = tower.eta();
           bool isHF = tower.ietaAbs() > 29;
              if(isHF && eta > 0){
                 etHFtowerSumPlus += tower.pt();
                 if(eta > 4)  etHFtowerSumPlusEta4 += tower.pt();
              }
              if(isHF && eta < 0){
                 etHFtowerSumMinus += tower.pt();
                 if(eta < -4)  etHFtowerSumMinusEta4 += tower.pt();
              }
              //if(isHF && (eta > 4 || eta < -4) ){
              //   etHFtowerSumEta4 += tower.pt();
              //}
        }

//	weight*=fVz->Eval(vz[1]);
   if(posDS && negDS && fabs(vz[1])<=15.0){
     hHFp4DS->Fill(etHFtowerSumPlusEta4,weight); 
     hHFpDS->Fill(etHFtowerSumPlus,weight); 
     hHFm4DS->Fill(etHFtowerSumMinusEta4,weight); 
     hHFmDS->Fill(etHFtowerSumMinus,weight); 
     hHF4DS->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight); 
     hMultDS->Fill(nTracks,weight);    
   }

   hvz->Fill(vz[1]);


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
AnaEff::beginJob()
{

  TH1D::SetDefaultSumw2();
  hvz = fs->make<TH1D>("vz","",300,-30.0,30.0);
  hpi0 = fs->make<TH1D>("pi0","",300,0.0,30.0);
  hpip = fs->make<TH1D>("pip","",300,0.0,30.0);

  hMultDS = fs->make<TH1D>("hMultDS","",500,0,500);
  hHFp4DS = fs->make<TH1D>("hHFp4DS","",20000,0,200);
  hHFpDS = fs->make<TH1D>("hHFpDS","",20000,0,200);
  hHFm4DS = fs->make<TH1D>("hHFm4DS","",20000,0,200);
  hHFmDS = fs->make<TH1D>("hHFmDS","",20000,0,200);
  hHF4DS = fs->make<TH1D>("hHF4DS","",20000,0,200);
  fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnaEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
AnaEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
AnaEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
AnaEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
AnaEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AnaEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnaEff);
