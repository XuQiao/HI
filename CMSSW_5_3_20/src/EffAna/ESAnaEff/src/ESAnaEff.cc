// -*- C++ -*-
//
// Package:    ESAnaEff
// Class:      ESAnaEff
// 
/**\class ESAnaEff ESAnaEff.cc EffAna/ESAnaEff/src/ESAnaEff.cc

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

//#include "EffAna/ESAnaEff/interface/Correction.C"
//#include "EffAna/ESAnaEff/interface/TrackCorrector2D.C"

//
// class declaration
//

class ESAnaEff : public edm::EDAnalyzer {
   public:
      explicit ESAnaEff(const edm::ParameterSet&);
      ~ESAnaEff();

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
   TH1D* hMultnonDS;
   TH1D* hMultDS;
   TH1D* hMultES;
   TH1D* hHFp4ES;
   TH1D* hHFpES;
   TH1D* hHFm4ES;
   TH1D* hHFmES;
   TH1D* hHF4ES;
   TF1* fVz;
/*
   TH1D* hNtrkEricMult;
   TH1D* hNtrkEricHFp4;
   TH1D* hNtrkEricHFp;
   TH1D* hNtrkEricHFm4;
   TH1D* hNtrkEricHFm;
   TH1D* hNtrkEricHF4;
	
   TH1D* hNtrkMult;
   TH1D* hNtrkHFp4;
   TH1D* hNtrkHFp;
   TH1D* hNtrkHFm4;
   TH1D* hNtrkHFm;
   TH1D* hNtrkHF4;

   TH1D* hNtrkbakMult;
   TH1D* hNtrkbakHFp4;
   TH1D* hNtrkbakHFp;
   TH1D* hNtrkbakHFm4;
   TH1D* hNtrkbakHFm;
   TH1D* hNtrkbakHF4;

   TH1D* hHFp4Mult;
   TH1D* hHFp4HFp4;
   TH1D* hHFp4HFp;
   TH1D* hHFp4HFm4;
   TH1D* hHFp4HFm;
   TH1D* hHFp4HF4;

   TH1D* hreHFp4Mult;
   TH1D* hreHFp4HFp4;
   TH1D* hreHFp4HFp;
   TH1D* hreHFp4HFm4;
   TH1D* hreHFp4HFm;
   TH1D* hreHFp4HF4;

   TH1D* hHFp4bakMult;
   TH1D* hHFp4bakHFp4;
   TH1D* hHFp4bakHFp;
   TH1D* hHFp4bakHFm4;
   TH1D* hHFp4bakHFm;
   TH1D* hHFp4bakHF4;

   TH1D* hreHFp4bakMult;
   TH1D* hreHFp4bakHFp4;
   TH1D* hreHFp4bakHFp;
   TH1D* hreHFp4bakHFm4;
   TH1D* hreHFp4bakHFm;
   TH1D* hreHFp4bakHF4;
*/
   TH1D* hpi0;
   TH1D* hpip;

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
ESAnaEff::ESAnaEff(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


ESAnaEff::~ESAnaEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ESAnaEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
//TrackCorrector2D corr_Eric;
//Corrector corr_Qiao("pPbHijing_EffCorr.root");

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
   //else hMultnonDS->Fill(mult);

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
       if(! (track.quality(reco::TrackBase::highPurity) !=1  || track.ptError()/track.pt()>0.1 || fabs(dz/dzsigma)>3.0 || fabs(dxy/dxysigma)>3.0)) {
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

 weight*=fVz->Eval(vz[1]);
   if(posDS && negDS && fabs(vz[1])<=15.0)
	hMultDS->Fill(nTracks,weight);
   if(!(posDS && negDS) && fabs(vz[1])<=15.0)
	hMultnonDS->Fill(nTracks,weight);
   if(fabs(vz[1])<=15.0){
     hHFp4ES->Fill(etHFtowerSumPlusEta4,weight); 
     hHFpES->Fill(etHFtowerSumPlus,weight); 
     hHFm4ES->Fill(etHFtowerSumMinusEta4,weight); 
     hHFmES->Fill(etHFtowerSumMinus,weight); 
     hHF4ES->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight); 
     hMultES->Fill(nTracks,weight); 
   /*
	weight = corr_Eric.getEventWeight(nTracks);
 hNtrkEricHFp4->Fill(etHFtowerSumPlusEta4,weight);
     hNtrkEricHFp->Fill(etHFtowerSumPlus,weight);
     hNtrkEricHFm4->Fill(etHFtowerSumMinusEta4,weight);
     hNtrkEricHFm->Fill(etHFtowerSumMinus,weight);
     hNtrkEricHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hNtrkEricMult->Fill(nTracks,weight);

	weight = corr_Qiao.getEventWeightNtrk(nTracks);
 hNtrkHFp4->Fill(etHFtowerSumPlusEta4,weight);
     hNtrkHFp->Fill(etHFtowerSumPlus,weight);
     hNtrkHFm4->Fill(etHFtowerSumMinusEta4,weight);
     hNtrkHFm->Fill(etHFtowerSumMinus,weight);
     hNtrkHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hNtrkMult->Fill(nTracks,weight);
	weight = corr_Qiao.getEventWeightNtrkbak(nTracks);
 hNtrkbakHFp4->Fill(etHFtowerSumPlusEta4,weight);
     hNtrkbakHFp->Fill(etHFtowerSumPlus,weight);
     hNtrkbakHFm4->Fill(etHFtowerSumMinusEta4,weight);
     hNtrkbakHFm->Fill(etHFtowerSumMinus,weight);
     hNtrkbakHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hNtrkbakMult->Fill(nTracks,weight);
	weight = corr_Qiao.getEventWeightHFPlus4(etHFtowerSumPlusEta4,kFALSE);
 hHFp4HFp4->Fill(etHFtowerSumPlusEta4,weight);
     hHFp4HFp->Fill(etHFtowerSumPlus,weight);
     hHFp4HFm4->Fill(etHFtowerSumMinusEta4,weight);
     hHFp4HFm->Fill(etHFtowerSumMinus,weight);
     hHFp4HF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hHFp4Mult->Fill(nTracks,weight);	
	weight = corr_Qiao.getEventWeightHFPlus4(etHFtowerSumPlusEta4,kTRUE);
 hreHFp4HFp4->Fill(etHFtowerSumPlusEta4,weight);
     hreHFp4HFp->Fill(etHFtowerSumPlus,weight);
     hreHFp4HFm4->Fill(etHFtowerSumMinusEta4,weight);
     hreHFp4HFm->Fill(etHFtowerSumMinus,weight);
     hreHFp4HF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hreHFp4Mult->Fill(nTracks,weight);
	weight = corr_Qiao.getEventWeightHFPlus4bak(etHFtowerSumPlusEta4,kFALSE);
 hHFp4bakHFp4->Fill(etHFtowerSumPlusEta4,weight);
     hHFp4bakHFp->Fill(etHFtowerSumPlus,weight);
     hHFp4bakHFm4->Fill(etHFtowerSumMinusEta4,weight);
     hHFp4bakHFm->Fill(etHFtowerSumMinus,weight);
     hHFp4bakHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hHFp4bakMult->Fill(nTracks,weight);
	weight = corr_Qiao.getEventWeightHFPlus4bak(etHFtowerSumPlusEta4,kTRUE);
 hreHFp4bakHFp4->Fill(etHFtowerSumPlusEta4,weight);
     hreHFp4bakHFp->Fill(etHFtowerSumPlus,weight);
     hreHFp4bakHFm4->Fill(etHFtowerSumMinusEta4,weight);
     hreHFp4bakHFm->Fill(etHFtowerSumMinus,weight);
     hreHFp4bakHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4,weight);
     hreHFp4bakMult->Fill(nTracks,weight);
*/

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
ESAnaEff::beginJob()
{

  TH1D::SetDefaultSumw2();
  hvz = fs->make<TH1D>("vz","",300,-30.0,30.0);
  hpi0 = fs->make<TH1D>("pi0","",300,0.0,30.0);
  hpip = fs->make<TH1D>("pip","",300,0.0,30.0);

  hMultDS = fs->make<TH1D>("hMultDS","",500,0,500);
  hMultnonDS = fs->make<TH1D>("hMultnonDS","",500,0,500);
  
  hMultES = fs->make<TH1D>("hMultES","",500,0,500);
  hHFp4ES = fs->make<TH1D>("hHFp4ES","",20000,0,200);
  hHFpES = fs->make<TH1D>("hHFpES","",20000,0,200);
  hHFm4ES = fs->make<TH1D>("hHFm4ES","",20000,0,200);
  hHFmES = fs->make<TH1D>("hHFmES","",20000,0,200);
  hHF4ES = fs->make<TH1D>("hHF4ES","",20000,0,200);
 fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);

/*
  hNtrkEricMult = fs->make<TH1D>("hNtrkEricMult","",500,0,500);
  hNtrkEricHFp4 = fs->make<TH1D>("hNtrkEricHFp4","",20000,0,200);
  hNtrkEricHFp = fs->make<TH1D>("hNtrkEricHFp","",20000,0,200);
  hNtrkEricHFm4 = fs->make<TH1D>("hNtrkEricHFm4","",20000,0,200);
  hNtrkEricHFm = fs->make<TH1D>("hNtrkEricHFm","",20000,0,200);
  hNtrkEricHF4 = fs->make<TH1D>("hNtrkEricHF4","",20000,0,200);


  hNtrkMult = fs->make<TH1D>("hNtrkMult","",500,0,500);
  hNtrkHFp4 = fs->make<TH1D>("hNtrkHFp4","",20000,0,200);
  hNtrkHFp = fs->make<TH1D>("hNtrkHFp","",20000,0,200);
  hNtrkHFm4 = fs->make<TH1D>("hNtrkHFm4","",20000,0,200);
  hNtrkHFm = fs->make<TH1D>("hNtrkHFm","",20000,0,200);
  hNtrkHF4 = fs->make<TH1D>("hNtrkHF4","",20000,0,200);

  hHFp4Mult = fs->make<TH1D>("hHFp4Mult","",500,0,500);
  hHFp4HFp4 = fs->make<TH1D>("hHFp4HFp4","",20000,0,200);
  hHFp4HFp = fs->make<TH1D>("hHFp4HFp","",20000,0,200);
  hHFp4HFm4 = fs->make<TH1D>("hHFp4HFm4","",20000,0,200);
  hHFp4HFm = fs->make<TH1D>("hHFp4HFm","",20000,0,200);
  hHFp4HF4 = fs->make<TH1D>("hHFp4HF4","",20000,0,200);

  hNtrkbakMult = fs->make<TH1D>("hNtrkbakMult","",500,0,500);
  hNtrkbakHFp4 = fs->make<TH1D>("hNtrkbakHFp4","",20000,0,200);
  hNtrkbakHFp = fs->make<TH1D>("hNtrkbakHFp","",20000,0,200);
  hNtrkbakHFm4 = fs->make<TH1D>("hNtrkbakHFm4","",20000,0,200);
  hNtrkbakHFm = fs->make<TH1D>("hNtrkbakHFm","",20000,0,200);
  hNtrkbakHF4 = fs->make<TH1D>("hNtrkbakHF4","",20000,0,200);

  hHFp4bakMult = fs->make<TH1D>("hHFp4bakMult","",500,0,500);
  hHFp4bakHFp4 = fs->make<TH1D>("hHFp4bakHFp4","",20000,0,200);
  hHFp4bakHFp = fs->make<TH1D>("hHFp4bakHFp","",20000,0,200);
  hHFp4bakHFm4 = fs->make<TH1D>("hHFp4bakHFm4","",20000,0,200);
  hHFp4bakHFm = fs->make<TH1D>("hHFp4bakHFm","",20000,0,200);
  hHFp4bakHF4 = fs->make<TH1D>("hHFp4bakHF4","",20000,0,200);

  hreHFp4Mult = fs->make<TH1D>("hreHFp4Mult","",500,0,500);
  hreHFp4HFp4 = fs->make<TH1D>("hreHFp4HFp4","",20000,0,200);
  hreHFp4HFp = fs->make<TH1D>("hreHFp4HFp","",20000,0,200);
  hreHFp4HFm4 = fs->make<TH1D>("hreHFp4HFm4","",20000,0,200);
  hreHFp4HFm = fs->make<TH1D>("hreHFp4HFm","",20000,0,200);
  hreHFp4HF4 = fs->make<TH1D>("hreHFp4HF4","",20000,0,200);

  hreHFp4bakMult = fs->make<TH1D>("hreHFp4bakMult","",500,0,500);
  hreHFp4bakHFp4 = fs->make<TH1D>("hreHFp4bakHFp4","",20000,0,200);
  hreHFp4bakHFp = fs->make<TH1D>("hreHFp4bakHFp","",20000,0,200);
  hreHFp4bakHFm4 = fs->make<TH1D>("hreHFp4bakHFm4","",20000,0,200);
  hreHFp4bakHFm = fs->make<TH1D>("hreHFp4bakHFm","",20000,0,200);
  hreHFp4bakHF4 = fs->make<TH1D>("hreHFp4bakHF4","",20000,0,200);
*/

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ESAnaEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ESAnaEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ESAnaEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ESAnaEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ESAnaEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ESAnaEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ESAnaEff);
