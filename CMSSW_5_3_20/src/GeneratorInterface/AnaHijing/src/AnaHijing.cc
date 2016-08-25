// -*- C++ -*-
//
// Package:    AnaHijing
// Class:      AnaHijing
// 

// system include files
#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <TH1D.h>
#include <TH2D.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// root include file
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

static const int MAXPARTICLES = 5000000;
static const int MAXVTX = 1000;
static const int ETABINS = 3; // Fix also in branch string

//
// class decleration
//

struct HIJINGEvent{

   int event;
   float b;
   float npart;
   float ncoll;
//   float nhard;
//   float phi0;

   int n[ETABINS];
   float ptav[ETABINS];

   int mult;
   float pt[MAXPARTICLES];
   float eta[MAXPARTICLES];
   float phi[MAXPARTICLES];
   int pdg[MAXPARTICLES];
   int chg[MAXPARTICLES];

//   float vx;
//   float vy;
//   float vz;
//   float vr;

};

class AnaHijing : public edm::EDAnalyzer {
   public:
      explicit AnaHijing(const edm::ParameterSet&);
      ~AnaHijing();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

   std::ofstream out_b;
   std::string fBFileName;

   std::ofstream out_n;
   std::string fNFileName;

   std::ofstream out_m;
   std::string fMFileName;

  
   TTree* hydjetTree_;
   HIJINGEvent hev_;

   TNtuple *nt;

   std::string output;           // Output filename
 
   bool doAnalysis_;
   bool printLists_;
   bool doCF_;
   bool doVertex_;
   double etaMax_;
   double ptMin_;
  edm::InputTag src_;
   edm::InputTag simVerticesTag_;

   edm::ESHandle < ParticleDataTable > pdt;
   edm::Service<TFileService> f;

   TH1D* hpt_1;
   TH1D* heta_1;
   TH1D* hphi_1;

	TH1D* hpt_2;
	TH1D* heta_2;
	TH1D* hphi_2;	

	TH1D* hInvMass;

	TH2D* hnpvsncp;
	TH2D* hnpvsnc;
	TH2D* hnpvsb;

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
AnaHijing::AnaHijing(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   fBFileName = iConfig.getUntrackedParameter<std::string>("output_b", "b_values.txt");
   fNFileName = iConfig.getUntrackedParameter<std::string>("output_n", "n_values.txt");
   fMFileName = iConfig.getUntrackedParameter<std::string>("output_m", "m_values.txt");
   doAnalysis_ = iConfig.getUntrackedParameter<bool>("doAnalysis", true);
   printLists_ = iConfig.getUntrackedParameter<bool>("printLists", true);
   doCF_ = iConfig.getUntrackedParameter<bool>("doMixed", false);
   doVertex_ = iConfig.getUntrackedParameter<bool>("doVertex", false);
   if (doVertex_) {
      simVerticesTag_ = iConfig.getParameter<edm::InputTag>("simVerticesTag");
   }
   etaMax_ = iConfig.getUntrackedParameter<double>("etaMax", 2);
   ptMin_ = iConfig.getUntrackedParameter<double>("ptMin", 0);
   src_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("generator"));
}


AnaHijing::~AnaHijing()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
AnaHijing::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace HepMC;
  
   iSetup.getData(pdt);

   hev_.event = iEvent.id().event();
/*   for(int ieta = 0; ieta < ETABINS; ++ieta){
      hev_.n[ieta] = 0;
      hev_.ptav[ieta] = 0;
   }*/
   hev_.mult = 0;
      
//   double phi0 = 0;
   double b = -1;
   int npart = -1;
   int ncoll = -1;
//   int nhard = -1;
//   double vx = -99;
//   double vy = -99;
//   double vz = -99;
//   double vr = -99;
   const GenEvent* evt;
//   const GenEvent* evt2;
  
//   int nmix = -1;
   int np = 0;
//   int sig = -1;
//   int src = -1;

   if(doCF_){

     Handle<CrossingFrame<HepMCProduct> > cf;
     iEvent.getByLabel(InputTag("mix","source"),cf);


     MixCollection<HepMCProduct> mix(cf.product());

  //   nmix = mix.size();

     cout<<"Mix Collection Size: "<<mix<<endl;

     MixCollection<HepMCProduct>::iterator mbegin = mix.begin();
     MixCollection<HepMCProduct>::iterator mend = mix.end();
     
     for(MixCollection<HepMCProduct>::iterator mixit = mbegin; mixit != mend; ++mixit){

       const GenEvent* subevt = (*mixit).GetEvent();
       int all = subevt->particles_size();
       np += all;

       /*
       HepMC::GenEvent::particle_const_iterator begin = subevt->particles_begin();
       HepMC::GenEvent::particle_const_iterator end = subevt->particles_end();
       for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
	 if((*it)->status() == 1){
           float pdg_id = (*it)->pdg_id();
           float eta = (*it)->momentum().eta();
           float pt = (*it)->momentum().perp();
           const ParticleData * part = pdt->particle(pdg_id );
           float charge = part->charge();
	   }
	 }
       }

       */
       
     }     

   }
   
      
   Handle<HepMCProduct> mc;
   iEvent.getByLabel(src_,mc);
   evt = mc->GetEvent();
      
      Handle<HepMCProduct> mc2;
      iEvent.getByLabel(src_,mc2);
//      evt2 = mc2->GetEvent();
   
   const HeavyIon* hi = evt->heavy_ion();
   if(hi){
      b = hi->impact_parameter();
      npart = hi->Npart_proj()+hi->Npart_targ();
      ncoll = hi->Ncoll();
//      nhard = hi->Ncoll_hard();
//      phi0 = hi->event_plane_angle();
      
      if(printLists_){
	 out_b<<b<<endl;
	 out_n<<npart<<endl;
      }
   }
   
//   src = evt->particles_size();
//   sig = evt2->particles_size();
   
   if(b<1.5){
   HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
   HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
   for(HepMC::GenEvent::particle_const_iterator it1 = begin; it1 != end; ++it1){ 
	for(HepMC::GenEvent::particle_const_iterator it2 = begin; it2 != end; ++it2){
	    if((*it1)->status() == 1 && (*it2)->status() == 1){
//if((*it)->status() != 1) cout<<(*it)->status()<<endl;
	 int pdg_id_1 = (*it1)->pdg_id();
	 int pdg_id_2 = (*it2)->pdg_id();
	
	if(pdg_id_1 == 211 && pdg_id_2 == -211){	//3122--lambda	2212--proton 211--pi+ -211--pi- 2112--neutron

	 //float eta_1 = (*it1)->momentum().eta();
	 //float phi_1 = (*it1)->momentum().phi();
	 //float pt_1 = (*it1)->momentum().perp();
	 float e_1 = (*it1)->momentum().e();
	 float px_1 = (*it1)->momentum().px();
	 float py_1 = (*it1)->momentum().py();
	 float pz_1 = (*it1)->momentum().pz();

	 //float eta_2 = (*it2)->momentum().eta();
	 //float phi_2 = (*it2)->momentum().phi();
	 //float pt_2 = (*it2)->momentum().perp();
	 float e_2 = (*it2)->momentum().e();
	 float px_2 = (*it2)->momentum().px();
	 float py_2 = (*it2)->momentum().py();
	 float pz_2 = (*it2)->momentum().pz();

	// const ParticleData * part_1 = pdt->particle(pdg_id_1 );
 	// const ParticleData * part_2 = pdt->particle(pdg_id_2 );

	 //int charge_1 = static_cast<int>(part_1->charge());
	 //int charge_2 = static_cast<int>(part_2->charge());
	 
//	 if(charge_1!=0&&pt_1>0.3&&fabs(eta_1)<4.0&&phi_1!=0){
//	 if(charge_2!=0&&pt_2>0.3&&fabs(eta_2)<4.0&&phi_2!=0){
	
	
//	 if(pt>0.3&&fabs(eta)<5.0&&phi!=0){
/*	 hev_.pt[hev_.mult] = pt;
	 hev_.eta[hev_.mult] = eta;
	 hev_.phi[hev_.mult] = phi;*/
	 //hev_.pdg[hev_.mult] = pdg_id;
	 //hev_.chg[hev_.mult] = charge;
/*       hpt_1->Fill(pt_1);
         heta_1->Fill(eta_1);
         hphi_1->Fill(phi_1);

         hpt_2->Fill(pt_2);
         heta_2->Fill(eta_2);
         hphi_2->Fill(phi_2);
*/
	/*
	 eta = fabs(eta);
	 int etabin = 0;
	 if(eta > 0.5) etabin = 1; 
	 if(eta > 1.) etabin = 2;
	 if(eta < 2.){
	    hev_.ptav[etabin] += pt;
	    ++(hev_.n[etabin]);
	 } */
	


	double invM2 = (e_1+e_2)*(e_1+e_2)-((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
	
	hInvMass->Fill(TMath::Sqrt(invM2));

	 ++(hev_.mult);

      	  }
	 }
	
    }
   }
   /*
   if(doVertex_){
      edm::Handle<edm::SimVertexContainer> simVertices;
      iEvent.getByLabel<edm::SimVertexContainer>(simVerticesTag_, simVertices);
      
      if (! simVertices.isValid() ) throw cms::Exception("FatalError") << "No vertices found\n";
      int inum = 0;
      
      edm::SimVertexContainer::const_iterator it=simVertices->begin();
      SimVertex vertex = (*it);
      cout<<" Vertex position "<< inum <<" " << vertex.position().rho()<<" "<<vertex.position().z()<<endl;
      vx = vertex.position().x();
      vy = vertex.position().y();
      vz = vertex.position().z();
      vr = vertex.position().rho();
   }
   
   for(int i = 0; i<3; ++i){
      hev_.ptav[i] = hev_.ptav[i]/hev_.n[i];
   }
   */
   hev_.b = b;
   hev_.npart = npart;
   hev_.ncoll = ncoll;
   //hev_.nhard= nhard;
   //hev_.phi0 = phi0;
   //hev_.vx = vx;
   //hev_.vy = vy;
   //hev_.vz = vz;
   //hev_.vr = vr;

   ///nt->Fill(nmix,np,src,sig);
	
//hnpvsncp->Fill(hev_.npart,hev_.mult);
//hnpvsnc->Fill(hev_.npart,hev_.ncoll);
//hnpvsb->Fill(hev_.npart,hev_.b);

  // hydjetTree_->Fill();
   //cout<<"one Event;"<<endl;
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
AnaHijing::beginJob()
{

   if(printLists_){
      out_b.open(fBFileName.c_str());
      if(out_b.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fBFileName;
      out_n.open(fNFileName.c_str());
      if(out_n.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fNFileName;
      out_m.open(fMFileName.c_str());
      if(out_m.good() == false)
	 throw cms::Exception("BadFile") << "Can\'t open file " << fMFileName;
   }   
   
   if(doAnalysis_){
 //     nt = f->make<TNtuple>("nt","Mixing Analysis","mix:np:src:sig");

   /*   hydjetTree_ = f->make<TTree>("hi","Tree of HIJING Events");
      hydjetTree_->Branch("event",&hev_.event,"event/I");
      hydjetTree_->Branch("b",&hev_.b,"b/F");
      hydjetTree_->Branch("npart",&hev_.npart,"npart/F");
      //hydjetTree_->Branch("ncoll",&hev_.ncoll,"ncoll/F");
      //hydjetTree_->Branch("nhard",&hev_.nhard,"nhard/F");
      //hydjetTree_->Branch("phi0",&hev_.phi0,"phi0/F");

      //hydjetTree_->Branch("n",hev_.n,"n[3]/I");
      //hydjetTree_->Branch("ptav",hev_.ptav,"ptav[3]/F");
      
      hydjetTree_->Branch("mult",&hev_.mult,"mult/I");
      hydjetTree_->Branch("pt",hev_.pt,"pt[mult]/F");
      hydjetTree_->Branch("eta",hev_.eta,"eta[mult]/F");
      hydjetTree_->Branch("phi",hev_.phi,"phi[mult]/F");
      //hydjetTree_->Branch("pdg",hev_.pdg,"pdg[mult]/I");
      //hydjetTree_->Branch("chg",hev_.chg,"chg[mult]/I");
*
      hydjetTree_->Branch("vx",&hev_.vx,"vx/F");
      hydjetTree_->Branch("vy",&hev_.vy,"vy/F");
      hydjetTree_->Branch("vz",&hev_.vz,"vz/F");
      hydjetTree_->Branch("vr",&hev_.vr,"vr/F");
*/

/*  
  hpt_1 = f->make<TH1D>("p_{T} of #pi^{+}","Tracks;p_{T}",200,0,50);
  heta_1 = f->make<TH1D>("#eta of #pi^{+}","Tracks;#eta",200,-4,4);
  hphi_1 = f->make<TH1D>("#phi of #pi^{+}","Tracks;#phi",200,-4,4);
	
//hnpvsncp = f->make<TH2D>("npart vs number of charged particles","Tracks;Npart;# of charged particles",100,0,100,500,0,500);
//hnpvsnc = f->make<TH2D>("npart vs ncoll","Tracks;Npart;Ncoll",100,0,100,100,0,100);
//hnpvsb = f->make<TH2D>("npart vs impart parameter","Tracks;Npart;b",100,0,100,200,0,20);


  hpt_2 = f->make<TH1D>("p_{T} of #pi^{-}","Tracks;p_{T}",200,0,50);
  heta_2 = f->make<TH1D>("#eta of #pi^{-}","Tracks;#eta",200,-4,4);
  hphi_2 = f->make<TH1D>("#phi of #pi{-}","Tracks;#phi",200,-4,4);

  hpt_1->Sumw2();
  heta_1->Sumw2();
  hphi_1->Sumw2();

  hpt_2->Sumw2();
  heta_2->Sumw2();
  hphi_2->Sumw2();
*/

hInvMass=f->make<TH1D>("Invariant Mass of K_{s}","Tracks;Invariant Mass",100,0.45,0.55);
//hInvMass->Sumw2();

   }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
AnaHijing::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(AnaHijing);
