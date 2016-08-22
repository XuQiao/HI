// -*- C++ -*-
//
// Package:    PythiaMC
// Class:      PythiaMC
// 
/**\class PythiaMC PythiaMC.cc dNchdeta/PythiaMC/plugins/PythiaMC.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qiao Xu
//         Created:  Sun, 24 Aug 2014 03:28:37 GMT
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1.h>
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
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

class PythiaMC : public edm::EDAnalyzer {
   public:
      explicit PythiaMC(const edm::ParameterSet&);
      ~PythiaMC();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      
   edm::Service<TFileService> fs;
   TH1D* hHFp4;
   TH1D* hHFp;
   TH1D* hHFm4;
   TH1D* hHFm;
   TH1D* hHF4;
   TH1D* hHF3;

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
PythiaMC::PythiaMC(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


PythiaMC::~PythiaMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PythiaMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;  
   using namespace reco;
   using namespace HepMC;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

     double etHFtowerSumPlus=0;
     double etHFtowerSumPlusEta4=0;
     double etHFtowerSumPlusEta3=0;
     double etHFtowerSumMinus=0;
     double etHFtowerSumMinusEta4=0;
     double etHFtowerSumMinusEta3=0;
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
                 if(eta > 3)  etHFtowerSumPlusEta3 += tower.pt();
              }
              if(isHF && eta < 0){
                 etHFtowerSumMinus += tower.pt();
                 if(eta < -4)  etHFtowerSumMinusEta4 += tower.pt();
                 if(eta < -3)  etHFtowerSumMinusEta3 += tower.pt();
              }
              //if(isHF && (eta > 4 || eta < -4) ){
              //   etHFtowerSumEta4 += tower.pt();
              //}
        }

     hHFp4->Fill(etHFtowerSumPlusEta4);
     hHFp->Fill(etHFtowerSumPlus);
     hHFm4->Fill(etHFtowerSumMinusEta4);
     hHFm->Fill(etHFtowerSumMinus);
     hHF4->Fill(etHFtowerSumPlusEta4 + etHFtowerSumMinusEta4);
     hHF3->Fill(etHFtowerSumPlusEta3 + etHFtowerSumMinusEta3);

}


// ------------ method called once each job just before starting event loop  ------------
void 
PythiaMC::beginJob()
{
  TH1D::SetDefaultSumw2();
  hHFp4 = fs->make<TH1D>("hHFp4","",2200,-20,200);
  hHFp = fs->make<TH1D>("hHFp","",2200,-20,200);
  hHFm4 = fs->make<TH1D>("hHFm4","",2200,-20,200);
  hHFm = fs->make<TH1D>("hHFm","",2200,-20,200);
  hHF4 = fs->make<TH1D>("hHF4","",2200,-20,200);
  hHF3 = fs->make<TH1D>("hHF3","",2200,-20,200);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PythiaMC::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
PythiaMC::beginRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a run  ------------

void 
PythiaMC::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------

void 
PythiaMC::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------

void 
PythiaMC::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PythiaMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PythiaMC);
