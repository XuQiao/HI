// -*- C++ -*-
//
// Package:    ExampleGenAnalyzer
// Class:      ExampleGenAnalyzer
// 
/**\class ExampleGenAnalyzer ExampleGenAnalyzer.cc Tutorial/ExampleGenAnalyzer/src/ExampleGenAnalyzer.cc

 Description: Makes some simple spectra histograms from genParticles and genJets

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Eric Appelt
//         Created:  Wed Sep 10 12:18:24 CDT 2014
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <TH1.h>
#include <TH2.h>

//
// class declaration
//

class ExampleGenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ExampleGenAnalyzer(const edm::ParameterSet&);
      ~ExampleGenAnalyzer();

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
      edm::InputTag genJetSrc_;
      edm::InputTag genParticleSrc_;

      double etaMin_;
      double etaMax_;
      std::vector<double> ptBins_;
      std::vector<double> etaBins_;

      bool applyPthatCut_;
      double pthatMin_;
      double pthatMax_;

      std::map<std::string,TH1F*> hist_;
      std::map<std::string,TH2F*> hist2D_;

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
ExampleGenAnalyzer::ExampleGenAnalyzer(const edm::ParameterSet& iConfig):
genJetSrc_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
genParticleSrc_(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
etaMin_(iConfig.getParameter<double>("etaMin")),
etaMax_(iConfig.getParameter<double>("etaMax")),
ptBins_(iConfig.getParameter<std::vector<double> >("ptBins")),
etaBins_(iConfig.getParameter<std::vector<double> >("etaBins")),
applyPthatCut_(iConfig.getParameter<bool>("applyPthatCut")),
pthatMin_(iConfig.getParameter<double>("pthatMin")),
pthatMax_(iConfig.getParameter<double>("pthatMax"))
{

   edm::Service<TFileService> fs;

   // charged particle histograms
   hist_["ch1d"] = fs->make<TH1F>("ch1d",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist2D_["ch2d"] = fs->make<TH2F>("ch2d",";p_{T};#eta",
                           ptBins_.size()-1, &ptBins_[0],
                           etaBins_.size()-1, &etaBins_[0]);

   // neutral pion histograms
   hist_["pizero1d"] = fs->make<TH1F>("pizero1d",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist2D_["pizero2d"] = fs->make<TH2F>("pizero2d",";p_{T};#eta",
                           ptBins_.size()-1, &ptBins_[0],
                           etaBins_.size()-1, &etaBins_[0]);
   
   // genjet histograms
   hist_["jet1d"] = fs->make<TH1F>("jet1d",";p_{T};counts",
                           ptBins_.size()-1, &ptBins_[0]);
   hist2D_["jet2d"] = fs->make<TH2F>("jet2d",";p_{T};#eta",
                           ptBins_.size()-1, &ptBins_[0],
                           etaBins_.size()-1, &etaBins_[0]);

   // keep track of the number of events analyzed
   hist_["events"] = fs->make<TH1F>("events",";;events",1,0.,2.);

}


ExampleGenAnalyzer::~ExampleGenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExampleGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Skip events not in the specified pt-hat range.
   //
   // Note: Pythia does not respect max pt-hat for MB process  
   // where the minimum is 0. In this case we need to 
   // remove events over the pt-hat maximum by hand.

   edm::Handle<GenEventInfoProduct> genEvtInfo;
   iEvent.getByLabel("generator", genEvtInfo);
   
   if( applyPthatCut_ && genEvtInfo->qScale() < pthatMin_ ) return;
   if( applyPthatCut_ && genEvtInfo->qScale() > pthatMax_ ) return;

   // increment the events histogram
   hist_["events"]->Fill(1);

   // access vector of reco::GenJet objects
   // stored in edm::Event
   Handle<reco::GenJetCollection> gcol;
   iEvent.getByLabel(genJetSrc_,gcol);


   // access vector of reco::GenParticle objects
   // stored in edm::Event
   Handle<reco::GenParticleCollection> pcol;
   iEvent.getByLabel(genParticleSrc_,pcol);

   // loop over GenJet objects in the event and 
   // fill spectra histograms
   for( const reco::GenJet & jet : *gcol )
   {
     if( jet.eta() >= etaMin_ && jet.eta() <= etaMax_ )
       hist_["jet1d"]->Fill(jet.pt());
     hist2D_["jet2d"]->Fill(jet.pt(),jet.eta());
   }

   // loop over particle objects in the event
   // and fill spectra histograms   
   for( const reco::GenParticle & particle : *pcol )
   {
     // charged long lived particles
     if( particle.status() == 1 || particle.charge() != 0  ) 
     {
       if( particle.eta() >= etaMin_ && particle.eta() <= etaMax_ )
         hist_["ch1d"]->Fill(particle.pt());
       hist2D_["ch2d"]->Fill(particle.pt(),particle.eta());       

     }

     // neutral pions
     if( particle.pdgId() == 111 )
     {
       if( particle.eta() >= etaMin_ && particle.eta() <= etaMax_ )
         hist_["pizero1d"]->Fill(particle.pt());
       hist2D_["pizero2d"]->Fill(particle.pt(),particle.eta());       
     }

   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
ExampleGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExampleGenAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ExampleGenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ExampleGenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ExampleGenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ExampleGenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExampleGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExampleGenAnalyzer);
