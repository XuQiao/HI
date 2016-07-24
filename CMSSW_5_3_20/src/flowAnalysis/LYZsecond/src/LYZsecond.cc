// -*- C++ -*-
//
// Package:    LYZsecond
// Class:      LYZsecond
// 
/**\class LYZsecond LYZsecond.cc Analysis/LYZsecond/src/LYZsecond.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qiao Xu
//         Created:  Wed Aug 20 14:20:34 CDT 2014
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

const double j01=2.404826;
const double ptbin[] = {0.3,0.4,0.7,1.2,1.7,2.2,2.7,3.4,4.4,5.4,6.0};
const int npt=sizeof(ptbin)/sizeof(double)-1;
const int nbin = 1;
const int ntheta = 5;
const int nstepr=200;


class LYZsecond : public edm::EDAnalyzer {
   public:
      explicit LYZsecond(const edm::ParameterSet&);
      ~LYZsecond();

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
       	int  m, n ,way;
       	double theta[ntheta];
	TVectorD* r0res[nbin];
	TVectorD Nevent, totmultall;
	TVectorD dDRe[nbin], dNRe[nbin][ntheta];
        TVectorD dDIm[nbin], dNIm[nbin][ntheta];
        TVectorD totmult[nbin];
        TVectorD totpt[nbin];
	TString res;
	
        edm::Service<TFileService> fs;
        edm::InputTag srcTrack;
	std::vector<int> trkbin;
        double vzcutmin;
        double vzcutmax;
        double etaMin;
        double etaMax;
        double ptMin;
        double ptMax;
	
	TH1D * hNtrk;

};

//
// constructors and destructor
//
LYZsecond::LYZsecond(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
        n=iConfig.getUntrackedParameter<int>("n");
        m=iConfig.getUntrackedParameter<int>("m");
        way=iConfig.getUntrackedParameter<int>("way");
        srcTrack=iConfig.getUntrackedParameter<edm::InputTag>("srcTrack");
        vzcutmin=iConfig.getUntrackedParameter<double>("vzcutmin");
        vzcutmax=iConfig.getUntrackedParameter<double>("vzcutmax");
        etaMin=iConfig.getUntrackedParameter<double>("etaMin");
        etaMax=iConfig.getUntrackedParameter<double>("etaMax");
        ptMin=iConfig.getUntrackedParameter<double>("ptMin");
        ptMax=iConfig.getUntrackedParameter<double>("ptMax");
	trkbin=iConfig.getUntrackedParameter<std::vector<int>>("trkbin");
	res=iConfig.getUntrackedParameter<std::string>("res");
}


LYZsecond::~LYZsecond()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LYZsecond::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   int mult=0; double phi[10000], pt[10000];// eta[10000];
   TComplex g0[ntheta], dDsum[ntheta], dNsum[ntheta][npt];

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
    int nTracks = 0;
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(srcTrack,tracks);

    for(unsigned it=0; it<tracks->size(); ++it){
      const reco::Track & track = (*tracks)[it];
      //if(pt>1&&pt<2) cout<<"eta= "<<eta<<"  pt="<<pt<<"  phi="<<phi<<" #####, track.size = "<<tracks->size()<<endl;
       double dz= track.dz(v1);
       double dzsigma = sqrt(track.dzError()*track.dzError()+vzError*vzError);
       double dxy= track.dxy(v1);
       double dxysigma = sqrt(track.dxyError()*track.dxyError()+vxError*vyError);
       if( track.quality(reco::TrackBase::highPurity) && track.pt()>0.4 && fabs(track.eta())<2.4 && track.ptError()/track.pt()<0.1 && fabs(dz/dzsigma)<3.0 && fabs(dxy/dxysigma)<3.0) {
         //if(track.pt()>0.4 && fabs(dz/dzsigma)<3.0 && fabs(dxy/dxysigma)<3.0) nTracks++;
         nTracks++;	//ntrkoffline
	}
       if( track.quality(reco::TrackBase::highPurity) && track.pt()>0.1 && fabs(track.eta())<2.4 && track.ptError()/track.pt()<0.1 && fabs(dz/dzsigma)<3.0 && fabs(dxy/dxysigma)<3.0) {
       //if( track.quality(reco::TrackBase::highPurity) && track.pt()>0.1 && fabs(track.eta())<2.4 && fabs(dz/dzsigma)<5.0 && fabs(dxy/dxysigma)<5.0) {  //looser
         //if(track.pt()>0.1 && track.ptError()/track.pt()<0.05 && fabs(dz/dzsigma)<2.0 && fabs(dxy/dxysigma)<2.0){
        //if(fabs(track.eta())<2.4 && track.pt()>0.1){
       		if(track.pt()>ptMin&&track.pt()<ptMax&&track.eta()>etaMin&&track.eta()<etaMax){
         		phi[mult]=track.phi(); pt[mult]=track.pt(); //eta[mult]=track.eta();
           		mult++;
        	}
	}
   } // tracks;

	for(int iptbin=0; iptbin<npt;iptbin++)
        	for(int itheta=0;itheta<ntheta;itheta++)
                        dNsum[itheta][iptbin]=0.;
                for(int itheta=0;itheta<ntheta;itheta++){
                        if(way==0)
                                g0[itheta]=1.;
                        else
                                g0[itheta]=0.;
                        dDsum[itheta]=0.;
                }
	
   	int xbin=-1;
   	for(int j=0;j<nbin;j++)
      		if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
        		xbin=j;
        if(xbin>=0 && xbin<nbin){
		hNtrk->Fill(nTracks);
		for(int imult=0;imult<mult;imult++){
                	int ipt=-1;
                	for(int iptbin=0;iptbin<npt; iptbin++){
                		if(pt[imult]>ptbin[iptbin] && pt[imult]<=ptbin[iptbin+1])
                        	ipt = iptbin;
                	}
                	if(ipt<0 || ipt==npt)   continue;
                	totpt[xbin][ipt]+=pt[imult];
                	totmult[xbin][ipt]++;
                	for(int itheta=0;itheta<ntheta;itheta++){
                		Double_t temp=TMath::Cos(n*(phi[imult]-theta[itheta]));
                        	TComplex temp1(1.,(*r0res[xbin])[itheta]*temp);
				
                        	if(way==0)
                        		g0[itheta]*=temp1;
                        	else
                        		g0[itheta]+=temp;
                       		dDsum[itheta]+=temp/temp1;
                        	dNsum[itheta][ipt]+=TMath::Cos(m*n*(phi[imult]-theta[itheta]))/temp1;
                        	}
       		totmultall[xbin]++;
        	}

        	for(int itheta=0;itheta<ntheta;itheta++){
        		dDRe[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Re();
                	dDIm[xbin][itheta]+=(g0[itheta]*dDsum[itheta]).Im();
                	for(int iptbin=0;iptbin<npt; iptbin++){
                		dNRe[xbin][itheta][iptbin]+=(g0[itheta]*dNsum[itheta][iptbin]).Re();
                        	dNIm[xbin][itheta][iptbin]+=(g0[itheta]*dNsum[itheta][iptbin]).Im();
                	}
        	}
        	Nevent[xbin]++;
        	}
	}

}


// ------------ method called once each job just before starting event loop  ------------
void 
LYZsecond::beginJob()
{
	TFile *fres = TFile::Open(res);
        for(int ibin=0;ibin<nbin;ibin++){
                r0res[ibin] = (TVectorD*)fres->Get(Form("D_%d/r0",ibin));
        }

        Nevent.ResizeTo(nbin);  totmultall.ResizeTo(nbin);
        Nevent.Zero();  totmultall.Zero();
	hNtrk = fs->make<TH1D>("hntrk","NTrackOffline Distribution", 800, 0., 800.);

        for(int itheta=0;itheta<ntheta;itheta++){
        	theta[itheta]=itheta*TMath::Pi()/ntheta/n;
        }

	for(int ibin=0; ibin<nbin; ibin++){
        	totmult[ibin].ResizeTo(npt);    totmult[ibin].Zero();
        	totpt[ibin].ResizeTo(npt);      totpt[ibin].Zero();
                dDRe[ibin].ResizeTo(ntheta);    dDRe[ibin].Zero();
                dDIm[ibin].ResizeTo(ntheta);    dDIm[ibin].Zero();
                for(int itheta=0;itheta<ntheta;itheta++){
                dNRe[ibin][itheta].ResizeTo(npt);       dNRe[ibin][itheta].Zero();
                dNIm[ibin][itheta].ResizeTo(npt);       dNIm[ibin][itheta].Zero();

                }
        }
        fres->Close();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LYZsecond::endJob() 
{
/*        ofstream fstrdNRe, fstrdNIm, fstrdDRe, fstrdDIm, fstrmult;
        fstrdNRe.open("dNRe.txt",ios::app);
        fstrdNIm.open("dNIm.txt",ios::app);
        fstrdDRe.open("dDRe.txt",ios::app);
        fstrdDIm.open("dDIm.txt",ios::app);
        fstrmult.open("mult_.txt",ios::app);
	 for(int ibin=0;ibin<nbin;ibin++){
                fstrmult<<Nevent[ibin]<<endl;
                for(int iptbin=0;iptbin<npt; iptbin++){
                        fstrmult<<(double)totmult[ibin][iptbin]/Nevent[ibin]<<"\t";
                        fstrmult<<(double)totpt[ibin][iptbin]/Nevent[ibin]<<"\t";
                }
                for(int itheta=0;itheta<ntheta;itheta++){
                        for(int iptbin=0;iptbin<npt; iptbin++){
                                fstrdNRe<<dNRe[ibin][itheta][iptbin]<<"\t";
                                fstrdNIm<<dNIm[ibin][itheta][iptbin]<<"\t";
                        }
                        fstrdNRe<<endl;
                        fstrdNIm<<endl;
                        fstrdDRe<<dDRe[ibin][itheta]<<"\t";
                        fstrdDIm<<dDIm[ibin][itheta]<<"\t";
                }
                fstrmult<<endl<<totmultall[ibin]<<endl;
        }
        fstrmult<<endl;
        fstrdNRe<<endl;
        fstrdNIm<<endl;
        fstrdDRe<<endl;
        fstrdDIm<<endl;
*/
	fs->cd();
        Nevent.Write("Nevent");
        totmultall.Write("totmultall");

        for(int ibin=0; ibin<nbin; ibin++){
                dDRe[ibin].Write(Form("dDRe_%d",ibin));
                dDIm[ibin].Write(Form("dDIm_%d",ibin));
                for(int itheta=0;itheta<ntheta;itheta++){
                        dNRe[ibin][itheta].Write(Form("dNRe_%d_%d",ibin,itheta));
                        dNIm[ibin][itheta].Write(Form("dNIm_%d_%d",ibin,itheta));
                }
                totmult[ibin].Write(Form("totmult_%d",ibin));
                totpt[ibin].Write(Form("totpt_%d",ibin));
        }

}

// ------------ method called when starting to processes a run  ------------
void 
LYZsecond::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
LYZsecond::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LYZsecond::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LYZsecond::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LYZsecond::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LYZsecond);
