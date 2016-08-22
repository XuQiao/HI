// -*- C++ -*-
//
// Package:    LYZfirst
// Class:      LYZfirst
// 
/**\class LYZfirst LYZfirst.cc Analysis/LYZfirst/src/LYZfirst.cc

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

const double j01=2.404826;
const int ntheta = 5;
const int nstepr=200;
const int nbin=1;

//
// class declaration
//

class LYZfirst : public edm::EDAnalyzer {
   public:
      explicit LYZfirst(const edm::ParameterSet&);
      ~LYZfirst();

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
        TVectorD Qx1, Qy1, Q2, Nevent, totmultall;
	double Qx2[nbin],Qy2[nbin];
	double theta[ntheta];
        double r[nbin][nstepr];
	TVectorD GRe[nbin][ntheta];
	TVectorD GIm[nbin][ntheta];
     	edm::Service<TFileService> fs;
	edm::InputTag srcTrack;
	double vzcutmin;
	double vzcutmax;
	double etaMin;
	double etaMax;
	double ptMin;
	double ptMax;
	std::vector<int> trkbin;
	TH1D* hMult[nbin];
	TH1D* hpt;
	TH1D* hptFull;
	TH1D* heta;
	TH1D* hetaFull;
	TH1D* hphi;
	TH1D* hphiFull;
	TH1D* hvz;
	TH1D * hNtrkBin;
	TH1D * hNtrk;
	TH2D * hNtrkvsMult;
	TH1D * hNtrkFull;
	TH1D * hNtrkMult[nbin];	
};

//
// constructors and destructor
//
LYZfirst::LYZfirst(const edm::ParameterSet& iConfig)

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
}


LYZfirst::~LYZfirst()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
LYZfirst::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

	TComplex g[ntheta][nstepr];
        int mult; double phi[10000], eta[10000], pt[10000];

	double Q[ntheta];
	double Qx=0,Qy=0;	mult=0;
	for(int itheta=0;itheta<ntheta;itheta++){
       			Q[itheta]=0;
                for(int ir=0; ir<nstepr; ir++)
                	if(way==0)
                        	g[itheta][ir]=1.;
                        else
                        	g[itheta][ir]=0.;
        }
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
      double eta_ = track.eta();
      double phi_ = track.phi();
      double pt_  = track.pt();
      hptFull->Fill(pt_);
      hetaFull->Fill(eta_);
      hphiFull->Fill(phi_);
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
      // if(track.pt()>0.1 && fabs(track.eta())<2.4  && fabs(dz/dzsigma)<3.0) {
       //if( track.quality(reco::TrackBase::highPurity) && track.pt()>0.1 && fabs(track.eta())<2.4 && fabs(dz/dzsigma)<5.0 && fabs(dxy/dxysigma)<5.0) {  //looser
         //if(track.pt()>0.1 && track.ptError()/track.pt()<0.05 && fabs(dz/dzsigma)<2.0 && fabs(dxy/dxysigma)<2.0){
        //if(fabs(track.eta())<2.4 && track.pt()>0.1){
       		if(track.pt()>ptMin&&track.pt()<ptMax&&track.eta()>etaMin&&track.eta()<etaMax){
         		phi[mult]=track.phi(); pt[mult]=track.pt(); eta[mult]=track.eta();
           		mult++;
        	}
	}
   } // tracks;

   	int xbin=-1;
   	for(int j=0;j<nbin;j++)
      		if(nTracks<trkbin[j]&&nTracks>=trkbin[j+1])
        		xbin=j;
   	hNtrkFull->Fill(nTracks);
   	hNtrkvsMult->Fill(nTracks,mult);
        if(xbin>=0 && xbin<nbin){
		hvz->Fill(vz);
   		hNtrk->Fill(nTracks);
   		hNtrkBin->Fill(xbin);
		hNtrkMult[xbin]->Fill(nTracks);
	for(int imult=0;imult<mult;imult++){
      		hpt->Fill(pt[imult]);
      		heta->Fill(eta[imult]);
      		hphi->Fill(phi[imult]);
		Qx+=1.*cos(n*phi[imult]);
		Qy+=1.*sin(n*phi[imult]);
		for(int itheta=0;itheta<ntheta;itheta++){
                        double temp=TMath::Cos(n*(phi[imult]-theta[itheta]));
			for(int ir=0; ir<nstepr; ir++)
	                           if(way==0)
					g[itheta][ir]*=TComplex(1.,r[xbin][ir]*temp);
                	          // else 
				//	g[itheta][0]+=temp;
                }
		totmultall[xbin]++;
	}
		Qx1[xbin]+=Qx;	Qy1[xbin]+=Qy;
		Qx2[xbin]+=Qx*Qx;	Qy2[xbin]+=Qy*Qy;
		for(int itheta=0;itheta<ntheta;itheta++){
			Q[itheta]=Qx*TMath::Cos(n*theta[itheta])+Qy*TMath::Sin(n*theta[itheta]);
			for(int ir=0; ir<nstepr; ir++){
				if(way!=0)	g[itheta][ir]=TComplex::Exp(TComplex(0,r[xbin][ir]*Q[itheta]));//g[itheta][0]
				GRe[xbin][itheta][ir]+=g[itheta][ir].Re();
				GIm[xbin][itheta][ir]+=g[itheta][ir].Im();
			}
		}

	Nevent[xbin]++;
	hMult[xbin]->Fill(mult);
	}	//ntrkoffline cut
	}	//vz cut
}


// ------------ method called once each job just before starting event loop  ------------
void 
LYZfirst::beginJob()
{
  double Vmax[nbin], eps[nbin];
for(int ibin=0; ibin<nbin ;ibin++){
        Vmax[ibin]=0.065*(trkbin[ibin]+30);
        eps[ibin]=0.00025*(trkbin[ibin]+30);
}

  TH1::SetDefaultSumw2();
  for(int ibin=0; ibin<nbin; ibin++){
   hMult[ibin] = fs->make<TH1D>(Form("hMult_%d",ibin),Form("Multiplicity distribution in %d bin",ibin),2000,0,2000);
   hNtrkMult[ibin] = fs->make<TH1D>(Form("hNtrkMult_%d",ibin),Form("NTrackOffline Distribution in %d bin",ibin), 500, 0., 500.);
  }
   hNtrkvsMult = fs->make<TH2D>(Form("hNtrkvsMult"),Form("NTrackOffline-Multiplicity Distribution"), 500, 0., 500., 2000,0.,2000.);

  hvz = fs->make<TH1D>("vz","track vz",200,-30,30);
  hptFull = fs->make<TH1D>("ptFull","track pt before cut",200,0,50);
  hpt = fs->make<TH1D>("pt","track pt after cut",200,0,50);
  hetaFull = fs->make<TH1D>("etaFull","track eta before cut",200,-4,4);
  heta = fs->make<TH1D>("eta","track eta after cut",200,-4,4);
  hphiFull = fs->make<TH1D>("phiFull","track phi before cut",200,-4,4);
  hphi = fs->make<TH1D>("phi","track phi after cut",200,-4,4);
	
	hNtrkBin= fs->make<TH1D>("hntrkbin","NTrackOffline Bin Distribution", 10, 0., 10.);
	hNtrk = fs->make<TH1D>("hntrk","NTrackOffline Distribution", 800, 0., 800.);
	hNtrkFull = fs->make<TH1D>("hntrkfull","NTrackOffline Distribution Full", 800, 0., 800.);
	Nevent.ResizeTo(nbin);	totmultall.ResizeTo(nbin);	Qx1.ResizeTo(nbin);	Qy1.ResizeTo(nbin);	Q2.ResizeTo(nbin);
	Nevent.Zero();	totmultall.Zero();	Qx1.Zero();	Qy1.Zero();	Q2.Zero();
	for(int ibin=0; ibin<nbin; ibin++){
        	Qx2[ibin]=0;  Qy2[ibin]=0;
	}

        	for(int itheta=0;itheta<ntheta;itheta++){
               		theta[itheta]=itheta*TMath::Pi()/ntheta/n;
			for(int ibin=0; ibin<nbin; ibin++){
				GRe[ibin][itheta].ResizeTo(nstepr); 	GRe[ibin][itheta].Zero();
				GIm[ibin][itheta].ResizeTo(nstepr);	GIm[ibin][itheta].Zero();
        	}
	}
	
	for(int ibin=0; ibin<nbin; ibin++){
        	for(int ir=0; ir<nstepr; ir++)
                	r[ibin][ir]=j01/(Vmax[ibin]-eps[ibin]*ir);
        //r[ir]=eps*(ir+1);

	}
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LYZfirst::endJob() 
{
        
/*	ofstream fstrGRe, fstrGIm, fstrQx1, fstrQy1, fstrQ2, fstrmult;
        fstrGRe.open("GRe.txt",std::ios::app);
        fstrGIm.open("GIm.txt",std::ios::app);
        fstrQx1.open("Qx1.txt",std::ios::app);
        fstrQy1.open("Qy1.txt",std::ios::app);
        fstrQ2.open("Q2.txt",std::ios::app);
        fstrmult.open("mult.txt",std::ios::app);
*/
	for(int ibin=0; ibin<nbin; ibin++){
		Q2[ibin]+=Qx2[ibin]+Qy2[ibin];}
/*		fstrmult<<Nevent[ibin]<<"\t";
        	fstrmult<<totmultall[ibin]<<std::endl;
                for(int ir=0; ir<nstepr; ir++){
                        for(int itheta=0;itheta<ntheta;itheta++){
                                fstrGRe<<GRe[ibin][itheta][ir]<<"\t";
                                fstrGIm<<GIm[ibin][itheta][ir]<<"\t";
			}
                        fstrGRe<<std::endl;
                        fstrGIm<<std::endl;
                }       
                fstrGRe<<std::endl<<std::endl;
                fstrGIm<<std::endl<<std::endl;
                fstrQx1<<Qx1[ibin]<<std::endl;
                fstrQy1<<Qy1[ibin]<<std::endl;
                fstrQ2<<Qx2[ibin]+Qy2[ibin]<<std::endl;
        }
	fstrmult.close();
        fstrGRe.close();
        fstrGIm.close();
        fstrQx1.close();
        fstrQy1.close();
        fstrQ2.close();
*/
	fs->cd();
	for(int ibin=0; ibin<nbin; ibin++){
        	for(int itheta=0;itheta<ntheta;itheta++){
			GRe[ibin][itheta].Write(Form("GRe_%d_%d",ibin,itheta));
			GIm[ibin][itheta].Write(Form("GIm_%d_%d",ibin,itheta));
		}
	}
	Nevent.Write("Nevent");
	totmultall.Write("totmultall");
	Qx1.Write("Qx1");
	Qy1.Write("Qy1");
	Q2.Write("Q2");
}

// ------------ method called when starting to processes a run  ------------
void 
LYZfirst::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
LYZfirst::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LYZfirst::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LYZfirst::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LYZfirst::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LYZfirst);
