#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include <TH1D.h>
#include <TNtuple.h>
#include "TChain.h"
#include <iostream>
#include <TString.h>
#include <TCut.h>
#include "hiForest.h"
#include "commonSetup.h"
#include "TrackCorrector.C"

using namespace std;

//define the kinematics cuts for the analysis

const double etacut = 2.0 ;
const double trackcut = 0.0;
const double tracketacut = 1.0; //make sure the tracks around jet cone inside eta range

int pthat; // //30 and 80 for pp; 30,50,80,120,170,200 for PbPb MC 
int ptmax; // = 170 ; = 9999 ;
TString coll = "PP"; //"PP2011"; 
TString TrigName = ""; //"Jet40" ; "Jet60" ; "Jet80"; "Jet100" ;"JetMB";
TString para ="Full"; //PYQUEN parameters setting: Wide or Full
const bool DoGenAna= kTRUE; //should be only be true for MC sample

//do tracking corrections
const bool doTrackCorrections = kTRUE; //for tracking efficiency correction
TString corrMet = "Hist" ; //Hist table from Frank, or Para from Krisztian
TrackingParam *trackCorrFromParam;

vector<TrackingCorrections*> trackCorrections;

//weight factor from vertex and centrality
double weight =1. ;

TString intputFile ;

//TString dataPath="/Users/ymao/group/CMS/hiForest";
TString dataPath ;

//if it is pp, no centrality bins, only one
const int nbin = 1 ;
const int centr[] ={0,100};
////for HI centrality bin
//const int nbin = 4 ;
//const int centr[] ={0,10,30,50,100};
//const int nbin = 6 ;
//const int centr[] ={0,5, 10, 30,50, 70, 90};

//const double pt[]={100., 120., 140., 160., 200., 300., 500.};
const int nptbin = 1 ;
const double pt[]={50., 500.};
//const int ntrkptbin = 6 ;
//const double trkpt[]={1., 2., 4., 8., 16., 32., 500.};

const int ntrkptbin = 4 ;
const double trkpt[]={1.,4., 16., 64., 500.};

const Double_t jetPtBin[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

//const int nPtBin = 27;
//double TrkBin[nPtBin+1] = {0, 0.5, 1, 1.203915, 1.449412, 1.74497, 2.100796, 2.529181, 3.04492, 3.665826, 4.413344, 5.313293, 6.396755, 7.701152, 9.271536, 11.16214, 13.43828, 16.17855, 19.47761, 23.44939, 28.23108, 33.98783, 40.91848, 49.26238, 59.30774, 71.40151, 85.96137, 103.4902}; 

//const Double_t TrkBin[]={0, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4., 5., 6., 8.,10.,12.,15.,20.,25.,30.,45.,60.,80.,100., 120.,150., 180.,300.,500.};
const Double_t TrkBin[]={0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
const Double_t TrkBin_Max[]={0,0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
//const Double_t TrkBin[]={0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 130}; 
const int nPtBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;
const int nPtBin_Max = sizeof(TrkBin_Max)/sizeof(Double_t)-1 ;

//const int netabin = 10 ;
//const double deta[]={-2.5, -2.0, -1.5, -1.0, -0.5, 0.0,0.5,1.0,1.5,2.0, 2.5} ;
const int netabin = 7 ;
const double deta[]={-1.8, -1.3, -0.8, -0.3, 0.3, 0.8,1.3,1.8} ;

const int ntrketabin = 7 ;
const double dtrketa[]={-1.8, -1.3, -0.8, -0.3, 0.3, 0.8,1.3,1.8} ;

class hist_class {
public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
    
    TString desc;
    bool IsMC ;
 /*   
    TH1F * Ntrack[nbin];
    TH1D * NEvents[nbin];
    TH1D * NevtCounter[nbin];
    TH1D * JetAbove40[nbin];
    TH1D * JetAbove60[nbin];
    TH1D * JetAbove75[nbin];
    TH1D * JetAbove95[nbin];
    TH1D * JetAbove120[nbin];
*/
    TH1F * jetpt[nbin];    
    TH2F * jetptEta[nbin];    
    TH2F * jetptphi[nbin]; //Added

    TH2F * jetptchMaxpt[nbin];	//Added
    TH2F * jetptchMaxpt_real[nbin];	//Added
    TH2F * jetptchMaxpt_fake[nbin];	//Added

    TH2F * jetptchMax[nbin];	//Added
    TH2F * jetptchMax_real[nbin];	//Added
    TH2F * jetptchMax_fake[nbin];	//Added
    
    TH2F * jetptneuMaxpt[nbin];	//Added
    TH2F * jetptneuMaxpt_real[nbin];	//Added
    TH2F * jetptneuMaxpt_fake[nbin];	//Added

    TH2F * jetptneuMaxr[nbin];   //Added
    TH2F * jetptneuMaxr_fake[nbin];      //Added
    TH2F * jetptneuMaxr_real[nbin];      //Added
   
    TH2F * jetptneuMax[nbin];	//Added
    TH2F * jetptneuMax_fake[nbin];	//Added
    TH2F * jetptneuMax_real[nbin];	//Added

    TH2F * jetptphoMaxpt[nbin]; //Added 
    TH2F * jetptphoMaxpt_fake[nbin];    //Added
    TH2F * jetptphoMaxpt_real[nbin];    //Added
    
    TH2F * jetptphoMax[nbin];   //Added
    TH2F * jetptphoMax_fake[nbin];      //Added
    TH2F * jetptphoMax_real[nbin];      //Added

    TH2F * jetptSumSum[nbin];	//Added
    TH2F * jetptSumSum_fake[nbin];	//Added
    TH2F * jetptSumSum_real[nbin];	//Added

    TH2F * jetptSumSumrawpt[nbin];   //Added
    TH2F * jetptSumSumrawpt_fake[nbin];      //Added
    TH2F * jetptSumSumrawpt_real[nbin];      //Added

    TH2F * jetptSumSumpt[nbin];   //Added
    TH2F * jetptSumSumpt_fake[nbin];      //Added
    TH2F * jetptSumSumpt_real[nbin];      //Added

/*    TH2F * trackpt[nbin];
    TH1F * trackrawpt[nbin];
    TH2F * trackphi[nbin];
    TH2F * tracketa[nbin];
    TH2F * trackpteff[nbin];
    

    TH1D * CenBin;
    TH1D * CenBinWt;
    TH1F * Vertex ;
    TH1F * VertexWt ;

    //For inclusive track analysis, outside jet loop
    TH1F * inctrkpt[nbin];
    TH1F * inctrkptM1P1[nbin];
    TH1F * incgenpartpt[nbin];
    TH1F * incgenmatchpt[nbin];



    TH2F * trackPtEtaBin[nbin][netabin]; //jet FF hist. in eta bins
    TH1F * inctrkptEtaBin[nbin][ntrketabin];
    TH1F * NtrkEtaBin[nbin][ntrketabin];
    TH1F * jetptEtaBin[nbin][netabin];
    TH1F * NjetsEtaBin[nbin][netabin];

    TH2F * genpartpt[nbin] ;
    TH2F * genpartphi[nbin] ;
    TH2F * genparteta[nbin] ;
    TH2F * bkgpartpt[nbin];
    
*/

    
};

hist_class::hist_class(TString the_desc)
{
    
    desc = the_desc;
    IsMC = kTRUE;
    for(int ibin = 0 ; ibin <nbin; ibin++){
/*
        NEvents[ibin] = new TH1D(Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), 100, 0, 2.);
        NEvents[ibin]->Sumw2();
        NevtCounter[ibin] = new TH1D(Form("NeventsCounter_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsCounter_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        NevtCounter[ibin]->Sumw2();
        JetAbove40[ibin] = new TH1D(Form("NeventsWithJetAbove40_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsWithJetAbove40_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        JetAbove60[ibin] = new TH1D(Form("NeventsWithJetAbove60_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsWithJetAbove60_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        JetAbove75[ibin] = new TH1D(Form("NeventsWithJetAbove75_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsWithJetAbove75_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        JetAbove95[ibin] = new TH1D(Form("NeventsWithJetAbove95_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsWithJetAbove95_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        JetAbove120[ibin] = new TH1D(Form("NeventsWithJetAbove120_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NeventsWithJetAbove120_%d-%d%%",centr[ibin],centr[ibin+1]), 1, 0, 1);
        Ntrack[ibin] = new TH1F(Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), 800, -1., 799);
*/
        jetpt[ibin] = new TH1F(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500.);
       jetpt[ibin]->Sumw2();

        jetptEta[ibin] = new TH2F(Form("jetptEta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptEta_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 100, -5., 5.);
        jetptEta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetptEta[ibin]->GetYaxis()->SetTitle("#eta");
        jetptEta[ibin]->Sumw2();

	jetptphi[ibin] = new TH2F(Form("jetptphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphi_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 20, -TMath::Pi(), TMath::Pi());  //Added
	jetptphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");	//Added
	jetptphi[ibin]->GetYaxis()->SetTitle("#phi");	//Added
	jetptphi[ibin]->Sumw2();	//Added
	
	jetptchMax[ibin] = new TH2F(Form("jetptchMax_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMax_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);	//Added
	jetptchMax[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");	//Added
	jetptchMax[ibin]->GetYaxis()->SetTitle("chargedMax (GeV/c)");	//Added
	jetptchMax[ibin]->Sumw2();	//Added

        jetptchMax_real[ibin] = new TH2F(Form("jetptchMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptchMax_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptchMax_real[ibin]->GetYaxis()->SetTitle("chargedMax (GeV/c)");     //Added
        jetptchMax_real[ibin]->Sumw2();      //Added

        jetptchMax_fake[ibin] = new TH2F(Form("jetptchMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptchMax_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptchMax_fake[ibin]->GetYaxis()->SetTitle("chargedMax (GeV/c)");     //Added
        jetptchMax_fake[ibin]->Sumw2();      //Added
	
	jetptchMaxpt[ibin] = new TH2F(Form("jetptchMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);    //Added
        jetptchMaxpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptchMaxpt[ibin]->GetYaxis()->SetTitle("chargedMax/p_{T}");   //Added
        jetptchMaxpt[ibin]->Sumw2();      //Added

	jetptchMaxpt_real[ibin] = new TH2F(Form("jetptchMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptchMaxpt_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptchMaxpt_real[ibin]->GetYaxis()->SetTitle("chargedMax/p_{T}");     //Added
        jetptchMaxpt_real[ibin]->Sumw2();      //Added

        jetptchMaxpt_fake[ibin] = new TH2F(Form("jetptchMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptchMaxpt_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptchMaxpt_fake[ibin]->GetYaxis()->SetTitle("chargedMax/p_{T}");     //Added
        jetptchMaxpt_fake[ibin]->Sumw2();      //Added	

//----------------------------------------------------------------------------------------------------------------

	jetptneuMaxpt[ibin] = new TH2F(Form("jetptneuMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);    //Added
        jetptneuMaxpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxpt[ibin]->GetYaxis()->SetTitle("neutralMax/p_{T}");   //Added
        jetptneuMaxpt[ibin]->Sumw2();      //Added

        jetptneuMaxpt_real[ibin] = new TH2F(Form("jetptneuMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptneuMaxpt_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxpt_real[ibin]->GetYaxis()->SetTitle("neutralMax/p_{T}");     //Added
        jetptneuMaxpt_real[ibin]->Sumw2();      //Added

        jetptneuMaxpt_fake[ibin] = new TH2F(Form("jetptneuMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptneuMaxpt_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxpt_fake[ibin]->GetYaxis()->SetTitle("neutralMax/p_{T}");     //Added
        jetptneuMaxpt_fake[ibin]->Sumw2();      //Added  

	jetptneuMax[ibin] = new TH2F(Form("jetptneuMax_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMax_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);    //Added
        jetptneuMax[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMax[ibin]->GetYaxis()->SetTitle("neutralMax");   //Added
        jetptneuMax[ibin]->Sumw2();      //Added

        jetptneuMax_real[ibin] = new TH2F(Form("jetptneuMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptneuMax_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMax_real[ibin]->GetYaxis()->SetTitle("neutralMax");     //Added
        jetptneuMax_real[ibin]->Sumw2();      //Added

        jetptneuMax_fake[ibin] = new TH2F(Form("jetptneuMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptneuMax_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMax_fake[ibin]->GetYaxis()->SetTitle("neutralMax");     //Added
        jetptneuMax_fake[ibin]->Sumw2();      //Added  

	jetptneuMaxr[ibin] = new TH2F(Form("jetptneuMaxr_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxr_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);    //Added
        jetptneuMaxr[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxr[ibin]->GetYaxis()->SetTitle("neutralMax/(neutralMax+chargedMax+photonMax)");     //Added
        jetptneuMaxr[ibin]->Sumw2();      //Added

        jetptneuMaxr_real[ibin] = new TH2F(Form("jetptneuMaxr_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxr_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptneuMaxr_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxr_real[ibin]->GetYaxis()->SetTitle("neutralMax/(neutralMax+chargedMax+photonMax)");     //Added
        jetptneuMaxr_real[ibin]->Sumw2();      //Added

	jetptneuMaxr_fake[ibin] = new TH2F(Form("jetptneuMaxr_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxr_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);     //Added
        jetptneuMaxr_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptneuMaxr_fake[ibin]->GetYaxis()->SetTitle("neutralMax/(neutralMax+chargedMax+photonMax)");     //Added
        jetptneuMaxr_fake[ibin]->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------------

        jetptphoMaxpt[ibin] = new TH2F(Form("jetptphoMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMaxpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);    //Added
        jetptphoMaxpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptphoMaxpt[ibin]->GetYaxis()->SetTitle("photonMax/p_{T}");   //Added
        jetptphoMaxpt[ibin]->Sumw2();      //Added

        jetptphoMaxpt_real[ibin] = new TH2F(Form("jetptphoMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptphoMaxpt_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptphoMaxpt_real[ibin]->GetYaxis()->SetTitle("photonMax/p_{T}");     //Added
        jetptphoMaxpt_real[ibin]->Sumw2();      //Added

        jetptphoMaxpt_fake[ibin] = new TH2F(Form("jetptphoMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 150, 0., 1.5);        //Added
        jetptphoMaxpt_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptphoMaxpt_fake[ibin]->GetYaxis()->SetTitle("photonMax/p_{T}");     //Added
        jetptphoMaxpt_fake[ibin]->Sumw2();      //Added  

        jetptphoMax[ibin] = new TH2F(Form("jetptphoMax_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMax_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);    //Added
        jetptphoMax[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptphoMax[ibin]->GetYaxis()->SetTitle("photonMax");   //Added
        jetptphoMax[ibin]->Sumw2();      //Added

        jetptphoMax_real[ibin] = new TH2F(Form("jetptphoMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptphoMax_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
	jetptphoMax_real[ibin]->GetYaxis()->SetTitle("photonMax");     //Added
        jetptphoMax_real[ibin]->Sumw2();      //Added

        jetptphoMax_fake[ibin] = new TH2F(Form("jetptphoMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., nPtBin_Max, TrkBin_Max);        //Added
        jetptphoMax_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptphoMax_fake[ibin]->GetYaxis()->SetTitle("photonMax");     //Added
        jetptphoMax_fake[ibin]->Sumw2();      //Added  

//---------------------------------------------------------------------------------------------------------

	jetptSumSum[ibin] = new TH2F(Form("jetptSumSum_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSum_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 100, 0., 1000);    //Added
        jetptSumSum[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSum[ibin]->GetYaxis()->SetTitle("chargedSum+photonSum+neutralSum+muSum+eSum");     //Added
        jetptSumSum[ibin]->Sumw2();      //Added

        jetptSumSum_real[ibin] = new TH2F(Form("jetptSumSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 100, 0., 1000);        //Added
        jetptSumSum_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSum_real[ibin]->GetYaxis()->SetTitle("chargedSum+photonSum+neutralSum+muSum+eSum");     //Added
        jetptSumSum_real[ibin]->Sumw2();      //Added

        jetptSumSum_fake[ibin] = new TH2F(Form("jetptSumSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 100, 0., 1000);      //Added
        jetptSumSum_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSum_fake[ibin]->GetYaxis()->SetTitle("chargedSum+photonSum+neutralSum+muSum+eSum");     //Added
        jetptSumSum_fake[ibin]->Sumw2();      //Added

	jetptSumSumrawpt[ibin] = new TH2F(Form("jetptSumSumrawpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumrawpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);    //Added
        jetptSumSumrawpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumrawpt[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt");     //Added
        jetptSumSumrawpt[ibin]->Sumw2();      //Added

        jetptSumSumrawpt_real[ibin] = new TH2F(Form("jetptSumSumrawpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumrawpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);        //Added
        jetptSumSumrawpt_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumrawpt_real[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt");     //Added
        jetptSumSumrawpt_real[ibin]->Sumw2();      //Added

        jetptSumSumrawpt_fake[ibin] = new TH2F(Form("jetptSumSumrawpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumrawpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);      //Added
        jetptSumSumrawpt_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumrawpt_fake[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt");     //Added
        jetptSumSumrawpt_fake[ibin]->Sumw2();      //Added	

        jetptSumSumpt[ibin] = new TH2F(Form("jetptSumSumpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumpt_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);    //Added
        jetptSumSumpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumpt[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/jetpt");     //Added        jetptSumSumrawpt[ibin]->Sumw2();      //Added
       
	jetptSumSumpt_real[ibin] = new TH2F(Form("jetptSumSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);        //Added
        jetptSumSumpt_real[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumpt_real[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/jetpt");     //Added
        jetptSumSumpt_real[ibin]->Sumw2();      //Added

        jetptSumSumpt_fake[ibin] = new TH2F(Form("jetptSumSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 50, 0., 500., 500, 0., 5);      //Added
        jetptSumSumpt_fake[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");  //Added
        jetptSumSumpt_fake[ibin]->GetYaxis()->SetTitle("(chargedSum+photonSum+neutralSum+muSum+eSum)/jetpt");     //Added
        jetptSumSumpt_fake[ibin]->Sumw2();      //Added

/*
     //   trackpt[ibin] = new TH2F(Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 200, 0., 200); 
        trackpt[ibin] = new TH2F(Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, nPtBin, TrkBin); 
        trackpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        trackpt[ibin]->GetYaxis()->SetTitle("p_{T}^{trk} (GeV/c)");   
        trackpt[ibin]->Sumw2();
        trackrawpt[ibin] = new TH1F(Form("trackrawpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackrawpt_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin);
        trackrawpt[ibin]->GetXaxis()->SetTitle("p_{T}^{trkRAW} (GeV/c)");
        trackrawpt[ibin]->Sumw2();
        trackphi[ibin] = new TH2F(Form("trackphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -5.05, 4.95); 
        trackphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        trackphi[ibin]->GetYaxis()->SetTitle("#phi^{trk}");   
        trackphi[ibin]->Sumw2();
        tracketa[ibin] = new TH2F(Form("tracketa_%d-%d%%",centr[ibin],centr[ibin+1]), Form("tracketa_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -5.05, 4.95); 
        tracketa[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        tracketa[ibin]->GetYaxis()->SetTitle("#eta^{trk}");   
        tracketa[ibin]->Sumw2();

        for(Int_t ieta = 0 ; ieta < netabin ; ieta++){
        jetptEtaBin[ibin][ieta] = new TH1F(Form("jetptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), Form("jetptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), 50, 0., 500.);
   //   else   jetpt[ibin] = new TH1F(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin);
        jetptEtaBin[ibin][ieta]->Sumw2();

        NjetsEtaBin[ibin][ieta] = new TH1F(Form("NjetsEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), Form("NjetsEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), 800, -1., 799);
       NjetsEtaBin[ibin][ieta]->Sumw2();

        trackPtEtaBin[ibin][ieta] = new TH2F(Form("trackptInJetEta%.f_%.f_Cen%d-%d%%",  deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), Form("trackptInJetEta%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), 500, 0., 500, nPtBin, TrkBin);
        trackPtEtaBin[ibin][ieta]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        trackPtEtaBin[ibin][ieta]->GetYaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
        trackPtEtaBin[ibin][ieta]->Sumw2();
        }

        trackpteff[ibin] = new TH2F(Form("jettrackpteff_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jettrackpteff_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin, 200, 0., 2.);
        trackpteff[ibin]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
        trackpteff[ibin]->Sumw2();


        //for inclusive track histgram
        inctrkpt[ibin] = new TH1F(Form("inclusivetrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("inclusivetrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin);
        inctrkpt[ibin]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
        inctrkpt[ibin]->Sumw2();
        inctrkptM1P1[ibin] = new TH1F(Form("inclusivetrackptInEtaM1P1_%d-%d%%",centr[ibin],centr[ibin+1]), Form("inclusivetrackptM1P1_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin);
        inctrkptM1P1[ibin]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
        inctrkptM1P1[ibin]->Sumw2();
         for(Int_t ieta = 0 ; ieta < ntrketabin ; ieta++){
         inctrkptEtaBin[ibin][ieta] = new TH1F(Form("IncTrkEtaBin%.f_%.f_Cen%d-%d%%",dtrketa[ieta]*10,dtrketa[ieta+1]*10, centr[ibin],centr[ibin+1]), Form("IncTrkEtaBin%.f_%.f_Cen%d-%d%%",dtrketa[ieta]*10,dtrketa[ieta+1]*10, centr[ibin],centr[ibin+1]), nPtBin, TrkBin);
        inctrkptEtaBin[ibin][ieta]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
        inctrkptEtaBin[ibin][ieta]->Sumw2();
 
        NtrkEtaBin[ibin][ieta] = new TH1F(Form("NtracksEtaBin%.f_%.f_Cen%d-%d%%",dtrketa[ieta]*10,dtrketa[ieta+1]*10, centr[ibin],centr[ibin+1]), Form("NtracksEtaBin%.f_%.f_Cen%d-%d%%",dtrketa[ieta]*10,dtrketa[ieta+1]*10, centr[ibin],centr[ibin+1]), 800, -1., 799);
       NtrkEtaBin[ibin][ieta]->Sumw2();

       }
        if(IsMC){
          if(DoGenAna){
          //  genpartpt[ibin] = new TH2F(Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 200, 0., 200);
            genpartpt[ibin] = new TH2F(Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., nPtBin, TrkBin);
            genpartpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genpartpt[ibin]->GetYaxis()->SetTitle("p_{T}^{part} (GeV/c)");   
            genpartpt[ibin]->Sumw2();
            genpartphi[ibin] = new TH2F(Form("genpartphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genpartphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 100, -5.05, 4.95);
            genpartphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genpartphi[ibin]->GetYaxis()->SetTitle("#phi^{part}");   
            genpartphi[ibin]->Sumw2();
            genparteta[ibin] = new TH2F(Form("genparteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genparteta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 100, -5.05, 4.95);
            genparteta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genparteta[ibin]->GetYaxis()->SetTitle("#eta^{part} ");   
            genparteta[ibin]->Sumw2();
          //  bkgpartpt[ibin] = new TH2F(Form("bkgpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 200, 0., 200); 
          
        //for inclusive track histgram
        incgenpartpt[ibin] = new TH1F(Form("inclusivegenpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("inclusivegenpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin);
        incgenpartpt[ibin]->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
        incgenpartpt[ibin]->Sumw2();
        incgenmatchpt[ibin] = new TH1F(Form("inclusivegenmatchpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("inclusivegenmatchpt_%d-%d%%",centr[ibin],centr[ibin+1]), nPtBin, TrkBin);   
       incgenmatchpt[ibin]->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
        incgenmatchpt[ibin]->Sumw2();

          } //generator level charge particles histogram
        } //MC histogram
*/
        }  //centrality bins loop
        
/*    CenBin = new TH1D((TString) (desc + "_Cent"), "", 100, 0, 100);
    CenBin->Sumw2();
    CenBinWt = new TH1D((TString) (desc + "_CentWeighted"), "", 100, 0, 100);
    CenBinWt->Sumw2();
    Vertex = new TH1F((TString) (desc + "_Vz"), "", 400, -20., 20.);
    Vertex->Sumw2();
    VertexWt = new TH1F((TString) (desc + "_VzWeighted"), "", 400, -20., 20.);
    VertexWt->Sumw2();
  */  
}

void hist_class::Delete()
{
    for(int ibin = 0 ; ibin <nbin; ibin++){
     /* 
	delete NEvents[ibin];
        delete NevtCounter[ibin];
        delete JetAbove40[ibin];
        delete JetAbove60[ibin];
        delete JetAbove75[ibin];
        delete JetAbove95[ibin];
        delete JetAbove120[ibin];
        delete Ntrack[ibin];  	commented */	

        delete jetpt[ibin];
        delete jetptEta[ibin];
 	delete jetptphi[ibin];  //Added
	delete jetptchMax[ibin];	//Added
	delete jetptchMax_fake[ibin];	//Added
	delete jetptchMax_real[ibin];	//Added
delete  jetptchMaxpt[ibin];       //Added
delete  jetptchMaxpt_real[ibin];       //Added
delete  jetptchMaxpt_fake[ibin];       //Added
delete  jetptneuMaxpt[ibin]; //Added
delete  jetptneuMaxpt_real[ibin];    //Added
delete  jetptneuMaxpt_fake[ibin];    //Added

delete    jetptneuMaxr[ibin];   //Added
delete    jetptneuMaxr_fake[ibin];      //Added
delete    jetptneuMaxr_real[ibin];      //Added

delete    jetptneuMax[ibin];   //Added
delete    jetptneuMax_fake[ibin];      //Added
delete    jetptneuMax_real[ibin];      //Added

delete   jetptphoMaxpt[ibin]; //Added 
delete    jetptphoMaxpt_fake[ibin];    //Added
delete    jetptphoMaxpt_real[ibin];    //Added

delete    jetptphoMax[ibin];   //Added
delete    jetptphoMax_fake[ibin];      //Added
delete    jetptphoMax_real[ibin];      //Added

delete    jetptSumSum[ibin];   //Added
delete    jetptSumSum_fake[ibin];      //Added
delete    jetptSumSum_real[ibin];      //Added

delete    jetptSumSumrawpt[ibin];   //Added
delete    jetptSumSumrawpt_fake[ibin];      //Added
delete    jetptSumSumrawpt_real[ibin];      //Added

delete    jetptSumSumpt[ibin];   //Added
delete    jetptSumSumpt_fake[ibin];      //Added
delete    jetptSumSumpt_real[ibin];      //Added
	


/*
        delete trackpt[ibin];
        delete trackrawpt[ibin];
        delete trackphi[ibin];
        delete tracketa[ibin];
       
        delete inctrkpt[ibin] ; 
        delete inctrkptM1P1[ibin] ; 
        delete trackpteff[ibin];
        for(Int_t ieta = 0 ; ieta < netabin ; ieta++){
        delete trackPtEtaBin[ibin][ieta]; 
        delete jetptEtaBin[ibin][ieta]; 
        delete NjetsEtaBin[ibin][ieta]; 
        }
       
        for(Int_t ieta = 0 ; ieta < ntrketabin ; ieta++)  {
           delete inctrkptEtaBin[ibin][ieta] ;
           delete NtrkEtaBin[ibin][ieta];  
            }
        if(IsMC){
          if(DoGenAna){
            delete genpartpt[ibin];
            delete genpartphi[ibin];
            delete genparteta[ibin];
            delete incgenpartpt[ibin];
            delete incgenmatchpt[ibin];
         }
        }	commented*/
        
    } //centrality loop
    //   for(int i=0; i<6;i++) delete ptbin[i];
 /*   delete CenBin;
    delete CenBinWt;
    delete Vertex;
    delete VertexWt;	commented*/

}

void hist_class::Write()
{
    TString out_name ;
    TString dataType ; 
    if(IsMC) dataType="MC" ;
    else dataType="DATA";
   TString anagen ;
    if(DoGenAna) anagen="GenCharge";
    else anagen="";
    if(doTrackCorrections){
          out_name=Form("%s%s_%sAkPu3PFJet%s2013%sIterTrkCorrtest%.fEtaCut%.f_%s",dataType.Data(),coll.Data(),TrigName.Data(), anagen.Data(), corrMet.Data(), trackcut,tracketacut*10, intputFile.Data());
      }
       else {
          out_name=Form("%s%s_%sAkPu3PFJet%sTrk%.fEtaCut%.f_%s",dataType.Data(),coll.Data(),TrigName.Data(), anagen.Data(),  trackcut,tracketacut*10,intputFile.Data());
    }
     
       TFile *out_file = new TFile(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/%s",out_name.Data()), "RECREATE");
    for(int ibin = 0 ; ibin <nbin; ibin++){
/*        NEvents[ibin]->Write();
        NevtCounter[ibin]->Write();
        JetAbove40[ibin]->Write();
        JetAbove60[ibin]->Write();
        JetAbove75[ibin]->Write();
        JetAbove95[ibin]->Write();
        JetAbove120[ibin]->Write();
        Ntrack[ibin]->Write();	commented*/
        jetpt[ibin]->Write();
         jetptEta[ibin]->Write();
	jetptphi[ibin]->Write();
	
jetptchMax[ibin]->Write();
jetptchMax_real[ibin]->Write();
jetptchMax_fake[ibin]->Write();

jetptchMaxpt[ibin]->Write();    //Added
jetptchMaxpt_real[ibin]->Write();       //Added
jetptchMaxpt_fake[ibin]->Write();       //Added

jetptneuMaxpt[ibin]->Write(); //Added
 jetptneuMaxpt_real[ibin]->Write();    //Added
 jetptneuMaxpt_fake[ibin]->Write();    //Added

 jetptneuMaxr[ibin]->Write();   //Added
 jetptneuMaxr_fake[ibin]->Write();      //Added
 jetptneuMaxr_real[ibin]->Write();      //Added

 jetptneuMax[ibin]->Write();   //Added
 jetptneuMax_fake[ibin]->Write();      //Added
 jetptneuMax_real[ibin]->Write();      //Added

 jetptphoMaxpt[ibin]->Write(); //Added 
 jetptphoMaxpt_fake[ibin]->Write();    //Added
 jetptphoMaxpt_real[ibin]->Write();    //Added

 jetptphoMax[ibin]->Write();   //Added
 jetptphoMax_fake[ibin]->Write();      //Added
 jetptphoMax_real[ibin]->Write();      //Added

 jetptSumSum[ibin]->Write();   //Added
 jetptSumSum_fake[ibin]->Write();      //Added
 jetptSumSum_real[ibin]->Write();      //Added

 jetptSumSumrawpt[ibin]->Write();   //Added
 jetptSumSumrawpt_fake[ibin]->Write();      //Added
 jetptSumSumrawpt_real[ibin]->Write();      //Added

 jetptSumSumpt[ibin]->Write();   //Added
 jetptSumSumpt_fake[ibin]->Write();      //Added
 jetptSumSumpt_real[ibin]->Write();      //Added

/*         trackpt[ibin]->Write();
         trackrawpt[ibin]->Write();
        trackphi[ibin]->Write();
        tracketa[ibin]->Write();
        inctrkpt[ibin]->Write();
        inctrkptM1P1[ibin]->Write();
        trackpteff[ibin]->Write();
        for(Int_t ieta = 0 ; ieta < netabin ; ieta++){
         trackPtEtaBin[ibin][ieta]->Write();
         jetptEtaBin[ibin][ieta]->Write();
         NjetsEtaBin[ibin][ieta]->Write();
        }
       for(Int_t ieta = 0 ; ieta < ntrketabin ; ieta++)  {
          inctrkptEtaBin[ibin][ieta]->Write();
          NtrkEtaBin[ibin][ieta]->Write();
        }
        if(IsMC){
         if(DoGenAna){ 
           genpartpt[ibin]->Write();
            genpartphi[ibin]->Write();
            genparteta[ibin]->Write();
            incgenpartpt[ibin]->Write();
            incgenmatchpt[ibin]->Write();
          }    
        }	commented*/
        } //centrality bins
    
/*    CenBin->Write();
    CenBinWt->Write();
    Vertex->Write();
    VertexWt->Write();	commented*/

    out_file->Close();
    cout <<"Output file: " <<Form("%s",out_name.Data()) <<endl ;
    
}


void anaAsyJetTrig()
{
   //for ak3PF jets, applied for pPb data only 
 //   TFile* fcrel3 = TFile::Open("/afs/cern.ch/user/d/dgulhan/public/Corrections/Casym_pPb_double_hcalbins_algo_ak3PF_pt100_140_jet80_alphahigh_20_phicut250.root", "readonly");
    //for akPu3PF jets and pPb data 
    //TFile* fcrel3 = TFile::Open("/afs/cern.ch/user/d/dgulhan/public/Corrections/Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root", "readonly");
    TFile* fcrel3 = TFile::Open("Casym_pPb_double_hcalbins_algo_akPu3PF_pt100_140_jet80_alphahigh_20_phicut250.root", "readonly");
    TH1D * C_rel=(TH1D*)fcrel3->Get("C_asym");
  //  fcrel3->Close();

    hist_class *my_hists = new hist_class("pfjet");
    cout <<my_hists->IsMC<<endl ;

  //this function is to correct for UE subtraction (we are using akPu3PF algorithm)
   TF1 * fUE ;
   TF1 * fgaus ;
   if(my_hists->IsMC==kTRUE){
    fUE = new TF1("fUE","[0]/pow(x,[1])");
  //   //for ak3PF jets
  //   fUE->SetParameters(0.4183,0.4244);
   //  //for akPu3PF jets
     fUE->SetParameters(1.052,0.5261);
    //for gauss smearing
    fgaus=new TF1("fgaus","gaus(0)",-20,20);
    fgaus->SetParameters(1,0,1);
    }
   else {
       fUE = new TF1("fUE","1-[0]/pow(x,[1])",20,300);
        //   //for ak3PF jets
  //   fUE->SetParameters(0.8648,0.8167);
   //  //for akPu3PF jets
     fUE->SetParameters(0.3015,0.8913);
    }

   TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
    fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);

    if(my_hists->IsMC==kTRUE){
      pthat=atoi(getenv("PTHAT")) ;
      ptmax=atoi(getenv("PTMAX")) ;
       cout <<"pthat = " <<pthat <<"  pthatmax =" <<ptmax <<endl ;
     }

    

    std::cout << "start working\n";
    if(my_hists->IsMC!=kTRUE){ 
      if(coll=="HI")
        dataPath="/net/hidsk0001/d00/scratch/yjlee/merge/pbpbDijet_v20" ;//mit PbPb data path
    else {
      if(TrigName=="Jet20")   
          dataPath="root://eoscms//eos/cms/store/group/phys_heavyions/krajczar/inbound/mnt/hadoop/cms/store/user/krajczar/pPb_Jet20_Full_v1" ;
   else if(TrigName=="Jet40" || TrigName=="Jet60")
          dataPath="root://eoscms//eos/cms/store/group/phys_heavyions/krajczar/inbound/mnt/hadoop/cms/store/user/krajczar/pPb_Jet40Jet60_Full_v1" ;
    else  
          dataPath="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/pPb2013/promptReco";
    }
  } //data Path
    else { //MC analysis
        if(coll=="HI") {
          if(pthat==50||pthat==80||pthat==100||pthat==170)
             dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/"); //MIT MC normial
           else 
               dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/"); //MIT MC normial
       }
       else if(coll=="PPb")
      //  dataPath=Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%d/HiForest_v77_merged01", pthat);
      //  dataPath=Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt%d/HiForest_v77_v2_merged01", pthat);
	dataPath=Form("/store/user/ymao/HiForest/MC/PYTHIA_HIJING/");	//ACCRE PATH
       else if(coll=="PbP")
        dataPath=Form("/mnt/hadoop/cms/store/user/dgulhan/Pbp/HP05/prod24/Hijing_Pythia_pt%d/HiForest_v84_merged02", pthat);
        else if(coll=="PP2011")
         dataPath= Form("/net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat%d", pthat); //lxplus path for pp
       else {       
       if(pthat==220)
        // dataPath= Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt%d/HiForest_v77_v2_merged02", pthat); //2013 pp tracking for 5.02TeV
	dataPath=Form("/store/user/ymao/HiForest/MC/PYTHIA_Z2_5TeV/");	//ACCRE PATH // ! 2013 pp 5.02TeV MC with pp tracking
      else 
     //  dataPath= Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Signal_Pythia_pt%d/HiForest_v77_v2_merged01", pthat); //2013 pp tracking for 5.02TeV
	dataPath=Form("/store/user/ymao/HiForest/MC/PYTHIA_Z2_5TeV/");	//ACCRE PATH // ! 2013 pp 5.02TeV MC with pp tracking
        }
    }
    if(my_hists->IsMC!=kTRUE){  //real data analysis
        if(coll=="HI")             
            intputFile="promptskim-hihighpt-hltjet80-pt90-v20.root" ; //full dataset
	else if(coll=="PP2011")  
             intputFile="hiForest2_pp_ppreco_415_90percent.root";  //! 2011 pp data rereco
        else if(coll=="PbPb")
           intputFile="PbPHiForest2_PbPbPAHighPtJet80_cent50-100_pprereco.root"; 	
        else if(coll=="PPb"){
         if(TrigName=="Jet20")
            intputFile="mergedJet20_KK.root" ;
         else    if(TrigName=="Jet40" || TrigName=="Jet60")
            intputFile="mergedJet40Jet60_KK.root" ;
         else    if(TrigName=="Jet80" || TrigName=="Jet100")
            intputFile="PA2013_HiForest_PromptReco_JSonPPb_forestv77.root" ;
         else
            intputFile="PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root" ;
        }
        else if(coll=="PbP")
         intputFile="PA2013_HiForest_PromptReco_JSonPbp_JECdb_forestv84.root" ;
        else 
         //  intputFile="PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root";  //! 2013 pp data with HI tracking
           intputFile="PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root";  //! 2013 pp data with pp tracking
    }
    else { //MC sample
        if(coll=="HI"){             
          if(pthat==50||pthat==80||pthat==100||pthat==170)
              intputFile=Form("Dijet%d_HydjetDrum_v27_mergedV1.root", pthat);
         else 
              intputFile=Form("Dijet%d_HydjetDrum_v28_mergedV1.root", pthat);
   }
      else if(coll=="PPb"){
      //  intputFile=Form("pt%d_HP04_prod16_v77_merged_forest_0.root", pthat);
      //  intputFile=Form("pt%d_HP04_hiforest77_hiSignal.root", pthat);
      	if(pthat==280)
		intputFile=Form("pt%d_HP04_prod25_v85_merged_forest_0_v2.root", pthat);	//ACCRE root files
	else
		intputFile=Form("pt%d_HP04_prod25_v85_merged_forest_0.root", pthat);	//ACCRE root files
	}
      else if(coll=="PbP")
        intputFile=Form("pt%d_HP05_prod24_v84_merged_forest_0.root", pthat);
      else if(coll=="PP2011")
          intputFile=Form("mergedFile.root");  // 2011 pp MC
        else 
         intputFile=Form("pt%d_HP04_hiforest77_hiSignal.root", pthat); // ! 2013 pp 5.02TeV MC with pp tracking 
    }
    
    TString inname=Form("%s/%s", dataPath.Data(),intputFile.Data());
    // Define the input file and HiForest
   HiForest * c ; 
   if(coll=="PP" || coll=="PP2011")
        c = new HiForest(inname,"forest",cPP);
     else if(coll=="PPb" || coll=="PbP")
        c = new HiForest(inname,"forest",cPPb);
     else 
        c = new HiForest(inname,"forest",cPbPb);
       c->doTrackCorrections = false;
    c->doTrackingSeparateLeadingSubleading = false;
    c->InitTree();
  //  cout << "start working222222\n";
  //  TFile *my_file=TFile::Open(Form("%s/%s", dataPath.Data(),intputFile.Data()));
    cout <<"Input file" << inname<<endl ;
    
    TrackCorrector corr("trackCorrections_HIN12017v4_HijingCombined.root");
   if(doTrackCorrections){
   corr.load("trkCorr_HIN12017");
   corr.setOption1(true);
   corr.setOption2(true);
  }

    Evts * offSel = &(c->evt); 

    Skims * my_skim = &(c->skim); 

    Hlts * trigSel = &(c->hlt); 
    //jet tree
//    if(coll=="HI") 
        Jets * my_ct = &(c->akPu3PF); 
 //   else 
 //   Jets * my_ct = &(c->ak3PF);
    //  Jets * jetthres =  &(c->icPu5); 
   //track tree
    Tracks * my_tr = &(c->track);
    
    //GenParticle tree
    GenParticles * my_GenPart = &(c->genparticle);

    int curr_bin = nbin-1 ;
    int Nevents = 0 ;
	int NFill_0_500=0;	int NFill_0_500_fake=0;	int NFill_0_500_real=0;	//Added
	int NFill_30_300=0;	int NFill_30_300_fake=0;	int NFill_30_300_real=0;	//Added
      Int_t Nevt_40_60 = 0 ;
      Int_t Nevt_60_75 = 0 ;
      Int_t Nevt_75_95 = 0 ;
      Int_t Nevt_95_120 = 0 ;
      Int_t Nevt_120 = 0 ;
    cout <<"Number of events ="<<c->GetEntries()<<endl ;

    for(int evi = 0; evi < c->GetEntries(); evi++) {
        c->GetEntry(evi);
        int noise_evt = my_skim->pHBHENoiseFilter ;
        //        int ecal_noise = my_skim->phiEcalRecHitSpikeFilter ;
        //        if(ecal_noise==0) continue ;
        
        double vz = offSel->vz ;
        int hiBin = offSel->hiBin ;
        weight = 1. ;      
 
         int pileup_Gplus ;
       if(my_hists->IsMC!=kTRUE){
            int evt_sel ;
            pileup_Gplus = my_skim->pVertexFilterCutGplus ;
           if(coll=="PbPb"|| coll=="HI"|| coll=="PP2011")  evt_sel = my_skim->pcollisionEventSelection ;
            else  evt_sel = my_skim->pPAcollisionEventSelectionPA;
         //   if(evt_sel==0) continue ;
            if(!evt_sel) continue ;
        }
        if(my_hists->IsMC!=kTRUE){
          //  if(noise_evt==0) continue ;
            if(!noise_evt) continue ;
             int jetTr2 ;
              if(coll=="HI"|| coll=="PbPb")
               jetTr2 = trigSel->HLT_HIJet80_v1 ;
                else if (coll=="PP2011")
                 jetTr2 = trigSel->HLT_Jet60_v1 ;
               else {
               //  jetTr2 = trigSel->HLT_PAJet80_NoJetID_v1 ;
                if(TrigName=="Jet20") jetTr2 = trigSel->HLT_PAJet20_NoJetID_v1  ;
                else  if(TrigName=="Jet40") jetTr2 = trigSel->HLT_PAJet40_NoJetID_v1  ;
                else  if(TrigName=="Jet60") jetTr2 = trigSel->HLT_PAJet60_NoJetID_v1  ;
                else  if(TrigName=="Jet80") jetTr2 = trigSel->HLT_PAJet80_NoJetID_v1  ;
                else  if(TrigName=="Jet100") jetTr2 = trigSel->HLT_PAJet100_NoJetID_v1  ;
                else  jetTr2 = trigSel->HLT_PAZeroBiasPixel_SingleTrack_v1 ;
               }
            if(!jetTr2) continue ;
            int run=offSel->run ;
             if( !(my_skim->phfPosFilter1 && my_skim->phfNegFilter1 && my_skim->pBeamScrapingFilter && my_skim->pprimaryvertexFilter)) continue ;
            if(coll=="PPb"){
            // if( my_skim->phfPosFilter1==0 || my_skim->phfNegFilter1==0 ||my_skim->pBeamScrapingFilter==0 || my_skim->pprimaryvertexFilter==0) continue ;
             if(!pileup_Gplus) continue ;
             if(run>211256) continue ;
             if(run<210676) continue ;  //remove the runs with old alignment
            }
           if(coll=="PbP"){
             if(pileup_Gplus==0) continue ;
            if(run<=211256) continue ;
           }
        }

        if(TMath::Abs(vz)>15.) continue ;
        if(my_hists->IsMC==kTRUE) weight*=fVz->Eval(vz);
     //run selection
    if(my_hists->IsMC!=kTRUE && coll=="PPb") {
       if(offSel->run<210676 ||offSel->run>211256) //211256: last pPb run (Pb goes to +eta)
         continue;
     }
      bool event_accepted = true;
/*      if(!(my_skim->pPAcollisionEventSelectionPA && my_skim->phfPosFilter1 && my_skim->phfNegFilter1
          && my_skim->pBeamScrapingFilter
          && my_skim->pprimaryvertexFilter
          && my_skim->pVertexFilterCutGplus
          && TMath::Abs(offSel->vz)<15.
          )
        ) event_accepted = false;
*/
      if(event_accepted == false)
         continue;

                
    /*    //if there is no jets or no PF candidates, skip the event
        if(my_ct->nref==0) continue ;
   */     //put the higher pthat cut
        if(my_hists->IsMC==kTRUE && my_ct->pthat>ptmax) continue ;
         if(my_ct->pthat>ptmax) cout <<"pthat =" <<my_ct->pthat <<endl ;
//       if(coll=="HI"|| coll=="PP2011") //commented
//        my_hists->CenBin->Fill(hiBin*2.5); //commented
//      else my_hists->CenBin->Fill(hiBin);  //commented

//        my_hists->Vertex->Fill(vz);	//commented

        //   cout <<"vz =" <<vz <<endl ;
        
         if(coll=="HI"){
            double centrality = hiBin*2.5 ;
            //   my_hists->CenBin->Fill(offSel->hiBin);
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(centrality >=centr[ibin] && centrality<centr[ibin+1]) curr_bin = ibin ;
            }
        }
        else {
            curr_bin=nbin-1 ;
         //    weight = 1. ;
        }
       //   weight = 1. ;

       //    cout << "  cent_bin:" <<curr_bin <<endl ;
        if(evi%10000==1)cout <<" coll = " <<coll <<" weight = " <<weight <<" evt = " <<evi <<endl ;
        
        //cout << "start working222222\n";
        
        
        //  cout << "still working222222\n";
//        if(my_hists->IsMC==kFALSE)my_hists->VertexWt->Fill(vz+0.4847, weight);	//commented
//        else  my_hists->VertexWt->Fill(vz, weight);	//commented
//         if(coll=="HI"|| coll=="PP2011")     
//           my_hists->CenBinWt->Fill(offSel->hiBin*2.5,weight);	//commented
//        else
//         my_hists->CenBinWt->Fill(offSel->hiBin,weight);	//commented

      //Jets for event weights; eta-pt cut removed
      int trackMult = 0;
      double leadingJet = -999.;
  for(int j4i = 0; j4i < my_ct->nref ; j4i++) {
    if( my_hists->IsMC==kTRUE && my_ct->subid[j4i] != 0) continue;
    if (my_ct->rawpt[j4i]<15) continue;
	    if(TMath::Abs(my_ct->jteta[j4i])>3.) continue ;
		 trackMult++;
        if (my_ct->jtpt[j4i]>leadingJet) {
            leadingJet = my_ct->jtpt[j4i];
        }
 }
       if(leadingJet<0.) continue ;
       
      // Don't analyze 0 multiplicity events; correction added later
      if(trackMult==0)
         continue;

//       my_hists->NEvents[curr_bin]->Fill(1, weight);	//commented
 
     //Tracks for event classification: same as for analysis
      bool jetAbove = false ;
      bool  jetAbove40 = false;
      bool  jetAbove60 = false;
      bool  jetAbove75 = false;
      bool  jetAbove95 = false;
      bool  jetAbove120 = false;
      if(leadingJet>40. ){ 
         jetAbove40 = true ;
       }
      if(leadingJet>60. ){                   
         jetAbove60 = true ;
       }
      if(leadingJet>75. ){                   
         jetAbove75 = true ;
       }
       if(leadingJet>95. ){
         jetAbove95 = true ;
       }
       if(leadingJet>120. ){
         jetAbove120 = true ;
       }
      if(!jetAbove40)                 Nevents++ ;
      if(jetAbove40 && !jetAbove60) Nevt_40_60++;
      if(jetAbove60 && !jetAbove75) Nevt_60_75++;
      if(jetAbove75 && !jetAbove95) Nevt_75_95++;
      if(jetAbove95 && !jetAbove120) Nevt_95_120++;
      if(jetAbove120)                  Nevt_120++;
    
      if(TrigName=="Jet20") jetAbove = jetAbove40 && !jetAbove60 ;
      else if(TrigName=="Jet40") jetAbove = jetAbove60 && !jetAbove75 ;
      else   if(TrigName=="Jet60") jetAbove = jetAbove75 && !jetAbove95  ;
      else   if(TrigName=="Jet80") jetAbove = jetAbove95 && !jetAbove120  ;
      else   if(TrigName=="Jet100") jetAbove = jetAbove120  ;
      else jetAbove = !jetAbove40 ;

   //   if(!jetAbove) continue ;
    
   // for inclusive jet analysis
   for(int j4i = 0; j4i < my_ct->nref ; j4i++) {
    double jetweight = 1; 
    double jet_pt= my_ct->jtpt[j4i];
    double jet_eta = my_ct->jteta[j4i];
    double jet_phi = my_ct->jtphi[j4i];	//Added
    double rawpt = my_ct->rawpt[j4i];
    double chargedMax = my_ct->chargedMax[j4i];	//Added
    double neutralMax = my_ct->neutralMax[j4i];
    double photonMax = my_ct->photonMax[j4i];
    double chargedSum = my_ct->chargedSum[j4i];
    double photonSum = my_ct->photonSum[j4i];
    double neutralSum = my_ct->neutralSum[j4i];
    double muSum = my_ct->muSum[j4i];
    double eSum = my_ct->eSum[j4i];
    double refpt = my_ct->refpt[j4i];	//Added
    if (my_ct->rawpt[j4i]<15) continue;
    int dEtaBin = -1. ;
    if( my_hists->IsMC==kTRUE && my_ct->subid[j4i] != 0) continue;	//commented
     //for jet kinematcis cuts
    // if(jet_pt>80. && TMath::Abs(jet_eta)<2.){
     if(TMath::Abs(jet_eta)<=3.){
       if( my_hists->IsMC!=kTRUE ) jetweight*=(fUE->Eval(jet_pt))*C_rel->GetBinContent(C_rel->FindBin(jet_eta));
       else
          jetweight*=((fUE->Eval(jet_pt))*fgaus->GetRandom()+1);
          my_hists->jetptEta[curr_bin]->Fill(jet_pt*jetweight, jet_eta, weight);
	  my_hists->jetptphi[curr_bin]->Fill(jet_pt*jetweight, jet_phi, weight);	//Added

	  my_hists->jetptchMax[curr_bin]->Fill(jet_pt*jetweight, chargedMax, weight);	//Added
          my_hists->jetptneuMax[curr_bin]->Fill(jet_pt*jetweight, neutralMax, weight);
          my_hists->jetptphoMax[curr_bin]->Fill(jet_pt*jetweight, photonMax, weight);

          my_hists->jetptchMaxpt[curr_bin]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);
          my_hists->jetptneuMaxpt[curr_bin]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);
          my_hists->jetptphoMaxpt[curr_bin]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);

	  my_hists->jetptneuMaxr[curr_bin]->Fill(jet_pt*jetweight, neutralMax/(neutralMax+chargedMax+photonMax),weight);
	 
	  my_hists->jetptSumSum[curr_bin]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
	  my_hists->jetptSumSumpt[curr_bin]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	  my_hists->jetptSumSumrawpt[curr_bin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt,weight);
		if(jet_pt*jetweight>0 && jet_pt*jetweight<500){	NFill_0_500++;}		//Added
		if(jet_pt*jetweight>30 && jet_pt*jetweight<300){ NFill_30_300++;}		//Added

	if(refpt<0){	//Added
	  my_hists->jetptchMax_fake[curr_bin]->Fill(jet_pt*jetweight, chargedMax, weight);	//Added
	  my_hists->jetptneuMax_fake[curr_bin]->Fill(jet_pt*jetweight, neutralMax, weight);	//Added
	  my_hists->jetptphoMax_fake[curr_bin]->Fill(jet_pt*jetweight, photonMax, weight);	//Added
	  my_hists->jetptchMaxpt_fake[curr_bin]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);	//Added
	  my_hists->jetptneuMaxpt_fake[curr_bin]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);	//Added
	  my_hists->jetptphoMaxpt_fake[curr_bin]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);	//Added
	  my_hists->jetptneuMaxr_fake[curr_bin]->Fill(jet_pt*jetweight, neutralMax/(neutralMax+chargedMax+photonMax),weight);
          
          my_hists->jetptSumSum_fake[curr_bin]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_fake[curr_bin]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_fake[curr_bin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt,weight);
		if(jet_pt*jetweight>0 && jet_pt*jetweight<500){	NFill_0_500_fake++;}		//Added
		if(jet_pt*jetweight>30 && jet_pt*jetweight<300){ NFill_30_300_fake++;}	
	}	//Added
	else{	//Added
	  my_hists->jetptchMax_real[curr_bin]->Fill(jet_pt*jetweight, chargedMax, weight);	//Added
 	  my_hists->jetptneuMax_real[curr_bin]->Fill(jet_pt*jetweight, neutralMax, weight);     //Added 
          my_hists->jetptphoMax_real[curr_bin]->Fill(jet_pt*jetweight, photonMax, weight);      //Added 
          my_hists->jetptchMaxpt_real[curr_bin]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);     //Added
          my_hists->jetptneuMaxpt_real[curr_bin]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);    //Added
          my_hists->jetptphoMaxpt_real[curr_bin]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);     //Added
          my_hists->jetptneuMaxr_real[curr_bin]->Fill(jet_pt*jetweight, neutralMax/(neutralMax+chargedMax+photonMax),weight);

          my_hists->jetptSumSum_real[curr_bin]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_real[curr_bin]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_real[curr_bin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/rawpt,weight);
		if(jet_pt*jetweight>0 && jet_pt*jetweight<500){	NFill_0_500_real++;}		//Added
		if(jet_pt*jetweight>30 && jet_pt*jetweight<300){ NFill_30_300_real++;}	
	}	//Added
       if(TMath::Abs(jet_eta+0.465)<=1.) my_hists->jetpt[curr_bin]->Fill(jet_pt*jetweight, weight);
         for(Int_t ieta = 0 ; ieta <netabin; ieta++){
		if(my_hists->IsMC!=kTRUE ){
            if(coll=="PPb"){
               if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
             }
            else if(coll=="PbP"){
              if((jet_eta-0.465)>deta[ieta]&&(jet_eta-0.465)<=deta[ieta+1]) dEtaBin = ieta ;
             }
           else {
             if((jet_eta)>deta[ieta]&&(jet_eta)<=deta[ieta+1]) dEtaBin = ieta ;
            } 
		}
		else 
		if((jet_eta)>deta[ieta]&&(jet_eta)<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
      
//  my_hists->jetptEtaBin[curr_bin][dEtaBin]->Fill(jet_pt*jetweight, weight);   //commented
//        my_hists->NjetsEtaBin[curr_bin][dEtaBin]->Fill(1);	//commented
     } 
    }// for jet kinematics cuts
   } //! jet loop

        //Leading Jets seach, for tracking efficiency
      double leadingJetPt = -1. ;
      Int_t leadingJetIndex = -1 ;
     for(int j = 0; j < my_ct->nref ; j++) {
         if (fabs(my_ct->jteta[j])>2.5) continue;
         if (my_ct->rawpt[j]<15) continue;
           if( my_hists->IsMC==kTRUE && my_ct->subid[j] != 0) continue;
         if (my_ct->jtpt[j]>leadingJetPt) {
            leadingJetPt = my_ct->jtpt[j];
            leadingJetIndex = j;
         }
      }
      if(leadingJetPt==-1) //for corrections
            leadingJetPt=10.;


     //for inclusive track analysis, without jet selection and requirement
       for(int itr = 0 ; itr < my_tr->nTrk ; itr++){
                            double tr_pt = my_tr->trkPt[itr];
                            double tr_phi = my_tr->trkPhi[itr];
                            double tr_eta = my_tr->trkEta[itr];
                            if(TMath::Abs(tr_eta)>2.4) continue ;
         //                    if(my_tr->trkPtError[itr]/my_tr->trkPt[itr]>=0.1 || TMath::Abs(my_tr->trkDz1[itr]/my_tr->trkDzError1[itr])>=3.0 ||TMath::Abs(my_tr->trkDxy1[itr]/my_tr->trkDxyError1[itr])>=3.0) continue ; //ridge cut for tracks                           
         if(!((my_tr->highPurity[itr])
             && (fabs(my_tr->trkDz1[itr]/my_tr->trkDzError1[itr])<3)
             && (fabs(my_tr->trkDxy1[itr]/my_tr->trkDxyError1[itr])<3)
             && (my_tr->trkPtError[itr]/my_tr->trkPt[itr]<0.1)
            ))
            continue;
                            Int_t TrkEtaBin = -1 ;
                       for(Int_t ieta = 0 ; ieta <ntrketabin; ieta++){
			if(my_hists->IsMC!=kTRUE){
                        if(coll=="PPb"){             
                           if((tr_eta+0.465)>dtrketa[ieta]&&(tr_eta+0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else if(coll=="PbP"){
                            if((tr_eta-0.465)>dtrketa[ieta]&&(tr_eta-0.465)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
                           }
                        else  
                      if((tr_eta)>dtrketa[ieta]&&(tr_eta)<=dtrketa[ieta+1]) TrkEtaBin = ieta ; 
                     } 
			else
			if((tr_eta)>dtrketa[ieta]&&(tr_eta)<=dtrketa[ieta+1]) TrkEtaBin = ieta ;
			}
			
                            if((my_tr->highPurity[itr])){
                                if(tr_pt<trackcut) continue ;
                            double trkweight=1. ;
                            if(doTrackCorrections){
                             //   if(corrMet=="Hist")trkweight = c->getTrackCorrection(itr);
                             //   else trkweight = c->getTrackCorrectionPara(itr);
                             trkweight = corr.getWeight(tr_pt,tr_eta,leadingJetPt);
                            }
 //                       my_hists->inctrkpt[curr_bin]->Fill(tr_pt, weight*trkweight); //commented
 //                       if(TMath::Abs(tr_eta+0.465)<=1.)	//commented
//                          my_hists->inctrkptM1P1[curr_bin]->Fill(tr_pt, weight*trkweight);	//commented
  //                      my_hists->Ntrack[curr_bin]->Fill(1);	//commented
//                        if(TrkEtaBin!=-1) //my_hists->inctrkptEtaBin[curr_bin][TrkEtaBin]->Fill(tr_pt, weight*trkweight);//commented
//                        if(TrkEtaBin!=-1) my_hists->NtrkEtaBin[curr_bin][TrkEtaBin]->Fill(1);	//commented
                      } //! high purity track cuts
              }  //! inclusive track loop
            if(my_hists->IsMC==kTRUE&&DoGenAna){
                          //using the sim track to calculate the tracking efficiency 
                    for(int ipart = 0 ; ipart < my_tr->nParticle ; ipart++){ //sim track loop 
                        double gen_pt = my_tr->pPt[ipart];
                        double gen_phi = my_tr->pPhi[ipart];
                        double gen_eta = my_tr->pEta[ipart];
                        if(gen_pt<trackcut)continue ;
                        if(TMath::Abs(gen_eta)>2.4)continue ;
                              if(my_tr->pNRec[ipart]>0&&((my_tr->mtrkQual[ipart]))) {
                                 if((my_tr->mtrkPtError[ipart]/my_tr->mtrkPt[ipart]<0.1 && TMath::Abs(my_tr->mtrkDz1[ipart]/my_tr->mtrkDzError1[ipart])<3.0 && TMath::Abs(my_tr->mtrkDxy1[ipart]/my_tr->mtrkDxyError1[ipart])<3.0)){
//                                my_hists->incgenmatchpt[curr_bin]->Fill(gen_pt, weight);	//commented
                                 } //tracking cut
                              } // matching hist
//                             my_hists->incgenpartpt[curr_bin]->Fill(gen_pt,weight);	//commented
                       } //! sim track loop
      } //only runs on MC
    }  ///event loop
    
/*       my_hists->NevtCounter[curr_bin]->SetBinContent(1, Nevents);
       my_hists->JetAbove40[curr_bin]->SetBinContent(1, Nevt_40_60);
       my_hists->JetAbove60[curr_bin]->SetBinContent(1, Nevt_60_75);
       my_hists->JetAbove75[curr_bin]->SetBinContent(1, Nevt_75_95);
       my_hists->JetAbove95[curr_bin]->SetBinContent(1, Nevt_95_120);
       my_hists->JetAbove120[curr_bin]->SetBinContent(1, Nevt_120);	commented */
    my_hists->Write();
    //   my_hists->Delete();
    //  delete my_hists;
    std::cout << "working done\n";
std::cout<<"NFill_0_500:"<<NFill_0_500<<'\t'<<"NFill_30_300:"<<NFill_30_300<<endl;
std::cout<<"NFill_0_500_fake:"<<NFill_0_500_fake<<'\t'<<"NFill_30_300_fake:"<<NFill_30_300_fake<<endl;
std::cout<<"NFill_0_500_real:"<<NFill_0_500_real<<'\t'<<"NFill_30_300_real:"<<NFill_30_300_real<<endl;
}




