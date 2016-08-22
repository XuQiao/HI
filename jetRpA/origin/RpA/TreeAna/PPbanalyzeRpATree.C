// this code is modified to use data only
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TCut.h"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/Correction.C"
using namespace std;
const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const int netabin = sizeof(deta)/sizeof(double)-1;

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PPb";
double weight;

class hist_class{
public:
    hist_class();
    void Delete();
    void Write();

    TH1F * jetpt[5];
    TH2F * jetptEta[5];
    TH2F * jetptphi[5];
    TH2F * jetptchMaxpt[5];      //Added
    TH2F * jetptchMax[5];    //Added
    TH2F * jetptneuMaxpt[5];     //Added
    TH2F * jetptneuMaxr[5];   //Added
    TH2F * jetptneuMax[5];       //Added
    TH2F * jetptphoMaxpt[5]; //Added 
    TH2F * jetptphoMax[5];   //Added
    TH2F * jetptchSum[5];    //Added
    TH2F * jetptneuSum[5];    //Added
    TH2F * jetptphoSum[5];    //Added
    TH2F * jetptchSumpt[5];    //Added
    TH2F * jetptneuSumpt[5];    //Added
    TH2F * jetptphoSumpt[5];    //Added
    TH2F* jetptchN[5];       //Added
    TH2F* jetptneuN[5];      //Added
    TH2F* jetptphoN[5];      //Added
    TH2F * jetptSumSum[5];   //Added
    TH2F * jetptSumSumrawpt[5];   //Added
    TH2F * jetptSumSumpt[5];   //Added

    TH2F * jetptchMaxpt_real[nbin];	//Added
    TH2F * jetptchMaxpt_fake[nbin];	//Added
    TH2F * jetptchMax_real[nbin];	//Added
    TH2F * jetptchMax_fake[nbin];	//Added
    TH2F * jetptneuMaxpt_real[nbin];	//Added
    TH2F * jetptneuMaxpt_fake[nbin];	//Added
    TH2F * jetptneuMaxr_fake[nbin];      //Added
    TH2F * jetptneuMaxr_real[nbin];      //Added
    TH2F * jetptneuMax_fake[nbin];	//Added
    TH2F * jetptneuMax_real[nbin];	//Added
    TH2F * jetptphoMaxpt_fake[nbin];    //Added
    TH2F * jetptphoMaxpt_real[nbin];    //Added
    TH2F * jetptphoMax_fake[nbin];      //Added
    TH2F * jetptphoMax_real[nbin];      //Added
    TH2F * jetptchSum_fake[nbin];	//Added
    TH2F * jetptchSum_real[nbin];	//Added
    TH2F * jetptneuSum_fake[nbin];       //Added
    TH2F * jetptneuSum_real[nbin];       //Added
    TH2F * jetptphoSum_fake[nbin];       //Added
    TH2F * jetptphoSum_real[nbin];       //Added
    TH2F * jetptchSumpt_fake[nbin];       //Added
    TH2F * jetptchSumpt_real[nbin];       //Added
    TH2F * jetptneuSumpt_fake[nbin];       //Added
    TH2F * jetptneuSumpt_real[nbin];       //Added
    TH2F * jetptphoSumpt_fake[nbin];       //Added
    TH2F * jetptphoSumpt_real[nbin];       //Added
    TH2F* jetptchN_fake[nbin];	
    TH2F* jetptchN_real[nbin];
    TH2F* jetptneuN_fake[nbin];
    TH2F* jetptneuN_real[nbin];
    TH2F* jetptphoN_fake[nbin];
    TH2F* jetptphoN_real[nbin];
    TH2F * jetptSumSum_fake[nbin];	//Added
    TH2F * jetptSumSum_real[nbin];	//Added
    TH2F * jetptSumSumrawpt_fake[nbin];      //Added
    TH2F * jetptSumSumrawpt_real[nbin];      //Added
    TH2F * jetptSumSumpt_fake[nbin];      //Added
    TH2F * jetptSumSumpt_real[nbin];      //Added

    TH1F* Pthat;
    TH1F* PthatAfXw;
    TH2F * rawptJES[5];

    TH2F * rawptJESEtaBin[5][7];
TH2F * jetptchMaxptEtaBin[nbin][netabin];	//Added
TH2F * jetptneuMaxptEtaBin[nbin][netabin];	//Added
TH2F * jetptneuMaxrEtaBin[nbin][netabin];	//Added
TH2F * jetptphoMaxptEtaBin[nbin][netabin];	//Added
TH2F * jetptchSumptEtaBin[nbin][netabin];	//Added
TH2F * jetptneuSumptEtaBin[nbin][netabin];	//Added
TH2F * jetptphoSumptEtaBin[nbin][netabin];	//Added
TH2F * jetptSumSumptEtaBin[nbin][netabin];	//Added
TH2F * jetptSumSumrawptEtaBin[nbin][netabin];		//Added
    TH1F * jetptEtaBin[5][7];
};

hist_class::hist_class()
{
  for(int itrig=0; itrig< 5; itrig++){
    jetpt[itrig] = new TH1F(Form("%sjetpt_0-100%%",TrigName[itrig].Data()),Form("%sjetpt_0-100%%",TrigName[itrig].Data()),1000,0.,1000.);
    jetpt[itrig] -> Sumw2(); 
    jetptEta[itrig] = new TH2F(Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),Form("%sjetptEta_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,100,-5.,5.);
    jetptEta[itrig] -> Sumw2(); 
    jetptphi[itrig] = new TH2F(Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),Form("%sjetptphi_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi[itrig] -> Sumw2(); 

       Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
        Pthat->Sumw2();
        PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
        PthatAfXw->Sumw2();

    jetptchMax[itrig] = new TH2F(Form("%sjetptchMax_0-100%%",TrigName[itrig].Data()), Form("jetptchMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);       //Added
        jetptchMax[itrig]->Sumw2();      //Added

        jetptchMaxpt[itrig] = new TH2F(Form("%sjetptchMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptchMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1);    //Added
    jetptchMaxpt[itrig]->Sumw2();      //Added

        jetptneuMaxpt[itrig] = new TH2F(Form("%sjetptneuMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptneuMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1.);    //Added
    jetptneuMaxpt[itrig]->Sumw2();      //Added

        jetptneuMax[itrig] = new TH2F(Form("%sjetptneuMax_0-100%%",TrigName[itrig].Data()), Form("jetptneuMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
    jetptneuMax[itrig]->Sumw2();      //Added

        jetptneuMaxr[itrig] = new TH2F(Form("%sjetptneuMaxr_0-100%%",TrigName[itrig].Data()), Form("jetptneuMaxr_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxr[itrig]->Sumw2();      //Added

       jetptphoMaxpt[itrig] = new TH2F(Form("%sjetptphoMaxpt_0-100%%",TrigName[itrig].Data()), Form("jetptphoMaxpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1.);    //Added
        jetptphoMaxpt[itrig]->Sumw2();      //Added

        jetptphoMax[itrig] = new TH2F(Form("%sjetptphoMax_0-100%%",TrigName[itrig].Data()), Form("jetptphoMax_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoMax[itrig]->Sumw2();      //Added

        jetptchSum[itrig] = new TH2F(Form("%sjetptchSum_0-100%%",TrigName[itrig].Data()), Form("jetptchSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptchSum[itrig]->Sumw2();      //Added

        jetptneuSum[itrig] = new TH2F(Form("%sjetptneuSum_0-100%%",TrigName[itrig].Data()), Form("jetptneuSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptneuSum[itrig]->Sumw2();      //Added

        jetptphoSum[itrig] = new TH2F(Form("%sjetptphoSum_0-100%%",TrigName[itrig].Data()), Form("jetptphoSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoSum[itrig]->Sumw2();      //Added

        jetptchSumpt[itrig] = new TH2F(Form("%sjetptchSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptchSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1.);    //Added
        jetptchSumpt[itrig]->Sumw2();      //Added

        jetptneuSumpt[itrig] = new TH2F(Form("%sjetptneuSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptneuSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1.);    //Added
        jetptneuSumpt[itrig]->Sumw2();      //Added

        jetptphoSumpt[itrig] = new TH2F(Form("%sjetptphoSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptphoSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 200, 0., 1.);    //Added
        jetptphoSumpt[itrig]->Sumw2();      //Added

        jetptchN[itrig] = new TH2F(Form("%sjetptchN_0-100%%",TrigName[itrig].Data()), Form("jetptchN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptchN[itrig]->Sumw2();

        jetptneuN[itrig] = new TH2F(Form("%sjetptneuN_0-100%%",TrigName[itrig].Data()), Form("jetptneuN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptneuN[itrig]->Sumw2();

        jetptphoN[itrig] = new TH2F(Form("%sjetptphoN_0-100%%",TrigName[itrig].Data()), Form("jetptphoN_0-100%%",TrigName[itrig].Data()),1000,0., 1000., 100, 0, 100);
        jetptphoN[itrig]->Sumw2();

        jetptSumSum[itrig] = new TH2F(Form("%sjetptSumSum_0-100%%",TrigName[itrig].Data()), Form("jetptSumSum_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 1000, 0., 1000);    //Added
        jetptSumSum[itrig]->Sumw2();      //Added

        jetptSumSumrawpt[itrig] = new TH2F(Form("%sjetptSumSumrawpt_0-100%%",TrigName[itrig].Data()), Form("jetptSumSumrawpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumrawpt[itrig]->Sumw2();      //Added
        jetptSumSumpt[itrig] = new TH2F(Form("%sjetptSumSumpt_0-100%%",TrigName[itrig].Data()), Form("jetptSumSumpt_0-100%%",TrigName[itrig].Data()), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumpt[itrig]->Sumw2();      //Added
x_real[ibin] = new TH2F(Form("jetptchMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added


        jetptchMax_real[ibin]->Sumw2();      //Added

        jetptchMax_fake[ibin] = new TH2F(Form("jetptchMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchMax_fake[ibin]->Sumw2();      //Added
	

	jetptchMaxpt_real[ibin] = new TH2F(Form("jetptchMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);        //Added
        jetptchMaxpt_real[ibin]->Sumw2();      //Added

        jetptchMaxpt_fake[ibin] = new TH2F(Form("jetptchMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);        //Added
        jetptchMaxpt_fake[ibin]->Sumw2();      //Added	

----------------------------------------------------------------------------------------------------------------

        jetptneuMaxpt_real[ibin] = new TH2F(Form("jetptneuMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptneuMaxpt_real[ibin]->Sumw2();      //Added

        jetptneuMaxpt_fake[ibin] = new TH2F(Form("jetptneuMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptneuMaxpt_fake[ibin]->Sumw2();      //Added  



        jetptneuMax_real[ibin] = new TH2F(Form("jetptneuMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_real[ibin]->Sumw2();      //Added

        jetptneuMax_fake[ibin] = new TH2F(Form("jetptneuMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_fake[ibin]->Sumw2();      //Added  



        jetptneuMaxr_real[ibin] = new TH2F(Form("jetptneuMaxr_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxr_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 1000, 0., 5.);        //Added
        jetptneuMaxr_real[ibin]->Sumw2();      //Added

		jetptneuMaxr_fake[ibin] = new TH2F(Form("jetptneuMaxr_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuMaxr_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 1000, 0., 5.);     //Added
        jetptneuMaxr_fake[ibin]->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------------

        jetptphoMaxpt_real[ibin] = new TH2F(Form("jetptphoMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMaxpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptphoMaxpt_real[ibin]->Sumw2();      //Added

        jetptphoMaxpt_fake[ibin] = new TH2F(Form("jetptphoMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMaxpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptphoMaxpt_fake[ibin]->Sumw2();      //Added  



        jetptphoMax_real[ibin] = new TH2F(Form("jetptphoMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMax_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_real[ibin]->Sumw2();      //Added

        jetptphoMax_fake[ibin] = new TH2F(Form("jetptphoMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoMax_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_fake[ibin]->Sumw2();      //Added  

//---------------------------------------------------------------------------------------------------------

        jetptchSum_real[ibin] = new TH2F(Form("jetptchSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchSum_real[ibin]->Sumw2();      //Added

        jetptchSum_fake[ibin] = new TH2F(Form("jetptchSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptchSum_fake[ibin]->Sumw2();      //Added

        jetptneuSum_real[ibin] = new TH2F(Form("jetptneuSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuSum_real[ibin]->Sumw2();      //Added

        jetptneuSum_fake[ibin] = new TH2F(Form("jetptneuSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptneuSum_fake[ibin]->Sumw2();      //Added


        jetptphoSum_real[ibin] = new TH2F(Form("jetptphoSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoSum_real[ibin]->Sumw2();      //Added

        jetptphoSum_fake[ibin] = new TH2F(Form("jetptphoSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptphoSum_fake[ibin]->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------

        jetptchSumpt_real[ibin] = new TH2F(Form("jetptchSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptchSumpt_real[ibin]->Sumw2();      //Added

        jetptchSumpt_fake[ibin] = new TH2F(Form("jetptchSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);      //Added
        jetptchSumpt_fake[ibin]->Sumw2();      //Added


        jetptneuSumpt_real[ibin] = new TH2F(Form("jetptneuSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptneuSumpt_real[ibin]->Sumw2();      //Added

        jetptneuSumpt_fake[ibin] = new TH2F(Form("jetptneuSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);      //Added
        jetptneuSumpt_fake[ibin]->Sumw2();      //Added


        jetptphoSumpt_real[ibin] = new TH2F(Form("jetptphoSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);        //Added
        jetptphoSumpt_real[ibin]->Sumw2();      //Added

        jetptphoSumpt_fake[ibin] = new TH2F(Form("jetptphoSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);      //Added
        jetptphoSumpt_fake[ibin]->Sumw2();      //Added

//--------------------------------------------------------------------------------------------------------
	
	jetptchN_fake[ibin] = new TH2F(Form("jetptchN_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchN_fake_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptchN_fake[ibin]->Sumw2();
	
	jetptchN_real[ibin] = new TH2F(Form("jetptchN_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptchN_real_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptchN_real[ibin]->Sumw2();

	jetptneuN_fake[ibin] = new TH2F(Form("jetptneuN_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuN_fake_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptneuN_fake[ibin]->Sumw2();

	jetptneuN_real[ibin] = new TH2F(Form("jetptneuN_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptneuN_real_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptneuN_real[ibin]->Sumw2();

	jetptphoN_fake[ibin] = new TH2F(Form("jetptphoN_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoN_fake_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptphoN_fake[ibin]->Sumw2();

	jetptphoN_real[ibin] = new TH2F(Form("jetptphoN_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptphoN_real_%d-%d%%",centr[ibin],centr[ibin+1]),1000,0., 1000., 100, 0, 100);
        jetptphoN_real[ibin]->Sumw2();

//---------------------------------------------------------------------------------------------------------

        jetptSumSum_real[ibin] = new TH2F(Form("jetptSumSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSum_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 1000, 0., 1000);        //Added
        jetptSumSum_real[ibin]->Sumw2();      //Added

        jetptSumSum_fake[ibin] = new TH2F(Form("jetptSumSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSum_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 1000, 0., 1000);      //Added
        jetptSumSum_fake[ibin]->Sumw2();      //Added

        jetptSumSumrawpt_real[ibin] = new TH2F(Form("jetptSumSumrawpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumrawpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumrawpt_real[ibin]->Sumw2();      //Added

        jetptSumSumrawpt_fake[ibin] = new TH2F(Form("jetptSumSumrawpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumrawpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumrawpt_fake[ibin]->Sumw2();      //Added	
       
		jetptSumSumpt_real[ibin] = new TH2F(Form("jetptSumSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumpt_real_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumpt_real[ibin]->Sumw2();      //Added

        jetptSumSumpt_fake[ibin] = new TH2F(Form("jetptSumSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetptSumSumpt_fake_%d-%d%%",centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumpt_fake[ibin]->Sumw2();      //Added

}

    rawptJES[itrig] = new TH2F(Form("%srawptJES_0-100%%",TrigName[itrig].Data()),Form("%srawptJES_0-100%%",TrigName[itrig].Data()),1000,0.,1000.,200,0.,2.);
    rawptJES[itrig] -> Sumw2(); 
    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[itrig][ieta] = new TH1F(Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[itrig][ieta]->Sumw2();

	jetptchMaxptEtaBin[ibin][ieta] = new TH2F(Form("jetptchMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptchMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10, centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchMaxptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptphoMaxptEtaBin[ibin][ieta] = new TH2F(Form("jetptphoMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptphoMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoMaxptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptneuMaxptEtaBin[ibin][ieta] = new TH2F(Form("jetptneuMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptneuMaxptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuMaxptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptchSumptEtaBin[ibin][ieta] = new TH2F(Form("jetptchSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptchSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchSumptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptneuSumptEtaBin[ibin][ieta] = new TH2F(Form("jetptneuSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptneuSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuSumptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptphoSumptEtaBin[ibin][ieta] = new TH2F(Form("jetptphoSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptphoSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoSumptEtaBin[ibin][ieta]->Sumw2();      //Added
        jetptSumSumptEtaBin[ibin][ieta] = new TH2F(Form("jetptSumSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptSumSumptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumptEtaBin[ibin][ieta]->Sumw2();        
	jetptSumSumrawptEtaBin[ibin][ieta] = new TH2F(Form("jetptSumSumrawptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptSumSumrawptEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumrawptEtaBin[ibin][ieta]->Sumw2();
	jetptneuMaxrEtaBin[ibin][ieta] = new TH2F(Form("jetptneuMaxrEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), Form("jetptneuMaxrEtaBin%.f_%.f_Cen%d-%d%%",deta[ieta]*10,deta[ieta+1]*10,centr[ibin],centr[ibin+1]), 1000, 0., 1000., 200, 0., 1.);    //Added
        jetptneuMaxrEtaBin[ibin][ieta]->Sumw2();

      jetptSumSumptEtaBin[itrig][ieta] = new TH2F(Form("%sjetptSumSumptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptSumSumptEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,2000, 0.,10);
      jetptSumSumptEtaBin[itrig][ieta]->Sumw2();
      jetptneuMaxrEtaBin[itrig][ieta] = new TH2F(Form("%sjetptneuMaxrEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%sjetptneuMaxrEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,1000, 0. ,5.);
      jetptneuMaxrEtaBin[itrig][ieta]->Sumw2();
      rawptJESEtaBin[itrig][ieta] = new TH2F(Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),Form("%srawptJESEtaBin%.f_%.f_Cen0-100%%",TrigName[itrig].Data(),deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[itrig][ieta]->Sumw2();
      }
    }
}

/*
void hist_class::Delete()
{
  for(int itrig=0; itrig < 5; itrig++){
    delete jetpt[itrig];
    delete jetptEta[itrig];
    delete jetptphi[itrig];
            delete jetptchMax[itrig];        //Added
        delete  jetptchMaxpt[itrig];       //Added
        delete  jetptneuMaxpt[itrig]; //Added
        delete    jetptneuMaxr[itrig];   //Added
        delete    jetptneuMax[itrig];   //Added
        delete    jetptphoMaxpt[itrig]; //Added 
        delete    jetptphoMax[itrig];   //Added
        delete    jetptchSum[itrig];
        delete    jetptneuSum[itrig];
        delete    jetptphoSum[itrig];
        delete    jetptchSumpt[itrig];
        delete    jetptneuSumpt[itrig];
        delete    jetptphoSumpt[itrig];
        delete    jetptchN[itrig];
        delete    jetptneuN[itrig];
        delete    jetptphoN[itrig];
        delete    jetptSumSum[itrig];   //Added
        delete    jetptSumSumrawpt[itrig];   //Added
        delete    jetptSumSumpt[itrig];   //Added
    delete rawptJES[itrig];
    for(int ieta=0; ieta < netabin; ieta++){
	 delete jetptEtaBin[itrig][ieta];
	 delete jetptSumSumptEtaBin[itrig][ieta];
	 delete jetptneuMaxrEtaBin[itrig][ieta];
	 delete rawptJESEtaBin[itrig][ieta];
	}
    }
}
*/

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "DATA";
  out_name=Form("AllTrigLumi%s%s%sDiJetMass_NoPtSel.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/scratch/xuq7/RpA/TreeAna/%s",out_name.Data()),"RECREATE");  

  for(int itrig=0; itrig < 5; itrig++){ 
    jetpt[itrig]->Write();
    jetptEta[itrig]->Write();
    jetptphi[itrig]->Write();
        jetptchMax[itrig]->Write();      //Added
jetptchMaxpt[itrig]->Write();       //Added
jetptneuMaxpt[itrig]->Write(); //Added
jetptneuMaxr[itrig]->Write();   //Added
jetptneuMax[itrig]->Write();   //Added
jetptphoMaxpt[itrig]->Write(); //Added 
jetptphoMax[itrig]->Write();   //Added
jetptchSum[itrig]->Write();
jetptneuSum[itrig]->Write();
jetptphoSum[itrig]->Write();
jetptchSumpt[itrig]->Write();
jetptneuSumpt[itrig]->Write();
jetptphoSumpt[itrig]->Write();
jetptchN[itrig]->Write();
jetptneuN[itrig]->Write();
jetptphoN[itrig]->Write();
jetptSumSum[itrig]->Write();   //Added
jetptSumSumrawpt[itrig]->Write();   //Added
jetptSumSumpt[itrig]->Write();   //Added

 jetptchMax_fake[ibin]->Write();	//Added
 jetptchMax_real[ibin]->Write();	//Added
  jetptchMaxpt_real[ibin]->Write();       //Added
  jetptchMaxpt_fake[ibin]->Write();       //Added
  jetptneuMaxpt_real[ibin]->Write();    //Added
  jetptneuMaxpt_fake[ibin]->Write();    //Added
    jetptneuMaxr_fake[ibin]->Write();      //Added
    jetptneuMaxr_real[ibin]->Write();      //Added
    jetptneuMax_fake[ibin]->Write();      //Added
    jetptneuMax_real[ibin]->Write();      //Added
    jetptphoMaxpt_fake[ibin]->Write();    //Added
    jetptphoMaxpt_real[ibin]->Write();    //Added
    jetptphoMax_fake[ibin]->Write();      //Added
    jetptphoMax_real[ibin]->Write();      //Added
 	  jetptchSum_fake[ibin]->Write();
 	  jetptchSum_real[ibin]->Write();
 	  jetptneuSum_fake[ibin]->Write();
 	  jetptneuSum_real[ibin]->Write();
 	  jetptphoSum_fake[ibin]->Write();
 	  jetptphoSum_real[ibin]->Write();
 	  jetptchSumpt_fake[ibin]->Write();
 	  jetptchSumpt_real[ibin]->Write();
 	  jetptneuSumpt_fake[ibin]->Write();
 	  jetptneuSumpt_real[ibin]->Write();
 	  jetptphoSumpt_fake[ibin]->Write();
 	  jetptphoSumpt_real[ibin]->Write();
 	  jetptchN_fake[ibin]->Write();
 	  jetptchN_real[ibin]->Write();
    jetptneuN_fake[ibin]->Write();
    jetptneuN_real[ibin]->Write();
    jetptphoN_fake[ibin]->Write();
    jetptphoN_real[ibin]->Write();
    jetptSumSum_fake[ibin]->Write();      //Added
    jetptSumSum_real[ibin]->Write();      //Added
    jetptSumSumrawpt_fake[ibin]->Write();      //Added
    jetptSumSumrawpt_real[ibin]->Write();      //Added
    jetptSumSumpt_fake[ibin]->Write();      //Added
    jetptSumSumpt_real[ibin]->Write();      //Added
    Pthat->Write();
    PthatAfXw->Write();
    rawptJES[itrig]->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[itrig][ieta]->Write(); 
	 jetptSumSumptEtaBin[itrig][ieta]->Write();
	 jetptneuMaxrEtaBin[itrig][ieta]->Write();
	 rawptJESEtaBin[itrig][ieta]->Write(); 
 	 }
	}
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void JetIDAllTriganalyzeRpATreeLumi()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing MC!"<<endl;

  TFile *f = new TFile("/scratch/xuq7/RpA/TreeAna/PPbMCKurtForest_ppReco_akPu3PF_QCDjetTrig_JetPt0noIPupperCut.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt,jteta,jtphi,rawpt;
	Int_t chargedN, neutralN, photonN;
Double_t chargedSum, photonSum, neutralSum, chargedMax, neutralMax,photonMax, eSum, muSum;
  Double_t weight;

  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("jtphi",&jtphi);
  nt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    nt->SetBranchAddress("chargedN",&chargedN);
    nt->SetBranchAddress("photonN",&photonN);
    nt->SetBranchAddress("neutralN",&neutralN);
    nt->SetBranchAddress("chargedMax",&chargedMax);
    nt->SetBranchAddress("photonMax",&photonMax);
    nt->SetBranchAddress("neutralMax",&neutralMax);
    nt->SetBranchAddress("chargedSum",&chargedSum);
    nt->SetBranchAddress("photonSum",&photonSum);
    nt->SetBranchAddress("neutralSum",&neutralSum);
    nt->SetBranchAddress("muSum",&muSum);
    nt->SetBranchAddress("eSum",&eSum);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
//    for(int i=0; i<50000; i++){
    nt->GetEntry(i);

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
 
    
    double jet_pt = jtpt;
    double ref_pt = refpt;
    double jet_eta = jteta;  
    double jet_phi = jtphi; 
    double raw_pt = rawpt;
  
      int dEtaBin = -1.;
	  my_hists->jetptchMax[itrig]->Fill(jet_pt*jetweight, chargedMax, weight);  //Added
    	  my_hists->jetptneuMax[itrig]->Fill(jet_pt*jetweight, neutralMax, weight);
    	  my_hists->jetptphoMax[itrig]->Fill(jet_pt*jetweight, photonMax, weight);
    	  my_hists->jetptchMaxpt[itrig]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);
    	  my_hists->jetptneuMaxpt[itrig]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);
    	  my_hists->jetptphoMaxpt[itrig]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr[itrig]->Fill(jet_pt*jetweight, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum[itrig]->Fill(jet_pt*jetweight, chargedSum, weight);
          my_hists->jetptneuSum[itrig]->Fill(jet_pt*jetweight, neutralSum, weight);
          my_hists->jetptphoSum[itrig]->Fill(jet_pt*jetweight, photonSum, weight);
          my_hists->jetptchSumpt[itrig]->Fill(jet_pt*jetweight, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt[itrig]->Fill(jet_pt*jetweight, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt[itrig]->Fill(jet_pt*jetweight, photonSum/jet_pt, weight);

          my_hists->jetptchN[itrig]->Fill(jet_pt*jetweight, chargedN, weight);
          my_hists->jetptneuN[itrig]->Fill(jet_pt*jetweight, neutralN, weight);
          my_hists->jetptphoN[itrig]->Fill(jet_pt*jetweight, photonN, weight);

          my_hists->jetptSumSum[itrig]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt[itrig]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt[itrig]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	my_hists->jetptEta[itrig]->Fill(jet_pt*jetweight,jet_eta,weight);
         my_hists->jetptphi[itrig]->Fill(jet_pt*jetweight,jet_phi,weight);
         my_hists->jetEtaphi[itrig]->Fill(jet_eta,jet_phi,weight);
	
	if(ref_pt<0){	//Added
	 my_hists->jetptchMax_fake[itrig]->Fill(jet_pt*jetweight, chargedMax, weight);  //Added
          my_hists->jetptneuMax_fake[itrig]->Fill(jet_pt*jetweight, neutralMax, weight);
          my_hists->jetptphoMax_fake[itrig]->Fill(jet_pt*jetweight, photonMax, weight);
          my_hists->jetptchMaxpt_fake[itrig]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);
          my_hists->jetptneuMaxpt_fake[itrig]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);
          my_hists->jetptphoMaxpt_fake[itrig]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr_fake[itrig]->Fill(jet_pt*jetweight, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum_fake[itrig]->Fill(jet_pt*jetweight, chargedSum, weight);
          my_hists->jetptneuSum_fake[itrig]->Fill(jet_pt*jetweight, neutralSum, weight);
          my_hists->jetptphoSum_fake[itrig]->Fill(jet_pt*jetweight, photonSum, weight);
          my_hists->jetptchSumpt_fake[itrig]->Fill(jet_pt*jetweight, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt_fake[itrig]->Fill(jet_pt*jetweight, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt_fake[itrig]->Fill(jet_pt*jetweight, photonSum/jet_pt, weight);

          my_hists->jetptchN_fake[itrig]->Fill(jet_pt*jetweight, chargedN, weight);
          my_hists->jetptneuN_fake[itrig]->Fill(jet_pt*jetweight, neutralN, weight);
          my_hists->jetptphoN_fake[itrig]->Fill(jet_pt*jetweight, photonN, weight);

          my_hists->jetptSumSum_fake[itrig]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_fake[itrig]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_fake[itrig]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	}
	else{
          my_hists->jetptchMax_real[itrig]->Fill(jet_pt*jetweight, chargedMax, weight);  //Added
          my_hists->jetptneuMax_real[itrig]->Fill(jet_pt*jetweight, neutralMax, weight);
          my_hists->jetptphoMax_real[itrig]->Fill(jet_pt*jetweight, photonMax, weight);
          my_hists->jetptchMaxpt_real[itrig]->Fill(jet_pt*jetweight, chargedMax/jet_pt, weight);
          my_hists->jetptneuMaxpt_real[itrig]->Fill(jet_pt*jetweight, neutralMax/jet_pt, weight);
          my_hists->jetptphoMaxpt_real[itrig]->Fill(jet_pt*jetweight, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr_real[itrig]->Fill(jet_pt*jetweight, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum_real[itrig]->Fill(jet_pt*jetweight, chargedSum, weight);
          my_hists->jetptneuSum_real[itrig]->Fill(jet_pt*jetweight, neutralSum, weight);
          my_hists->jetptphoSum_real[itrig]->Fill(jet_pt*jetweight, photonSum, weight);
          my_hists->jetptchSumpt_real[itrig]->Fill(jet_pt*jetweight, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt_real[itrig]->Fill(jet_pt*jetweight, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt_real[itrig]->Fill(jet_pt*jetweight, photonSum/jet_pt, weight);
        
          my_hists->jetptchN_real[itrig]->Fill(jet_pt*jetweight, chargedN, weight);
          my_hists->jetptneuN_real[itrig]->Fill(jet_pt*jetweight, neutralN, weight);
          my_hists->jetptphoN_real[itrig]->Fill(jet_pt*jetweight, photonN, weight);
      
          my_hists->jetptSumSum_real[itrig]->Fill(jet_pt*jetweight, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_real[itrig]->Fill(jet_pt*jetweight, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_real[itrig]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);	
	}	

      if(TMath::Abs(jet_eta+0.465)<=1.) {
		my_hists->jetpt[itrig]->Fill(jet_pt*jetweight,weight);
        	my_hists->rawptJES[itrig]->Fill(raw_pt,jet_pt/raw_pt,weight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,weight);
	 my_hists->jetptSumSumptEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	 my_hists->jetptneuMaxrEtaBin[itrig][dEtaBin]->Fill(jet_pt*jetweight,neutralMax/TMath::Max(chargedSum,neutralSum),weight);
	 my_hists->rawptJESEtaBin[itrig][dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	}
    } //loop over itrig
  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  









  











