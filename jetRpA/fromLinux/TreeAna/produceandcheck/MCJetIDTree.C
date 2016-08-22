// this code is modified to use MC only
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
using namespace std;
//const double deta[] = {-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2,2.2};
const double deta[] = {-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0};
const int netabin = sizeof(deta)/sizeof(double)-1;

TString algo="akPu3PF";//"akPu3PF"
TString coll = "PbP";
TString wtype = "";
    
void MCJetIDTree(){
    TH1F * jetpt;
    TH1F * jetpt_real;
    TH1F * jetpt_fake;
    TH2F * jetptchMaxpt;      //Added
    TH2F * jetptchMax;    //Added
    TH2F * jetptneuMaxpt;     //Added
    TH2F * jetptneuMaxr;   //Added
    TH2F * jetptneuMax;       //Added
    TH2F * jetptphoMaxpt; //Added 
    TH2F * jetptphoMax;   //Added
    TH2F * jetptchSum;    //Added
    TH2F * jetptneuSum;    //Added
    TH2F * jetptphoSum;    //Added
    TH2F * jetptchSumpt;    //Added
    TH2F * jetptneuSumpt;    //Added
    TH2F * jetptphoSumpt;    //Added
    TH2F * jetptchN;       //Added
    TH2F * jetptneuN;      //Added
    TH2F * jetptphoN;      //Added
    TH2F * jetptSumSum;   //Added
    TH2F * jetptSumSumrawpt;   //Added
    TH2F * jetptSumSumpt;   //Added

    TH2F * jetptchMaxpt_real;	//Added
    TH2F * jetptchMaxpt_fake;	//Added
    TH2F * jetptchMax_real;	//Added
    TH2F * jetptchMax_fake;	//Added
    TH2F * jetptneuMaxpt_real;	//Added
    TH2F * jetptneuMaxpt_fake;	//Added
    TH2F * jetptneuMaxr_fake;      //Added
    TH2F * jetptneuMaxr_real;      //Added
    TH2F * jetptneuMax_fake;	//Added
    TH2F * jetptneuMax_real;	//Added
    TH2F * jetptphoMaxpt_fake;    //Added
    TH2F * jetptphoMaxpt_real;    //Added
    TH2F * jetptphoMax_fake;      //Added
    TH2F * jetptphoMax_real;      //Added
    TH2F * jetptchSum_fake;	//Added
    TH2F * jetptchSum_real;	//Added
    TH2F * jetptneuSum_fake;       //Added
    TH2F * jetptneuSum_real;       //Added
    TH2F * jetptphoSum_fake;       //Added
    TH2F * jetptphoSum_real;       //Added
    TH2F * jetptchSumpt_fake;       //Added
    TH2F * jetptchSumpt_real;       //Added
    TH2F * jetptneuSumpt_fake;       //Added
    TH2F * jetptneuSumpt_real;       //Added
    TH2F * jetptphoSumpt_fake;       //Added
    TH2F * jetptphoSumpt_real;       //Added
    TH2F* jetptchN_fake;	
    TH2F* jetptchN_real;
    TH2F* jetptneuN_fake;
    TH2F* jetptneuN_real;
    TH2F* jetptphoN_fake;
    TH2F* jetptphoN_real;
    TH2F * jetptSumSum_fake;	//Added
    TH2F * jetptSumSum_real;	//Added
    TH2F * jetptSumSumrawpt_fake;      //Added
    TH2F * jetptSumSumrawpt_real;      //Added
    TH2F * jetptSumSumpt_fake;      //Added
    TH2F * jetptSumSumpt_real;      //Added

TH1F * jetptEtaBin[7];
TH1F * jetptEtaBin_real[7];
TH1F * jetptEtaBin_fake[7];
TH2F * jetptchMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxrEtaBin[7];	//Added
TH2F * jetptphoMaxptEtaBin[7];	//Added
TH2F * jetptchSumptEtaBin[7];	//Added
TH2F * jetptneuSumptEtaBin[7];	//Added
TH2F * jetptphoSumptEtaBin[7];	//Added
TH2F * jetptSumSumptEtaBin[7];	//Added
TH2F * jetptSumSumptEtaBin_fake[7];	//Added
TH2F * jetptSumSumptEtaBin_real[7];	//Added
TH2F * jetptSumSumrawptEtaBin[7];		//Added

    jetpt = new TH1F(Form("jetptEtaBin-10_10"), Form("jetptEtaBin-10_10"), 1000, 0., 1000.);       //Added
    jetpt->Sumw2(); 
    jetpt_real = new TH1F(Form("jetptEtaBin_real-10_10"), Form("jetptEtaBin_real-10_10"), 1000, 0., 1000.);       //Added
    jetpt_real->Sumw2(); 
    jetpt_fake = new TH1F(Form("jetptEtaBin_fake-10_10"), Form("jetptEtaBin_fake-10_10"), 1000, 0., 1000.);       //Added
    jetpt_fake->Sumw2(); 

    jetptchMax = new TH2F(Form("jetptchMaxEtaBin-10_10"), Form("jetptchMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);       //Added
        jetptchMax->Sumw2();      //Added

        jetptchMaxpt = new TH2F(Form("jetptchMaxptEtaBin-10_10"), Form("jetptchMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptchMaxpt->Sumw2();      //Added

        jetptneuMaxpt = new TH2F(Form("jetptneuMaxptEtaBin-10_10"), Form("jetptneuMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptneuMaxpt->Sumw2();      //Added

        jetptneuMax = new TH2F(Form("jetptneuMaxEtaBin-10_10"), Form("jetptneuMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
    jetptneuMax->Sumw2();      //Added

        jetptneuMaxr = new TH2F(Form("jetptneuMaxrEtaBin-10_10"), Form("jetptneuMaxrEtaBin-10_10"), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxr->Sumw2();      //Added

       jetptphoMaxpt = new TH2F(Form("jetptphoMaxptEtaBin-10_10"), Form("jetptphoMaxptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoMaxpt->Sumw2();      //Added

        jetptphoMax = new TH2F(Form("jetptphoMaxEtaBin-10_10"), Form("jetptphoMaxEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoMax->Sumw2();      //Added

        jetptchSum = new TH2F(Form("jetptchSumEtaBin-10_10"), Form("jetptchSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptchSum->Sumw2();      //Added

        jetptneuSum = new TH2F(Form("jetptneuSumEtaBin-10_10"), Form("jetptneuSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptneuSum->Sumw2();      //Added

        jetptphoSum = new TH2F(Form("jetptphoSumEtaBin-10_10"), Form("jetptphoSumEtaBin-10_10"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoSum->Sumw2();      //Added

        jetptchSumpt = new TH2F(Form("jetptchSumptEtaBin-10_10"), Form("jetptchSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptchSumpt->Sumw2();      //Added

        jetptneuSumpt = new TH2F(Form("jetptneuSumptEtaBin-10_10"), Form("jetptneuSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptneuSumpt->Sumw2();      //Added

        jetptphoSumpt = new TH2F(Form("jetptphoSumptEtaBin-10_10"), Form("jetptphoSumptEtaBin-10_10"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoSumpt->Sumw2();      //Added

        jetptchN = new TH2F(Form("jetptchNEtaBin-10_10"), Form("jetptchNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptchN->Sumw2();

        jetptneuN = new TH2F(Form("jetptneuNEtaBin-10_10"), Form("jetptneuNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptneuN->Sumw2();

        jetptphoN = new TH2F(Form("jetptphoNEtaBin-10_10"), Form("jetptphoNEtaBin-10_10"),1000,0., 1000., 100, 0, 100);
        jetptphoN->Sumw2();

        jetptSumSum = new TH2F(Form("jetptSumSumEtaBin-10_10"), Form("jetptSumSumEtaBin-10_10"), 1000, 0., 1000., 1000, 0., 1000);    //Added
        jetptSumSum->Sumw2();      //Added

        jetptSumSumrawpt = new TH2F(Form("jetptSumSumrawptEtaBin-10_10"), Form("jetptSumSumrawptEtaBin-10_10"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumrawpt->Sumw2();      //Added
        jetptSumSumpt = new TH2F(Form("jetptSumSumptEtaBin-10_10"), Form("jetptSumSumptEtaBin-10_10"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumpt->Sumw2();      //Added
		
		jetptchMax_real = new TH2F(Form("jetptchMaxEtaBin_real-10_10"), Form("jetptchMaxEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchMax_real->Sumw2();      //Added

        jetptchMax_fake = new TH2F(Form("jetptchMaxEtaBin_fake-10_10"), Form("jetptchMaxEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchMax_fake->Sumw2();      //Added
	

	jetptchMaxpt_real = new TH2F(Form("jetptchMaxptEtaBin_real-10_10"), Form("jetptchMaxptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchMaxpt_real->Sumw2();      //Added

        jetptchMaxpt_fake = new TH2F(Form("jetptchMaxptEtaBin_fake-10_10"), Form("jetptchMaxptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchMaxpt_fake->Sumw2();      //Added	

//----------------------------------------------------------------------------------------------------------------

        jetptneuMaxpt_real = new TH2F(Form("jetptneuMaxptEtaBin_real-10_10"), Form("jetptneuMaxptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuMaxpt_real->Sumw2();      //Added

        jetptneuMaxpt_fake = new TH2F(Form("jetptneuMaxptEtaBin_fake-10_10"), Form("jetptneuMaxptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuMaxpt_fake->Sumw2();      //Added  



        jetptneuMax_real = new TH2F(Form("jetptneuMaxEtaBin_real-10_10"), Form("jetptneuMaxEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_real->Sumw2();      //Added

        jetptneuMax_fake = new TH2F(Form("jetptneuMaxEtaBin_fake-10_10"), Form("jetptneuMaxEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_fake->Sumw2();      //Added  



        jetptneuMaxr_real = new TH2F(Form("jetptneuMaxrEtaBin_real-10_10"), Form("jetptneuMaxrEtaBin_real-10_10"), 1000, 0., 1000., 1000, 0., 5.);        //Added
        jetptneuMaxr_real->Sumw2();      //Added

		jetptneuMaxr_fake = new TH2F(Form("jetptneuMaxrEtaBin_fake-10_10"), Form("jetptneuMaxrEtaBin_fake-10_10"), 1000, 0., 1000., 1000, 0., 5.);     //Added
        jetptneuMaxr_fake->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------------

        jetptphoMaxpt_real = new TH2F(Form("jetptphoMaxptEtaBin_real-10_10"), Form("jetptphoMaxptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoMaxpt_real->Sumw2();      //Added

        jetptphoMaxpt_fake = new TH2F(Form("jetptphoMaxptEtaBin_fake-10_10"), Form("jetptphoMaxptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoMaxpt_fake->Sumw2();      //Added  



        jetptphoMax_real = new TH2F(Form("jetptphoMaxEtaBin_real-10_10"), Form("jetptphoMaxEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_real->Sumw2();      //Added

        jetptphoMax_fake = new TH2F(Form("jetptphoMaxEtaBin_fake-10_10"), Form("jetptphoMaxEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_fake->Sumw2();      //Added  

//---------------------------------------------------------------------------------------------------------

        jetptchSum_real = new TH2F(Form("jetptchSumEtaBin_real-10_10"), Form("jetptchSumEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchSum_real->Sumw2();      //Added

        jetptchSum_fake = new TH2F(Form("jetptchSumEtaBin_fake-10_10"), Form("jetptchSumEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptchSum_fake->Sumw2();      //Added

        jetptneuSum_real = new TH2F(Form("jetptneuSumEtaBin_real-10_10"), Form("jetptneuSumEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuSum_real->Sumw2();      //Added

        jetptneuSum_fake = new TH2F(Form("jetptneuSumEtaBin_fake-10_10"), Form("jetptneuSumEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptneuSum_fake->Sumw2();      //Added


        jetptphoSum_real = new TH2F(Form("jetptphoSumEtaBin_real-10_10"), Form("jetptphoSumEtaBin_real-10_10"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoSum_real->Sumw2();      //Added

        jetptphoSum_fake = new TH2F(Form("jetptphoSumEtaBin_fake-10_10"), Form("jetptphoSumEtaBin_fake-10_10"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptphoSum_fake->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------

        jetptchSumpt_real = new TH2F(Form("jetptchSumptEtaBin_real-10_10"), Form("jetptchSumptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchSumpt_real->Sumw2();      //Added

        jetptchSumpt_fake = new TH2F(Form("jetptchSumptEtaBin_fake-10_10"), Form("jetptchSumptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptchSumpt_fake->Sumw2();      //Added


        jetptneuSumpt_real = new TH2F(Form("jetptneuSumptEtaBin_real-10_10"), Form("jetptneuSumptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuSumpt_real->Sumw2();      //Added

        jetptneuSumpt_fake = new TH2F(Form("jetptneuSumptEtaBin_fake-10_10"), Form("jetptneuSumptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptneuSumpt_fake->Sumw2();      //Added


        jetptphoSumpt_real = new TH2F(Form("jetptphoSumptEtaBin_real-10_10"), Form("jetptphoSumptEtaBin_real-10_10"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoSumpt_real->Sumw2();      //Added

        jetptphoSumpt_fake = new TH2F(Form("jetptphoSumptEtaBin_fake-10_10"), Form("jetptphoSumptEtaBin_fake-10_10"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptphoSumpt_fake->Sumw2();      //Added

//--------------------------------------------------------------------------------------------------------
	
	jetptchN_fake = new TH2F(Form("jetptchNEtaBin_fake-10_10"), Form("jetptchNEtaBin_fake-10_10"),1000,0., 1000., 100, 0, 100);
        jetptchN_fake->Sumw2();
	
	jetptchN_real = new TH2F(Form("jetptchNEtaBin_real-10_10"), Form("jetptchNEtaBin_real-10_10"),1000,0., 1000., 100, 0, 100);
        jetptchN_real->Sumw2();

	jetptneuN_fake = new TH2F(Form("jetptneuNEtaBin_fake-10_10"), Form("jetptneuNEtaBin_fake-10_10"),1000,0., 1000., 100, 0, 100);
        jetptneuN_fake->Sumw2();

	jetptneuN_real = new TH2F(Form("jetptneuNEtaBin_real-10_10"), Form("jetptneuNEtaBin_real-10_10"),1000,0., 1000., 100, 0, 100);
        jetptneuN_real->Sumw2();

	jetptphoN_fake = new TH2F(Form("jetptphoNEtaBin_fake-10_10"), Form("jetptphoNEtaBin_fake-10_10"),1000,0., 1000., 100, 0, 100);
        jetptphoN_fake->Sumw2();

	jetptphoN_real = new TH2F(Form("jetptphoNEtaBin_real-10_10"), Form("jetptphoNEtaBin_real-10_10"),1000,0., 1000., 100, 0, 100);
        jetptphoN_real->Sumw2();

//---------------------------------------------------------------------------------------------------------

        jetptSumSum_real = new TH2F(Form("jetptSumSumEtaBin_real-10_10"), Form("jetptSumSumEtaBin_real-10_10"), 1000, 0., 1000., 1000, 0., 1000);        //Added
        jetptSumSum_real->Sumw2();      //Added

        jetptSumSum_fake = new TH2F(Form("jetptSumSumEtaBin_fake-10_10"), Form("jetptSumSumEtaBin_fake-10_10"), 1000, 0., 1000., 1000, 0., 1000);      //Added
        jetptSumSum_fake->Sumw2();      //Added

        jetptSumSumrawpt_real = new TH2F(Form("jetptSumSumrawptEtaBin_real-10_10"), Form("jetptSumSumrawptEtaBin_real-10_10"), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumrawpt_real->Sumw2();      //Added

        jetptSumSumrawpt_fake = new TH2F(Form("jetptSumSumrawptEtaBin_fake-10_10"), Form("jetptSumSumrawptEtaBin_fake-10_10"), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumrawpt_fake->Sumw2();      //Added	
       
		jetptSumSumpt_real = new TH2F(Form("jetptSumSumptEtaBin_real-10_10"), Form("jetptSumSumptEtaBin_real-10_10"), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumpt_real->Sumw2();      //Added

        jetptSumSumpt_fake = new TH2F(Form("jetptSumSumptEtaBin_fake-10_10"), Form("jetptSumSumptEtaBin_fake-10_10"), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumpt_fake->Sumw2();      //Added

    for(int ieta=0; ieta < netabin; ieta++){

	jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000.);    //Added
        jetptEtaBin[ieta]->Sumw2();      //Added
	jetptEtaBin_real[ieta] = new TH1F(Form("jetptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000.);    //Added
        jetptEtaBin_real[ieta]->Sumw2();      //Added
	jetptEtaBin_fake[ieta] = new TH1F(Form("jetptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000.);    //Added
        jetptEtaBin_fake[ieta]->Sumw2();      //Added
	jetptchMaxptEtaBin[ieta] = new TH2F(Form("jetptchMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptchMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptphoMaxptEtaBin[ieta] = new TH2F(Form("jetptphoMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptphoMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptneuMaxptEtaBin[ieta] = new TH2F(Form("jetptneuMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuMaxptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuMaxptEtaBin[ieta]->Sumw2();      //Added
        jetptchSumptEtaBin[ieta] = new TH2F(Form("jetptchSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptchSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptchSumptEtaBin[ieta]->Sumw2();      //Added
        jetptneuSumptEtaBin[ieta] = new TH2F(Form("jetptneuSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptneuSumptEtaBin[ieta]->Sumw2();      //Added
        jetptphoSumptEtaBin[ieta] = new TH2F(Form("jetptphoSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptphoSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 200, 0., 1);    //Added
        jetptphoSumptEtaBin[ieta]->Sumw2();      //Added
        jetptSumSumptEtaBin[ieta] = new TH2F(Form("jetptSumSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumptEtaBin[ieta]->Sumw2();        
        jetptSumSumptEtaBin_fake[ieta] = new TH2F(Form("jetptSumSumptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumptEtaBin_fake%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumptEtaBin_fake[ieta]->Sumw2();        
        jetptSumSumptEtaBin_real[ieta] = new TH2F(Form("jetptSumSumptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumptEtaBin_real%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumptEtaBin_real[ieta]->Sumw2();        
	jetptSumSumrawptEtaBin[ieta] = new TH2F(Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumrawptEtaBin[ieta]->Sumw2();
	jetptneuMaxrEtaBin[ieta] = new TH2F(Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxrEtaBin[ieta]->Sumw2();

      }

  std::cout<<"start working\n";

  cout<<"Analyzing MC!"<<endl;

if(coll=="PbP")
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PbPMCKurtForestVzWeight_ppReco_akPu3PF_QCDjetQiaoIDCut_JetPt0noIPupperCut.root");  
if(coll=="PPb")	
  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/PPbMCKurtForest_ppReco_akPu3PF_QCDjetQiaoIDCut_JetPt0noIPupperCut.root");  
  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt,jteta,rawpt,refpt;
	Int_t t_chargedN, t_neutralN, t_photonN;
Double_t t_chargedSum, t_photonSum, t_neutralSum, t_chargedMax, t_neutralMax, t_photonMax, t_eSum, t_muSum;
Int_t nref, subid;
  Double_t weight,vz;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("subid",&subid);
    nt->SetBranchAddress("chargedN",&t_chargedN);
    nt->SetBranchAddress("photonN",&t_photonN);
    nt->SetBranchAddress("neutralN",&t_neutralN);
    nt->SetBranchAddress("chargedMax",&t_chargedMax);
    nt->SetBranchAddress("photonMax",&t_photonMax);
    nt->SetBranchAddress("neutralMax",&t_neutralMax);
    nt->SetBranchAddress("chargedSum",&t_chargedSum);
    nt->SetBranchAddress("photonSum",&t_photonSum);
    nt->SetBranchAddress("neutralSum",&t_neutralSum);
    nt->SetBranchAddress("muSum",&t_muSum);
    nt->SetBranchAddress("eSum",&t_eSum);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
	if(TMath::Abs(vz)>15) continue;
        
    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
    double jet_pt = jtpt;
    double jet_eta = jteta;  
    double ref_pt = refpt;  
    int sub_id = subid;  
    double raw_pt = rawpt;

		double chargedMax = t_chargedMax;
		double chargedSum = t_chargedSum;
		double neutralMax = t_neutralMax;
		double neutralSum = t_neutralSum;
		double photonMax = t_photonMax;
		double photonSum = t_photonSum;
		int chargedN = t_chargedN;
		int neutralN = t_neutralN;
		int photonN = t_photonN;
		double muSum = t_muSum;
		double eSum = t_eSum;
	 if(raw_pt<22) continue;
      int dEtaBin = -1;
	if(coll=="PbP")	jet_eta=jet_eta-0.465;
	if(coll=="PPb")	jet_eta=jet_eta+0.465;
     if(TMath::Abs(jet_eta)<=1.) {
	  jetpt->Fill(jet_pt,weight);
	  jetptchMax->Fill(jet_pt, chargedMax, weight);  //Added
    	  jetptneuMax->Fill(jet_pt, neutralMax, weight);
    	  jetptphoMax->Fill(jet_pt, photonMax, weight);
    	  jetptchMaxpt->Fill(jet_pt, chargedMax/jet_pt, weight);
    	  jetptneuMaxpt->Fill(jet_pt, neutralMax/jet_pt, weight);
    	  jetptphoMaxpt->Fill(jet_pt, photonMax/jet_pt, weight);
          jetptneuMaxr->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          jetptchSum->Fill(jet_pt, chargedSum, weight);
          jetptneuSum->Fill(jet_pt, neutralSum, weight);
          jetptphoSum->Fill(jet_pt, photonSum, weight);
          jetptchSumpt->Fill(jet_pt, chargedSum/jet_pt, weight);
          jetptneuSumpt->Fill(jet_pt, neutralSum/jet_pt, weight);
          jetptphoSumpt->Fill(jet_pt, photonSum/jet_pt, weight);

          jetptchN->Fill(jet_pt, chargedN, weight);
          jetptneuN->Fill(jet_pt, neutralN, weight);
          jetptphoN->Fill(jet_pt, photonN, weight);

          jetptSumSum->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          jetptSumSumpt->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          jetptSumSumrawpt->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	
	if(ref_pt<0){	//Added
	  jetpt_fake->Fill(jet_pt,weight);
	  jetptchMax_fake->Fill(jet_pt, chargedMax, weight);  //Added
          jetptneuMax_fake->Fill(jet_pt, neutralMax, weight);
          jetptphoMax_fake->Fill(jet_pt, photonMax, weight);
          jetptchMaxpt_fake->Fill(jet_pt, chargedMax/jet_pt, weight);
          jetptneuMaxpt_fake->Fill(jet_pt, neutralMax/jet_pt, weight);
          jetptphoMaxpt_fake->Fill(jet_pt, photonMax/jet_pt, weight);
          jetptneuMaxr_fake->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          jetptchSum_fake->Fill(jet_pt, chargedSum, weight);
          jetptneuSum_fake->Fill(jet_pt, neutralSum, weight);
          jetptphoSum_fake->Fill(jet_pt, photonSum, weight);
          jetptchSumpt_fake->Fill(jet_pt, chargedSum/jet_pt, weight);
          jetptneuSumpt_fake->Fill(jet_pt, neutralSum/jet_pt, weight);
          jetptphoSumpt_fake->Fill(jet_pt, photonSum/jet_pt, weight);

          jetptchN_fake->Fill(jet_pt, chargedN, weight);
          jetptneuN_fake->Fill(jet_pt, neutralN, weight);
          jetptphoN_fake->Fill(jet_pt, photonN, weight);

          jetptSumSum_fake->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          jetptSumSumpt_fake->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          jetptSumSumrawpt_fake->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	}
	else{
	  jetpt_real->Fill(jet_pt,weight);
          jetptchMax_real->Fill(jet_pt, chargedMax, weight);  //Added
          jetptneuMax_real->Fill(jet_pt, neutralMax, weight);
          jetptphoMax_real->Fill(jet_pt, photonMax, weight);
          jetptchMaxpt_real->Fill(jet_pt, chargedMax/jet_pt, weight);
          jetptneuMaxpt_real->Fill(jet_pt, neutralMax/jet_pt, weight);
          jetptphoMaxpt_real->Fill(jet_pt, photonMax/jet_pt, weight);
          jetptneuMaxr_real->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          jetptchSum_real->Fill(jet_pt, chargedSum, weight);
          jetptneuSum_real->Fill(jet_pt, neutralSum, weight);
          jetptphoSum_real->Fill(jet_pt, photonSum, weight);
          jetptchSumpt_real->Fill(jet_pt, chargedSum/jet_pt, weight);
          jetptneuSumpt_real->Fill(jet_pt, neutralSum/jet_pt, weight);
          jetptphoSumpt_real->Fill(jet_pt, photonSum/jet_pt, weight);
        
          jetptchN_real->Fill(jet_pt, chargedN, weight);
          jetptneuN_real->Fill(jet_pt, neutralN, weight);
          jetptphoN_real->Fill(jet_pt, photonN, weight);
      
          jetptSumSum_real->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          jetptSumSumpt_real->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          jetptSumSumrawpt_real->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);	
	}	
	}
	
      for(int ieta =0; ieta<netabin;ieta++){
          if(jet_eta>deta[ieta]&&jet_eta<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	if(ref_pt<0){
	 jetptEtaBin_fake[dEtaBin]->Fill(jet_pt,weight);
	 jetptSumSumptEtaBin_fake[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	}
	else{
	 jetptEtaBin_real[dEtaBin]->Fill(jet_pt,weight);
	 jetptSumSumptEtaBin_real[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	}
	 jetptchMaxptEtaBin[dEtaBin]->Fill(jet_pt,chargedMax/jet_pt,weight);     //Added
	 jetptneuMaxptEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/jet_pt,weight);     //Added
	 jetptphoMaxptEtaBin[dEtaBin]->Fill(jet_pt,photonMax/jet_pt,weight);     //Added
	 jetptchSumptEtaBin[dEtaBin]->Fill(jet_pt,chargedSum/jet_pt,weight);     //Added
	 jetptneuSumptEtaBin[dEtaBin]->Fill(jet_pt,neutralSum/jet_pt,weight);     //Added
	 jetptphoSumptEtaBin[dEtaBin]->Fill(jet_pt,photonSum/jet_pt,weight);     //Added
	 jetptSumSumrawptEtaBin[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	 jetptSumSumptEtaBin[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	 jetptneuMaxrEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/TMath::Max(chargedSum,neutralSum),weight);
	}
	
}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s%sJetID.root",dataType.Data(),coll.Data(),algo.Data(),wtype.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

jetpt->Write();       //Added
jetpt_real->Write();       //Added
jetpt_fake->Write();       //Added
jetptchMaxpt->Write();       //Added
jetptneuMaxpt->Write(); //Added
jetptneuMaxr->Write();   //Added
jetptneuMax->Write();   //Added
jetptphoMaxpt->Write(); //Added 
jetptphoMax->Write();   //Added
jetptchSum->Write();
jetptneuSum->Write();
jetptphoSum->Write();
jetptchSumpt->Write();
jetptneuSumpt->Write();
jetptphoSumpt->Write();
jetptchN->Write();
jetptneuN->Write();
jetptphoN->Write();
jetptSumSum->Write();   //Added
jetptSumSumrawpt->Write();   //Added
jetptSumSumpt->Write();   //Added

 jetptchMax_fake->Write();	//Added
 jetptchMax_real->Write();	//Added
  jetptchMaxpt_real->Write();       //Added
  jetptchMaxpt_fake->Write();       //Added
  jetptneuMaxpt_real->Write();    //Added
  jetptneuMaxpt_fake->Write();    //Added
    jetptneuMaxr_fake->Write();      //Added
    jetptneuMaxr_real->Write();      //Added
    jetptneuMax_fake->Write();      //Added
    jetptneuMax_real->Write();      //Added
    jetptphoMaxpt_fake->Write();    //Added
    jetptphoMaxpt_real->Write();    //Added
    jetptphoMax_fake->Write();      //Added
    jetptphoMax_real->Write();      //Added
 	  jetptchSum_fake->Write();
 	  jetptchSum_real->Write();
 	  jetptneuSum_fake->Write();
 	  jetptneuSum_real->Write();
 	  jetptphoSum_fake->Write();
 	  jetptphoSum_real->Write();
 	  jetptchSumpt_fake->Write();
 	  jetptchSumpt_real->Write();
 	  jetptneuSumpt_fake->Write();
 	  jetptneuSumpt_real->Write();
 	  jetptphoSumpt_fake->Write();
 	  jetptphoSumpt_real->Write();
 	  jetptchN_fake->Write();
 	  jetptchN_real->Write();
    jetptneuN_fake->Write();
    jetptneuN_real->Write();
    jetptphoN_fake->Write();
    jetptphoN_real->Write();
    jetptSumSum_fake->Write();      //Added
    jetptSumSum_real->Write();      //Added
    jetptSumSumrawpt_fake->Write();      //Added
    jetptSumSumrawpt_real->Write();      //Added
    jetptSumSumpt_fake->Write();      //Added
    jetptSumSumpt_real->Write();      //Added
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write();     //Added
	 jetptEtaBin_real[ieta]->Write();     //Added
	 jetptEtaBin_fake[ieta]->Write();     //Added
	 jetptchMaxptEtaBin[ieta]->Write();     //Added
	 jetptneuMaxptEtaBin[ieta]->Write();    //Added
	 jetptneuMaxrEtaBin[ieta]->Write();     //Added
	 jetptphoMaxptEtaBin[ieta]->Write();    //Added
	 jetptchSumptEtaBin[ieta]->Write();     //Added
  	 jetptneuSumptEtaBin[ieta]->Write();    //Added
	 jetptphoSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumptEtaBin_fake[ieta]->Write();    //Added
	 jetptSumSumptEtaBin_real[ieta]->Write();    //Added
	 jetptSumSumrawptEtaBin[ieta]->Write();         //Added
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
