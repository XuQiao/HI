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
TString coll = "PPb";
    
void PPbanalyzeRpAskimTree(){
    TH1F * refjetpt;
    TH1F * genjetpt;
    TH1F * jetpt;
    TH2F * jetptEta;
    TH2F * jetptphi;
    TH2F * jetEtaphi;
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
    TH2F* jetptchN;       //Added
    TH2F* jetptneuN;      //Added
    TH2F* jetptphoN;      //Added
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

    TH2F * rawptJES;
    TH2F * refptJES;
    TH1F * Vz;
    TH1F * VzW;
    TH1F * Cent;
    TH1F * CentW;
    TH1F * Pthat;
    TH1F * PthatAfXw;

    TH2F * rawptJESEtaBin[7];
    TH1F * refjetptEtaBin[7];
    TH1F * genjetptEtaBin[7];
    TH2F * refptJESEtaBin[7];
TH2F * jetptchMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxptEtaBin[7];	//Added
TH2F * jetptneuMaxrEtaBin[7];	//Added
TH2F * jetptphoMaxptEtaBin[7];	//Added
TH2F * jetptchSumptEtaBin[7];	//Added
TH2F * jetptneuSumptEtaBin[7];	//Added
TH2F * jetptphoSumptEtaBin[7];	//Added
TH2F * jetptSumSumptEtaBin[7];	//Added
TH2F * jetptSumSumrawptEtaBin[7];		//Added
    TH1F * jetptEtaBin[7];

    refjetpt = new TH1F(Form("refjetpt"),Form("refjetpt"),1000,0.,1000.);
    refjetpt -> Sumw2(); 
    genjetpt = new TH1F(Form("genjetpt"),Form("genjetpt"),1000,0.,1000.);
    genjetpt -> Sumw2(); 
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 

    jetptchMax = new TH2F(Form("jetptchMax"), Form("jetptchMax"), 1000, 0., 1000., 3000, 0., 300.);       //Added
        jetptchMax->Sumw2();      //Added

        jetptchMaxpt = new TH2F(Form("jetptchMaxpt"), Form("jetptchMaxpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptchMaxpt->Sumw2();      //Added

        jetptneuMaxpt = new TH2F(Form("jetptneuMaxpt"), Form("jetptneuMaxpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
    jetptneuMaxpt->Sumw2();      //Added

        jetptneuMax = new TH2F(Form("jetptneuMax"), Form("jetptneuMax"), 1000, 0., 1000., 3000, 0., 300.);    //Added
    jetptneuMax->Sumw2();      //Added

        jetptneuMaxr = new TH2F(Form("jetptneuMaxr"), Form("jetptneuMaxr"), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxr->Sumw2();      //Added

       jetptphoMaxpt = new TH2F(Form("jetptphoMaxpt"), Form("jetptphoMaxpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoMaxpt->Sumw2();      //Added

        jetptphoMax = new TH2F(Form("jetptphoMax"), Form("jetptphoMax"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoMax->Sumw2();      //Added

        jetptchSum = new TH2F(Form("jetptchSum"), Form("jetptchSum"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptchSum->Sumw2();      //Added

        jetptneuSum = new TH2F(Form("jetptneuSum"), Form("jetptneuSum"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptneuSum->Sumw2();      //Added

        jetptphoSum = new TH2F(Form("jetptphoSum"), Form("jetptphoSum"), 1000, 0., 1000., 3000, 0., 300.);    //Added
        jetptphoSum->Sumw2();      //Added

        jetptchSumpt = new TH2F(Form("jetptchSumpt"), Form("jetptchSumpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptchSumpt->Sumw2();      //Added

        jetptneuSumpt = new TH2F(Form("jetptneuSumpt"), Form("jetptneuSumpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptneuSumpt->Sumw2();      //Added

        jetptphoSumpt = new TH2F(Form("jetptphoSumpt"), Form("jetptphoSumpt"), 1000, 0., 1000., 200, 0., 2.);    //Added
        jetptphoSumpt->Sumw2();      //Added

        jetptchN = new TH2F(Form("jetptchN"), Form("jetptchN"),1000,0., 1000., 100, 0, 100);
        jetptchN->Sumw2();

        jetptneuN = new TH2F(Form("jetptneuN"), Form("jetptneuN"),1000,0., 1000., 100, 0, 100);
        jetptneuN->Sumw2();

        jetptphoN = new TH2F(Form("jetptphoN"), Form("jetptphoN"),1000,0., 1000., 100, 0, 100);
        jetptphoN->Sumw2();

        jetptSumSum = new TH2F(Form("jetptSumSum"), Form("jetptSumSum"), 1000, 0., 1000., 1000, 0., 1000);    //Added
        jetptSumSum->Sumw2();      //Added

        jetptSumSumrawpt = new TH2F(Form("jetptSumSumrawpt"), Form("jetptSumSumrawpt"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumrawpt->Sumw2();      //Added
        jetptSumSumpt = new TH2F(Form("jetptSumSumpt"), Form("jetptSumSumpt"), 1000, 0., 1000., 2000, 0., 10);    //Added
        jetptSumSumpt->Sumw2();      //Added
		
		jetptchMax_real = new TH2F(Form("jetptchMax_real"), Form("jetptchMax_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchMax_real->Sumw2();      //Added

        jetptchMax_fake = new TH2F(Form("jetptchMax_fake"), Form("jetptchMax_fake"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchMax_fake->Sumw2();      //Added
	

	jetptchMaxpt_real = new TH2F(Form("jetptchMaxpt_real"), Form("jetptchMaxpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchMaxpt_real->Sumw2();      //Added

        jetptchMaxpt_fake = new TH2F(Form("jetptchMaxpt_fake"), Form("jetptchMaxpt_fake"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchMaxpt_fake->Sumw2();      //Added	

//----------------------------------------------------------------------------------------------------------------

        jetptneuMaxpt_real = new TH2F(Form("jetptneuMaxpt_real"), Form("jetptneuMaxpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuMaxpt_real->Sumw2();      //Added

        jetptneuMaxpt_fake = new TH2F(Form("jetptneuMaxpt_fake"), Form("jetptneuMaxpt_fake"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuMaxpt_fake->Sumw2();      //Added  



        jetptneuMax_real = new TH2F(Form("jetptneuMax_real"), Form("jetptneuMax_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_real->Sumw2();      //Added

        jetptneuMax_fake = new TH2F(Form("jetptneuMax_fake"), Form("jetptneuMax_fake"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuMax_fake->Sumw2();      //Added  



        jetptneuMaxr_real = new TH2F(Form("jetptneuMaxr_real"), Form("jetptneuMaxr_real"), 1000, 0., 1000., 1000, 0., 5.);        //Added
        jetptneuMaxr_real->Sumw2();      //Added

		jetptneuMaxr_fake = new TH2F(Form("jetptneuMaxr_fake"), Form("jetptneuMaxr_fake"), 1000, 0., 1000., 1000, 0., 5.);     //Added
        jetptneuMaxr_fake->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------------

        jetptphoMaxpt_real = new TH2F(Form("jetptphoMaxpt_real"), Form("jetptphoMaxpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoMaxpt_real->Sumw2();      //Added

        jetptphoMaxpt_fake = new TH2F(Form("jetptphoMaxpt_fake"), Form("jetptphoMaxpt_fake"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoMaxpt_fake->Sumw2();      //Added  



        jetptphoMax_real = new TH2F(Form("jetptphoMax_real"), Form("jetptphoMax_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_real->Sumw2();      //Added

        jetptphoMax_fake = new TH2F(Form("jetptphoMax_fake"), Form("jetptphoMax_fake"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoMax_fake->Sumw2();      //Added  

//---------------------------------------------------------------------------------------------------------

        jetptchSum_real = new TH2F(Form("jetptchSum_real"), Form("jetptchSum_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptchSum_real->Sumw2();      //Added

        jetptchSum_fake = new TH2F(Form("jetptchSum_fake"), Form("jetptchSum_fake"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptchSum_fake->Sumw2();      //Added

        jetptneuSum_real = new TH2F(Form("jetptneuSum_real"), Form("jetptneuSum_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptneuSum_real->Sumw2();      //Added

        jetptneuSum_fake = new TH2F(Form("jetptneuSum_fake"), Form("jetptneuSum_fake"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptneuSum_fake->Sumw2();      //Added


        jetptphoSum_real = new TH2F(Form("jetptphoSum_real"), Form("jetptphoSum_real"), 1000, 0., 1000., 3000, 0., 300.);        //Added
        jetptphoSum_real->Sumw2();      //Added

        jetptphoSum_fake = new TH2F(Form("jetptphoSum_fake"), Form("jetptphoSum_fake"), 1000, 0., 1000., 3000, 0., 300.);      //Added
        jetptphoSum_fake->Sumw2();      //Added

//---------------------------------------------------------------------------------------------------

        jetptchSumpt_real = new TH2F(Form("jetptchSumpt_real"), Form("jetptchSumpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptchSumpt_real->Sumw2();      //Added

        jetptchSumpt_fake = new TH2F(Form("jetptchSumpt_fake"), Form("jetptchSumpt_fake"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptchSumpt_fake->Sumw2();      //Added


        jetptneuSumpt_real = new TH2F(Form("jetptneuSumpt_real"), Form("jetptneuSumpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptneuSumpt_real->Sumw2();      //Added

        jetptneuSumpt_fake = new TH2F(Form("jetptneuSumpt_fake"), Form("jetptneuSumpt_fake"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptneuSumpt_fake->Sumw2();      //Added


        jetptphoSumpt_real = new TH2F(Form("jetptphoSumpt_real"), Form("jetptphoSumpt_real"), 1000, 0., 1000., 200, 0., 2.);        //Added
        jetptphoSumpt_real->Sumw2();      //Added

        jetptphoSumpt_fake = new TH2F(Form("jetptphoSumpt_fake"), Form("jetptphoSumpt_fake"), 1000, 0., 1000., 200, 0., 2.);      //Added
        jetptphoSumpt_fake->Sumw2();      //Added

//--------------------------------------------------------------------------------------------------------
	
	jetptchN_fake = new TH2F(Form("jetptchN_fake"), Form("jetptchN_fake"),1000,0., 1000., 100, 0, 100);
        jetptchN_fake->Sumw2();
	
	jetptchN_real = new TH2F(Form("jetptchN_real"), Form("jetptchN_real"),1000,0., 1000., 100, 0, 100);
        jetptchN_real->Sumw2();

	jetptneuN_fake = new TH2F(Form("jetptneuN_fake"), Form("jetptneuN_fake"),1000,0., 1000., 100, 0, 100);
        jetptneuN_fake->Sumw2();

	jetptneuN_real = new TH2F(Form("jetptneuN_real"), Form("jetptneuN_real"),1000,0., 1000., 100, 0, 100);
        jetptneuN_real->Sumw2();

	jetptphoN_fake = new TH2F(Form("jetptphoN_fake"), Form("jetptphoN_fake"),1000,0., 1000., 100, 0, 100);
        jetptphoN_fake->Sumw2();

	jetptphoN_real = new TH2F(Form("jetptphoN_real"), Form("jetptphoN_real"),1000,0., 1000., 100, 0, 100);
        jetptphoN_real->Sumw2();

//---------------------------------------------------------------------------------------------------------

        jetptSumSum_real = new TH2F(Form("jetptSumSum_real"), Form("jetptSumSum_real"), 1000, 0., 1000., 1000, 0., 1000);        //Added
        jetptSumSum_real->Sumw2();      //Added

        jetptSumSum_fake = new TH2F(Form("jetptSumSum_fake"), Form("jetptSumSum_fake"), 1000, 0., 1000., 1000, 0., 1000);      //Added
        jetptSumSum_fake->Sumw2();      //Added

        jetptSumSumrawpt_real = new TH2F(Form("jetptSumSumrawpt_real"), Form("jetptSumSumrawpt_real"), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumrawpt_real->Sumw2();      //Added

        jetptSumSumrawpt_fake = new TH2F(Form("jetptSumSumrawpt_fake"), Form("jetptSumSumrawpt_fake"), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumrawpt_fake->Sumw2();      //Added	
       
		jetptSumSumpt_real = new TH2F(Form("jetptSumSumpt_real"), Form("jetptSumSumpt_real"), 1000, 0., 1000., 2000, 0., 10);        //Added
        jetptSumSumpt_real->Sumw2();      //Added

        jetptSumSumpt_fake = new TH2F(Form("jetptSumSumpt_fake"), Form("jetptSumSumpt_fake"), 1000, 0., 1000., 2000, 0., 10);      //Added
        jetptSumSumpt_fake->Sumw2();      //Added

    rawptJES = new TH2F(Form("rawptJES"),Form("rawptJES"),1000,0.,1000.,200,0.,2.);
    rawptJES -> Sumw2(); 
    refptJES = new TH2F(Form("refptJES"),Form("refptJES"),1000,0.,1000.,200,0.,2.);
    refptJES -> Sumw2(); 
    Vz = new TH1F(Form("Vz"),Form("Vz"),600,-30.,30.);
    Vz-> Sumw2(); 
    VzW = new TH1F(Form("VzW"),Form("VzW"),600,-30.,30.);
    VzW-> Sumw2(); 
    Cent = new TH1F(Form("Cent"),Form("Cent"),100,0.,100.);
    Cent-> Sumw2(); 
    CentW = new TH1F(Form("CentW"),Form("CentW"),100,0.,100.);
    CentW-> Sumw2(); 
    Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
    Pthat-> Sumw2(); 
    PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
    PthatAfXw-> Sumw2(); 

    for(int ieta=0; ieta < netabin; ieta++){
      jetptEtaBin[ieta] = new TH1F(Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("jetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      jetptEtaBin[ieta]->Sumw2();

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
	jetptSumSumrawptEtaBin[ieta] = new TH2F(Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptSumSumrawptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 2000, 0., 10.);    //Added
        jetptSumSumrawptEtaBin[ieta]->Sumw2();
	jetptneuMaxrEtaBin[ieta] = new TH2F(Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), Form("jetptneuMaxrEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10), 1000, 0., 1000., 1000, 0., 5.);    //Added
        jetptneuMaxrEtaBin[ieta]->Sumw2();

      refjetptEtaBin[ieta] = new TH1F(Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      refjetptEtaBin[ieta]->Sumw2();
      genjetptEtaBin[ieta] = new TH1F(Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("genjetptEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.);
      genjetptEtaBin[ieta]->Sumw2();
      refptJESEtaBin[ieta] = new TH2F(Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("refptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      refptJESEtaBin[ieta]->Sumw2();
      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }

  std::cout<<"start working\n";

  cout<<"Analyzing MC!"<<endl;

  TFile *f = new TFile("/home/xuq/Documents/HI/RpA/output/skimTree/MCPPbakPu3PFskimfile0_8.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

  Float_t jtpt[100],jteta[100],geneta[100],jtphi[100],rawpt[100],refpt[100],refeta[100],genpt[100];
	Int_t t_chargedN[100], t_neutralN[100], t_photonN[100],subid[100];
Float_t t_chargedSum[100], t_photonSum[100], t_neutralSum[100], t_chargedMax[100], t_neutralMax[100], t_photonMax[100], t_eSum[100], t_muSum[100];
Int_t nref,ngen,hiBin;
  Float_t weight,xSecWeight,vzWeight,centWeight;
  Float_t vz,pthat;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("ngen",&ngen);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("refpt",refpt);
  nt->SetBranchAddress("refeta",refeta);
  nt->SetBranchAddress("genpt",genpt);
  nt->SetBranchAddress("subid",subid);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("geneta",geneta);
  nt->SetBranchAddress("jtphi",jtphi);
    nt->SetBranchAddress("chargedN",t_chargedN);
    nt->SetBranchAddress("photonN",t_photonN);
    nt->SetBranchAddress("neutralN",t_neutralN);
    nt->SetBranchAddress("chargedMax",t_chargedMax);
    nt->SetBranchAddress("photonMax",t_photonMax);
    nt->SetBranchAddress("neutralMax",t_neutralMax);
    nt->SetBranchAddress("chargedSum",t_chargedSum);
    nt->SetBranchAddress("photonSum",t_photonSum);
    nt->SetBranchAddress("neutralSum",t_neutralSum);
    nt->SetBranchAddress("muSum",t_muSum);
    nt->SetBranchAddress("eSum",t_eSum);
    nt->SetBranchAddress("pthat",&pthat);
    nt->SetBranchAddress("xSecWeight",&xSecWeight);
    nt->SetBranchAddress("vzWeight",&vzWeight);
    nt->SetBranchAddress("centWeight",&centWeight);
    nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();
 
  for(int i=0; i<nentries; i++){
    nt->GetEntry(i);
	if(TMath::Abs(vz)>15) continue;
        
	Vz->Fill(vz,xSecWeight);
	VzW->Fill(vz,xSecWeight*vzWeight);
	Cent->Fill(hiBin,xSecWeight*vzWeight);
	CentW->Fill(hiBin,weight);
	Pthat->Fill(pthat);
	PthatAfXw->Fill(pthat,xSecWeight);

    if(i%10000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double ref_pt = refpt[j4i];
    double ref_eta = refeta[j4i];
    double jet_eta = jteta[j4i];  
    double jet_phi = jtphi[j4i]; 
    double raw_pt = rawpt[j4i];
    int sub_id = subid[j4i];

		double chargedMax = t_chargedMax[j4i];
		double chargedSum = t_chargedSum[j4i];
		double neutralMax = t_neutralMax[j4i];
		double neutralSum = t_neutralSum[j4i];
		double photonMax = t_photonMax[j4i];
		double photonSum = t_photonSum[j4i];
		int chargedN = t_chargedN[j4i];
		int neutralN = t_neutralN[j4i];
		int photonN = t_photonN[j4i];
		double muSum = t_muSum[j4i];
		double eSum = t_eSum[j4i];
	 
	if(TMath::Abs(ref_eta+0.465) <= 1.) {
        refjetpt->Fill(ref_pt, weight);
       }
    int dEtaBin = -1. ;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if((ref_eta+0.465)>deta[ieta]&&(ref_eta+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        refjetptEtaBin[dEtaBin]->Fill(ref_pt, weight);
       }

	        if(rawpt<22) continue;
        if(sub_id==0 && TMath::Abs(jet_eta)<=3){
                jetptEta->Fill(jet_pt,jet_eta,weight);
                jetptphi->Fill(jet_pt,jet_phi,weight);
                jetEtaphi->Fill(jet_eta,jet_phi,weight);
        } 

      int dEtaBin = -1.;
     if(TMath::Abs(jet_eta+0.465)<=1.) {
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
		
		jetpt->Fill(jet_pt,weight);
		refptJES->Fill(ref_pt,jet_pt/ref_pt,weight);
        	rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	 jetptchMaxptEtaBin[dEtaBin]->Fill(jet_pt,chargedMax/jet_pt,weight);     //Added
	 jetptneuMaxptEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/jet_pt,weight);     //Added
	 jetptphoMaxptEtaBin[dEtaBin]->Fill(jet_pt,photonMax/jet_pt,weight);     //Added
	 jetptchSumptEtaBin[dEtaBin]->Fill(jet_pt,chargedSum/jet_pt,weight);     //Added
	 jetptneuSumptEtaBin[dEtaBin]->Fill(jet_pt,neutralSum/jet_pt,weight);     //Added
	 jetptphoSumptEtaBin[dEtaBin]->Fill(jet_pt,photonSum/jet_pt,weight);     //Added
	 jetptSumSumptEtaBin[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	 jetptneuMaxrEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/TMath::Max(chargedSum,neutralSum),weight);
	 rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	 refptJESEtaBin[dEtaBin]->Fill(ref_pt,jet_pt/ref_pt,weight);
	}
      } //loop over jet
	
	for(int j5i = 0; j5i < ngen ; j5i++) {
    double gen_pt=genpt[j5i];
    double gen_eta=geneta[j5i];
        if(TMath::Abs(gen_eta+0.465)<=1.) {
          genjetpt->Fill(gen_pt, weight);
       }
     int    dEtaBin=-1;
        for(Int_t ieta = 0 ; ieta <netabin; ieta++){
               if((gen_eta+0.465)>deta[ieta]&&(gen_eta+0.465)<=deta[ieta+1]) dEtaBin = ieta ;
         } //assign the eta bin for jets
      if(dEtaBin!=-1){
        genjetptEtaBin[dEtaBin]->Fill(gen_pt, weight);
        }
	} //end of gen jet loop
}	//loop over tree

  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s_useskim.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/home/xuq/Documents/HI/RpA/output/TreeAna/%s",out_name.Data()),"RECREATE");  

	
    refjetpt->Write();
    genjetpt->Write();
    jetpt->Write();
    jetptEta->Write();
    jetptphi->Write();
    jetEtaphi->Write();
    jetptchMax->Write();      //Added
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
    rawptJES->Write();
    refptJES->Write();
    Vz->Write();
    VzW->Write();
    Cent->Write();
    CentW->Write();
    Pthat->Write();
    PthatAfXw->Write();
    for(int ieta=0; ieta<netabin; ieta++){
	 jetptEtaBin[ieta]->Write(); 
	 jetptchMaxptEtaBin[ieta]->Write();     //Added
	 jetptneuMaxptEtaBin[ieta]->Write();    //Added
	 jetptneuMaxrEtaBin[ieta]->Write();     //Added
	 jetptphoMaxptEtaBin[ieta]->Write();    //Added
	 jetptchSumptEtaBin[ieta]->Write();     //Added
  	 jetptneuSumptEtaBin[ieta]->Write();    //Added
	 jetptphoSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumptEtaBin[ieta]->Write();    //Added
	 jetptSumSumrawptEtaBin[ieta]->Write();         //Added
	 genjetptEtaBin[ieta]->Write(); 
	 refjetptEtaBin[ieta]->Write(); 
	 refptJESEtaBin[ieta]->Write(); 
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;

  cout<<"working done\n";
 
} 
