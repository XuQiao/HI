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

    TH1F* Pthat;
    TH1F* PthatAfXw;
    TH2F * rawptJES;

    TH2F * rawptJESEtaBin[7];
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
};

hist_class::hist_class()
{
    jetpt = new TH1F(Form("jetpt"),Form("jetpt"),1000,0.,1000.);
    jetpt -> Sumw2(); 
    jetptEta = new TH2F(Form("jetptEta"),Form("jetptEta"),1000,0.,1000.,100,-5.,5.);
    jetptEta -> Sumw2(); 
    jetptphi = new TH2F(Form("jetptphi"),Form("jetptphi"),1000,0.,1000.,200, -TMath::Pi(), TMath::Pi());
    jetptphi -> Sumw2(); 
    jetEtaphi = new TH2F(Form("jetEtaphi"),Form("jetEtaphi"),100,-5.,5.,200, -TMath::Pi(), TMath::Pi());
    jetEtaphi -> Sumw2(); 

       Pthat = new TH1F(Form("Pthat"),Form("Pthat"),1000,0.,1000.);
        Pthat->Sumw2();
        PthatAfXw = new TH1F(Form("PthatAfXw"),Form("PthatAfXw"),1000,0.,1000.);
        PthatAfXw->Sumw2();

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

      rawptJESEtaBin[ieta] = new TH2F(Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),Form("rawptJESEtaBin%.f_%.f",deta[ieta]*10,deta[ieta+1]*10),1000,0.,1000.,200,0.,2.);
      rawptJESEtaBin[ieta]->Sumw2();
      }
}

void hist_class::Write()
{
  TString dataType;
  TString out_name;
  dataType = "MC";
  out_name=Form("%s%s%s.root",dataType.Data(),coll.Data(),algo.Data());

  TFile *out_file = new TFile(Form("/scratch/xuq7/RpA/TreeAna/%s",out_name.Data()),"RECREATE");  

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
    Pthat->Write();
    PthatAfXw->Write();
    rawptJES->Write();
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
	 rawptJESEtaBin[ieta]->Write(); 
 	 }
  out_file->Close();
  cout<<"Output file: "<<Form("%s",out_name.Data())<<endl;
}

void PPbanalyzeRpATree()
{
  std::cout<<"start working\n";
  hist_class *my_hists = new hist_class();

  cout<<"Analyzing MC!"<<endl;

  TFile *f = new TFile("/scratch/xuq7/RpA/TreeAna/PPbMCKurtForest_ppReco_akPu3PF_QCDjetTrig_JetPt0noIPupperCut.root");  
  
  TTree *nt = (TTree*)f->Get("nt");

  Double_t jtpt,jteta,jtphi,rawpt,refpt;
	Int_t chargedN, neutralN, photonN,subid;
Double_t chargedSum, photonSum, neutralSum, chargedMax, neutralMax,photonMax, eSum, muSum;
  Double_t weight;

  nt->SetBranchAddress("rawpt",&rawpt);
  nt->SetBranchAddress("jtpt",&jtpt);
  nt->SetBranchAddress("refpt",&refpt);
  nt->SetBranchAddress("subid",&subid);
  nt->SetBranchAddress("jteta",&jteta);
  nt->SetBranchAddress("jtphi",&jtphi);
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

	if(rawpt<22) continue;
	if(subid==0){
		my_hists->jetptEta->Fill(jet_pt,jet_etaweight);
		my_hists->jetptphi->Fill(jet_pt,jet_phi,weight);
		my_hists->jetEtaphi->Fill(jet_eta,jet_phi,weight);
	}
     if(TMath::Abs(jet_eta+0.465)<=1.) {
	  my_hists->jetptchMax->Fill(jet_pt, chargedMax, weight);  //Added
    	  my_hists->jetptneuMax->Fill(jet_pt, neutralMax, weight);
    	  my_hists->jetptphoMax->Fill(jet_pt, photonMax, weight);
    	  my_hists->jetptchMaxpt->Fill(jet_pt, chargedMax/jet_pt, weight);
    	  my_hists->jetptneuMaxpt->Fill(jet_pt, neutralMax/jet_pt, weight);
    	  my_hists->jetptphoMaxpt->Fill(jet_pt, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum->Fill(jet_pt, chargedSum, weight);
          my_hists->jetptneuSum->Fill(jet_pt, neutralSum, weight);
          my_hists->jetptphoSum->Fill(jet_pt, photonSum, weight);
          my_hists->jetptchSumpt->Fill(jet_pt, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt->Fill(jet_pt, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt->Fill(jet_pt, photonSum/jet_pt, weight);

          my_hists->jetptchN->Fill(jet_pt, chargedN, weight);
          my_hists->jetptneuN->Fill(jet_pt, neutralN, weight);
          my_hists->jetptphoN->Fill(jet_pt, photonN, weight);

          my_hists->jetptSumSum->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	
	if(ref_pt<0){	//Added
	 my_hists->jetptchMax_fake->Fill(jet_pt, chargedMax, weight);  //Added
          my_hists->jetptneuMax_fake->Fill(jet_pt, neutralMax, weight);
          my_hists->jetptphoMax_fake->Fill(jet_pt, photonMax, weight);
          my_hists->jetptchMaxpt_fake->Fill(jet_pt, chargedMax/jet_pt, weight);
          my_hists->jetptneuMaxpt_fake->Fill(jet_pt, neutralMax/jet_pt, weight);
          my_hists->jetptphoMaxpt_fake->Fill(jet_pt, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr_fake->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum_fake->Fill(jet_pt, chargedSum, weight);
          my_hists->jetptneuSum_fake->Fill(jet_pt, neutralSum, weight);
          my_hists->jetptphoSum_fake->Fill(jet_pt, photonSum, weight);
          my_hists->jetptchSumpt_fake->Fill(jet_pt, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt_fake->Fill(jet_pt, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt_fake->Fill(jet_pt, photonSum/jet_pt, weight);

          my_hists->jetptchN_fake->Fill(jet_pt, chargedN, weight);
          my_hists->jetptneuN_fake->Fill(jet_pt, neutralN, weight);
          my_hists->jetptphoN_fake->Fill(jet_pt, photonN, weight);

          my_hists->jetptSumSum_fake->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_fake->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_fake->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);
	}
	else{
          my_hists->jetptchMax_real->Fill(jet_pt, chargedMax, weight);  //Added
          my_hists->jetptneuMax_real->Fill(jet_pt, neutralMax, weight);
          my_hists->jetptphoMax_real->Fill(jet_pt, photonMax, weight);
          my_hists->jetptchMaxpt_real->Fill(jet_pt, chargedMax/jet_pt, weight);
          my_hists->jetptneuMaxpt_real->Fill(jet_pt, neutralMax/jet_pt, weight);
          my_hists->jetptphoMaxpt_real->Fill(jet_pt, photonMax/jet_pt, weight);
          my_hists->jetptneuMaxr_real->Fill(jet_pt, neutralMax/TMath::Max(chargedSum,neutralSum) ,weight);
          my_hists->jetptchSum_real->Fill(jet_pt, chargedSum, weight);
          my_hists->jetptneuSum_real->Fill(jet_pt, neutralSum, weight);
          my_hists->jetptphoSum_real->Fill(jet_pt, photonSum, weight);
          my_hists->jetptchSumpt_real->Fill(jet_pt, chargedSum/jet_pt, weight);
          my_hists->jetptneuSumpt_real->Fill(jet_pt, neutralSum/jet_pt, weight);
          my_hists->jetptphoSumpt_real->Fill(jet_pt, photonSum/jet_pt, weight);
        
          my_hists->jetptchN_real->Fill(jet_pt, chargedN, weight);
          my_hists->jetptneuN_real->Fill(jet_pt, neutralN, weight);
          my_hists->jetptphoN_real->Fill(jet_pt, photonN, weight);
      
          my_hists->jetptSumSum_real->Fill(jet_pt, chargedSum+photonSum+neutralSum+muSum+eSum,weight);
          my_hists->jetptSumSumpt_real->Fill(jet_pt, (chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
          my_hists->jetptSumSumrawpt_real->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/raw_pt,weight);	
	}	
		
		my_hists->jetpt->Fill(jet_pt,weight);
        	my_hists->rawptJES->Fill(raw_pt,jet_pt/raw_pt,weight);
	}

      for(int ieta =0; ieta<netabin;ieta++){
          if((jet_eta+0.465)>deta[ieta]&&(jet_eta+0.465)<=deta[ieta+1]) dEtaBin=ieta;
      }//assign the eta bin for jets
      if(dEtaBin!=-1){
	 my_hists->jetptEtaBin[dEtaBin]->Fill(jet_pt,weight);
	 my_hists->jetptchMaxptEtaBin[dEtaBin]->Fill(jet_pt,chargedMax/jet_pt,weight);     //Added
	 my_hists->jetptneuMaxptEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/jet_pt,weight);     //Added
	 my_hists->jetptphoMaxptEtaBin[dEtaBin]->Fill(jet_pt,photonMax/jet_pt,weight);     //Added
	 my_hists->jetptchSumptEtaBin[dEtaBin]->Fill(jet_pt,chargedSum/jet_pt,weight);     //Added
	 my_hists->jetptneuSumptEtaBin[dEtaBin]->Fill(jet_pt,neutralSum/jet_pt,weight);     //Added
	 my_hists->jetptphoSumptEtaBin[dEtaBin]->Fill(jet_pt,photonSum/jet_pt,weight);     //Added
	 my_hists->jetptSumSumptEtaBin[dEtaBin]->Fill(jet_pt,(chargedSum+photonSum+neutralSum+muSum+eSum)/jet_pt,weight);
	 my_hists->jetptneuMaxrEtaBin[dEtaBin]->Fill(jet_pt,neutralMax/TMath::Max(chargedSum,neutralSum),weight);
	 my_hists->rawptJESEtaBin[dEtaBin]->Fill(raw_pt,jet_pt/raw_pt,weight);
	}
  }//loop over tree

  my_hists->Write();
  cout<<"working done\n";
 
}  
