// this code is modified to run 8 eta bins and save all the combined spectra into a root file
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TF1.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include <TMath.h>
#include <TSystem.h>
#include "THStack.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include <iostream>
#include <iostream.h>
#include <stdlib.h>
#include <fstream.h>
#include <iomanip>
using namespace std;

bool SavePlot = kFALSE ;
const Double_t TrkBin[]={0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
const int nTrkBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;
const int trigRes[]={20,55,75,95,120,999};
//const int trigRes[]={20,};
//const double trigPres[]={2238, 40.8, 3.89, 1.41, 1.00};
const int trigThres[]={20, 40, 60, 80, 100};
//const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const double deta[]={-2.0, -1.5, -1.0, -0.5, 0.5, 1.0,1.5,2.0};
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const float jetptcut = 30. ;
const float trkcut = 0.;
enum Display_t {kJet, kTrack} ;

void CombinePreScale()
{

    double jetptbin[] = {27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};   
    int nbins = sizeof(jetptbin)/sizeof(double)-1; 

    TString kDir="/home/xuq/Documents/HI/RpA/output/TreeAna";
    TString algo ="akPu3PF"; //"AkPu3PF" ;
    TString residual = "NoResidual"; //"NoResidual" ;
    TString coll = "PPb" ; // or "PbP" ;
    TString class = "HFsumEta4Bin1" ; // "" for inclusive ;
    TString effTab = "TrkEffHIN12017v5TrkCorr2DCut" ;  //"HistIterTrkCorrtestFilterCut"; //  "Trk" ; "HistIterTrkCorrtest" ;
    const int Nfile = 5 ;
    TFile * f ;
//    TString name[]={"jetptEta","jetptphi","jetEtaphi","jetptchMax","jetptchSum","jetptneuMax","jetptneuSum","jetptphoMax","jetptphoSum","jetptchMaxpt","jetptchSumpt","jetptneuMaxpt","jetptneuSumpt","jetptphoMaxpt","jetptphoSumpt","jetptSumSumpt","jetptSumSumrawpt","jetptneuMaxr","jetptchN","jetptneuN","jetptphoN"};
    TString name[]={"jetptEta","jetptphi","jetEtaphi"};
 //   TString etaname[]={"jetpt","rawptJES","jetptchMaxpt","jetptneuMaxr","jetptSumSumpt"};
    TString etaname[]={"jetpt","rawptJES"};
    int nname=sizeof(name)/sizeof(TString);
    int netaname=sizeof(etaname)/sizeof(TString);
    f = TFile::Open(Form("%s/DATA%s%s_0_15_Univ_woCom.root",kDir.Data(),coll.Data(),algo.Data()),"readonly");
    TFile * outf = new TFile(Form("%s/DATA%s%s_0_15_Univ.root", kDir.Data(), coll.Data(),algo.Data()), "RECREATE");
//    f =  TFile::Open("AllTrigLumiDATAPPbakPu3PFDiJetMass.root", "readonly");
for(int m=0;m<nname;m++){
    TH1 *hCombined;
    TH1 *hTrkPt[5];
  
                hTrkPt[0] = (TH1*)f->Get(Form("Jet20%s_0-100%%",name[m].Data()));
                hTrkPt[1] = (TH1*)f->Get(Form("Jet40%s_0-100%%",name[m].Data()));
                hTrkPt[2] = (TH1*)f->Get(Form("Jet60%s_0-100%%",name[m].Data()));
                hTrkPt[3] = (TH1*)f->Get(Form("Jet80%s_0-100%%",name[m].Data()));
                hTrkPt[4] = (TH1*)f->Get(Form("Jet100%s_0-100%%",name[m].Data()));
    for( int i = 0; i<5; i++)
    {
        if(i==0){ 
	hCombined = (TH1*)hTrkPt[i]->Clone(Form("%s",name[m].Data()));
	hCombined->SetTitle("");
	}
        else hCombined->Add(hTrkPt[i]);
    }
    hCombined->Write();
}
 //    double etaMin[8] = {-1.,-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2};
//    double etaMax[8] = { 1.,-1.2,-0.7,-0.3, 0.3,0.7,1.2,2.2};
double etaMin[8] = {-1, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0,1.5};
double etaMax[8] = {1, -1.5, -1.0, -0.5, 0.5, 1.0,1.5,2.0};
TH1* hTrk[8][5];
TH1* hCombinedEtaBin[8];
for(int m=0;m<netaname;m++){
	    for(int ieta=0; ieta<8; ieta++){
            if(ieta==0){
                hTrk[ieta][0] = (TH1*)f->Get(Form("Jet20%s_0-100%%",etaname[m].Data()));
                hTrk[ieta][1] = (TH1*)f->Get(Form("Jet40%s_0-100%%",etaname[m].Data()));
                hTrk[ieta][2] = (TH1*)f->Get(Form("Jet60%s_0-100%%",etaname[m].Data()));
                hTrk[ieta][3] = (TH1*)f->Get(Form("Jet80%s_0-100%%",etaname[m].Data()));
                hTrk[ieta][4] = (TH1*)f->Get(Form("Jet100%s_0-100%%",etaname[m].Data()));
            }
            else {
                hTrk[ieta][0] = (TH1*)f->Get(Form("Jet20%sEtaBin%.f_%.f_Cen0-100%%",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][1] = (TH1*)f->Get(Form("Jet40%sEtaBin%.f_%.f_Cen0-100%%",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][2] = (TH1*)f->Get(Form("Jet60%sEtaBin%.f_%.f_Cen0-100%%",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][3] = (TH1*)f->Get(Form("Jet80%sEtaBin%.f_%.f_Cen0-100%%",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][4] = (TH1*)f->Get(Form("Jet100%sEtaBin%.f_%.f_Cen0-100%%",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
                }

 for(int i = 0; i<5; i++)
    {
        if(i==0) {
	hCombinedEtaBin[ieta] = (TH1*)hTrk[ieta][i]->Clone(Form("%sEtaBin%.f_%.f",etaname[m].Data(), etaMin[ieta]*10, etaMax[ieta]*10));
	hCombinedEtaBin[ieta]->SetTitle("");
	}
        else hCombinedEtaBin[ieta]->Add(hTrk[ieta][i]);
    }
    hCombinedEtaBin[ieta]->Write();
}
}

    outf->Close();//write into a root file
    
}
