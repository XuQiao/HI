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
const double deta[]={-2.2, -1.2, -0.7, -0.3, 0.3, 0.7,1.2,2.2} ;
const int netabin = sizeof(deta)/sizeof(Double_t)-1 ;
const float jetptcut = 30. ;
const float trkcut = 0.;
enum Display_t {kJet, kTrack} ;

void AllEtaCombineJetTrigUsingPreScale()
{

    double etaMin[8] = {-1.,-2.2,-1.2,-0.7,-0.3,0.3,0.7,1.2};
    double etaMax[8] = { 1.,-1.2,-0.7,-0.3, 0.3,0.7,1.2,2.2};
 
    double etaMinRev[8]; 
    double etaMaxRev[8];
  
    double jetptbin[] = {27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};   
    int nbins = sizeof(jetptbin)/sizeof(double)-1; 

    for(int j=0; j<8; j++){
      etaMinRev[j] = -etaMax[j];
      etaMaxRev[j] = -etaMin[j];//swap eta convention
      }  

    TString algo ="akPu3PF"; //"AkPu3PF" ;
    TString residual = "NoResidual"; //"NoResidual" ;
    TString coll = "PPbNoGplus" ; // or "PbP" ;
    TString class = "HFsumEta4Bin1" ; // "" for inclusive ;
    TString effTab = "TrkEffHIN12017v5TrkCorr2DCut" ;  //"HistIterTrkCorrtestFilterCut"; //  "Trk" ; "HistIterTrkCorrtest" ;
    const int Nfile = 5 ;
    TFile * f ;
    
    TString plotsdir = "/home/nih/CMSSW_5_3_8_HI_patch2/src/NewTree/output/CombinedSpectra/plots";
    const char * kHomeDir = "/home/nih/CMSSW_5_3_8_HI_patch2/src/NewTree/output/CombinedSpectra" ; 
    f = TFile::Open("JESAllTrigLumiDATAPPbakPu3PFDiJetMass.root","readonly");
//    f =  TFile::Open("AllTrigLumiDATAPPbakPu3PFDiJetMass.root", "readonly");

    TH1F *hTrk[8][5];
    TH1F *hCombined[8];
    TH1F *hTrkPt[8][5];
  
    TString outname ;
    outname ="RebinJESAllEtaLumiCombinedSpectraAllEtaBin.root";
//    outname ="AllEtaLumiCombinedSpectraAllEtaBin.root";    
  
    TCanvas *c1[8];   
    TH1F *hDum[8];  
    TLatex *tex[8];
    TLegend *leg[8];    

    TH1F *hCombinedRebin[8];

    for(int ieta=0; ieta<8; ieta++){
            if(TMath::Abs(etaMin[ieta])==1.){
                hTrk[ieta][0] = (TH1F*)f->Get("Jet20jetpt_0-100%");
                hTrk[ieta][1] = (TH1F*)f->Get("Jet40jetpt_0-100%");
                hTrk[ieta][2] = (TH1F*)f->Get("Jet60jetpt_0-100%");
                hTrk[ieta][3] = (TH1F*)f->Get("Jet80jetpt_0-100%");
                hTrk[ieta][4] = (TH1F*)f->Get("Jet100jetpt_0-100%");
            }
            else {
                hTrk[ieta][0] = (TH1F*)f->Get(Form("Jet20jetptEtaBin%.f_%.f_Cen0-100%%", etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][1] = (TH1F*)f->Get(Form("Jet40jetptEtaBin%.f_%.f_Cen0-100%%", etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][2] = (TH1F*)f->Get(Form("Jet60jetptEtaBin%.f_%.f_Cen0-100%%", etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][3] = (TH1F*)f->Get(Form("Jet80jetptEtaBin%.f_%.f_Cen0-100%%", etaMin[ieta]*10, etaMax[ieta]*10));
                hTrk[ieta][4] = (TH1F*)f->Get(Form("Jet100jetptEtaBin%.f_%.f_Cen0-100%%", etaMin[ieta]*10, etaMax[ieta]*10));
                }
    
     // create 1D histos
         for( int i = 0; i<5; i++)
         {
             hTrkPt[ieta][i] = (TH1F*)hTrk[ieta][i]->Clone(Form("jetPt%d%d_%dInEtaBin%.f_%.f",i,trigRes[i],trigRes[i+1], etaMin[ieta]*10, etaMax[ieta]*10));
             for(int ibin = 1; ibin < hTrkPt[ieta][i]->GetNbinsX(); ibin++){
               if(ibin >= hTrkPt[ieta][i]->FindBin(trigRes[i]) && ibin < hTrkPt[ieta][i]->FindBin(trigRes[i+1])){
                hTrkPt[ieta][i]->SetBinContent(ibin, hTrk[ieta][i]->GetBinContent(ibin));
                hTrkPt[ieta][i]->SetBinError(ibin, hTrk[ieta][i]->GetBinError(ibin));
                }
               else {
                 hTrkPt[ieta][i]->SetBinContent(ibin, 0.);
                 hTrkPt[ieta][i]->SetBinError(ibin, 0.);
                 }
               }
         }
    
    // add together spectra
    for( int i = 0; i<5; i++)
    {
        if(i==0) hCombined[ieta] = (TH1F*)hTrkPt[ieta][i]->Clone(Form("CombinedSpectraInEtaBin%.f_%.f", etaMin[ieta]*10, etaMax[ieta]*10));
        else hCombined[ieta]->Add(hTrkPt[ieta][i]);
    }
    //    hCombined->Scale(2.);
    //    convertToInvYield(hCombined);
    ////////////////////////////////////////////////////////////
 

    TFile * outf = new TFile(Form("%s/%s", kHomeDir, outname.Data()), "UPDATE");
    hCombinedRebin[ieta]=(TH1F*)hCombined[ieta]->Rebin(nbins,Form("hCombinedRebin_%d",ieta),jetptbin);
    hCombinedRebin[ieta]->Write();
//    hCombined[ieta]->Write();
    outf->Close();//write into a root file
    
    
   hCombined[ieta]->SetMarkerStyle(20);

   c1[ieta] = new TCanvas(Form("c1a_%d",ieta),Form("c1_%d",ieta),25*ieta,25*ieta,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c1[ieta]->Range(0,0,1,1);
   c1[ieta]->SetFillColor(0);
   c1[ieta]->SetBorderMode(0);
   c1[ieta]->SetBorderSize(2);
   c1[ieta]->SetTickx(1);
   c1[ieta]->SetTicky(1);
   c1[ieta]->SetLeftMargin(0.13);
   c1[ieta]->SetRightMargin(0.06);
   c1[ieta]->SetTopMargin(0.05);
   c1[ieta]->SetBottomMargin(0.16);
   c1[ieta]->SetFrameFillStyle(0);
   c1[ieta]->SetFrameBorderMode(0);

   c1[ieta]->cd();
   c1[ieta]->SetLogy();

   hDum[ieta] = new TH1F(Form("hhdum_%d",ieta),"",30,20.,320.);
   hDum[ieta]->SetMinimum(1.);
   hDum[ieta]->SetMaximum(1.e9);
   hDum[ieta]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum[ieta]->GetXaxis()->SetTitleSize(0.035);
   hDum[ieta]->GetXaxis()->SetTitleOffset(1.5);
   hDum[ieta]->GetXaxis()->SetLabelSize(0.035);
   hDum[ieta]->GetYaxis()->SetTitle("dN_{jet}/dp_{T}");
   hDum[ieta]->GetYaxis()->SetTitleSize(0.035);
   hDum[ieta]->GetYaxis()->SetLabelSize(0.035);
   hDum[ieta]->GetYaxis()->SetTitleOffset(1.5);
   hDum[ieta]->GetXaxis()->CenterTitle();
   hDum[ieta]->GetYaxis()->CenterTitle();

   hDum[ieta]->Draw();

   hCombined[ieta]->SetMarkerStyle(24);
   hCombined[ieta]->SetMarkerColor(1);
   hCombined[ieta]->SetMarkerSize(1.5);
   hCombined[ieta]->SetLineColor(1);  
   hCombined[ieta]->DrawCopy("same");

   for( int i = 0; i<5; i++){
      hTrk[ieta][i]->SetMarkerStyle(29);
      hTrk[ieta][i]->SetMarkerColor(1+i);
      hTrk[ieta][i]->SetLineColor(1+i);
      hTrk[ieta][i]->DrawCopy("same");
      }

      
    if(TMath::Abs(etaMin[ieta])==TMath::Abs(etaMax[ieta]))
      tex[ieta] = new TLatex(0.18,0.80,Form("| #eta_{CM} | < %.1f", etaMax[ieta]));
    else
      tex[ieta] = new TLatex(0.18,0.80,Form("%.1f < #eta_{CM} < %.1f", etaMinRev[ieta], etaMaxRev[ieta]));
    tex[ieta]->SetNDC();
    tex[ieta]->SetTextSize(0.04);
    tex[ieta]->Draw("same");


    leg[ieta] = new TLegend(0.45,0.60,0.9,0.9);
    leg[ieta]->AddEntry(hCombined[ieta],"pPb Combined:","lp");
    for( int i = 0; i<5; i++) leg[ieta]->AddEntry(hTrk[ieta][i],Form("Jet%d *Effective Prescale (to Jet100)",trigThres[i]),"lp");
    leg[ieta]->SetFillColor(kWhite);
    leg[ieta]->Draw("same");

   // if(!SavePlot)c1[ieta]->Print(Form("CombinedSpectraEtaBin%.f_%.f.pdf",10*etaMin[ieta],10*etaMax[ieta]));
   //   if(!SavePlot)c1[ieta]->Print(Form("JESCombinedSpectraEtaBin%.f_%.f.pdf",10*etaMin[ieta],10*etaMax[ieta]));
   }//end for eta loop 

}

TH1 * normalizeByBinWidth(TH1 *histo) {
    for(int i = 1; i <= histo->GetNbinsX(); i++) {
        float content = histo->GetBinContent(i);
        float error = histo->GetBinError(i);
        histo->SetBinContent(i,content/histo->GetBinWidth(i));
        histo->SetBinError(i,error/histo->GetBinWidth(i));
    }
    return histo ;
}
TH1 * convertToInvYield(TH1 *hist) {
    for(int i = 1; i<=hist->GetNbinsX(); i++) {
        double content = hist->GetBinContent(i);
        double pt = hist->GetBinCenter(i);
        double error = hist->GetBinError(i);
        
        double new_content = content/(2.*TMath::Pi()*pt);
        double new_error = error/(2.*TMath::Pi()*pt);
        
        hist->SetBinContent(i,new_content);
        hist->SetBinError(i,new_error);
    }
    return hist ;
}
