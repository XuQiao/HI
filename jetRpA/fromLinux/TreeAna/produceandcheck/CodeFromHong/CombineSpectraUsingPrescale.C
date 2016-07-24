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
const int trigRes[]={20,53,75,95,120,999};
const int trigThres[]={20, 40, 60, 80, 100};
const float jetptcut = 30. ;
const float trkcut = 0.;
enum Display_t {kJet, kTrack} ;

void CombineSpectraUsingPrescale()
{

    double etaMin[8] = {-1.,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5};
    double etaMax[8] = { 1.,-1.5,-1.0,-0.5, 0.5,1.0,1.5,2.0};
 
    double etaMinRev[8]; 
    double etaMaxRev[8];

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
    
    TString PlotsDir = ".";
    TString InputDir = "/home/xuq/Documents/HI/RpA/output/TreeAna";
 
   f =  TFile::Open(Form("%s/DATAPPbakPu3PFskimUniv.root",InputDir.Data()), "readonly");//without JetID
//   f =  TFile::Open(Form("%s/DATAPPbakPu3PFJetID_skimUniv_woCom.root",InputDir.Data()), "readonly");//with JetID
ofstream fstr("AfterJetIDbins.txt");


    TH1F *hTrk[8][5];
    TH1F *hCombined[8];
    TH1F *hTrkPt[8][5];
    TH1F *hRatio[8][5];  
    TH1F *hCombinedRebin[8];
    TH1F *hCombinedBeforeRebin[8];   
    TH1F *hTrkRebin[8];   
    TH1F *hTrkBeforeRebin[8];
    double prevbins[25];
    for(int j =0; j<25; j++) prevbins[j] = 20.+2*j;
 
    double finalbins[4] = {75.,90.,115.,150.};  

    double ptbins[29];
    for(int j=0; j<29; j++){
      if(j<25) ptbins[j]=prevbins[j];
      else ptbins[j] = finalbins[j-25];
    }
   
    TCanvas *c1[8];  
    TCanvas *c2[8]; 
    TH1F *hDum[8];  
    TLatex *tex[8];
    TLatex *tex2[8];
    TLatex *tex3[8];
    TLatex *tex4[8];
    TLegend *leg1[8];    
    TLegend *leg2[8]; 
    TLegend *leg3[8];
    TLegend *leg4[8];

    for(int ieta=0; ieta<1; ieta++){
for(int itrig=0;itrig<5;itrig++){
            if(TMath::Abs(etaMin[ieta])==1.){
                hTrk[ieta][itrig] = (TH1F*)f->Get(Form("Jet%djetpt",trigThres[itrig]));
               // hTrk[ieta][itrig] = (TH1F*)f->Get(Form("Jet%djetpt",trigThres[itrig]));

            }
            else {
                hTrk[ieta][itrig] = (TH1F*)f->Get(Form("Jet%djetptEtaBin%d_%d", trigThres[itrig],etaMin[ieta]*10, etaMax[ieta]*10));
              //  hTrk[ieta][itrig] = (TH1F*)f->Get(Form("Jet%djetptEtaBin%d_%d", trigThres[itrig],etaMin[ieta]*10, etaMax[ieta]*10));
                }
    }
     // create 1D histos
         for( int i = 0; i<5; i++)
         {
             hTrkPt[ieta][i] = (TH1F*)hTrk[ieta][i]->Clone(Form("jetPt%d%d_%dInEtaBin%.f_%.f",i,trigRes[i],trigRes[i+1], etaMin[ieta]*10, etaMax[ieta]*10));
         /*    for(int ibin = 1; ibin < hTrkPt[ieta][i]->GetNbinsX(); ibin++){
               if(ibin >= hTrkPt[ieta][i]->FindBin(trigRes[i]) && ibin < hTrkPt[ieta][i]->FindBin(trigRes[i+1])){
                 hTrkPt[ieta][i]->SetBinContent(ibin, hTrk[ieta][i]->GetBinContent(ibin));
                hTrkPt[ieta][i]->SetBinError(ibin, hTrk[ieta][i]->GetBinError(ibin));
                }
               else {
                 hTrkPt[ieta][i]->SetBinContent(ibin, 0.);
                 hTrkPt[ieta][i]->SetBinError(ibin, 0.);
                 }
               }*/
         }
    
    // add together spectra
/*    for( int i = 0; i<5; i++)
    {
        if(i==0) hCombined[ieta] = (TH1F*)hTrkPt[ieta][i]->Clone(Form("CombinedSpectraInEtaBin%.f_%.f", etaMin[ieta]*10, etaMax[ieta]*10));
        else hCombined[ieta]->Add(hTrkPt[ieta][i]);
    }*/
 hCombined[ieta] = (TH1F*)f->Get("jetpt");
    //    hCombined->Scale(2.);
    //    convertToInvYield(hCombined);
    ////////////////////////////////////////////////////////////
   for(int ibin=0;ibin<hCombined[ieta]->GetNbinsX();ibin++){
if(hCombined[ieta]->GetBinCenter(ibin)>20 && hCombined[ieta]->GetBinCenter(ibin)<150){
	fstr<<hCombined[ieta]->GetBinCenter(ibin)<<": "<<hCombined[ieta]->GetBinContent(ibin)<<endl;
}
}
   hCombined[ieta]->SetMarkerStyle(20);

   c1[ieta] = new TCanvas(Form("c1a_%d",ieta),Form("c1_%d",ieta),99,26,599,599);

   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c1[ieta]->Range(-38.59656,-12.02472,644.9871,2.486126);
   c1[ieta]->SetFillColor(0);
   c1[ieta]->SetBorderMode(0);
   c1[ieta]->SetBorderSize(2);
   c1[ieta]->SetTickx(1);
   c1[ieta]->SetTicky(1);
   c1[ieta]->SetLeftMargin(0.1252174);
   c1[ieta]->SetRightMargin(0.06434783);
   c1[ieta]->SetTopMargin(0.08166969);
   c1[ieta]->SetBottomMargin(0.1288566);
   c1[ieta]->SetFrameFillStyle(0);
   c1[ieta]->SetFrameBorderMode(0);

   c1[ieta]->cd();
   c1[ieta]->cd()->SetLogy();

   hDum[ieta] = new TH1F(Form("hhdum_%d",ieta),"",13,20.,150.);
   hDum[ieta]->SetMinimum(1.e3);
   hDum[ieta]->SetMaximum(1.e8);
   hDum[ieta]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum[ieta]->GetXaxis()->SetLimits(20,150);
   hDum[ieta]->GetXaxis()->SetTitleSize(0.035);
   hDum[ieta]->GetXaxis()->SetLabelSize(0.035);
   hDum[ieta]->GetXaxis()->SetTitleOffset(1.5);
   hDum[ieta]->GetYaxis()->SetTitle("dN_{jet}/dp_{T}");
   hDum[ieta]->GetYaxis()->SetTitleSize(0.035);
   hDum[ieta]->GetYaxis()->SetLabelSize(0.035);
   hDum[ieta]->GetYaxis()->SetTitleOffset(1.5);
   hDum[ieta]->GetXaxis()->CenterTitle();
   hDum[ieta]->GetYaxis()->CenterTitle();

   hDum[ieta]->Draw();
      
   hCombinedBeforeRebin[ieta] =(TH1F*) hCombined[ieta]->Clone(hCombined[ieta]->GetName());
   hCombinedRebin[ieta] =(TH1F*) hCombinedBeforeRebin[ieta]->Rebin(29,Form("hCombinedRebin_%d",ieta),ptbins);
  

   hTrkBeforeRebin[ieta] =(TH1F*) hTrk[ieta][0]->Clone(hTrk[ieta][0]->GetName());
   hTrkRebin[ieta] =(TH1F*) hTrkBeforeRebin[ieta]->Rebin(29,Form("hTrkRebin_%d",ieta),ptbins);

   hCombined[ieta]->SetMarkerStyle(24);
   hCombined[ieta]->SetMarkerColor(1);
   hCombined[ieta]->SetMarkerSize(1.5);
   hCombined[ieta]->SetLineColor(1); 
   hCombined[ieta]->Rebin(2);
   normalizeByBinWidth(hCombined[ieta]);
   hCombined[ieta]->DrawCopy("same");

   for( int i = 0; i<5; i++){
      hTrk[ieta][i]->SetMarkerStyle(20);
      hTrk[ieta][i]->SetMarkerColor(1+i);
      hTrk[ieta][i]->SetLineColor(1+i);
      hTrk[ieta][i]->Rebin(2);
      normalizeByBinWidth(hTrk[ieta][i]);
      hTrk[ieta][i]->DrawCopy("same");
      }
  
  
      
    if(TMath::Abs(etaMin[ieta])==TMath::Abs(etaMax[ieta]))
      tex[ieta] = new TLatex(0.50,0.59,Form("| #eta_{CM} | < %.1f", etaMax[ieta]));
    else
      tex[ieta] = new TLatex(0.60,0.60,Form("%.1f < #eta_{CM} < %.1f", etaMinRev[ieta], etaMaxRev[ieta]));
    tex[ieta]->SetTextFont(43);
    tex[ieta]->SetTextSize(20);
    tex[ieta]->SetNDC();
    tex[ieta]->Draw("same");


    leg1[ieta] = new TLegend(0.45,0.72,0.7,0.88);
    leg1[ieta]->SetBorderSize(0);
    leg1[ieta]->SetTextFont(43);
    leg1[ieta]->SetTextSize(20);
    leg1[ieta]->AddEntry(hCombined[ieta],"Combined","p");
    for( int i = 0; i<2; i++) leg1[ieta]->AddEntry(hTrk[ieta][i],Form("Jet%d",trigThres[i]),"p");
    leg1[ieta]->SetFillColor(kWhite);
    leg1[ieta]->Draw("same");

    leg2[ieta] = new TLegend(0.72,0.72,0.82,0.88);
    leg2[ieta]->SetBorderSize(0);
    leg2[ieta]->SetTextFont(43);
    leg2[ieta]->SetTextSize(20);
    for( int i = 2; i<5; i++) leg2[ieta]->AddEntry(hTrk[ieta][i],Form("Jet%d",trigThres[i]),"p");
    leg2[ieta]->SetFillColor(kWhite);
    leg2[ieta]->Draw("same");

    drawText("CMS Preliminary, pPb #sqrt{s_{NN}} = 5.02 TeV",0.2,0.95,20);
    drawCMS(0.50,0.70,35);

     tex2[ieta] = new TLatex(0.20,0.17,"Anti-k_{T} Particle Flow PU Jets: R=0.3");
     tex2[ieta]->SetTextFont(43);
     tex2[ieta]->SetTextSize(20);
     tex2[ieta]->SetNDC();
     tex2[ieta]->Draw("same");

    //do the ratio plot   
    c2[ieta] = new TCanvas(Form("c2a_%d",ieta),Form("c2_%d",ieta),99,26,599,599);

   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c2[ieta]->Range(-38.59656,-12.02472,644.9871,2.486126);
   c2[ieta]->SetFillColor(0);
   c2[ieta]->SetBorderMode(0);
   c2[ieta]->SetBorderSize(2);
   c2[ieta]->SetTickx(1);
   c2[ieta]->SetTicky(1);
   c2[ieta]->SetLeftMargin(0.1252174);
   c2[ieta]->SetRightMargin(0.06434783);
   c2[ieta]->SetTopMargin(0.08166969);
   c2[ieta]->SetBottomMargin(0.1288566);
   c2[ieta]->SetFrameFillStyle(0);
   c2[ieta]->SetFrameBorderMode(0);

   c2[ieta]->cd();
    hDum2[ieta] = new TH1F(Form("hhdum2_%d",ieta),"",13,20.,150.);
    hDum2[ieta]->SetMinimum(0.1);
    hDum2[ieta]->SetMaximum(5);
    hDum2[ieta]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hDum2[ieta]->GetXaxis()->SetLimits(20,150);
    hDum2[ieta]->GetXaxis()->SetTitleSize(0.035);
    hDum2[ieta]->GetXaxis()->SetTitleOffset(1.5);
    hDum2[ieta]->GetXaxis()->SetLabelSize(0.035);
    hDum2[ieta]->GetYaxis()->SetTitle("Combined/SingleJetTrigger");
    hDum2[ieta]->GetYaxis()->SetTitleSize(0.035);
    hDum2[ieta]->GetYaxis()->SetLabelSize(0.035);
    hDum2[ieta]->GetYaxis()->SetTitleOffset(1.5);
    hDum2[ieta]->GetXaxis()->CenterTitle();
    hDum2[ieta]->GetYaxis()->CenterTitle();
    hDum2[ieta]->Draw();      

    hRatio[ieta][0] = (TH1F*)hCombinedRebin[ieta]->Clone(hCombined[ieta]->GetName());
    hRatio[ieta][0] ->Divide(hTrkRebin[ieta]); 
    hRatio[ieta][0]->SetMarkerStyle(20);
    hRatio[ieta][0]->SetMarkerSize(1.5);
    hRatio[ieta][0]->SetMarkerColor(kBlack);
    hRatio[ieta][0]->SetLineColor(kBlack);
    hRatio[ieta][0]->DrawCopy("same");

    for( int i = 1; i<5; i++) {
      hRatio[ieta][i] =(TH1F*)hCombined[ieta]->Clone(hCombined[ieta]->GetName());
      hRatio[ieta][i] -> Divide( hTrk[ieta][i]); 
      hRatio[ieta][i]->SetMarkerStyle(20);
      hRatio[ieta][i]->SetMarkerColor(1+i);
      hRatio[ieta][i]->SetLineColor(1+i);
     if(i>0)    hRatio[ieta][i] -> DrawCopy("same"); 
      } 

     tex3[ieta] = new TLatex(0.20,0.17,"Anti-k_{T} Particle Flow PU Jets: R=0.3");
     tex3[ieta]->SetTextFont(43);
     tex3[ieta]->SetTextSize(20);
     tex3[ieta]->SetNDC();
     tex3[ieta]->Draw("same");

    leg3[ieta] = new TLegend(0.62,0.72,0.72,0.88);
    leg3[ieta]->SetBorderSize(0);
    leg3[ieta]->SetTextFont(43);
    leg3[ieta]->SetTextSize(20);
    for( int i = 0; i<2; i++) leg3[ieta]->AddEntry(hRatio[ieta][i],Form("Jet%d",trigThres[i]),"p");
    leg3[ieta]->SetFillColor(kWhite);
    leg3[ieta]->Draw("same");

    leg4[ieta] = new TLegend(0.75,0.72,0.85,0.88);
    leg4[ieta]->SetBorderSize(0);
    leg4[ieta]->SetTextFont(43);
    leg4[ieta]->SetTextSize(20);
    for( int i = 2; i<5; i++) leg4[ieta]->AddEntry(hRatio[ieta][i],Form("Jet%d",trigThres[i]),"p");
    leg4[ieta]->SetFillColor(kWhite);
    leg4[ieta]->Draw("same");

    drawText("CMS Preliminary, pPb #sqrt{s_{NN}} = 5.02 TeV",0.2,0.95,20);
    drawCMS(0.65,0.70,35);
    tex4[ieta] = new TLatex(0.65,0.59,Form("| #eta_{CM} | < %.1f", etaMax[ieta]));
    tex4[ieta]->SetTextFont(43);
    tex4[ieta]->SetTextSize(20);
    tex4[ieta]->SetNDC();
    tex4[ieta]->Draw("same");

    if(SavePlot)
    {
       c1[ieta]->SaveAs(Form("%s/JetID2CombineSpectraUsingPrescaleEtaBin_-10_10.pdf",PlotsDir.Data()));
       c2[ieta]->SaveAs(Form("%s/JetID2RatioCombineSpectraUsingPrescaleEtaBin-10_10.pdf",PlotsDir.Data()));
    } 
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

void drawText(const char *text, float xp, float yp, int textSize=20){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(textSize); 
  tex->SetTextColor(kBlack); 
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw(); 
}


void drawCMS(float px, float py, float nLumi) {
    TLatex *lumi = new TLatex(px-0.02,py-0.05,Form("#intL dt = %.0f nb^{-1}",nLumi));
    lumi->SetTextFont(43);
    lumi->SetTextSize(20);
    lumi->SetNDC();
    lumi->Draw();
}
