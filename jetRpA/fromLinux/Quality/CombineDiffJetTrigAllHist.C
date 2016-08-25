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

bool SavePlot = kFALSE;
//bool SavePlot = kTRUE ;
const bool SaveFile = kFALSE ;
const Double_t TrkBin[]={0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 9.6, 12.0, 14.4, 19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8, 73.6, 86.4, 103.6, 120.8, 138., 155.2, 172.4, 189.6, 206.8};
const int nTrkBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;
const int trigRes[]={0, 40, 60, 75, 95, 120, 999};
const int trigThres[]={20, 40, 60, 80, 100};
enum Display_t {kJet, kTrack} ;

void CombineDiffJetTrigAllHist(Display_t dis = kJet)
{

  

    if(dis==kTrack){
        double etaMin = -1.0;
        double etaMax = 1.0 ;
    }
    else {
        double etaMin = -1.0;
        double etaMax = 1.0 ;
    }
  double etaMinRev = -etaMax;
  double etaMaxRev = -etaMin;
  
   TString algo ="ak3PF"; //"AkPu3PF" ;
//   TString residual = "NoResidual"; //"NoResidual" ;
//   TString shiftType = "hong_2minus_";
TString JetID="";
    const int Nfile = 6 ;
    TFile * f[Nfile] ;
    
    TString plotsdir = "/home/xuq7";
    const char * kHomeDir = "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/JetID/ak3PF/NoResidual";

  f[0] =  TFile::Open(Form("%s/DATAPPb_MB%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedMB_pPb_SingleTrack_Full_UsingKKForest_v1.root", kHomeDir,algo.Data(),JetID.Data()), "readonly");
  f[1] =  TFile::Open(Form("%s/DATAPPb_Jet20%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet20_pPb_Jet20_Full_UsingKKForest_v1.root", kHomeDir, algo.Data(), JetID.Data()), "readonly");
  f[2] =  TFile::Open(Form("%s/DATAPPb_Jet40%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root", kHomeDir, algo.Data(),JetID.Data()), "readonly");
  f[3] =  TFile::Open(Form("%s/DATAPPb_Jet60%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_mergedJet40Jet60_pPb_Jet40Jet60_Full_UsingKKForest_v1.root", kHomeDir, algo.Data(), JetID.Data()), "readonly");
  f[4] =  TFile::Open(Form("%s/DATAPPb_Jet80%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root", kHomeDir, algo.Data(), JetID.Data()), "readonly");
  f[5] =  TFile::Open(Form("%s/DATAPPb_Jet100%s%sDiJetMassJetPt30TrkEffHIN12017v5TrkCorr2DCut0EtaBin7_HFsumEta4Bin1_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root", kHomeDir, algo.Data(), JetID.Data()), "readonly");

//  f[3] =  TFile::Open(Form("%s/%sDATAPPb_Jet60%sJet2013HistIterTrkCorrtest0EtaCut10_mergedJet40Jet60_KK.root", kHomeDir,shiftType.Data(), algo.Data()), "readonly");

 //  f[4] =  TFile::Open(Form("%s/%sDATAPPb_Jet80%sJet2013HistIterTrkCorrtest0EtaCut10_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root", kHomeDir,shiftType.Data(), algo.Data()), "readonly");
    
 //  f[5] =  TFile::Open(Form("%s/%sDATAPPb_Jet100%sJet2013HistIterTrkCorrtest0EtaCut10_PA2013_HiForest_PromptReco_JSonPPb_forestv77.root", kHomeDir,shiftType.Data(), algo.Data()), "readonly");

  // array goes as 0=MB, 1=Trk12, 2=Trk20, 3=Trk30
        TH1D * hnevt[6];
        TH1D * hnevtNorm[5];
        
        hnevt[0] = (TH1D *) f[0]->Get("Nevents_0-100%");
        hnevt[1] = (TH1D *) f[1]->Get("NeventsWithJetAbove40_0-100%");
        hnevt[2] = (TH1D *) f[2]->Get("NeventsWithJetAbove60_0-100%");
        hnevt[3] = (TH1D *) f[3]->Get("NeventsWithJetAbove75_0-100%");
        hnevt[4] = (TH1D *) f[4]->Get("NeventsWithJetAbove95_0-100%");
        hnevt[5] = (TH1D *) f[5]->Get("NeventsWithJetAbove120_0-100%");
        
        hnevtNorm[0] = (TH1D *) f[0]->Get("NeventsWithJetAbove40_0-100%");
        hnevtNorm[1] = (TH1D *) f[1]->Get("NeventsWithJetAbove60_0-100%");
        hnevtNorm[2] = (TH1D *) f[2]->Get("NeventsWithJetAbove75_0-100%");
        hnevtNorm[3] = (TH1D *) f[3]->Get("NeventsWithJetAbove95_0-100%");
        hnevtNorm[4] = (TH1D *) f[4]->Get("NeventsWithJetAbove120_0-100%");
    
    // Calculate Nevt for each sample
    
    double nevt[6];
    
    // Adding the nevent histogram
    nevt[0] = hnevt[0]->Integral();
    nevt[1] = nevt[0] * hnevt[1]->Integral()/hnevtNorm[0]->Integral();
    nevt[2] = nevt[1] * hnevt[2]->Integral()/hnevtNorm[1]->Integral();
    nevt[3] = nevt[2] * hnevt[3]->Integral()/hnevtNorm[2]->Integral();
    nevt[4] = nevt[3] * hnevt[4]->Integral()/hnevtNorm[3]->Integral();
    nevt[5] = nevt[4] * hnevt[5]->Integral()/hnevtNorm[4]->Integral();

for(int ifile=0;ifile<Nfile;ifile++){
if(ifile!=Nfile-1){
cout<<ifile<<":"<<hnevt[ifile]->Integral()<<'\t'<<hnevtNorm[ifile]->Integral()<<'\t'<<nevt[ifile]<<endl;
cout<<ifile<<":"<<hnevt[ifile]->GetEntries()<<'\t'<<hnevtNorm[ifile]->GetEntries()<<'\t'<<nevt[ifile]<<endl;
}
else{
cout<<ifile<<":"<<hnevt[ifile]->Integral()<<'\t'<<nevt[ifile]<<endl;
cout<<ifile<<":"<<hnevt[ifile]->GetEntries()<<'\t'<<nevt[ifile]<<endl;
}
}


    TList * list[Nfile] ;
    int N ;
   for(int ifile=0; ifile<Nfile; ifile++){
        
        list[ifile] =(TList*)f[ifile]->GetListOfKeys();
        if(ifile==0){
            N = list[ifile]->GetEntries();
        }
    }
    const int Nhist = N ;
    TH1 * h[Nfile][Nhist];
    TH1 * htot[Nhist];
    for(int ifile=0; ifile<Nfile; ifile++){
        for(int ihist = 0 ; ihist < Nhist ; ihist++){
            h[ifile][ihist] = (TH1*)f[ifile]->Get(list[ifile]->At(ihist)->GetName()) ;
            //  h[ifile][ihist]->Sumw2();
            h[ifile][ihist]->Scale(1./nevt[ifile]);
        }
    }
    
    TString outname ;

    outname =Form("%sCombinedJetTrig%sJetAllHist.root",JetID.Data(), algo.Data());
    
    TFile * outf = new TFile(Form("%s/%s", kHomeDir, outname.Data()), "RECREATE");
   
    for(int ihist = 0 ; ihist < Nhist ; ihist++){
        
        for(int ifile=0; ifile<Nfile; ifile++){
            if(ifile==0)
                htot[ihist] = (TH1*)h[ifile][ihist]->Clone(Form("%s",list[ifile]->At(ihist)->GetName()));
            else
                htot[ihist]->Add(h[ifile][ihist]);
            
        }
        
        if(htot[ihist]->GetEntries()>0){
            htot[ihist]->Write();
        }
    }
    outf->Close() ;
    
TString JetID_Pro="SumSumpt";
    
  TH2F *hTrk1[6]; TProfile *hTrk[6];
    switch(dis){
        case kJet:
            if(TMath::Abs(etaMin)==1.){
/*                hTrk[0] = (TH1F*)f[0]->Get("jetpt_0-100%");
                hTrk[1] = (TH1F*)f[1]->Get("jetpt_0-100%");
                hTrk[2] = (TH1F*)f[2]->Get("jetpt_0-100%");
                hTrk[3] = (TH1F*)f[3]->Get("jetpt_0-100%");
                hTrk[4] = (TH1F*)f[4]->Get("jetpt_0-100%");
                hTrk[5] = (TH1F*)f[5]->Get("jetpt_0-100%");*/
	for(int i=0;i<6;i++){
	hTrk1[i]=(TH2F*)f[i]->Get(Form("jetpt%s_0-100%%",JetID_Pro.Data()));
	hTrk[i]=(TProfile*)hTrk1[i]->ProfileX(Form("%d",i));}
            }
            else {
                hTrk[0] = (TH1F*)f[0]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[1] = (TH1F*)f[1]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[2] = (TH1F*)f[2]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[3] = (TH1F*)f[3]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[4] = (TH1F*)f[4]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[5] = (TH1F*)f[5]->Get(Form("jetptEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
            }
            
            break;
        case kTrack :
            if(TMath::Abs(etaMin)==TMath::Abs(etaMax)){
                hTrk[0] = (TH1F*)f[0]->Get("inclusivetrackptInEtaM1P1_0-100%");
                hTrk[1] = (TH1F*)f[1]->Get("inclusivetrackptInEtaM1P1_0-100%");
                hTrk[2] = (TH1F*)f[2]->Get("inclusivetrackptInEtaM1P1_0-100%");
                hTrk[3] = (TH1F*)f[3]->Get("inclusivetrackptInEtaM1P1_0-100%");
                hTrk[4] = (TH1F*)f[4]->Get("inclusivetrackptInEtaM1P1_0-100%");
                hTrk[5] = (TH1F*)f[5]->Get("inclusivetrackptInEtaM1P1_0-100%");
            }
            else {
                hTrk[0] = (TH1F*)f[0]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[1] = (TH1F*)f[1]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[2] = (TH1F*)f[2]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[3] = (TH1F*)f[3]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[4] = (TH1F*)f[4]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));
                hTrk[5] = (TH1F*)f[5]->Get(Form("IncTrkEtaBin%.f_%.f_Cen0-100%%", etaMin*10, etaMax*10));

            }
            break ;
    }
 
  
    
    for( int i = 0; i<6; i++){
    // normalizeByBinWidth(hTrk[i]);
//        hTrk[i]->Scale(1./(etaMax-etaMin));
    }
/////////////////////////////////////////////////////////////////////////////
    
  // Select Jet-Et Range


  // create 1D histos
  TProfile * hTrkPt[6];
  for( int i = 0; i<6; i++)
  {
     hTrkPt[i] = (TProfile*)hTrk[i]->Clone(Form("jetPt%d_%dInEtaBin%.f_%.f",trigRes[i],trigRes[i+1], etaMin*10, etaMax*10));
  //    hTrkPt[i]->Scale(1./nevt[i]);
  }

    // add together spectra
    TH2F * hCombined ; TProfile *hCombinedPro;
    for( int i = 0; i<6; i++)
    {
       if(i==0) hCombined = (TH2F*)hTrk1[i]->Clone(Form("CombinedSpectraInEtaBin%.f_%.f", etaMin*10, etaMax*10));
       else hCombined->Add(hTrk1[i]);
    }
hCombinedPro=(TProfile*)hCombined->ProfileX();
//    hCombined->Scale(2.);//
//    convertToInvYield(hCombined);//
  ////////////////////////////////////////////////////////////
    
  hCombinedPro->SetMarkerStyle(20);

   TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);

   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   //c1->SetLogy();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.13);
   c1->SetRightMargin(0.06);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);

   //c1->Divide(1,2,0,0);
   c1->cd();
//   c1->SetLogy();
  if(dis==kTrack) c1->SetLogx();
  // c1->SetLogy();
  // c1->SetLogx();

   if(dis==kTrack) TH1F * hDum = new TH1F("hhdum","",10,1.,210.);
    else TH1F * hDum = new TH1F("hhdum","",10,20.,320.);
    if(dis==kTrack) {
        hDum->SetMinimum(1e-12);
       hDum->SetMaximum(100);
    }
    else {
        hDum->SetMinimum(0.8);
       hDum->SetMaximum(1.1);
    }
   hDum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   hDum->GetXaxis()->SetTitleSize(0.035);
   hDum->GetXaxis()->SetTitleOffset(1.5);
   hDum->GetXaxis()->SetLabelSize(0.035);
  //  hDum->GetYaxis()->SetTitle("1/N_{evt} dN_{ch}/dp_{T}");//just normalize by binwidth
    hDum->GetYaxis()->SetTitle(JetID_Pro);
//   hDum->GetYaxis()->SetTitle("1/N_{evt} d^{2}N_{ch}/d_p{T}d#eta");
//    hhdum->GetYaxis()->SetTitle("1/N_{ev} 1/(2#pip_{T}) d^{2}N/dp_{T}d#eta");//invariant yield
   hDum->GetYaxis()->SetTitleSize(0.035);
   hDum->GetYaxis()->SetLabelSize(0.035);
   hDum->GetYaxis()->SetTitleOffset(1.5);
    hhdum->GetXaxis()->CenterTitle();
    hhdum->GetYaxis()->CenterTitle();

  hDum->Draw();
    hCombinedPro->SetMarkerStyle(24);
    hCombinedPro->SetMarkerColor(1);
    hCombinedPro->SetMarkerSize(1.5);
    hCombinedPro->SetLineColor(1);
  
	const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
//hCombinedPro=(TProfile*)hCombinedPro->Rebin(Nbin_pt,"hCombinedPro",binbound_pt);
hCombinedPro->DrawCopy("same");
    for( int i = 0; i<6; i++){
	hTrkPt[i]=(TProfile*)hTrkPt[i]->Rebin(Nbin_pt,hTrkPt[i]->GetName(),binbound_pt);
        hTrkPt[i]->SetMarkerStyle(20+i);//29 set by yaxian
        hTrkPt[i]->SetMarkerColor(1+i);
        hTrkPt[i]->SetLineColor(1+i);
	if(i==4){
        hTrkPt[i]->SetMarkerColor(7);
        hTrkPt[i]->SetLineColor(7);
}
        hTrkPt[i]->DrawCopy("same");
        }
/*
//---do interplate check------
  TFile checkfile("CombinedJetTrigAkPu3PFJetAllHist.root");
  TH1F *hcheck = (TH1F*)checkfile.Get("jetpt_0-100%;1"); 
  hcheck->SetMarkerStyle(29);
  hcheck->SetMarkerSize(2);
  hcheck->SetMarkerColor(kRed);
  hcheck->Scale(1./10);
  hcheck->DrawCopy("same"); 
//--------------------------------
*/  

    

    
    if(TMath::Abs(etaMin)==TMath::Abs(etaMax))
      TLatex * tex = new TLatex(60.5,1.0e-9,Form("| #eta_{CM} | < %.1f", etaMax));
    else
        TLatex * tex = new TLatex(60.5,1.0e-9,Form("%.1f < #eta_{CM} < %.1f", etaMin, etaMax));
 //   TLatex * tex = new TLatex(1.5,1.0e-7,Form("| #eta_{lab} |< 1.0"));
  //TLatex * tex = new TLatex(2.6,1.05,Form(|#eta_{COM}| < 1.0");
  tex->SetTextSize(0.04);
  tex->Draw();


  TLegend * leg = new TLegend(0.5,0.23,0.9,0.5);
    leg->SetBorderSize(0);
    leg->AddEntry(hCombinedPro,"pPb Combined","lp");
//  leg->AddEntry(hCombined_CombpPbPbp_TrackTriggered_EtaCM_P13_P18_copy1,"12-017","lp");
   leg->AddEntry(hTrkPt[0],"Minimum Bias","lp");
   for( int i = 1; i<6; i++) leg->AddEntry(hTrkPt[i],Form("Jet%d",trigThres[i-1]),"lp");
 // leg->AddEntry(hTrkPtRev[0],"Corrected Tracks Pbp","lp");
  leg->SetFillColor(kWhite);
  leg->Draw();

    if(SavePlot)c1->Print(Form("%s/%sCombinedJetTrig%s.gif",plotsdir.Data(),algo.Data(),JetID_Pro.Data()));
    
//----------do the ratio for interplate check-------
/*
TCanvas *c2 = new TCanvas("c2a", "c2",0,0,600,600);
c2->cd();
hCombined->Clone("hratio");
hratio->Divide(hcheck);
hratio->DrawCopy();
*/



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
