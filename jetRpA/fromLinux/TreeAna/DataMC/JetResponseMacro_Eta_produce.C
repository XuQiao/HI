#include <TROOT.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <TVector3.h> 
#include <stdio.h>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TCut.h>

#define PI 3.14159;
const double leadingjetcut = 20. ;
bool SavePlot = kFALSE ;
bool DoGausFit = kFALSE ;
enum Algo_t {kPF1, kPF2, kPF3, kPF4, kPF5, kPF6, kCalo1, kCalo2, kCalo3, kCalo4, kCalo5, kCalo6} ;

//const char* AlgoNames[] = {"ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF","ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo"};
//const Int_t nBin = 12 ;
//const char* AlgoNames[] = {"ak2PF","ak3PF","ak4PF","ak5PF","ak2Calo","ak3Calo","ak4Calo","ak5Calo"};
//const Int_t nBin = 8 ;
const char* AlgoNames[] = {"ak3PF","ak4PF","ak5PF","ak3Calo","ak4Calo","ak5Calo", "akPu3PF","akPu4PF","akPu5PF","akPu3Calo","akPu4Calo","akPu5Calo"};
//const char* AlgoNames[] = {"ak3PF","ak4PF","ak5PF","ak3Calo","ak4Calo","ak5Calo"};
const Int_t nBin = 12 ;

const Double_t jetPtBin[] = {10, 20, 30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400, 500};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

double etabins[] ={-3., -2.5, -2., -1.5, -1.2,-1.0,-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2., 2.5, 3.};
const int hetab  = sizeof(etabins)/sizeof(Double_t) - 1;

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns, 
			  const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0., 
			  const Float_t leftMargin=0.2, 
			  const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05);

void JetResponseMacro_Eta(){

    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

//  gROOT->ProcessLine(".x rootlogon.C");
  gROOT->ProcessLine(".L format1Dhisto.C");
  gROOT->ProcessLine(".L formatTGraph.C");
//  gROOT->ProcessLine(".x betterColors.C");
  //TFile *FileA = TFile::Open(Form("correctedFileForBtag.root"));
//    TFile *FileA = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/pPb/JEC/corrected_dijet_pp_mergedpthatbins_2013MCV3.root"));
 //   TFile *FileA = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/pPb/JEC/dijet_Corrpp_mergedpthatbins_2013MC.root"));
 //   TFile *FileA = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JEC/jet_CorrPbp502_mergedpthatbins_2013MC.root"));
    TFile *FileA = TFile::Open(Form("/afs/cern.ch/work/y/ymao/public/RpA/ForQiao/jet_CorrPbp502_Combined2013_2014MC.root"));
//    TFile *FileA = TFile::Open(Form("/Volumes/WORKDISK/CMS/jet_Corrpp276HITrack_mergedpthatbins_2013MC.root"));

  TString outname = "JetResponse_Plots_woWeight.root"; 
  TFile* outf = new TFile(outname,"recreate");
    TString plotsdir = "./";
/*
    TString canv_name = "c1";
    const Double_t kw = 1300;
    const Double_t kh = 480;
    c1 = new TCanvas(canv_name," ",10,10,kw,kh);
  //  makeMultiPanelCanvas(c1,nBin/4,2,0.0,0.0,0.1,0.2,0.05);
    makeMultiPanelCanvas(c1,nBin/3,1,0.0,0.0,0.1,0.2,0.05);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    
    TLegend *t1=new TLegend(0.2,0.75,0.8,0.92);
    //   TLegend *t1=new TLegend(0.25,0.6,0.8,0.92);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(15);
    TLegend *t2=new TLegend(0.20,0.45,0.35,0.6);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(17);

    TH1F * dummy = new TH1F("dummy", "dummy", 500, 0., 500.);
    dummy->SetAxisRange(0., 500., "X") ;
    dummy->GetXaxis()->SetTitle("p_{T}^{ref} (GeV/c)");

//    const char* AlgoNames[] = {"ak3PF"};
//  TCanvas *c1 = new TCanvas("c1","c1",1100,700);
  //makeMultiPanelCanvas(c1,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c1B = new TCanvas("c1B","c1B",1100,700);
  //makeMultiPanelCanvas(c1B,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c1C = new TCanvas("c1C","c1C",1100,700);
  //makeMultiPanelCanvas(c1C,3,2,0.0,0.0,0.2,0.15,0.07);
//  const int cCanDiv = nBin/2;
  const int cCanDiv = 2;
//  makeMultiPanelCanvas(c1,cCanDiv,2,0.0,0.0,0.2,0.15,0.07);
  int nCanvDiv = 2;
  const int CnXbins = hetab;
    int nXbins = CnXbins;
*/
	int nXbins = hetab;
  const int npt=4;
  TH2D* iHistoCorrEta[2][4];
  const int ptbinlow[npt]={50,90,160,240};
  const int ptbinhigh[npt]={70,120,200,260};
/*
  TH2D* iHistoRefEta[nBin];
  TH2D* iHistoJtPt[nBin];
  TH2D* iHistoRawPt[nBin];
  TH1D* h1[nBin][CnXbins];
    
  TF1* Gauss6[nBin][CnXbins];
  double mean[nBin][CnXbins];
  double sigma[nBin][CnXbins];
  double meanErr[nBin][CnXbins];
  double sigmaErr[nBin][CnXbins];
  double xPoint[nBin][CnXbins];
  double xPointErr[nBin][CnXbins];
  TLegend*  leg[npt];

    for(Int_t nj= 0 ; nj <nBin ; nj++){
        for(Int_t ipt = 0 ; ipt < nXbins ; ipt++){
            mean[nj][ipt] = 0.;
            sigma[nj][ipt]= 0.;
            meanErr[nj][ipt] = 0.;
            sigmaErr[nj][ipt] =0.; 
            xPoint[nj][ipt] = 0.;
            xPointErr[nj][ipt] =0.;
           
        }
    }
*/
 //   Int_t i = 2 ;
    for(int ipt =0; ipt<npt; ipt++)
    for(int i =0; i<1; i++)
  //  for(int i =1; i<nCanvDiv-1; i++)
    {
//    c1->cd(i+1);
    
    TTree* t      = (TTree*)FileA->Get(Form("%sJetAnalyzer/t", AlgoNames[i]));
        TTree* tt2      = (TTree*)FileA->Get(Form("%sJetAnalyzer/t", AlgoNames[i+6]));
      
      leg[i]= new TLegend(0.22,0.63,0.45,0.90);//top right
      leg[i]->SetFillColor(0);
      leg[i]->SetTextSize(0.05);
      leg[i]->SetBorderSize(0);
      
        iHistoCorrEta[i][ipt] = new TH2D(Form("%s_CorrEtaInpt%d",AlgoNames[i],ipt),Form("%s_CorrEtaInpt%d",AlgoNames[i],ipt),nXbins,etabins,500,0,5);
        iHistoCorrEta[i+1][ipt] = new TH2D(Form("%s_CorrEtaInpt%d",AlgoNames[i+6],ipt),Form("%s_CorrEtaInpt%d",AlgoNames[i+6],ipt),nXbins,etabins,500,0,5);
        cout<<ipt<<endl;
     //   t->Draw(Form("corrpt/refpt:jteta>>%s_CorrEta",AlgoNames[i]),Form("corrpt>0 && fabs(jteta)<2.0"),"goff");
     //   t->Draw(Form("corrpt/refpt:jteta>>%s_CorrEta",AlgoNames[i]),"weight*(refpt>50)","goff");
 //       t->Draw(Form("corrpt/refpt:jteta>>%s_CorrEta",AlgoNames[i]),Form("(rawpt>15 && refpt>30 && fabs(jteta)<3.0)"),"goff");
       t->Draw(Form("corrpt/refpt:refeta>>%s_CorrEtaInpt%d",AlgoNames[i],ipt),Form("corrpt>0 && refpt>%d && refpt<=%d && fabs(jteta)<3.0 && rawpt>22.",ptbinlow[ipt],ptbinhigh[ipt]),"goff");
        tt2->Draw(Form("corrpt/refpt:refeta>>%s_CorrEtaInpt%d",AlgoNames[i+6],ipt),Form("corrpt>0 && refpt>%d && refpt<=%d && fabs(jteta)<3.0 && rawpt>22.",ptbinlow[ipt],ptbinhigh[ipt]),"same");

	outf->cd();
	iHistoCorrEta[i][ipt]->Write();
	iHistoCorrEta[i+1][ipt]->Write();
 /*       iHistoCorrEta[i]->Draw("colz");
        iHistoCorrEta[i+6]->Draw("same colz");
        leg[i]->AddEntry(iHistoCorrEta[i],Form("%s |#eta|<3.0",AlgoNames[i]),"");
        leg[i]->AddEntry(iHistoCorrEta[i+6],Form("%s |#eta|<3.0",AlgoNames[i+6]),"");
        leg[i]->AddEntry(iHistoCorrEta[i],Form("Corrected / ref p_{T}"),"");
        format1Dhisto(*iHistoCorrEta[i],250,-1,2,20,2,1,"Jet #eta","Corrected p_{T}/ ref p_{T}");//Red      
        format1Dhisto(*iHistoCorrEta[i+6],250,-1,2,20,2,1,"Jet #eta","Corrected p_{T}/ ref p_{T}");//Red
        leg[i]->Draw();
        for(int iX=iHistoCorrEta[i]->GetXaxis()->GetFirst(); iX<iHistoCorrEta[i]->GetXaxis()->GetLast(); iX++ )
        {
            //cout<<"bin ( "<<iX<<" ) --limits ["<<iHistoCorrEta[i]->GetXaxis()->GetBinLowEdge(iX)<<" , "<<iHistoCorrEta[i]->GetXaxis()->GetBinUpEdge(iX)<<" ] "<<endl;
            h1[i][iX] = new TH1D(Form("h1_%d_%d_py",i,iX),Form("h1_%d_%d_py",i,iX),500,0,5);
            h1[i][iX] =(TH1D*)iHistoCorrEta[i]->ProjectionY(Form("h1_%d_%d_py",i,iX),iX,iX,"e");
            Gauss6[i][iX] = new TF1(Form("F6_c%d_d%d",i,iX),"gaus",0,5);
            Gauss6[i][iX]->SetParLimits(1,0.7,1.5);
            h1[i][iX]->Fit(Form("F6_c%d_d%d",i,iX),"RQ0WL", "", 0.7, 1.5);
            
            h1[i+6][iX] = new TH1D(Form("h1_%d_%d_py",i+6,iX),Form("h1_%d_%d_py",i+6,iX),500,0,5);
            h1[i+6][iX] =(TH1D*)iHistoCorrEta[i+6]->ProjectionY(Form("h1_%d_%d_py",i+6,iX),iX,iX,"e");
            Gauss6[i+6][iX] = new TF1(Form("F6_c%d_d%d",i+6,iX),"gaus",0,5);
            Gauss6[i+6][iX]->SetParLimits(1,0.7,1.5);
            h1[i+6][iX]->Fit(Form("F6_c%d_d%d",i+6,iX),"RQ0WL", "NO", 0.7, 1.5);
     //   h1[i][iX]->Fit(Gauss6[i][iX],"0Q");
        if(DoGausFit){
	  mean[i][iX] = Gauss6[i][iX]->GetParameter(1);
	  sigma[i][iX] = Gauss6[i][iX]->GetParameter(2);
	  meanErr[i][iX] = Gauss6[i][iX]->GetParError(1);
	  sigmaErr[i][iX] = Gauss6[i][iX]->GetParError(2);
            mean[i+6][iX] = Gauss6[i+6][iX]->GetParameter(1);
            sigma[i+6][iX] = Gauss6[i+6][iX]->GetParameter(2);
            meanErr[i+6][iX] = Gauss6[i+6][iX]->GetParError(1);
            sigmaErr[i+6][iX] = Gauss6[i+6][iX]->GetParError(2);
        }
        else {
            mean[i][iX] = h1[i][iX]->GetMean();
            sigma[i][iX]= h1[i][iX]->GetRMS();
            meanErr[i][iX] = h1[i][iX]->GetMeanError();
            sigmaErr[i][iX] =h1[i][iX]->GetRMSError();  
            mean[i+6][iX] = h1[i+6][iX]->GetMean();
            sigma[i+6][iX]= h1[i+6][iX]->GetRMS();
            meanErr[i+6][iX] = h1[i+6][iX]->GetMeanError();
            sigmaErr[i+6][iX] =h1[i+6][iX]->GetRMSError();
        }
	  xPoint[i][iX] =iHistoCorrEta[i]->GetXaxis()->GetBinCenter(iX);
	  xPointErr[i][iX] = 0;
            xPoint[i+6][iX] =iHistoCorrEta[i+6]->GetXaxis()->GetBinCenter(iX);
            xPointErr[i+6][iX] = 0;
	  
        
        if(i==4 && iX==10)//Just to test 
	    {
            cout<<"bin ( "<<iX<<" ) --limits ["<<iHistoCorrEta[i]->GetXaxis()->GetBinLowEdge(iX)<<" , "<<iHistoCorrEta[i]->GetXaxis()->GetBinUpEdge(iX)<<" ] "<<endl;
	      TCanvas *c11 = new TCanvas("c11","c11",500,500);
	      c11->cd();	    
	      h1[i][iX]->Draw();
	      Gauss6[i][iX]->Draw("same");
	    }	  
	}

    }    
    if(SavePlot){
        c1->Print(Form("%s/Pbp502CorrectedEtaOverRefJet.pdf", plotsdir.Data()));
    }
    */
}
/*
    c1->Update();
    
  TLegend*  leg2B[cCanDiv];
  TLegend*  leg2A[cCanDiv];
  TGraphErrors *CorrPt_Mean[nBin];
    TCanvas *c2 = new TCanvas("c2","c2",10,10,kw,kh);
 //   makeMultiPanelCanvas(c2,nBin/4,2,0.0,0.0,0.1,0.2,0.05);
    makeMultiPanelCanvas(c2,nBin/4,1,0.0,0.0,0.1,0.2,0.05);
    
    //makeMultiPanelCanvas(c2,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c2B = new TCanvas("c2B","c2B",1100,700);
  //makeMultiPanelCanvas(c2B,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c2C = new TCanvas("c2C","c2C",1100,700);
  //makeMultiPanelCanvas(c2C,3,2,0.0,0.0,0.2,0.15,0.07);

  TGraphErrors *CorrPt_Sigma[nBin];
    TCanvas *c3 = new TCanvas("c3","c3",10,10,kw,kh);
  //  makeMultiPanelCanvas(c3,nBin/4,2,0.0,0.0,0.1,0.2,0.05);
    makeMultiPanelCanvas(c3,nBin/4,1,0.0,0.0,0.1,0.2,0.05);
    //makeMultiPanelCanvas(c3,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c3B = new TCanvas("c3B","c3B",1100,700);
  //makeMultiPanelCanvas(c3B,3,2,0.0,0.0,0.2,0.15,0.07);
  //TCanvas *c3C = new TCanvas("c3C","c3C",1100,700);
  //makeMultiPanelCanvas(c3C,3,2,0.0,0.0,0.2,0.15,0.07);

  for(int i2 =0; i2<nCanvDiv; i2++)
    {
      leg2A[i2]= new TLegend(0.42,0.65,0.65,0.9);//top right
      leg2A[i2]->SetFillColor(0);
      leg2A[i2]->SetTextSize(0.05);
      leg2A[i2]->SetBorderSize(0);
      
      leg2B[i2]= new TLegend(0.42,0.65,0.65,0.9);//top right
      leg2B[i2]->SetFillColor(0);
      leg2B[i2]->SetTextSize(0.05);
      leg2B[i2]->SetBorderSize(0);

      double xLoc[CnXbins];
      double y2Loc[CnXbins];
      double y3Loc[CnXbins];
      double xLocErr[CnXbins];
      double y2LocErr[CnXbins];
      double y3LocErr[CnXbins];
        double PuxLoc[CnXbins];
        double Puy2Loc[CnXbins];
        double Puy3Loc[CnXbins];
        double PuxLocErr[CnXbins];
        double Puy2LocErr[CnXbins];
        double Puy3LocErr[CnXbins];

        for(Int_t ipt = 0 ; ipt < nXbins ; ipt++){
                xLoc[ipt] = 0.;
                xLocErr[ipt] = 0.;
                y2Loc[ipt] = 0.;
                y2LocErr[ipt] = 0.;
                y3Loc[ipt] = 0.;
                y3LocErr[ipt] = 0.;
            PuxLoc[ipt] = 0.;
            PuxLocErr[ipt] = 0.;
            Puy2Loc[ipt] = 0.;
            Puy2LocErr[ipt] = 0.;
            Puy3Loc[ipt] = 0.;
            Puy3LocErr[ipt] = 0.;
            
            }
       for(int iBin=0; iBin<nXbins; iBin++)
	{	  
	  xLoc[iBin] = xPoint[i2][iBin];
	  xLocErr[iBin] = xPointErr[i2][iBin];
	  y2Loc[iBin] = mean[i2][iBin];
	  y2LocErr[iBin] = meanErr[i2][iBin];
	  y3Loc[iBin] = sigma[i2][iBin];
	  y3LocErr[iBin] = sigmaErr[i2][iBin];
        PuxLoc[iBin] = xPoint[i2+6][iBin];
        PuxLocErr[iBin] = xPointErr[i2+6][iBin];
        Puy2Loc[iBin] = mean[i2+6][iBin];
        Puy2LocErr[iBin] = meanErr[i2+6][iBin];
        Puy3Loc[iBin] = sigma[i2+6][iBin];
        Puy3LocErr[iBin] = sigmaErr[i2+6][iBin];
	}

//     if(i2<=5)c2->cd(i2+1);
//     if(i2>5 && i2<12)c2B->cd(i2-5);
//     if(i2>=12)c2C->cd(i2-12);
      c2->cd(i2+1);
      CorrPt_Mean[i2] = new TGraphErrors(nXbins,xLoc,y2Loc,xLocErr,y2LocErr);
        CorrPt_Mean[i2+6] = new TGraphErrors(nXbins,PuxLoc,Puy2Loc,PuxLocErr,Puy2LocErr);
      formatTGraph(*CorrPt_Mean[i2],250,-1,1,20,1,1,"Jet #eta","#mu ");
    formatTGraph(*CorrPt_Mean[i2+6],250,-1,1,20,1,1,"p_{T}^{ref} (GeV/c)","#mu ");
      leg2A[i2]->AddEntry(CorrPt_Mean[i2],Form("%s |#eta|<3.0",AlgoNames[i2]),"PL");
        leg2A[i2]->AddEntry(CorrPt_Mean[i2+6],Form("%s |#eta|<3.0",AlgoNames[i2+6]),"PL");
    leg2A[i2]->AddEntry(CorrPt_Mean[i2],Form("p_{T}^{ref} > 50 GeV/c"),"");
      leg2A[i2]->AddEntry(CorrPt_Mean[i2],Form("p_{T}^{corrected}/p_{T}^{ref}"),"");
     CorrPt_Mean[i2]->GetXaxis()->SetRangeUser(-3., 3.);
      CorrPt_Mean[i2]->SetMaximum(1.288);
      CorrPt_Mean[i2]->SetMinimum(0.84);
        CorrPt_Mean[i2]->SetMarkerStyle(20);
        CorrPt_Mean[i2]->SetMarkerColor(1);
        CorrPt_Mean[i2]->SetLineColor(1);
        CorrPt_Mean[i2]->SetMarkerSize(1.3);
        CorrPt_Mean[i2]->Draw("AP");
        CorrPt_Mean[i2+6]->SetMarkerStyle(24);
        CorrPt_Mean[i2+6]->SetMarkerColor(2);
        CorrPt_Mean[i2+6]->SetLineColor(2);
        CorrPt_Mean[i2+6]->SetMarkerSize(1.3);
        CorrPt_Mean[i2+6]->Draw("psame");
      leg2A[i2]->Draw();
        regSun(-4.,1.,4.,1.,1, 1);

//     if(i2<=5)c3->cd(i2+1);
//     if(i2>5 && i2<12)c3B->cd(i2-5);
//     if(i2>=12)c3C->cd(i2-12);
      c3->cd(i2+1);
      CorrPt_Sigma[i2] = new TGraphErrors(nXbins,xLoc,y3Loc,xLocErr,y3LocErr);
        CorrPt_Sigma[i2+6] = new TGraphErrors(nXbins,PuxLoc,Puy3Loc,PuxLocErr,Puy3LocErr);
      formatTGraph(*CorrPt_Sigma[i2],250,-1,1,20,1,1,"Jet #eta","#sigma ");
        formatTGraph(*CorrPt_Sigma[i2+6],250,-1,1,20,1,1,"p_{T}^{ref} (GeV/c)","#sigma ");
      leg2B[i2]->AddEntry(CorrPt_Sigma[i2],Form("%s |#eta|<3.0",AlgoNames[i2]),"PL");
    leg2B[i2]->AddEntry(CorrPt_Sigma[i2+6],Form("%s |#eta|<3.0",AlgoNames[i2+6]),"PL");
        leg2B[i2]->AddEntry(CorrPt_Mean[i2],Form("p_{T}^{ref} > 50 GeV/c"),"");
      leg2B[i2]->AddEntry(CorrPt_Sigma[i2],Form("p_{T}^{corrected}/p_{T}^{ref}"),"");
      CorrPt_Sigma[i2]->GetXaxis()->SetRangeUser(-3., 3.);
        CorrPt_Sigma[i2]->SetMaximum(0.43);
        CorrPt_Sigma[i2]->SetMinimum(0.0001);
        CorrPt_Sigma[i2]->SetMarkerStyle(20);
        CorrPt_Sigma[i2]->SetMarkerColor(1);
        CorrPt_Sigma[i2]->SetLineColor(1);
        CorrPt_Sigma[i2]->SetMarkerSize(1.3);
        CorrPt_Sigma[i2]->Draw("AP");
        CorrPt_Sigma[i2+6]->SetMarkerStyle(24);
        CorrPt_Sigma[i2+6]->SetMarkerColor(2);
        CorrPt_Sigma[i2+6]->SetLineColor(2);
        CorrPt_Sigma[i2+6]->SetMarkerSize(1.3);
        CorrPt_Sigma[i2+6]->Draw("psame");
      leg2B[i2]->Draw();
        regSun(-3.,1.,3.,1.,1, 1);

    }
    if(SavePlot){
        c2->Print(Form("%s/Pbp502CorrectedEtaOverRefJetMean.pdf", plotsdir.Data()));
        c3->Print(Form("%s/Pbp502CorrectedEtaOverRefJetSigma.pdf",plotsdir.Data())) ; 
    }
    c2->Update();
    c3->Update();
*/
}



///////////////////////////////////////////////////////////////////////
//         TOOL BOX
//////////////////////////////////////////////////////////////////////
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();

   TPad* pad[columns][rows];
   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth =
   (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
   (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
   (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }
   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
            Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}
void drawText(const char *text, float xp, float yp){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(17);
  //tex->SetTextSize(0.05);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
void regSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(3);
    t1->SetLineColor(color);
    t1->Draw();
}
