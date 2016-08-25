#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include "iostream"
#include "fstream"
#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"

void ComparePPbPbPtemp(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------
int isMatch=0;
bool Save = 0;
TFile *ftempPPb = new TFile("ftempPPb.root","ReadOnly");
TFile *ftempPbP = new TFile("ftempPbP.root","ReadOnly");
TFile *fjetidcorrPPb = new TFile("fjetidcorrPPb.root","ReadOnly");
TFile *fjetidcorrPbP = new TFile("fjetidcorrPbP.root","ReadOnly");

//--------------------------------------------------------SetUp canvas---------------------------------------------------
c1 = new TCanvas("c1"," ",1200,600);
c2 = new TCanvas("c2"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
makeMultiPanelCanvas(c2,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame1=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame1,2.0,3.0);
hFrame1->SetTitle("");
TLegend *leg1=new TLegend(0.36,0.62,0.60,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.055);
TLatex T1;
T1.SetNDC();
T1.SetTextSize(0.065);
T1.SetTextFont(42);
TLine *l = new TLine(hFrame1->GetXaxis()->GetXmin(),1,hFrame1->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);

for(int i=0;i<Neta;i++){
	if(canvas[i]==0){
		double ybase=0.15;	double xbase=0.28;
	}
	else if(canvas[i]>4){
		double ybase=0.30;	double xbase=0.05;
	}
	else{
		double ybase=0.15;	double xbase=0.05;
	}
if(canvas[i]!=4){
	TH1F *h1PPb= (TH1F*)ftempPPb->Get(Form("ratio_hPPb_JetIDcutvsnocut_pt_real_%s_%d",etabinnameswap[i].Data(),isMatch));
	TH1F *h2PPb= (TH1F*)ftempPPb->Get(Form("ratio_hPPb_JetIDcut_pt_realvsInc_%s_%d",etabinnameswap[i].Data(),isMatch));
	TH1F *h1PbP= (TH1F*)ftempPbP->Get(Form("ratio_hPbP_JetIDcutvsnocut_pt_real_%s_%d",etabinname[i].Data(),isMatch));
	TH1F *h2PbP= (TH1F*)ftempPbP->Get(Form("ratio_hPbP_JetIDcut_pt_realvsInc_%s_%d",etabinname[i].Data(),isMatch));
	TH1F *hcorrPPb= (TH1F*)fjetidcorrPPb->Get(Form("jetidcorr_%s",etabinnameswap[i].Data()));
	TH1F *hcorrPbP= (TH1F*)fjetidcorrPbP->Get(Form("jetidcorr_%s",etabinname[i].Data()));
        TH1F* hratio = (TH1F*)hcorrPPb->Clone(Form("corr_%d",i));
        hratio->Divide(hcorrPbP);
    h1PPb->SetMarkerSize(1.3);
    h1PbP->SetMarkerSize(1.3);
    h2PPb->SetMarkerSize(1.3);
    h2PbP->SetMarkerSize(1.3);
    hcorrPPb->SetMarkerSize(1.3);
    hcorrPbP->SetMarkerSize(1.3);
    h1PPb->SetMarkerColor(1);
    h1PPb->SetLineColor(1);
    h1PPb->SetMarkerStyle(20);
    h1PbP->SetMarkerColor(2);
    h1PbP->SetLineColor(2);
    h1PbP->SetMarkerStyle(24);
    h2PPb->SetMarkerColor(4);
    h2PPb->SetLineColor(4);
    h2PPb->SetMarkerStyle(33);
    h2PbP->SetMarkerColor(6);
    h2PbP->SetLineColor(6);
    h2PbP->SetMarkerStyle(27);
    hcorrPPb->SetMarkerStyle(20);
    hcorrPPb->SetMarkerColor(1);
    hcorrPPb->SetLineColor(1);
    hcorrPbP->SetMarkerStyle(24);
    hcorrPbP->SetMarkerColor(2);
    hcorrPbP->SetLineColor(2);
}

c1->cd(canvas[i]+1)->SetGridx();
    if(canvas[i]==0  || canvas[i]==4){
        hFrame1->GetYaxis()->SetTitle("Ratio");
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
        hFrame1->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
	
   hFrame1->GetXaxis()->SetNdivisions(510);
   hFrame1->GetXaxis()->SetRangeUser(32,597);
   hFrame1->SetMinimum(0.902);
   hFrame1->SetMaximum(1.012);
   hFrame1->DrawCopy();

if(canvas[i]!=4){
h1PPb->DrawCopy("Psame");
h2PPb->DrawCopy("Psame");
h1PbP->DrawCopy("Psame");
h2PbP->DrawCopy("Psame");
T1.SetTextSize(0.065);
T1.DrawLatex(0.4,ybase,etastring[i]);
l->Draw("same");
}
if(canvas[i]==4){
    T1.DrawLatex(0.4,0.8,"cut 3:");
    leg1->AddEntry(h1PPb,"cut eff, PPb","lp");
    leg1->AddEntry(h1PbP,"cut eff, PbP","lp");
    leg1->AddEntry(h2PPb,"cut cont, PPb","lp");
    leg1->AddEntry(h2PbP,"cut cont, PbP","lp");
    leg1->Draw("same");
}


c2->cd(canvas[i]+1)->SetGridx();
    if(canvas[i]==0  || canvas[i]==4){
        hFrame1->GetYaxis()->SetTitle("Ratio");
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
        hFrame1->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
	
   hFrame1->GetXaxis()->SetNdivisions(510);
   hFrame1->GetXaxis()->SetRangeUser(32,597);
   hFrame1->SetMinimum(0.902);
   hFrame1->SetMaximum(1.042);
   hFrame1->DrawCopy();

if(canvas[i]!=4){
//hcorrPPb->DrawCopy("Psame");
//hcorrPbP->DrawCopy("Psame");
hratio->Draw("Psame");
T1.SetTextSize(0.065);
T1.DrawLatex(0.4,ybase,etastring[i]);
l->Draw("same");
}
if(canvas[i]==4){
    T1.DrawLatex(0.4,0.8,"cut 3:");
    TLegend *leg2 = (TLegend*)leg1;
    leg2->Clear();
    leg2->AddEntry(hcorrPPb,"Corr factor, PPb","lp");
    leg2->AddEntry(hcorrPbP,"Corr factor, PbP","lp");
    leg2->Draw("same");
}
}

if(Save){
if(isMatch){
c1->Print(Form("pic/ComparePPbPbPtemp.png"));
c1->Print(Form("pic/ComparePPbPbPcorr.png"));
}
else{
c1->Print(Form("pic/ComparePPbPbPtemp.png"));
c1->Print(Form("pic/ComparePPbPbPcorr.png"));
}
}


}
