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

void Comparejetidclosure(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------
int isMatch=0;
bool Save = 1;
//--------------------------------------------------------SetUp canvas---------------------------------------------------
c1 = new TCanvas("c1"," ",1200,600);
c2 = new TCanvas("c2"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
makeMultiPanelCanvas(c2,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame1=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame1,2.0,3.0);
hFrame1->SetTitle("");
TLegend *leg1=new TLegend(0.26,0.52,0.60,0.85);
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
if(i==7){
	TH1F *hjet_pt_real= (TH1F*)fMCPPb->Get(Form("jetpt_real"));
	TH1F *hjet_pt_jetidclosure= (TH1F*)fMCPPb->Get(Form("jetpt_jetidclosure"));
}
else{
	TH1F *hjet_pt_real= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_real",etabinnameswap[i].Data()));
	TH1F *hjet_pt_jetidclosure= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_jetidclosure",etabinnameswap[i].Data()));
}
        hjet_pt_real = (TH1F*)hjet_pt_real->Rebin(Nbin_pt,hjet_pt_real->GetName(),binbound_pt);
        normalizeByBinWidth(hjet_pt_real);
        hjet_pt_jetidclosure = (TH1F*)hjet_pt_jetidclosure->Rebin(Nbin_pt,hjet_pt_jetidclosure->GetName(),binbound_pt);
        normalizeByBinWidth(hjet_pt_jetidclosure);
        hratio = (TH1F*)hjet_pt_jetidclosure->Clone(Form("jetidclosure_%d",i));
        hratio->Divide(hjet_pt_real);
    hjet_pt_real->SetMarkerSize(1.3);
    hjet_pt_real->SetMarkerColor(1);
    hjet_pt_real->SetLineColor(1);
    hjet_pt_real->SetMarkerStyle(20);
    hjet_pt_jetidclosure->SetMarkerSize(1.3);
    hjet_pt_jetidclosure->SetMarkerColor(2);
    hjet_pt_jetidclosure->SetLineColor(2);
    hjet_pt_jetidclosure->SetMarkerStyle(24);
    hratio->SetMarkerStyle(20);
    hratio->SetMarkerColor(1);
    hratio->SetLineColor(1);
    getRidYError(hratio);

c1->cd(canvas[i]+1)->SetGridx();
c1->cd(canvas[i]+1)->SetLogy();
    if(canvas[i]==0  || canvas[i]==4){
        hFrame1->GetYaxis()->SetTitle("jet Yield");
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
   hFrame1->SetMinimum(1e-13);
   hFrame1->SetMaximum(1e-2);
   hFrame1->DrawCopy();

if(canvas[i]!=4){
hjet_pt_real->DrawCopy("Psame");
hjet_pt_jetidclosure->DrawCopy("Psame");
l->Draw("same");
}
if(canvas[i]==4){
    leg1->AddEntry(hjet_pt_real,"Real before cut","lp");
    leg1->AddEntry(hjet_pt_jetidclosure,"Inclusive after cut & corr","lp");
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
//hratio->DrawCopy("Psame");
//hcorrPbP->DrawCopy("Psame");
hratio->Draw("Psame");
T1.SetTextSize(0.065);
T1.DrawLatex(0.4,ybase,etastring[i]);
l->Draw("same");
}
if(canvas[i]==4){
    T1.DrawLatex(0.4,0.8,"cut 3 closure");
//    TLegend *leg2 = (TLegend*)leg1;
//    leg2->Clear();
//    leg2->AddEntry(hratio,"Ratio","lp");
//    leg2->Draw("same");
}
}

if(Save){
if(isMatch){
c1->Print(Form("pic/JetIDclosure_cut1.png"));
c1->Print(Form("pic/JetIDclosure_cut1.pdf"));
c2->Print(Form("pic/JetIDclosureratio_cut1.png"));
c2->Print(Form("pic/JetIDclosureratio_cut1.pdf"));
}
else{
c1->Print(Form("pic/JetIDclosure.png"));
c1->Print(Form("pic/JetIDclosure.pdf"));
c2->Print(Form("pic/JetIDclosureratio.png"));
c2->Print(Form("pic/JetIDclosureratio.pdf"));
}
}


}
