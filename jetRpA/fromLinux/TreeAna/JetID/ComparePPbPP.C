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
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

//TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
normalizeByBinWidth(rehisto);
}
if(binbound[Nbin]!=1000)
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
else 
rehisto->GetXaxis()->SetRangeUser(30,600);
fixedFontHist(rehisto,1.2,1.3);

return rehisto;
}

void ComparePPbPP(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
//gStyle->SetLabelFont(70);
/*for(int iptbin=0;iptbin<=500;iptbin++)
binbound_pt[iptbin]=iptbin*2;
int Nbin_pt=500;*/
//------------------------------------------------------------Get Histograms---------------------------------------------
const int Neta=8;
const TString etabinnameswap[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
TString etabinname[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
int ieta=7;


TH1F *hPP_pt= (TH1F*)fMCPPxSecJetID->Get(Form("jetptEtaBin%s",etabinname[ieta].Data()));
TH1F *hPP_pt_real= (TH1F*)fMCPPxSecJetID->Get(Form("jetptEtaBin_real%s",etabinname[ieta].Data()));
TH1F *hPP_pt_fake= (TH1F*)fMCPPxSecJetID->Get(Form("jetptEtaBin_fake%s",etabinname[ieta].Data()));
TH1F *hPPb_pt= (TH1F*)fMCPPbxSecJetID->Get(Form("jetptEtaBin%s",etabinnameswap[ieta].Data()));
TH1F *hPPb_pt_real= (TH1F*)fMCPPbxSecJetID->Get(Form("jetptEtaBin_real%s",etabinnameswap[ieta].Data()));
TH1F *hPPb_pt_fake= (TH1F*)fMCPPbxSecJetID->Get(Form("jetptEtaBin_fake%s",etabinnameswap[ieta].Data()));
TH1F *hPbP_pt= (TH1F*)fMCPbPJetIDOld->Get(Form("jetptEtaBin%s",etabinname[ieta].Data()));
TH1F *hPbP_pt_real= (TH1F*)fMCPbPJetIDOld->Get(Form("jetptEtaBin_real%s",etabinname[ieta].Data()));
TH1F *hPbP_pt_fake= (TH1F*)fMCPbPJetIDOld->Get(Form("jetptEtaBin_fake%s",etabinname[ieta].Data()));
TH1F *hhij_pt= (TH1F*)fhijing->Get(Form("jetptEtaBin%s",etabinname[ieta].Data()));
TH1F *hhij_pt_real= (TH1F*)fhijing->Get(Form("jetptEtaBin_real%s",etabinname[ieta].Data()));
TH1F *hhij_pt_fake= (TH1F*)fhijing->Get(Form("jetptEtaBin_fake%s",etabinname[ieta].Data()));

const double binbound_pt_coarse[]={0,30,70,150,200,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;

//---------------------------------------------------Rebin and Format histograms-----------------------------------------
TH1D* rehisto_hPP_pt=(TH1D*)myRebin(hPP_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_pt_real=(TH1D*)myRebin(hPP_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_pt_fake=(TH1D*)myRebin(hPP_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt=(TH1D*)myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_real=(TH1D*)myRebin(hPPb_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_fake=(TH1D*)myRebin(hPPb_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPbP_pt=(TH1D*)myRebin(hPbP_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPbP_pt_real=(TH1D*)myRebin(hPbP_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPbP_pt_fake=(TH1D*)myRebin(hPbP_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hhij_pt=(TH1D*)myRebin(hhij_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hhij_pt_real=(TH1D*)myRebin(hhij_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hhij_pt_fake=(TH1D*)myRebin(hhij_pt_fake,Nbin_pt,binbound_pt);


/*cout<<"low to high: "<<"PP real"<<"\t\t"<<"PP fake"<<"\t\t"<<"PP Total jets"<<"\t\t"<<"fake fraction"<<endl;
for(int i=1;i<rehisto_hProfPP_fake->GetXaxis()->GetNbins();i++){
if(rehisto_hPPb_pt->GetBinContent(i)!=0 && rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)>110 && rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)<190){
cout<<rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)<<" to "<<rehisto_hPP_pt_real->GetXaxis()->GetBinUpEdge(i)<<" : ";
//cout<<rehisto_hPPb_pt_real->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt_fake->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt_fake->GetBinContent(i)/rehisto_hPPb_pt->GetBinContent(i)<<endl;
cout<<rehisto_hPP_pt_real->GetBinContent(i)<<"\t\t"<<rehisto_hPP_pt_fake->GetBinContent(i)<<"\t\t"<<rehisto_hPP_pt->GetBinContent(i)<<"\t\t"<<rehisto_hPP_pt_fake->GetBinContent(i)/rehisto_hPP_pt->GetBinContent(i)<<endl;
}
}*/
rehisto_hPP_pt->SetMarkerStyle(24);
rehisto_hPP_pt->SetMarkerSize(1.2);
rehisto_hPP_pt->SetMarkerColor(2);
rehisto_hPP_pt->SetLineColor(2);
rehisto_hPP_pt_fake->SetMarkerStyle(24);
rehisto_hPP_pt_fake->SetMarkerSize(1.2);
rehisto_hPP_pt_fake->SetMarkerColor(2);
rehisto_hPP_pt_fake->SetLineColor(2);
rehisto_hPP_pt_real->SetMarkerStyle(24);
rehisto_hPP_pt_real->SetMarkerSize(1.2);
rehisto_hPP_pt_real->SetMarkerColor(2);
rehisto_hPP_pt_real->SetLineColor(2);


rehisto_hPPb_pt->SetMarkerStyle(20);
rehisto_hPPb_pt->SetMarkerSize(1.2);
rehisto_hPPb_pt->SetMarkerColor(1);
rehisto_hPPb_pt->SetLineColor(1);
rehisto_hPPb_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_pt_fake->SetMarkerSize(1.2);
rehisto_hPPb_pt_fake->SetMarkerColor(1);
rehisto_hPPb_pt_fake->SetLineColor(1);
rehisto_hPPb_pt_real->SetMarkerStyle(20);
rehisto_hPPb_pt_real->SetMarkerSize(1.2);
rehisto_hPPb_pt_real->SetMarkerColor(1);
rehisto_hPPb_pt_real->SetLineColor(1);

rehisto_hPbP_pt->SetMarkerStyle(29);
rehisto_hPbP_pt->SetMarkerSize(1.2);
rehisto_hPbP_pt->SetMarkerColor(4);
rehisto_hPbP_pt->SetLineColor(4);
rehisto_hPbP_pt_fake->SetMarkerStyle(29);
rehisto_hPbP_pt_fake->SetMarkerSize(1.2);
rehisto_hPbP_pt_fake->SetMarkerColor(4);
rehisto_hPbP_pt_fake->SetLineColor(4);
rehisto_hPbP_pt_real->SetMarkerStyle(29);
rehisto_hPbP_pt_real->SetMarkerSize(1.2);
rehisto_hPbP_pt_real->SetMarkerColor(4);
rehisto_hPbP_pt_real->SetLineColor(4);

rehisto_hhij_pt->SetMarkerStyle(29);
rehisto_hhij_pt->SetMarkerSize(1.2);
rehisto_hhij_pt->SetMarkerColor(4);
rehisto_hhij_pt->SetLineColor(4);
rehisto_hhij_pt_fake->SetMarkerStyle(29);
rehisto_hhij_pt_fake->SetMarkerSize(1.2);
rehisto_hhij_pt_fake->SetMarkerColor(4);
rehisto_hhij_pt_fake->SetLineColor(4);
rehisto_hhij_pt_real->SetMarkerStyle(29);
rehisto_hhij_pt_real->SetMarkerSize(1.2);
rehisto_hhij_pt_real->SetMarkerColor(4);
rehisto_hhij_pt_real->SetLineColor(4);


TH1F* hFrame = new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetRangeUser(32,599);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
fixedFontHist(hFrame,1.3,1.3);
hFrame->GetXaxis()->SetLabelSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.06);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetYaxis()->SetTitleSize(0.06);
//--------------------------------------------------------pT distribution------------------------------------------------

TCanvas* c1 = new TCanvas("c1"," ",800,600);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.18,0.15);
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);

c1->cd(1)->SetLogy();
hFrame->GetYaxis()->SetRangeUser(5e-13,1e-2);
hFrame->GetYaxis()->SetTitle("d#sigma/dp_{T}");
hFrame->DrawCopy();
rehisto_hPPb_pt->Draw("E1same");
rehisto_hPbP_pt->Draw("E1same");
rehisto_hPP_pt->Draw("E1same");
ratio_pt=(TH1D*)rehisto_hPPb_pt->Clone("ratio_pt");
ratio_pt->Divide(rehisto_hPP_pt);
ratio_pt_Pbp=(TH1D*)rehisto_hPbP_pt->Clone("ratio_pt_Pbp");
ratio_pt_Pbp->Divide(rehisto_hPP_pt);
TLegend *leg1=new TLegend(0.6,0.68,0.8,0.82);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(rehisto_hPP_pt,"PYTHIA","lpf");
leg1->AddEntry(rehisto_hPPb_pt,"PYTHIA+HIJING","lp");
leg1->AddEntry(rehisto_hPbP_pt,"HIJING+PYTHIA","lp");
leg1->Draw("same");
T1->DrawLatex(0.7,0.4,"Inclusive");

c1->cd(2)->SetLogy();
hFrame->DrawCopy();
rehisto_hPPb_pt_real->Draw("E1same");
rehisto_hPbP_pt_real->Draw("E1same");
rehisto_hPP_pt_real->Draw("E1same");
ratio_pt_real=(TH1D*)rehisto_hPPb_pt_real->Clone("ratio_pt_real");
ratio_pt_real_Pbp=(TH1D*)rehisto_hPbP_pt_real->Clone("ratio_pt_real_Pbp");
ratio_pt_real->Divide(rehisto_hPP_pt_real);
ratio_pt_real_Pbp->Divide(rehisto_hPP_pt_real);
T1->DrawLatex(0.7,0.4,"Real");
T1->DrawLatex(0.7,0.6,etastring[ieta]);

c1->cd(3)->SetLogy();
hFrame->GetYaxis()->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(5e-16,1e-2);
hFrame->DrawCopy("Y+");
rehisto_hPPb_pt_fake->Draw("E1same");
rehisto_hPbP_pt_fake->Draw("E1same");
rehisto_hPP_pt_fake->Draw("E1same");
ratio_pt_fake=(TH1D*)rehisto_hPPb_pt_fake->Clone("ratio_pt_fake");
ratio_pt_fake_Pbp=(TH1D*)rehisto_hPbP_pt_fake->Clone("ratio_pt_real_Pbp");
ratio_pt_fake->Divide(rehisto_hPP_pt_fake);
ratio_pt_fake_Pbp->Divide(rehisto_hPP_pt_fake);
T1->DrawLatex(0.7,0.4,"Fake");

TLine *l = new TLine(30,1,600,1);
l->SetLineStyle(2);
c1->cd(4);
hFrame->GetYaxis()->SetRangeUser(0,1.78);
hFrame->GetYaxis()->SetTitle("ratio");
hFrame->DrawCopy();
T1->DrawLatex(0.3,0.32,"PYTHIA+HIJING/PYTHIA");
ratio_pt->Draw("same");
ratio_pt_Pbp->Draw("same");
l->Draw("same");

c1->cd(5);
hFrame->DrawCopy();
ratio_pt_real->Draw("same");
ratio_pt_real_Pbp->Draw("same");
l->Draw("same");

c1->cd(6);
hFrame->GetYaxis()->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(0,1.78);
hFrame->DrawCopy("Y+");
ratio_pt_fake->Draw("same");
ratio_pt_fake_Pbp->Draw("same");
l->Draw("same");


//c1->Print(Form("pic/%s/jetpt_PPbvsPP.png",JetIDName.Data()));

}
