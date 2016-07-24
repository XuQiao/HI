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
#include "root_setting.h"
#include "file.h"

bool isMC=kTRUE;
void WriteBins(TH1D *histo){
cout<<"points:"<<endl;
for(int i=1;i<histo->GetXaxis()->GetNbins()+1;i++){
cout<<i<<"th bin, from "<<histo->GetBinLowEdge(i)<<" to "<<histo->GetBinLowEdge(i)+histo->GetBinWidth(i)<<" content: "<<histo->GetBinContent(i)<<endl;
}
}

TH1D* myRebin(TH1D* histo,int Nbin, double binbound[]){

TH1D* rehisto=(TH1D*)histo->Clone(histo->GetName());
TH1D* rehisto=(TH1D*)histo->Rebin(Nbin,histo->GetName(),binbound);
if(binbound[Nbin]==1000)
rehisto->GetXaxis()->SetRangeUser(binbound[10],binbound[Nbin]);
else
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
for(int i=1;i<Nbin+1;i++){
rehisto->SetBinContent(i,rehisto->GetBinContent(i)/rehisto->GetBinWidth(i));
rehisto->SetBinError(i,rehisto->GetBinError(i)/rehisto->GetBinWidth(i));
}
if(!isMC)
rehisto->Scale(1/rehisto->Integral());
fixedFontHist(rehisto);
return rehisto;
}

void ptphi_test(){
using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);

TH2F *hdata2D= (TH2F*)fdata->Get("jetptphi_0-100%");
TH2F *hPP2D= (TH2F*)fPP->Get("jetptphi_0-100%");
TH2F *hPPb2D= (TH2F*)fPPb->Get("jetptphi_0-100%");
double xrange_phi[2]={-TMath::Pi(),TMath::Pi()};
const int Nbin_phi=20;
double binbound_phi[Nbin_phi+1];
for(int i=0;i<Nbin_phi+1;i++)
binbound_phi[i]=-TMath::Pi()+i*TMath::Pi()/10;

TH1D* hdata_phi=hdata2D->ProjectionY("hdata_phi",hdata2D->GetXaxis()->FindBin(xrange_pt[0]),hdata2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hdata_pt=hdata2D->ProjectionX("hdata_pt",hdata2D->GetYaxis()->FindBin(xrange_phi[0]),hdata2D->GetYaxis()->FindBin(xrange_phi[1])-1);
TH1D* hPPb_phi=hPPb2D->ProjectionY("hPPb_phi",hPPb2D->GetXaxis()->FindBin(xrange_pt[0]),hPPb2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb_pt=hPPb2D->ProjectionX("hPPb_pt",hPPb2D->GetYaxis()->FindBin(xrange_phi[0]),hPPb2D->GetYaxis()->FindBin(xrange_phi[1])-1);
TH1D* hPP_phi=hPP2D->ProjectionY("hPP_phi",hPP2D->GetXaxis()->FindBin(xrange_pt[0]),hPP2D->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPP_pt=hPP2D->ProjectionX("hPP_pt",hPP2D->GetYaxis()->FindBin(xrange_phi[0]),hPP2D->GetYaxis()->FindBin(xrange_phi[1])-1);

TH1D* rehisto_hdata_pt=myRebin(hdata_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPP_pt=myRebin(hPP_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt=myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hdata_phi=myRebin(hdata_phi,Nbin_phi,binbound_phi);
TH1D* rehisto_hPP_phi=myRebin(hPP_phi,Nbin_phi,binbound_phi);
TH1D* rehisto_hPPb_phi=myRebin(hPPb_phi,Nbin_phi,binbound_phi);


rehisto_hdata_pt->SetMarkerStyle(20);
rehisto_hdata_pt->SetMarkerSize(1);
rehisto_hdata_pt->SetMarkerColor(1);
if(!isMC){
rehisto_hPP_pt->SetFillStyle(3004);
rehisto_hPP_pt->SetMarkerStyle(0);
rehisto_hPP_pt->SetFillColor(2);
rehisto_hPP_pt->SetMarkerColor(2);
rehisto_hPP_pt->SetLineColor(2);
}
else{
rehisto_hPP_pt->SetMarkerStyle(24);
rehisto_hPP_pt->SetMarkerSize(1.2);
}
rehisto_hPP_pt->SetTitle(Form("#phi: %.3f--%.3f",xrange_phi[0],xrange_phi[1]));
rehisto_hPP_pt->GetYaxis()->SetTitle("Event Fraction");//#frac{d^{2}#sigma}{dp_{T}d#phi}");
rehisto_hPP_pt->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
rehisto_hPP_pt->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
rehisto_hPPb_pt->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
rehisto_hPPb_pt->SetMarkerStyle(29);
rehisto_hPPb_pt->SetMarkerSize(1.2);
rehisto_hPPb_pt->SetMarkerColor(4);
rehisto_hPPb_pt->SetLineColor(4);

if(!isMC){
rehisto_hdata_phi->SetMarkerStyle(20);
rehisto_hdata_phi->SetMarkerSize(1);
rehisto_hdata_phi->SetMarkerColor(1);
rehisto_hPP_phi->SetFillStyle(3004);
rehisto_hPP_phi->SetMarkerStyle(0);
rehisto_hPP_phi->SetFillColor(2);
rehisto_hPP_phi->SetMarkerColor(2);
rehisto_hPP_phi->SetLineColor(2);
}
else{
rehisto_hPP_phi->SetMarkerStyle(24);
rehisto_hPP_phi->SetMarkerSize(1.2);
}
rehisto_hPP_phi->SetTitle(Form("p_{T}: %3.f--%3.f [GeV/c]",xrange_pt[0],xrange_pt[1]));
rehisto_hPP_phi->GetYaxis()->SetTitle("Event Fraction");//#frac{d^{2}#sigma}{dp_{T}d#phi}");
rehisto_hPP_phi->GetXaxis()->SetTitle("#phi");
rehisto_hPP_phi->GetXaxis()->SetRangeUser(xrange_phi[0],xrange_phi[1]);
rehisto_hPPb_phi->GetXaxis()->SetRangeUser(xrange_phi[0],xrange_phi[1]);
rehisto_hPPb_phi->SetMarkerStyle(29);
rehisto_hPPb_phi->SetMarkerSize(1.2);
rehisto_hPPb_phi->SetMarkerColor(4);
rehisto_hPPb_phi->SetLineColor(4);
if(isMC){
rehisto_hPP_phi->SetMaximum(1e-4);
rehisto_hPP_phi->SetMinimum(0);
}
else{
rehisto_hPP_phi->SetMaximum(0.1);
rehisto_hPP_phi->SetMinimum(0.01);
}
c1 = new TCanvas("c1"," ",600,600);
c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c1,1,2,0.03,0.03,0.08,0.08,0.03);
makeMultiPanelCanvas(c2,1,2,0.03,0.03,0.08,0.08,0.03);
c1->cd(1)->SetLogy();
rehisto_hPP_pt->Draw("HIST E1");
if(!isMC)
rehisto_hdata_pt->Draw("E1same");
rehisto_hPPb_pt->Draw("E1same");
TLegend *leg1=new TLegend(0.7,0.8,0.9,0.9);
leg1->SetTextSize(0.03);
if(!isMC)
leg1->AddEntry(rehisto_hdata_pt,"Data","lp");
leg1->AddEntry(rehisto_hPP_pt,"PPMC","lpf");
leg1->AddEntry(rehisto_hPPb_pt,"PPbMC","lp");
leg1->Draw("same");
c1->cd(2);
if(isMC){
TH1D* ratio_pt=(TH1D*)rehisto_hPPb_pt->Clone("ratio_pt");
ratio_pt->Divide(rehisto_hPP_pt);
ratio_pt->SetTitle(Form("#phi: %.3f--%.3f",xrange_phi[0],xrange_phi[1]));
ratio_pt->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratio_pt->SetMarkerStyle(20);
ratio_pt->SetMarkerSize(1);
ratio_pt->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
ratio_pt->SetMinimum(0.8);
ratio_pt->SetMaximum(1.2);
ratio_pt->Draw("E1");
TLine *l1 =new TLine(xrange_pt[0],1,xrange_pt[1],1);
l1->SetLineStyle(2);
l1->SetLineColor(kBlack);
l1->Draw("same");
}
else{
TH1D* ratioppb_pt=(TH1D*)rehisto_hdata_pt->Clone();
TH1D* ratiopp_pt=(TH1D*)rehisto_hdata_pt->Clone();
ratioppb_pt->Divide(rehisto_hPPb_pt);
ratiopp_pt->Divide(rehisto_hPP_pt);
ratioppb_pt->SetTitle(Form("p_{T}: %.3f--%.3f [GeV/c]",xrange_phi[0],xrange_phi[1]));
ratioppb_pt->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratioppb_pt->SetMarkerStyle(20);
ratioppb_pt->SetMarkerColor(4);
ratiopp_pt->SetMarkerStyle(29);
ratiopp_pt->SetMarkerColor(2);
ratioppb_pt->SetMarkerSize(1.2);
ratiopp_pt->SetMarkerSize(1.2);
ratioppb_pt->GetXaxis()->SetRangeUser(xrange_pt[0],xrange_pt[1]);
ratioppb_pt->SetMinimum(0.5);
ratioppb_pt->SetMaximum(1.5);
ratioppb_pt->Draw();
ratiopp_pt->Draw("same");
TLine *l2 =new TLine(xrange_pt[0],1,xrange_pt[1],1);
l2->SetLineStyle(2);
l2->SetLineColor(kBlack);
l2->Draw("same");
}

//c1->Print(Form("pic/test/jetptphi_ProX%.3f-%.3f_%3.f-%.3f_PPbvsPP.png",xrange_phi[0],xrange_phi[1],xrange_pt[0],xrange_pt[1]));
c2->cd(1);
rehisto_hPP_phi->Draw("HIST E1");
if(!isMC)
rehisto_hdata_phi->Draw("E1same");
rehisto_hPPb_phi->Draw("E1same");
TLegend *leg2=new TLegend(0.7,0.8,0.9,0.9);
leg2->SetTextSize(0.03);
if(!isMC)
leg2->AddEntry(rehisto_hdata_phi,"Data","lp");
leg2->AddEntry(rehisto_hPP_phi,"PPMC","lpf");
leg2->AddEntry(rehisto_hPPb_phi,"PPbMC","lp");
leg2->Draw("same");
c2->cd(2);
if(isMC){
TH1D* ratio_phi=(TH1D*)rehisto_hPPb_phi->Clone();
ratio_phi->Divide(rehisto_hPP_phi);
ratio_phi->SetTitle(Form("p_{T}: %.3f--%.3f [GeV/c]",xrange_pt[0],xrange_pt[1]));
ratio_phi->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratio_phi->SetMarkerStyle(20);
ratio_phi->SetMarkerSize(1);
ratio_phi->GetXaxis()->SetRangeUser(xrange_phi[0],xrange_phi[1]);
ratio_phi->SetMinimum(0.8);
ratio_phi->SetMaximum(1.2);
ratio_phi->Draw();
}
else{
TH1D* ratioppb_phi=(TH1D*)rehisto_hdata_phi->Clone();
TH1D* ratiopp_phi=(TH1D*)rehisto_hdata_phi->Clone();
ratioppb_phi->Divide(rehisto_hPPb_phi);
ratiopp_phi->Divide(rehisto_hPP_phi);
ratioppb_phi->SetTitle(Form("p_{T}: %.3f--%.3f [GeV/c]",xrange_pt[0],xrange_pt[1]));
ratioppb_phi->GetYaxis()->SetTitle("#frac{PPb}{PP} ratio");
ratioppb_phi->SetMarkerStyle(20);
ratioppb_phi->SetMarkerColor(4);
ratiopp_phi->SetMarkerStyle(30);
ratiopp_phi->SetMarkerColor(2);
ratioppb_phi->SetMarkerSize(1.2);
ratiopp_phi->SetMarkerSize(1.2);
ratioppb_phi->GetXaxis()->SetRangeUser(xrange_phi[0],xrange_phi[1]);
ratioppb_phi->SetMinimum(0.8);
ratioppb_phi->SetMaximum(1.2);
ratioppb_phi->Draw();
ratiopp_phi->Draw("same");
TLine *l2 =new TLine(xrange_phi[0],1,xrange_phi[1],1);
l2->SetLineStyle(2);
l2->SetLineColor(kBlack);
l2->Draw("same");
//c2->Print(Form("pic/test/ratio_jetptphi_ProY%.3f-%.3f_%.3f-%.3f_PPbvsPP.png",xrange_phi[0],xrange_phi[1],xrange_pt[0],xrange_pt[1]));
}
}
