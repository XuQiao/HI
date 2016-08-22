#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
void CompareDataMC(){
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(kTRUE);
gStyle->SetErrorX(0);
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 24, 30, 28,27,20
};

const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12",""};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TH2F* hPPb[N];	 TProfile* hProPPb[N]; TH1D* hPPb1D[N];
TH2F* hdata[N];	 TProfile* hProdata[N]; TH1D* hdata1D[N];
TH1D* hRatio[N]; 

c1 = new TCanvas("c1"," ",1200,700);
c2 = new TCanvas("c2","c2",600,600);
makeMultiPanelCanvas(c2,1,1,-0.12,0,0.14,0.12,0.03);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.06);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
TLegend* leg1=new TLegend(0.6,0.70,0.78,0.92);
TLegend* leg2=new TLegend(0.6,0.65,0.78,0.95);
leg1->SetTextSize(0.06);
leg2->SetTextSize(0.04);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetRangeUser(10,599);
hFrame->SetTitle("");
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.06);
T->SetTextColor(1);
T->SetTextFont(42);
TLine *l=new TLine(10,1,599,1);
l->SetLineStyle(2);
l->SetLineColor(1);
for(int i=0;i<N;i++){
if(i==N-1){
hPPb[i]=(TH2F*)fPPb->Get("rawptJES_0-100%");
hdata[i]=(TH2F*)fdata->Get("rawptJESCombinedSpectraInEtaBin-10_10");
}
else{
hPPb[i]=(TH2F*)fPPb->Get(Form("rawptJESEtaBin%s_Cen0-100%%",etabinname[i].Data()));
hdata[i]=(TH2F*)fdata->Get(Form("rawptJESCombinedSpectraInEtaBin%s",etabinname[i].Data()));
}
hPPb[i]->SetName(Form("hPPb_%d",i));
hdata[i]->SetName(Form("hdata_%d",i));
hProPPb[i]=(TProfile*)hPPb[i]->ProfileX();
hProdata[i]=(TProfile*)hdata[i]->ProfileX();
hProPPb[i]=(TProfile*)hProPPb[i]->Rebin(Nbin_pt,hProPPb[i]->GetName(),binbound_pt);
hProdata[i]=(TProfile*)hProdata[i]->Rebin(Nbin_pt,hProdata[i]->GetName(),binbound_pt);
hPPb1D[i]=(TH1D*)hProPPb[i]->ProjectionX();
hdata1D[i]=(TH1D*)hProdata[i]->ProjectionX();
cout<<hPPb1D[i]->GetBinContent(10)<<endl;
cout<<hdata1D[i]->GetBinContent(10)<<endl;
hRatio[i]=(TH1D*)hdata1D[i]->Clone(Form("hRatio_%d",i));
hRatio[i]->Divide(hPPb1D[i]);
//cout<<hRatio[i]->GetBinContent(10)<<endl;
}
c2->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("JES Data/MC");
hFrame->GetYaxis()->SetRangeUser(0.99,1.01);
fixedFontHist(hFrame,1.2,1.6);
hFrame->DrawCopy();


for(int i=0;i<N;i++){
c1->cd(i+1);
if(i==0  || i==4)
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{raw}");
else
hFrame->GetYaxis()->SetTitle("");
if(i>=4)
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("");
fixedFontHist(hFrame,2.2,2.7);
hFrame->GetYaxis()->SetRangeUser(1.04,1.29);
hFrame->DrawCopy();
hProPPb[i]->SetMarkerStyle(24);
hProdata[i]->SetMarkerStyle(20);
hProPPb[i]->SetMarkerColor(2);
hProdata[i]->SetMarkerColor(1);
hProPPb[i]->SetLineColor(2);
hProdata[i]->SetLineColor(1);
hProPPb[i]->SetMarkerSize(1.2);
hProdata[i]->SetMarkerSize(1.2);
if(i==0){
leg1->AddEntry(hProPPb[i],Form("PYTHIA+HIJING"),"lp");
leg1->AddEntry(hProdata[i],Form("Data"),"lp");
leg1->Draw("same");
}
hProdata[i]->DrawCopy("same");
hProPPb[i]->DrawCopy("same");
T->DrawLatex(0.3,0.82,etastring[i]);
c2->cd(1);
hRatio[i]->SetMarkerStyle(markerCode[i]);
hRatio[i]->SetMarkerColor(colorCode[i]);
hRatio[i]->SetLineColor(colorCode[i]);
hRatio[i]->SetMarkerSize(1.2);
hRatio[i]->Draw("same");
leg2->AddEntry(hRatio[i],etastring[i],"lp");
}
for(int i=0;i<N-1;i++) 	hRatio[i]->Draw("same");

leg2->Draw("same");
l->Draw("same");
c1->Print(Form("DataMCrawptJES_etabin.pdf"));
c2->Print(Form("DataMCrawptJESratio_etabin.pdf"));

}
