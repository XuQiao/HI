#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"
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
const TString etabinname[N]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etastring[N]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
TH2F* hMC[N];	 TProfile* hProMC[N]; TH1D* hMC1D[N];
TH2F* hdata[N];	 TProfile* hProdata[N]; TH1D* hdata1D[N];
TH1D* hRatio[N]; 
TString coll="PbP";

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
if(coll=="PPb"){
TFile *fMC=new TFile(MCPPbfile);
TFile *fData=new TFile(DataPPbfile);
}
else if(coll=="PbP"){
TFile *fMC=new TFile(MCPbPfileOld);
TFile *fData=new TFile(DataPbPfileOld);
}
for(int i=0;i<N;i++){
if(i==N-1){
hMC[i]=(TH2F*)fMC->Get(Form("rawptJES"));
hdata[i]=(TH2F*)fData->Get(Form("rawptJES"));
}
else{
hMC[i]=(TH2F*)fMC->Get(Form("rawptJESEtaBin%s",etabinname[i].Data()));
hdata[i]=(TH2F*)fData->Get(Form("rawptJESEtaBin%s",etabinname[i].Data()));
}
hMC[i]->SetName(Form("hMC_%d",i));
hdata[i]->SetName(Form("hdata_%d",i));
hProMC[i]=(TProfile*)hMC[i]->ProfileX();
hProdata[i]=(TProfile*)hdata[i]->ProfileX();
hProMC[i]=(TProfile*)hProMC[i]->Rebin(Nbin_pt,hProMC[i]->GetName(),binbound_pt);
hProdata[i]=(TProfile*)hProdata[i]->Rebin(Nbin_pt,hProdata[i]->GetName(),binbound_pt);
hMC1D[i]=(TH1D*)hProMC[i]->ProjectionX();
hdata1D[i]=(TH1D*)hProdata[i]->ProjectionX();
cout<<hMC1D[i]->GetBinContent(10)<<endl;
cout<<hdata1D[i]->GetBinContent(10)<<endl;
hRatio[i]=(TH1D*)hdata1D[i]->Clone(Form("hRatio_%d",i));
hRatio[i]->Divide(hMC1D[i]);
//cout<<hRatio[i]->GetBinContent(10)<<endl;
}
c2->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("JES Data/MC");
hFrame->GetYaxis()->SetRangeUser(0.95,1.05);
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
hProMC[i]->SetMarkerStyle(24);
hProdata[i]->SetMarkerStyle(20);
hProMC[i]->SetMarkerColor(2);
hProdata[i]->SetMarkerColor(1);
hProMC[i]->SetLineColor(2);
hProdata[i]->SetLineColor(1);
hProMC[i]->SetMarkerSize(1.2);
hProdata[i]->SetMarkerSize(1.2);
if(i==0){
leg1->AddEntry(hProMC[i],Form("PYTHIA+HIJING"),"lp");
leg1->AddEntry(hProdata[i],Form("Data"),"lp");
leg1->Draw("same");
}
hProdata[i]->DrawCopy("same");
hProMC[i]->DrawCopy("same");
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
c1->Print(Form("%sDataMCrawptJES_etabin.pdf",coll.Data()));
c2->Print(Form("%sDataMCrawptJESratio_etabin.pdf",coll.Data()));

}
