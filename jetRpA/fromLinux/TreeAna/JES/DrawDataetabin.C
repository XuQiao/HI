#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"
void DrawDataetabin(){
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(kFALSE);
gStyle->SetErrorX(0);
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

const int N=8;
TString coll="PPb";
const TString etabinname[N]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etabinnameswap[N]={"-20_-15","-15_-10","-10_-5","-5_5","5_10","10_15","15_20","-10_10"};
//const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
//ofstream fstr("JESsysetabin.txt");
TH2F* hrawptJES[N];	 TProfile* hProrawptJES[N]; 
c1 = new TCanvas("c1","c1",600,600);
    c1->SetTickx(1);
    c1->SetTicky(1);
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.14,0.12,0.03);
if(coll=="PPb"){
TFile *fData=new TFile(DataPPbfile);
}
else if(coll=="PbP"){
TFile *fData=new TFile(DataPbPfileOld);
}

TLegend *leg1=new TLegend(0.72,0.55,0.88,0.9);
TLegend *leg2=new TLegend(0.72,0.55,0.88,0.9);
leg1->SetTextSize(0.035);
leg2->SetTextSize(0.035);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(24.1,599.9);
hFrame->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(1.0,1.35);
fixedFontHist(hFrame,1.2,1.8);
TLatex *T1=new TLatex;
T1->SetNDC();
T1->SetTextAlign(10);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
TLine *l=new TLine(31,1,599.9,1);
l->SetLineStyle(2);
l->SetLineColor(1);
for(int i=0;i<N;i++){
if(i==N-1){
hrawptJES[i]=(TH2F*)fData->Get("rawptJES");
}
else{
if(coll=="PPb")
hrawptJES[i]=(TH2F*)fData->Get(Form("rawptJESEtaBin%s",etabinname[i].Data()));
else if(coll=="PbP")
hrawptJES[i]=(TH2F*)fData->Get(Form("rawptJESEtaBin%s",etabinnameswap[i].Data()));
}
hrawptJES[i]->SetName(Form("hrawptJES_%d",i));
hProrawptJES[i]=(TProfile*)hrawptJES[i]->ProfileX();
hProrawptJES[i]=(TProfile*)hProrawptJES[i]->Rebin(Nbin_pt,hProrawptJES[i]->GetName(),binbound_pt);
}

c1->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{raw}");
hFrame->DrawCopy();
for(int i=0;i<N;i++){
hProrawptJES[i]->SetMarkerStyle(markerCode[i]);
hProrawptJES[i]->SetMarkerColor(colorCode[i]);
hProrawptJES[i]->SetLineColor(colorCode[i]);
hProrawptJES[i]->SetMarkerSize(1.4);
leg1->AddEntry(hProrawptJES[i],etastring[i],"lp");
hProrawptJES[i]->DrawCopy("same");
}
for(int i=0;i<N-1;i++)	if(i!=3) hProrawptJES[i]->DrawCopy("same");
if(coll=="PPb"){
drawCMS(0.22,0.90,20.7);
T1->DrawLatex(0.22,0.85,"Pb going positive Data");}
else if(coll=="PbP"){
drawCMS(0.22,0.90,14);
T1->DrawLatex(0.22,0.85,"Proton going positive Data");}
T1->DrawLatex(0.28,0.15,"Anti-k_{T} Particle Flow PU Jets: R=0.3");
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->Draw("same");
//drawCMS(0.4,0.9);
c1->Print(Form("Data%srawptJES_etabin.pdf",coll.Data()));

}
