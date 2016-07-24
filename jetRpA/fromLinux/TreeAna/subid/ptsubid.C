#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
void ptsubid(){
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.12);

const int N=8;
const double binbound_pt[] = { 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
const int Nbin_pt=39;

const TString etabinname[N]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};

TFile *fPPb = new TFile("/home/xuq/Documents/HI/RpA/output/TreeAna/subidcut.root");

const int Npoint=1000;
TH1F* ratio[N];
TGraphAsymmErrors *g[N];
TGraphAsymmErrors *g0[N];
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.6);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: subidcut/nocut");
hFrame->GetXaxis()->SetLimits(30,300);
hFrame->GetYaxis()->SetRangeUser(0.6,1.1);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.20,0.78,0.52);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
for(int i=0;i<N;i++){
if(i==N-1){
TH1F* histo1=(TH1F*)fPPb->Get(Form("jetpt"));
TH1F* histo2=(TH1F*)fPPb->Get(Form("subidjetpt"));
}
else{
TH1F* histo1=(TH1F*)fPPb->Get(Form("jetptEtaBin%s",etabinname[i].Data()));
TH1F* histo2=(TH1F*)fPPb->Get(Form("subidjetptEtaBin%s",etabinname[i].Data()));
}
//TH1D* histo2=(TH1D*)histo22F->ProjectionY("histo2",histo22F->GetXaxis()->FindBin(binbound_pt[5]),histo22F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]),"e");
//TH1D* histo1=(TH1D*)histo12F->ProjectionX("histo1");
//TH1D* histo2=(TH1D*)histo22F->ProjectionX("histo2");

histo1=(TH1F*)histo1->Rebin(Nbin_pt,histo1->GetName(),binbound_pt);
histo2=(TH1F*)histo2->Rebin(Nbin_pt,histo2->GetName(),binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
//rehisto1->Scale(1/rehisto1->Integral());
//rehisto2->Scale(1/rehisto2->Integral());

TLatex *T1=new TLatex(0.5,0.90,"PYTHIA+HIJING");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
/*
ratio[i]=(TH1F*)histo2->Clone("ratio");
ratio[i]->Divide(histo1);
ratio[i]->SetMarkerColor(colorCode[i]);
ratio[i]->SetMarkerStyle(markerCode[i]);
ratio[i]->SetMarkerSize(1.2);
ratio[i]->SetLineColor(colorCode[i]);
ratio[i]->Draw("E1same");
leg1->AddEntry(ratio[i],etastring[i],"lp");
*/

g[i]=new TGraphAsymmErrors();
g[i]->BayesDivide(histo2,histo1);
for(int ipoint=0;ipoint<g[i]->GetN();ipoint++){
x[ipoint]=g[i]->GetX()[ipoint];
y[ipoint]=g[i]->GetY()[ipoint];
yerrhigh[ipoint]=g[i]->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g[i]->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0[i]=new TGraphAsymmErrors(g[i]->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
g0[i]->SetMarkerColor(colorCode[i]);
g0[i]->SetMarkerStyle(markerCode[i]);
g0[i]->SetMarkerSize(1.2);
g0[i]->SetLineColor(colorCode[i]);
g0[i]->Draw("Psame");
leg1->AddEntry(g0[i],etastring[i],"lp");
TLine *l =new TLine(30,1,300,1);
l->SetLineStyle(2);
l->SetLineColor(1);
}
for(int i=0;i<N-1;i++)	g0[i]->Draw("Psame");
//for(int i=0;i<N-1;i++)	ratio[i]->Draw("Psame");
leg1->Draw("same");
l->Draw("same");
T1->Draw("same");
c1->Print("subidsysPPb.pdf");

}
