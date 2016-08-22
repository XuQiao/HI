#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

void ptrawptcut(){
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};
const double pPbLumi = 15.78 ;

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
    c1->SetBottomMargin(0.1);
TString gen="PPb";
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
TGraphAsymmErrors *g[2];
TGraphAsymmErrors *g0[2];
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
double rawptcut=22;
TH1F* ratio[2];
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1,1.4);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle(Form("Ratio: rawpt>=%.f(GeV/c)/nocut",rawptcut));
hFrame->GetXaxis()->SetLimits(1,100);
hFrame->GetYaxis()->SetRangeUser(0,1.1);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.65,0.78,0.75);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);

TLatex *T1=new TLatex(0.25,0.24,"Anti-k_{T} Particle Flow Jets R=0.3 with Pile-Up subtraction");
TLine *l =new TLine(1,1,100,1);
TLine *l0 =new TLine(29,0,29,1.1);

TString datatype;

for(int j=0;j<2;j++){
if(j==0){
TH2F* histo2F=(TH2F*)fPPb->Get(Form("ptrawpt_0-100%%"));
int colorst=2;
int markerst=24;
datatype="PYTHIA+HIJING";
}
else if(j==1){
TH2F* histo2F=(TH2F*)fdata->Get(Form("ptrawpt_0-100%%"));
int colorst=1;
int markerst=20;
datatype="Data";
}

TH1F* histo1=(TH1F*)histo2F->ProjectionX("histo1",histo2F->GetYaxis()->FindBin(0),-1);
TH1F* histo2=(TH1F*)histo2F->ProjectionX("histo2",histo2F->GetYaxis()->FindBin(rawptcut),-1);
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.03);
T1->SetTextColor(1);
T1->SetTextFont(42);

g[j]=new TGraphAsymmErrors();
g[j]->BayesDivide(histo2,histo1);
for(int ipoint=0;ipoint<g[j]->GetN();ipoint++){
x[ipoint]=g[j]->GetX()[ipoint];
y[ipoint]=g[j]->GetY()[ipoint];
yerrhigh[ipoint]=g[j]->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g[j]->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}

g0[j]=new TGraphAsymmErrors(g[j]->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
g0[j]->SetMarkerColor(colorst);
g0[j]->SetMarkerStyle(markerst);
g0[j]->SetMarkerSize(1.1);
g0[j]->SetLineColor(colorst);
g0[j]->Draw("Psame");
//leg1->AddEntry(ratio,"ratio","lp");
leg1->AddEntry(g0[j],Form("%s",datatype.Data()),"lp");
}
leg1->Draw("same");
l->SetLineStyle(2);
l0->SetLineStyle(2);
l->SetLineColor(1);
l0->SetLineColor(4);
l->Draw("same");
l0->Draw("same");
T1->Draw("same");
c1->Print("pic/centWeight/rawptcutsysPPbvsData.png");
c1->Print("pic/centWeight/rawptcutsysPPbvsData.pdf");

}
