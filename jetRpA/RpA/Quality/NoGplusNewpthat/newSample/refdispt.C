#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

const double pPbLumi = 15.78 ;

TH1* handlehisto(TH1* histo){
histo->SetMarkerColor(colorCode[nColor-1]);
histo->SetMarkerStyle(markerCode[nColor-1]);
histo->SetMarkerSize(1.2);
histo->SetLineColor(colorCode[nColor-1]);
return histo;
}

void refdispt(){

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
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
TGraphAsymmErrors *g;
TGraphAsymmErrors *g0;
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
double jetptcut=29;

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.6);
hFrame->SetTitle("");
//hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(0,100);
hFrame->GetYaxis()->SetRangeUser(0,0.14);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.65,0.78,0.92);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);
TH2F* histo2F=(TH2F*)fPPb->Get(Form("refptpt_0-100%%"));

TH1F* histo1=(TH1F*)histo2F->ProjectionX("histo1",histo2F->GetYaxis()->FindBin(0),-1);
TH1F* histo2=(TH1F*)histo2F->ProjectionX("histo2",histo2F->GetYaxis()->FindBin(jetptcut),-1);
TH1F* histo3=(TH1F*)histo2F->ProjectionX("histo3",histo2F->GetYaxis()->FindBin(30),histo2F->GetYaxis()->FindBin(33));
TH1F* histo4=(TH1F*)histo2F->ProjectionY("histo4",histo2F->GetXaxis()->FindBin(18),histo2F->GetXaxis()->FindBin(22));
/*
histo1=(TH1F*)histo1->Rebin(Nbin_pt,histo1->GetName(),binbound_pt);
histo2=(TH1F*)histo2->Rebin(Nbin_pt,histo2->GetName(),binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
rehisto1->Scale(1/rehisto1->Integral());
rehisto2->Scale(1/rehisto2->Integral());
*/
histo3->Scale(1./histo3->Integral());
histo4->Scale(1./histo4->Integral());
//cout<<(double)histo3->Integral(histo3->GetXaxis()->FindBin(0),histo4->GetYaxis()->FindBin(20-1e-4))/histo3->Integral()<<endl;
cout<<(double)histo4->Integral(histo4->GetXaxis()->FindBin(0),histo4->GetXaxis()->FindBin(31-1e-4))/histo4->Integral()<<endl;
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->DrawLatex(0.35,0.88,"p_{T}^{hat}=15(GeV/c)    PYTHIA+HIJING");
T1->DrawLatex(0.55,0.4,"30<=p_{T}^{jet}<33 (GeV/c)");

ratio=(TH1F*)histo2->Clone("ratio");
ratio->Divide(histo1);
handlehisto(ratio);
handlehisto(histo3);
handlehisto(histo4);
histo3->GetXaxis()->SetRangeUser(0,100);
histo3->GetYaxis()->SetRangeUser(0,0.1);
//histo4->Fit("gaus","","",0,100);
//cout<<"gaus->Integral(0,30)/gaus->Integral(0,200)="<<(double)gaus->Integral(0,30)/gaus->Integral(0,200)<<endl;

//ratio->Draw("E1same");

g=new TGraphAsymmErrors();
g->BayesDivide(histo2,histo1);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
g0->SetMarkerColor(colorCode[nColor-1]);
g0->SetMarkerStyle(markerCode[nColor-1]);
g0->SetMarkerSize(1.2);
g0->SetLineColor(colorCode[nColor-1]);
//g0->Draw("Psame");
//leg1->AddEntry(ratio,"ratio","lp");
leg1->AddEntry(g0,"ratio","lp");
histo3->Draw("same");
//histo4->Draw("same");
TLine *l =new TLine(0,1,100,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T1->Draw("same");
gPad->RedrawAxis();
TString outname="refpt";
c1->Print(Form("%sPPb.png",outname.Data()));
c1->Print(Form("%sPPb.pdf",outname.Data()));

}
