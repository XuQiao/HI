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

void recopthatdis(){

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
int pts=18,ptl=22;
const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12",""};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.6);
hFrame->SetTitle("");
//hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
//hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetYaxis()->SetTitle("Event Probability");
hFrame->GetXaxis()->SetLimits(10,60);
hFrame->GetYaxis()->SetRangeUser(0,6e-3);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.5,0.78,0.6);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->DrawLatex(0.2,0.88,"p_{T}^{hat}=15(GeV/c)      PYTHIA+HIJING");
T1->DrawLatex(0.2,0.83,"CMS simulation");
T1->DrawLatex(0.55,0.25,Form("%d<=p_{T}^{ref}<%d (GeV/c)",pts,ptl));

TH2F* histo2F[N]; TH1F* histo[N];	TH2F* histo2Fnocut;	TH1F* histonocut;
for(int i=0;i<1;i++){
if(i==0)
histo2F[i]=(TH2F*)fPPb->Get(Form("refptpt_0-100%%"));
else
histo2F[i]=(TH2F*)fPPb->Get(Form("refptptEtaBin%s_Cen0-100%%",etabinname[i].Data()));
histo2Fnocut=(TH2F*)fPPb->Get(Form("refptptnocut_0-100%%"));

histo[i]=(TH1F*)histo2F[i]->ProjectionY(Form("histo_%d",i),histo2F[i]->GetXaxis()->FindBin(pts),histo2F[i]->GetXaxis()->FindBin(ptl));
histonocut=(TH1F*)histo2Fnocut->ProjectionY(Form("histonocut"),histo2Fnocut->GetXaxis()->FindBin(pts),histo2Fnocut->GetXaxis()->FindBin(ptl));
//histo[i]->Scale(1./histo[i]->Integral());
histo[i]->SetMarkerColor(colorCode[7]);
histonocut->SetMarkerColor(colorCode[6]);
histo[i]->SetLineColor(colorCode[7]);
histonocut->SetLineColor(colorCode[6]);
histo[i]->SetMarkerSize(1.2);
histonocut->SetMarkerSize(1.2);
histo[i]->SetMarkerStyle(markerCode[7]);
histonocut->SetMarkerStyle(markerCode[6]);
histo[i]->Draw("same");
histonocut->Draw("same");
//leg1->AddEntry(histo[i],"|#eta| < 3","lp");
leg1->AddEntry(histo[i],"rawpt cut >22","lp");
leg1->AddEntry(histonocut,"no cut","lp");
}
leg1->Draw("same");
//for(int i=0;i<1;i++) if(i!=3) histo[i]->Draw("same"); 
T1->Draw("same");
gPad->RedrawAxis();
TString outname="recopthatdis";
c1->Print(Form("%sPPb.png",outname.Data()));
c1->Print(Form("%sPPb.pdf",outname.Data()));

}
