#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

void resolution(){

double binbound_pt[101];
for(int ipt=0;ipt<=100;ipt++)
binbound_pt[ipt]=ipt*10;
int Nbin_pt=100;
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

TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.6);
hFrame->SetTitle("");
//hFrame->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
hFrame->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
hFrame->GetYaxis()->SetTitle(Form("Resolution"));
//hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetXaxis()->SetLimits(30,600);
hFrame->GetYaxis()->SetRangeUser(0,0.25);
hFrame->DrawCopy();
TLegend* leg1=new TLegend(0.6,0.65,0.78,0.92);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);


TH2F* histoPPb2F=(TH2F*)fPPb->Get(Form("refptpt_0-100%%"));
TH2F* histoPP2F=(TH2F*)fPP->Get(Form("refptpt_0-100%%"));
TLatex *T1=new TLatex();
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextFont(42);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
TGraph *gPPb; TGraph *gPP; TGraph *gRel;
int N=histoPPb2F->GetNbinsX();
double sigPPb[1000],sigPP[1000],refpt[1000];
double meanPPb[1000],meanPP[1000],rel[1000];
for(int i=0;i<Nbin_pt;i++){
refpt[i]=histoPPb2F->GetXaxis()->GetBinLowEdge(histoPPb2F->GetXaxis()->FindBin(binbound_pt[i]));
TH1F* histoPPb=(TH1F*)histoPPb2F->ProjectionY("histoPPb",histoPPb2F->GetXaxis()->FindBin(binbound_pt[i]),histoPPb2F->GetXaxis()->FindBin(binbound_pt[i+1]));
normalizeByBinWidth(histoPPb);
TH1F* histoPP=(TH1F*)histoPP2F->ProjectionY("histoPP",histoPP2F->GetXaxis()->FindBin(binbound_pt[i]),histoPP2F->GetXaxis()->FindBin(binbound_pt[i+1]));
normalizeByBinWidth(histoPP);
sigPPb[i] = histoPPb->GetRMS();
meanPPb[i] = histoPPb->GetMean();
sigPPb[i] /= meanPPb[i];

sigPP[i] = histoPP->GetRMS();
meanPP[i] = histoPP->GetMean();
sigPP[i] /= meanPP[i];

rel[i] = TMath::Sqrt(TMath::Abs(sigPPb[i]**2-sigPP[i]**2));
}
gPPb=new TGraph(N,refpt,sigPPb);
gPP=new TGraph(N,refpt,sigPP);
gRel=new TGraph(N,refpt,rel);

gPPb->SetMarkerColor(1);
gPP->SetMarkerColor(2);
gRel->SetMarkerColor(4);
gPPb->SetMarkerStyle(20);
gPP->SetMarkerStyle(24);
gRel->SetMarkerStyle(30);
gPPb->SetMarkerSize(1.2);
gPP->SetMarkerSize(1.2);
gRel->SetMarkerSize(1.2);
gPPb->SetLineColor(1);
gPPb->SetLineColor(2);
gRel->SetLineColor(4);
//gPPb->Draw("Psame");
//gPP->Draw("Psame");
gRel->Draw("Psame");
leg1->AddEntry(gPPb,"PYTHIA+HIJING","lp");
leg1->AddEntry(gPP,"PYTHIA","lp");
//leg1->Draw("same");
TLine *l =new TLine(30,1,600,1);
TLine *l0 =new TLine(45,0,45,1.1);
l->SetLineStyle(2);
l0->SetLineStyle(2);
l->SetLineColor(1);
l0->SetLineColor(4);
l->Draw("same");
//l0->Draw("same");
T1->DrawLatex(0.2,0.8,"|#eta_{CM}| < 1");
TString outname="resolution";
c1->Print(Form("%s.png",outname.Data()));
c1->Print(Form("%s.pdf",outname.Data()));
}

