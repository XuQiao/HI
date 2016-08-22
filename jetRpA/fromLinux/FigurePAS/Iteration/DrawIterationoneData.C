#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

static const int nColor = 9;

static const int colorCode[nColor] = {
    1, 2, 4, 46,6, 7,8,  1,4
};
static const int markerCode[nColor] = {
    20, 34, 33, 25,27, 28,24 , 33,25
};

const double pPbLumi = 35.0 ; //excluded the old alignment run for pPb
void DrawIterationoneData(){

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
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

const int Neta=8;
int ieta=7;
const int binbound_pt[] = { 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
const int Nbin_pt=40;

const TString etabinname[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const double etabin[Neta]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};
f = TFile::Open(Form("/home/xuq/Documents/HI/RpA/output/UnfoldedOutputs/PPb_UnfoPriorGen_akPu3PFKurtForestEtaWeightNoIDCut_jtpt20_EtaBin%s_Inc_v13.root",etabinname[ieta].Data()));

const int N=9;
TH1D* hReco_cent[N];
TH1D* hReco_cent_Cl[N];
TH1D* hReco0=(TH1D*)f->Get("hReco0");
TH1D* hReco_cent0=(TH1D*)f->Get("hReco_cent0");
TH1D* hUnfoldedJeCsys_cent0=(TH1D*)f->Get("UnfoldedJeCSys_cent0");
for(int i=0;i<N;i++){
//int i=0;
hReco_cent[i]=(TH1D*)f->Get(Form("hRecoRAA_IterSys%d_cent0",i+2));
hReco_cent_Cl[i]= (TH1D*)hReco_cent[i]->Clone(Form("hRecoRAA_IterSys%d_cent0_Cl",i+2));
hReco_cent_Cl[i]->Divide(hReco_cent0);
hReco_cent_Cl[i]->SetMarkerStyle(markerCode[i]);
hReco_cent_Cl[i]->SetMarkerColor(colorCode[i]);
hReco_cent_Cl[i]->SetLineColor(colorCode[i]);
hReco_cent_Cl[i]->SetMarkerSize(1.2);
}
TH1D* hUnfoldedJeCsys_cent0_Cl=(TH1D*)hUnfoldedJeCsys_cent0->Clone("UnfoldedJeCsys_cent0_Cl");
hUnfoldedJeCsys_cent0_Cl->Divide(hReco_cent0);
hUnfoldedJeCsys_cent0_Cl->SetMarkerStyle(20);
hUnfoldedJeCsys_cent0_Cl->SetMarkerColor(1);
hUnfoldedJeCsys_cent0_Cl->SetMarkerSize(1);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(47.5,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.6,1.5);
hFrame->GetYaxis()->SetTitle("Ratio (Unfolded/Nominal)");
fixedFontHist(hFrame,1.2,1.6);
hFrame->DrawCopy();

TLegend* leg1=new TLegend(0.45,0.74,0.62,0.88);
TLegend* leg2=new TLegend(0.65,0.74,0.82,0.88);
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.035);
leg2->SetTextSize(0.035);
for(int i=0;i<N;i++){
//i=0;
hReco_cent_Cl[i]->DrawCopy("same");
if(i<5)
leg1->AddEntry(hReco_cent_Cl[i],Form("Iteration %d",i+2),"lp");
else
leg2->AddEntry(hReco_cent_Cl[i],Form("Iteration %d",i+2),"lp");
}
//hUnfoldedJeCsys_cent0_Cl->DrawCopy("same");
//leg2->AddEntry(hUnfoldedJeCsys_cent0_Cl,"UnfoldedJeCsys","lp");
leg1->Draw("same");
leg2->Draw("same");

drawCMS(0.2,0.88,pPbLumi);

TLatex *com0 = new TLatex(0.17,0.88,"PYTHIA+HIJING");
com0->SetTextFont(43);
com0->SetTextSize(20);
com0->SetNDC();
//com0->Draw();

TLatex *com2 = new TLatex(0.63,0.28,etastring[ieta]);
com2->SetTextFont(43);
com2->SetTextSize(20);
com2->SetNDC();
com2->Draw();
TLatex *com3 = new TLatex(0.30,0.20,"Anti-k_{T} Particle Flow Jets: R=0.3");
com3->SetTextFont(43);
com3->SetTextSize(20);
com3->SetNDC();
com3->Draw();


TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");

c1->Print(Form("Iteration_oneData_%s.pdf",etabinname[ieta].Data()));

}
