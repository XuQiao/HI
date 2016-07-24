#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};
const double pPbLumi = 20.7 ;
void DrawUnfoldingData(){

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);   
    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.12);
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
//c1 = new TCanvas("c1"," ",400,700);
//makeMultiPanelCanvas(c1,1,2,-0.14,0,0.2,0.14,0.03);
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

const int N=8;
//const TString etabinname[N]={"-22_-12","-12_-7","-7_-3","-3_3","3_7","7_12","12_22","-10_10"};
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(45,450);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.5,1.5);
hFrame->GetYaxis()->SetTitle("no unfolding/Baysian");
fixedFontHist(hFrame,1.2,1.4);
hFrame->DrawCopy();
TFile *f[N];
TH1D* hBay[N];
TH1F* hMeas[N];
TH1F* hBay_Cl[N];
TLegend* leg2=new TLegend(0.68,0.20,0.84,0.45);
ofstream fstr[N];
for(int i=0;i<N;i++){
fstr[i]<<setprecision(4)<<fixed;
fstr[i].open(Form("jetUnfsysData%s.txt",etabinname[i].Data()));
//TFile *f = TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/PPb_UnfoPriorGen_akPu3PFNoResidual_MC_jtpt30_EtaBin-10_10_Inc_v2.root");
//TFile *f = TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/PPb_UnfoPriorGen_akPu3PFNoResidual_MC_jtpt30_EtaBin-10_10_Inc_v2.root");
f[i] = TFile::Open(Form("/scratch/xuq7/RpA/JetTrig/DataMC/FromYX/PPb_UnfoPriorGen_akPu3PFKurtMC_jtpt20_EtaBin%s_Inc_v7.root",etabinname[i].Data()));
//TFile *f = TFile::Open(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/ForHong/PPb_UnfoPriorGenRecoEta_ak3PFKurtMC_jtpt30_EtaBin%s_Inc_v6.root",etabinname[i].Data()));
hBay[i] = (TH1D*)f[i]->Get("hReco0");
hMeas[i] = (TH1F*)f[i]->Get("hMeas0");
//TH1F* hBin = (TH1F*)f->Get("hRecoBinByBin0");

/*hGen_PPb=(TH1F*)hGen_PPb->Rebin(Nbin_pt,hGen_PPb->GetName(),binbound_pt);
hBay=(TH1F*)hBay->Rebin(Nbin_pt,hBay->GetName(),binbound_pt);
hMeas=(TH1F*)hMeas->Rebin(Nbin_pt,hMeas->GetName(),binbound_pt);
hBin=(TH1F*)hBin->Rebin(Nbin_pt,hBin->GetName(),binbound_pt);
*/

//hBin_Cl = (TH1F*)hBin->Clone("hBin_Cl");

//c1->cd(1)->SetLogy();
/*hMeas->SetMarkerStyle(20);
hMeas->SetLineColor(2);
hMeas->SetMarkerColor(2);
hMeas->SetMarkerSize(1.2);

hBin->SetMarkerStyle(25);
hBin->SetLineColor(1);
hBin->SetMarkerColor(1);
hBin->SetMarkerSize(1.2);
*/
hBay_Cl[i] = (TH1F*)hMeas[i]->Clone(Form("hBay_Cl_%d",i));
hBay_Cl[i]->SetMarkerStyle(markerCode[i]);
hBay_Cl[i]->SetLineColor(colorCode[i]);
hBay_Cl[i]->SetMarkerColor(colorCode[i]);
hBay_Cl[i]->SetMarkerSize(1.2);

/*
hMeas->DrawCopy("same");
hBin->DrawCopy("same");
hBay->DrawCopy("same");*/
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.035);
/*
TLatex *com0 = new TLatex(0.57,0.88,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.045);
com0->SetNDC();
com0->Draw();

TLatex *com1 = new TLatex(0.65,0.4,"CMS simulation");
com1->SetTextFont(63);
com1->SetTextSize(17);
com1->SetNDC();
com1->Draw();
*/
TLatex *com2 = new TLatex(0.65,0.32,etastring[i]);
com2->SetTextFont(63);
com2->SetTextSize(17);
com2->SetNDC();
//com2->Draw();
TLatex *com3 = new TLatex(0.45,0.84,"Anti-k_{T} Particle Flow PU Jets R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();

TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
//c1->cd(2);
//hBin_Cl->Divide(hMeas);
hBay_Cl[i]->Divide(hBay[i]);
//hBin_Cl->DrawCopy("same");
hBay_Cl[i]->DrawCopy("same");
for(int ibin=1;ibin<=hBay_Cl[i]->GetNbinsX();ibin++){
if(hBay_Cl[i]->GetBinCenter(ibin)>50 && hBay_Cl[i]->GetBinCenter(ibin)<600){
fstr[i]<<hBay_Cl[i]->GetBinCenter(ibin)<<'\t';
fstr[i]<<100*TMath::Abs((hBay_Cl[i]->GetBinContent(ibin)-1))<<endl;
}
}
leg2->AddEntry(hBay_Cl[i],etastring[i],"lp");
}
for(int i=0;i<N-1;i++) 	hBay_Cl[i]->DrawCopy("same");
//com2->Draw();
com3->Draw();
drawCMS(0.2,0.88,pPbLumi);
l->Draw("same");
leg2->Draw("same");
c1->Print(Form("Unfolding_Data_Etabin.pdf"));
}

