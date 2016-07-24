#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

void FitptTruth(){

	const int N=8;
    TCanvas *c1[N];
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);   

    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TFile *f; TF1 *fnlo[N];
ofstream fstr("refptpar.txt");
TH1D *hBay[N]; TH1F *hGen_PPb[N]; TH1F* hGen_PPb_Cl[N];
TPaveStats *stats;
TLegend* leg=new TLegend(0.6,0.5,0.78,0.55);
f = TFile::Open(Form("/scratch/xuq7/RpA/FromYX/MCPPbakPu3PFJetSpectraCombineTriggerJetWeight.root"));
for(int i=0;i<N;i++){
hGen_PPb[i]=new TH1F(Form("hGen_PPb_%d",i),"",1000,0,1000);
fstr<<etastring[i]<<endl;
c1[i] = new TCanvas(Form("c%da",i), Form("c%d",i),0,0,600,600);
 c1[i]->Range(0,0,1,1);
    c1[i]->SetFillColor(0);
    c1[i]->SetBorderMode(0);
    c1[i]->SetBorderSize(2);
    c1[i]->SetTickx(1);
    c1[i]->SetTicky(1);
    c1[i]->SetLeftMargin(0.16);
    c1[i]->SetRightMargin(0.06);
    c1[i]->SetTopMargin(0.05);
    c1[i]->SetBottomMargin(0.1);
    c1[i]->SetFrameFillStyle(0);
    c1[i]->SetFrameBorderMode(0);
c1[i]->cd()->SetLogy();
//hFrame->DrawCopy();
//f[i] = TFile::Open(Form("/scratch/xuq7/RpA/FromYX/PPb_UnfoPriorGen_akPu3PFKurtMC_MC_jtpt20_EtaBin%s_Inc_v5.root",etabinname[i].Data()));

fnlo[i] = new TF1(Form("fus%d",i),"[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/4000.,[3])", 20., 1000.);
//fnlo[i] = new TF1(Form("fus%d",i),"[0]*exp([1]/x)*pow(x,[2])", 20., 1000.);
fnlo[i]->SetParameters(4.657e6,-41.96,-5.566,2.922,-1.908);
//fnlo[i]->FixParameter(0,1);
fnlo[i]->SetLineStyle(1);
fnlo[i]->SetLineColor(2);
//hGen_PPb[i] = (TH1F*)f[i]->Get("hGen0");
hGen_PPb_Cl[i] = (TH1F*)f->Get(Form("jetptEtaBin%s",etabinname[i].Data()));
//hGen_PPb[i]=(TH1F*)hGen_PPb_Cl[i]->Clone();
hGen_PPb[i]->GetXaxis()->SetRangeUser(0,500);
for(int ibin=0;ibin<=hGen_PPb_Cl[i]->GetNbinsX();ibin++){
hGen_PPb[i]->SetBinContent(ibin,hGen_PPb_Cl[i]->GetBinContent(ibin));
hGen_PPb[i]->SetBinError(ibin,hGen_PPb_Cl[i]->GetBinError(ibin));
}
hGen_PPb[i]->GetXaxis()->SetTitle("p_{T}^{ref} (GeV/c)");
hGen_PPb[i]->GetYaxis()->SetRangeUser(5e-13,1e-2);
hGen_PPb[i]->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}}");
//hGen_PPb[i]->Scale(1./etabin[i]);
hGen_PPb[i]->SetMarkerStyle(markerCode[i]);
hGen_PPb[i]->SetLineColor(0);
hGen_PPb[i]->SetMarkerColor(0);
hGen_PPb[i]->SetMarkerSize(1);
hGen_PPb_Cl[i]->SetLineColor(1);
hGen_PPb_Cl[i]->SetMarkerColor(1);
hGen_PPb_Cl[i]->SetMarkerStyle(20);
hGen_PPb_Cl[i]->SetMarkerSize(0.9);
hGen_PPb[i]->Draw();
fixedFontHist(hGen_PPb[i],1.2,1.5);
hGen_PPb_Cl[i]->Fit(fnlo[i],"","same",30,600);
fstr<<"NDF="<<fnlo[i]->GetNDF()<<'\t'<<"Chi2="<<fnlo[i]->GetChisquare()<<endl;
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.035);
leg->AddEntry(hGen_PPb[i],etastring[i],"lp");
//fnlo[i]->Draw("same");
TPaveText *fitpar = new TPaveText(0.54,0.6,0.92,0.9,"NDC");
fitpar->SetTextFont(72);
fitpar->SetTextAlign(11);
fitpar->SetTextSize(0.028);
fitpar->SetFillColor(0);
fitpar->SetBorderSize(1);
fitpar->AddText(Form("#chi^{2}/Ndf = %.1f/%d",fnlo[i]->GetChisquare(),fnlo[i]->GetNDF()));
fitpar->AddText(Form("p0 = %.3f #pm %.3f",fnlo[i]->GetParameter(0),fnlo[i]->GetParError(0)));
fitpar->AddText(Form("p1 = %.3f #pm %.3f",fnlo[i]->GetParameter(1),fnlo[i]->GetParError(1)));
fitpar->AddText(Form("p2 = %.3f #pm %.3f",fnlo[i]->GetParameter(2),fnlo[i]->GetParError(2)));
fitpar->AddText(Form("p3 = %.3f #pm %.3f",fnlo[i]->GetParameter(3),fnlo[i]->GetParError(3)));
fitpar->AddText(Form("p4 = %.3f #pm %.3f",fnlo[i]->GetParameter(4),fnlo[i]->GetParError(4)));
fitpar->Draw("same");

TLatex *com0 = new TLatex(0.25,0.88,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
TLatex *com2 = new TLatex(0.20,0.32,"CMS simulation");
com2->SetTextFont(63);
com2->SetTextSize(17);
com2->SetNDC();
TLatex *com3 = new TLatex(0.20,0.24,"Anti-k_{T} Particle Flow PU Jets R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();

fstr<<"par0="<<fnlo[i]->GetParameter(0)<<'\t';
fstr<<"par1="<<fnlo[i]->GetParameter(1)<<'\t';
fstr<<"par2="<<fnlo[i]->GetParameter(2)<<'\t';
fstr<<"par3="<<fnlo[i]->GetParameter(3)<<'\t';
fstr<<"par4="<<fnlo[i]->GetParameter(4)<<endl;
fstr<<"par0 Error="<<fnlo[i]->GetParError(0)<<'\t';
fstr<<"par1 Error="<<fnlo[i]->GetParError(1)<<'\t';
fstr<<"par2 Error="<<fnlo[i]->GetParError(2)<<'\t';
fstr<<"par3 Error="<<fnlo[i]->GetParError(3)<<'\t';
fstr<<"par4 Error="<<fnlo[i]->GetParError(4)<<endl;

fstr<<"ftrue->SetParameters("<<fnlo[i]->GetParameter(0)<<","<<fnlo[i]->GetParameter(1)<<","<<fnlo[i]->GetParameter(2)<<","<<fnlo[i]->GetParameter(3)<<","<<fnlo[i]->GetParameter(4)<<")";

leg->Draw("same");
com0->Draw("same");
com2->Draw("same");
com3->Draw("same");
c1[i]->Print(Form("TruthFit_%s.pdf",etabinname[i].Data()));
leg->Clear();
}
}

