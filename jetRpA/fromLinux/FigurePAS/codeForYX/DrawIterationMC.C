#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};

double getMax(double arr[],int N){
double maxv=0;
for(int i=0;i<N;i++){
maxv=TMath::Max(maxv,arr[i]);
}
return maxv;
}


const double pPbLumi = 15.78 ; //excluded the old alignment run for pPb
void DrawIterationMC(){

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
    c1->SetBottomMargin(0.11);
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
//int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
const int Nbin_pt=30;
const int Neta=8;
//const TString etabinname[Neta]={"-22_-12","-12_-7","-7_-3","-3_3","3_7","7_12","12_22","-10_10"};
const TString etabinname[Neta]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[Neta]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
//const TString etastring[Neta]={"-2.2<#eta<-1.2","-1.2<#eta<-0.7","-0.7<#eta<-0.3","-0.3<#eta<0.3","0.3<#eta<0.7","0.7<#eta<1.2","1.2<#eta<2.2","-1.0<#eta<1.0"};
const TString etastring[Neta]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(50,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.9,1.1);
hFrame->GetYaxis()->SetTitle("Ratio (Unfolded/Nominal)");
fixedFontHist(hFrame,1.1,1.4);
hFrame->DrawCopy();
TFile *f[Neta];
TH1D *hBay[Neta]; TH1F *hGen_PPb[Neta]; TH1F* hBay_Cl[Neta];
TLegend* leg=new TLegend(0.6,0.65,0.78,0.92);
const int N=9;
TH1D* hReco_cent[Neta][N];
TH1D* hReco_cent0[Neta];
TH1D* IterMax[Neta];
TH1D* hReco_cent_Cl[Neta][N];
double IterMaxvalue[Neta][Nbin_pt][N];
ofstream fstr[Neta];
for(int i=0;i<Neta;i++){
fstr[i].open(Form("jetsysIter%s.txt",etabinname[i].Data()));
//f[i] = TFile::Open(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/ForHong/PPb_UnfoPriorGen_ak3PFKurtMC_jtpt30_EtaBin%s_Inc_v6.root",etabinname[i].Data()));
f[i] = TFile::Open(Form("/scratch/xuq7/RpA/JetTrig/DataMC/FromYX/PPb_UnfoPriorGen_akPu3PFKurtMC_MC_jtpt30_EtaBin%s_Inc_v3.root",etabinname[i].Data()));
//TH1D* hReco0=(TH1D*)f->Get("hReco0");
hReco_cent0[i]=(TH1D*)f[i]->Get("hReco_cent0");
//TH1D* hUnfoldedJeCsys_cent0=(TH1D*)f->Get("UnfoldedJeCSys_cent0");
for(int j=0;j<N;j++){
//int i=0;
hReco_cent[i][j]=(TH1D*)f[i]->Get(Form("hRecoRAA_IterSys%d_cent0",j+2));
hReco_cent_Cl[i][j]= (TH1D*)hReco_cent[i][j]->Clone(Form("hRecoRAA_IterSys%d_cent0_Cl_%d",j+2,i));
hReco_cent_Cl[i][j]->Divide(hReco_cent0[i]);
}
	IterMax[i]=(TH1D*)hReco_cent_Cl[i][0]->Clone(Form("IterMax_%d",i));
	for(int ibin=1;ibin<=IterMax[i]->GetNbinsX();ibin++){
		for(int j=0;j<N;j++)
		IterMaxvalue[i][ibin][j]=hReco_cent_Cl[i][j]->GetBinContent(ibin);
		IterMax[i]->SetBinContent(ibin,getMax(IterMaxvalue[i][ibin],N));
		IterMax[i]->SetBinError(ibin,1e-10);
	}
	for(int ibin=1;ibin<=IterMax[i]->GetNbinsX();ibin++){
	fstr[i]<<IterMax[i]->GetXaxis()->GetBinLowEdge(ibin)<<"to"<<IterMax[i]->GetXaxis()->GetBinUpEdge(ibin)<<'\t';
		fstr[i]<<IterMax[i]->GetBinContent(ibin)<<endl;
}
}	

TLegend* leg1=new TLegend(0.6,0.60,0.78,0.88);
//TLegend* leg2=new TLegend(0.65,0.70,0.82,0.85);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.035);
for(int i=0;i<Neta;i++){
IterMax[i]->SetMarkerStyle(markerCode[i]);
IterMax[i]->SetMarkerColor(colorCode[i]);
IterMax[i]->SetLineColor(colorCode[i]);
IterMax[i]->SetMarkerSize(1.2);
IterMax[i]->DrawCopy("same");
leg1->AddEntry(IterMax[i],etastring[i],"lp");
}
for(int i=0;i<Neta-1;i++) IterMax[i]->DrawCopy("same");
//hUnfoldedJeCsys_cent0_Cl->DrawCopy("same");
//leg2->AddEntry(hUnfoldedJeCsys_cent0_Cl,"UnfoldedJeCsys","lp");
leg1->Draw("same");

//drawCMS(0.2,0.88,pPbLumi);
TLatex *com0 = new TLatex(0.2,0.88,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
com0->Draw();
TLatex *com3 = new TLatex(0.45,0.24,"Anti-k_{T} Particle Flow Jets R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();


TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");
c1->Print("Iteration_MC_Etabin.pdf");

}
