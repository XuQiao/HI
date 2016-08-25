#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 21, 30, 28,27,20
};
const double pPbLumi = 15.78 ;

void DrawRatio(){
    TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    TCanvas *c2 = new TCanvas("c2a", "c2",0,0,600,600);
    TCanvas *c3 = new TCanvas("c3a", "c3",0,0,600,600);
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
    c2->SetLeftMargin(0.13);
    c2->SetRightMargin(0.06);
    c2->SetTopMargin(0.05);
    c2->SetBottomMargin(0.12);
    c2->SetFrameFillStyle(0);
    c2->SetFrameBorderMode(0);
    c3->SetLeftMargin(0.13);
    c3->SetRightMargin(0.06);
    c3->SetTopMargin(0.05);
    c3->SetBottomMargin(0.12);
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 362, 429, 692, 1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;

const int N=8;
//const TString etabinname[N]={"-22_-12","-12_-7","-7_-3","-3_3","3_7","7_12","12_22","-10_10"};
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12","-10_10"};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TFile *f_GenReco[N];
TFile *f_Gen[N];
TH1D *hBay_GenReco[N];
TH1D *hBay_GenReco_Cl[N];
TH1D *hBay_GenReco_Cl1[N];
TH1D *hBay_GenReco_Cl2[N];
TH1D *hBay_Gen[N];
ofstream fstr[N];
for(int i=0;i<N;i++){
fstr[i].open(Form("jetaxissys%s.txt",etabinname[i].Data()));
fstr[i]<<setprecision(4)<<fixed;
f_GenReco[i] = TFile::Open(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/ForHong/PPb_UnfoPriorGenRecoEta_ak3PFKurtMC_jtpt30_EtaBin%s_Inc_v6.root",etabinname[i].Data()));
f_Gen[i] = TFile::Open(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/DataMC/FromYX/ForHong/PPb_UnfoPriorGen_ak3PFKurtMC_jtpt30_EtaBin%s_Inc_v6.root",etabinname[i].Data()));
hBay_Gen[i] = (TH1D*)f_Gen[i]->Get("hReco0");hBay_Gen[i]->Scale(1./etabin[i]);
hBay_GenReco[i] = (TH1D*)f_GenReco[i]->Get("hReco0");hBay_GenReco[i]->Scale(1./etabin[i]);
}

//-------------------------
c1->cd();
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(47,692);
//hFrame->GetXaxis()->SetLimits(50,600);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetRangeUser(0.9,1.1);
hFrame->GetYaxis()->SetTitle("jet axis uncertainty");
fixedFontHist(hFrame,1.2,1.4);
hFrame->DrawCopy();

TLegend* leg1=new TLegend(0.6,0.65,0.78,0.92);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.03);

TLegend* leg2=new TLegend(0.5,0.8,0.68,0.92);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetTextSize(0.03);

TLegend* leg3=new TLegend(0.5,0.68,0.68,0.92);
leg3->SetBorderSize(0);
leg3->SetFillColor(0);
leg3->SetTextSize(0.03);

for(int i=0;i<N;i++){
hBay_GenReco_Cl[i]=(TH1D*)hBay_GenReco[i]->Clone(Form("GenRecoRatio_%d",i));
hBay_GenReco_Cl[i]->Divide(hBay_Gen[i]);
//fstr[i]<<etastring[i]<<endl<<"Data jet axis"<<endl;
for(int ibin=1;ibin<=hBay_GenReco_Cl[i]->GetNbinsX();ibin++){
//hBay_GenReco_Cl[i]->SetBinContent(ibin,hBay_GenReco_Cl[i]->GetBinContent(ibin));
//fstr[i]<<hBay_GenReco_Cl[i]->GetXaxis()->GetBinLowEdge(ibin)<<"to"<<hBay_GenReco_Cl[i]->GetXaxis()->GetBinUpEdge(ibin)<<'\t';
if(hBay_GenReco_Cl[i]->GetBinCenter(ibin)>50 && hBay_GenReco_Cl[i]->GetBinCenter(ibin)<600){
fstr[i]<<hBay_GenReco_Cl[i]->GetBinCenter(ibin)<<'\t';
double value=TMath::Abs((hBay_GenReco_Cl[i]->GetBinContent(ibin)-1))*100;
fstr[i]<<value<<endl;
}
hBay_GenReco_Cl[i]->SetBinError(ibin,1e-10);
}
hBay_GenReco_Cl[i]->SetMarkerColor(colorCode[i]);
hBay_GenReco_Cl[i]->SetLineColor(colorCode[i]);
hBay_GenReco_Cl[i]->SetMarkerStyle(markerCode[i]);
hBay_GenReco_Cl[i]->SetMarkerSize(1.2);
hBay_GenReco_Cl1[i]=(TH1D*)hBay_GenReco_Cl[i]->Clone(Form("GenRecodoubleRatio_%d",i));
hBay_GenReco_Cl2[i]=(TH1D*)hBay_GenReco_Cl[i]->Clone(Form("GenRecoFMRatio_%d",i));
hBay_GenReco_Cl[i]->DrawCopy("same");
leg1->AddEntry(hBay_GenReco_Cl[i],etastring[i],"lp");
}
for(int i=0;i<N-1;i++)	hBay_GenReco_Cl[i]->DrawCopy("same");
leg1->Draw("same");
drawCMS(0.2,0.85,pPbLumi);

/*TLatex *com0 = new TLatex(0.57,0.90,"PYTHIA+HIJING");
com0->SetTextFont(42);
com0->SetTextSize(0.04);
com0->SetNDC();
com0->Draw();

TLatex *com1 = new TLatex(0.65,0.4,"CMS simulation");
com1->SetTextFont(63);
com1->SetTextSize(17);
com1->SetNDC();
com1->Draw();*/
TLatex *com3 = new TLatex(0.45,0.24,"Anti-k_{T} Particle Flow Jets R=0.3");
com3->SetTextFont(63);
com3->SetTextSize(17);
com3->SetNDC();
com3->Draw();

TLine *l=new TLine(hFrame->GetXaxis()->GetXmin(),1,hFrame->GetXaxis()->GetXmax(),1);
l->SetLineStyle(2);
l->Draw("same");

c1->Print("RecoGensys.pdf");

//----------------------------------------Double Ratio
c2->cd();
const int fweta[3]={2,1,0};
const int bweta[3]={4,5,6};
const int meta[2]={3,7};
hFrame->GetYaxis()->SetTitle("jet axis uncertainty Forward/Backward");
hFrame->DrawCopy();

for(int k=0;k<3;k++){
hBay_GenReco_Cl1[fweta[k]]->Divide(hBay_GenReco_Cl1[bweta[k]]);
hBay_GenReco_Cl1[fweta[k]]->SetMarkerColor(colorCode[fweta[2-k]]);
hBay_GenReco_Cl1[fweta[k]]->SetLineColor(colorCode[fweta[2-k]]);
hBay_GenReco_Cl1[fweta[k]]->SetMarkerStyle(markerCode[fweta[2-k]]);
hBay_GenReco_Cl1[fweta[k]]->DrawCopy("same");
leg2->AddEntry(hBay_GenReco_Cl1[fweta[k]],Form("%s/%s",etastring[fweta[k]].Data(),etastring[bweta[k]].Data()),"lp");
}

leg2->Draw("same");
drawCMS(0.2,0.85,pPbLumi);
com3->Draw();
l->Draw("same");
c2->Print("RecoGendoubleratioFBsys.pdf");

c3->cd();
hFrame->GetYaxis()->SetTitle("jet axis uncertainty Double side/Middle");
hFrame->DrawCopy();

for(int m=0;m<N;m++){
if(m!=meta[1]){
hBay_GenReco_Cl2[m]->Divide(hBay_GenReco_Cl2[meta[1]]);
hBay_GenReco_Cl2[m]->SetMarkerColor(colorCode[m]);
hBay_GenReco_Cl2[m]->SetLineColor(colorCode[m]);
hBay_GenReco_Cl2[m]->SetMarkerStyle(markerCode[m]);
hBay_GenReco_Cl2[m]->DrawCopy("same");
leg3->AddEntry(hBay_GenReco_Cl2[m],Form("%s/%s",etastring[m].Data(),etastring[meta[1]].Data()),"lp");
}}

leg3->Draw("same");
drawCMS(0.2,0.85,pPbLumi);
com3->Draw();
l->Draw("same");
c3->Print("RecoGendoubleratioFMsys.pdf");
}

