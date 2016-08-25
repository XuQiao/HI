#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void overlayJetIDcut(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const int Npoint=1000;
TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
double binbound_JetID[];
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={30+1e-4,600-1e-4};
TString JetID;

int ilist=14;
TString JetIDName = JetIDNameList[ilist];

double cut;
if(ilist== 6){      JetID = "chargedMax/p_{T}";	cut=0.3;}
else if(ilist== 7){      JetID = "chargedSum/p_{T}";	cut=0.6;}
else if(ilist== 8){     JetID = "neutralMax/p_{T}";	cut=0.08;}
else if(ilist== 9){     JetID = "neutralSum/p_{T}";	cut=0.15;}
//else if(ilist== 10){     JetID = "photonMax/p_{T}";	}
else if(ilist== 11){     JetID = "photonSum/p_{T}";	cut=0.3;}
else if(ilist== 12){     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";cut =1.01;}
else if(ilist== 13){     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}";cut =1.13;}
else if(ilist== 14){     JetID = "neutralMax/Max(neutralSum,chargedSum)";cut =0.975;}
else{   exit();}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}


//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};


int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TString filename="/scratch/xuq7/RpA/TreeAna/Datacombined.root";

TString histoname=Form("jetpt%sCombinedSpectra",JetIDName.Data());
//TString histoname=Form("jetptEtaCombinedSpectra");

TFile *file=TFile::Open(filename);

TH2F* histo2F=(TH2F*)file->Get(histoname);
TH1D* histo1=(TH1D*)histo2F->ProjectionX("histo1",histo2F->GetYaxis()->FindBin(binbound_JetID[0]),histo2F->GetYaxis()->FindBin(cut));
//TH1D* histo1=(TH1D*)histo2F->ProjectionX("histo1");
TH1D* histo2=(TH1D*)histo2F->ProjectionX("histo2");
histo1=(TH1D*)histo1->Rebin(Nbin_pt,"histo1",binbound_pt);
histo2=(TH1D*)histo2->Rebin(Nbin_pt,"histo2",binbound_pt);
normalizeByBinWidth(histo1);
normalizeByBinWidth(histo2);
histo1->SetMarkerStyle(24);
histo1->SetMarkerSize(1.4);
histo1->SetMarkerColor(1);
histo1->GetYaxis()->SetTitle("Number of Events");
histo2->SetMarkerStyle(20);
histo2->SetMarkerSize(1.2);
histo2->SetMarkerColor(2);
fixedFontHist(histo1,1.8,2.0);
fixedFontHist(histo2,1.8,2.0);
cout<<histo1->GetBinCenter(9)<<'\t'<<histo1->GetBinContent(9)<<endl;
TGraphAsymmErrors *g;
TGraphAsymmErrors *g0;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];

TH1F* hFrame=new TH1F("","",1000,0,1000);
c1 = new TCanvas("c1","",600,1000);
makeMultiPanelCanvas(c1,1,2,0.03,0.03,0.1,0.12,0.03);
c1->cd(1)->SetLogy();
fixedFontHist(hFrame,1.2,1.4);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("Number of Events");
hFrame->GetXaxis()->SetLimits(25,600);
hFrame->GetYaxis()->SetRangeUser(2e-1,5e7);
hFrame->DrawCopy();

histo1->Draw("E1same");
histo2->Draw("E1same");
TLegend *leg1=new TLegend(0.60,0.75,0.85,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.06);
leg1->AddEntry(histo1,"After cut","lp");
leg1->AddEntry(histo2,"Before cut","lp");
leg1->Draw("same");
TLatex *T=new TLatex(0.25,0.30,Form("Cut: %s<%.2f",JetID.Data(),cut));
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.04);
T->SetTextColor(1);
T->SetTextFont(42);

c1->cd(2);
TH1D* ratio=(TH1D*)histo1->Clone("ratio");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle("Ratio: After cut/Before cut");
hFrame->GetYaxis()->SetRangeUser(0.7,1.09);
hFrame->DrawCopy();
ratio->Divide(histo2);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.4);
ratio->SetMarkerColor(1);
ratio->SetMinimum(ratio->GetMinimum()/1.2);
ratio->SetMaximum(1.09);
ratio->SetTitle("");
//ratio->Draw("E1");

g=new TGraphAsymmErrors();
g->BayesDivide(histo1,histo2);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
g0->SetMarkerColor(1);
g0->SetMarkerStyle(20);
g0->SetMarkerSize(1.2);
g0->SetLineColor(1);
g0->Draw("Psame");
TLine *l =new TLine(30,1,600,1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
T->Draw("same");

c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/NoGplusNewpthat/newSample/pic/JetIDcut/akPu3PF/%s/overlay_JetIDcut.png",JetIDName.Data()));
}
