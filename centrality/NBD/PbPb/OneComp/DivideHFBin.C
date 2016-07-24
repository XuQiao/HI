#include <iomanip>
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
void DivideHFBin(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
double goal;
bool IsMC=kFALSE, Isweight=kTRUE;
if(IsMC){
if(Isweight){
TString dir="/home/xuq/Documents/HI/RpA/GlobalEvent/CentrDep/RecoHI/VRw";
TString inputfile="Closure.root";
TString inputfileDS="pPbHijing_EffCorr.root";}
else{
TString dir="/home/xuq/Documents/HI/RpA/GlobalEvent/CentrDep/RecoHI/NoVRw";
TString inputfile="Closure_v1.root";
TString inputfileDS="pPbHijing_EffCorr.root";
}
}
else{
TString dir="/home/xuq/Documents/HI/RpA/output/GlobalEvent/CentrDep/";
if(!Isweight)
TString inputfile="MBCen_noweight_RecoHI_v1.root";
else
TString inputfile="MBCen_weight_RecoHI.root";
}
ofstream fstr;
//fstr<<setprecision(2)<<fixed<<endl;
TFile *f=TFile::Open(Form("%s/%s",dir.Data(),inputfile.Data()));
if(IsMC)
TFile *fDS=TFile::Open(Form("%s/%s",dir.Data(),inputfileDS.Data()));
fstr<<dir<<"/"<<inputfile<<endl;
if(IsMC){
TString Var[6]={"HFp4","HFp","HFm4","HFm","HF4","Mult"};
TString corr[9]={"ES","DS","NtrkEric","Ntrk","Ntrkbak","HFp4","HFp4bak","reHFp4","reHFp4bak"};
TH1F* hHF[6][9]; 	TH1F* hHF_tr[6];
}
else{
TString Var[4]={"HF4","HFPlus","HFPlus4","Ntrk"};
TString corr[7]={"Wt","_QiaoNtrkWt","_QiaoNtrkbakWt","_QiaoHFPlus4bakWt","_QiaoHFPlus4Wt","_QiaoReHFPlus4Wt","_QiaoReHFPlus4bakWt"};
TH1F* hHF[4][7];	TH1F* hHF_tr[4];
}

if(IsMC){
if(Isweight)
fstr.open(Form("HFbin_Hijing_weight_bak.txt"));
else
fstr.open(Form("HFbin_Hijing_noweight_bak.txt"));
}
else{
if(Isweight)
fstr.open(Form("HFbin_Data_weight.txt"));
else
fstr.open(Form("HFbin_Data_noweight.txt"));
}

int ncorr=sizeof(corr)/sizeof(TString)-1;
const int N=5;
double cent[N]={0.1,0.2,0.3,0.6,0.9};
double HFcentrbin[N];
double HFbin[]={0,0.23,0.47,0.67,0.87,1.05,1.24,1.43,1.62,1.82,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.77,4.07,4.38,4.72,5.06,5.41,5.79,6.17,6.57,6.98,7.41,7.85,8.3,8.76,9.24,9.74,10.25,10.79,11.33,11.91,12.51,13.14,13.81,14.51,15.27,16.09,17,17.99,19.14,20.48,22.14,24.35,27.93,90.29};
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
TLine *l[N];
//for(int k=0;k<nVar;k++){
int k=2;
//for(int j=0;j<=ncorr;j++){
int j=4;
if(IsMC){
if(j==0)
hHF[k][j]=(TH1F*)f->Get(Form("demo/h%s%s",Var[k].Data(),corr[j].Data()));
else if(j==1)
hHF[k][j]=(TH1F*)fDS->Get(Form("h%s%s",Var[k].Data(),corr[j].Data()));
else
hHF[k][j]=(TH1F*)f->Get(Form("demo/h%s%s",corr[j].Data(),Var[k].Data()));
hHF_tr[k]=(TH1F*)f->Get(Form("demo/h%s%s",Var[k].Data(),corr[0].Data()));
}
else{
if(j==2){
if(k==1)
corr[j]="_QiaoNtrkbalWt";
else
corr[j]="_QiaoNtrkbakWt";
}
hHF[k][j]=(TH1F*)f->Get(Form("pfjet_%s%s",Var[k].Data(),corr[j].Data()));
hHF_tr[k]=(TH1F*)f->Get(Form("pfjet_%s",Var[k].Data()));
}
fstr<<hHF[k][j]->GetName()<<" Integral="<<hHF_tr[k]->Integral()/hHF[k][j]->Integral()<<endl;
for(int icent=0;icent<N;icent++){
//fstr<<"start from 90-100%:"<<endl;
for(int ibin=1;ibin<hHF[k][j]->GetNbinsX();ibin++){
goal=cent[icent]*hHF[k][j]->Integral();
if(hHF[k][j]->Integral(ibin,hHF[k][j]->GetNbinsX())<=goal){
HFcentrbin[icent]=hHF[k][j]->GetXaxis()->GetBinLowEdge(ibin);
fstr<<icent<<"th: "<<hHF[k][j]->GetXaxis()->GetBinLowEdge(ibin)<<'\t'<<cent[icent]<<endl;
break;
}
}
l[icent]=new TLine(HFcentrbin[icent],0,HFcentrbin[icent],6e-2);
l[icent]->SetLineStyle(2);
l[icent]->SetLineColor(4);
l[icent]->SetLineWidth(2);
}
//}
//}
hHF[k][j]=(TH1F*)hHF[k][j]->Rebin(nHFbin,hHF[k][j]->GetName(),HFbin);
normalizeByBinWidth(hHF[k][j]);
hHF[k][j]->Scale((double)1.0/hHF[k][j]->Integral());
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.05,0.05,0.03);
c1->cd()->SetLogy();
//c1->cd()->SetLogx();
TH1F* hFrame=new TH1F("","",1000,0,100);
fixedFontHist(hFrame,1.3,1.2);
hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} #eta>4");
hFrame->GetXaxis()->SetRangeUser(1e-1,30);
hFrame->GetYaxis()->SetTitle("Event Fraction");
hFrame->GetYaxis()->SetRangeUser(1e-3,1e-1);
hFrame->Draw();
TLatex *T1=new TLatex(0.25,0.85,"CMS preliminary pPb #sqrt{s} = 5.02TeV");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);

TLatex T2;
T2.SetNDC();
T2.SetTextAngle(90);
T2.SetTextSize(0.03);
T2.SetTextColor(1);
T2.SetTextFont(42);
T2.DrawLatex(0.27,0.30,"90-100%");
T2.DrawLatex(0.54,0.30,"60-90%");
T2.DrawLatex(0.71,0.30,"30-60%");
T2.DrawLatex(0.79,0.30,"20-30%");
T2.DrawLatex(0.83,0.30,"10-20%");
T2.DrawLatex(0.88,0.30,"0-10%");


hHF[k][j]->SetMarkerStyle(20);
hHF[k][j]->SetMarkerColor(1);
hHF[k][j]->SetLineColor(1);
hHF[k][j]->SetMarkerSize(1);
hHF[k][j]->Draw("same");
T1->Draw("same");
for(int icent=0;icent<N;icent++){
l[icent]->Draw("same");
}
c1->Print("HFbin.png");
}
