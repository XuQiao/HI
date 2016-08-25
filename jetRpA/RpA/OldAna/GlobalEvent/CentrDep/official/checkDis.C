#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void checkDis(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
bool Isweight=kFALSE;
bool IsMC=kTRUE;
if(IsMC){
if(Isweight){
TString dir="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/official/output/weight";
TString inputfile="Closure_Hijing_DS_bak.root";}
else{
TString dir="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/official/output/noweight";
TString inputfile="Closure_Hijing_DS_bak.root";}
}
else{
TString dir="/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/CentrDep";
if(!Isweight)
TString inputfile="MBCen_noweight_bak.root";
else
TString inputfile="MBCen_weight_bak_2.root";
}
TFile *f=TFile::Open(Form("%s/%s",dir.Data(),inputfile.Data()));
cout<<dir<<"/"<<inputfile<<endl;
if(IsMC){
TString Var[4]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4","Ntrk"};
TString corr[6]={"_tr","_DS","_corr_Eric","_corr_QiaoNtrk","_corr_QiaoHFPlus4","_corr_QiaoReHFPlus4"};
TH1F* hHF[4][6];
TH1F* rehHF[4][6];
}
else{
TString Var[4]={"HF4","HFPlus","HFPlus4","Ntrk"};
TString corr[5]={"","Wt","_QiaoNtrkWt","_QiaoHFPlus4Wt","_QiaoReHFPlus4Wt"};
TH1F* hHF[4][5];
}
int nVar=sizeof(Var)/sizeof(TString)-1;
int ncorr=sizeof(corr)/sizeof(TString)-1;
for(int k=0;k<=nVar;k++){
for(int j=0;j<=ncorr;j++){
if(IsMC)
hHF[k][j]=(TH1F*)f->Get(Form("h%s%s",Var[k].Data(),corr[j].Data()));
else
hHF[k][j]=(TH1F*)f->Get(Form("pfjet_%s%s",Var[k].Data(),corr[j].Data()));
}
}

TLatex Tl;
Tl.SetNDC();
Tl.SetTextSize(0.04);
Tl.SetTextAlign(12);
int kVar=3;
int jcorr1=1; int jcorr2=2;
//double HFbin[]={0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,8,10,12,15,20,30,40,60,80,100};
//double HFbin[]={0,0.1,0.24,0.36,0.47,0.57,0.67,0.77,0.86,0.96,1.05,1.15,1.24,1.33,1.43,1.53,1.62,1.72,1.83,1.93,2.03,2.13,2.24,2.36,2.47,2.59,2.7,2.83,2.95,3.08,3.21,3.34,3.48,3.62,3.76,3.91,4.06,4.22,4.38,4.54,4.71,4.88,5.05,5.23,5.42,5.6,5.79,5.98,6.18,6.38,6.58,6.79,6.99,7.2,7.41,7.62,7.84,8.07,8.3,8.53,8.76,9,9.24,9.49,9.73,9.99,10.24,10.5,10.77,11.04,11.31,11.6,11.9,12.19,12.49,12.81,13.12,13.45,13.79,14.14,14.5,14.87,15.26,15.67,16.08,16.51,16.97,17.45,17.96,18.51,19.11,19.75,20.44,21.21,22.12,23.14,24.34,25.84,27.91,31.24,74.73}; 1% centrality
double HFbin[]={0,0.24,0.47,0.67,0.86,1.05,1.24,1.43,1.62,1.83,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.76,4.06,4.38,4.71,5.05,5.42,5.79,6.18,6.58,6.99,7.41,7.84,8.3,8.76,9.24,9.73,10.24,10.77,11.31,11.9,12.49,13.12,13.79,14.5,15.26,16.08,16.97,17.96,19.11,20.44,22.12,24.34,27.91,74.73};
int nHFbin=sizeof(HFbin)/sizeof(double)-1;
c1 = new TCanvas("c1"," ",600,800);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.12,0.12,0.03);
TH1F* hFrame=new TH1F("hFrame","",10000,0,100);
fixedFontHist(hFrame,2.2,2);
TLegend *leg = new TLegend(0.3,0.2,0.5,0.3);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.04);
if(kVar!=3){
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
rehHF[kVar][jcorr1]=(TH1F*)hHF[kVar][jcorr1]->Rebin(nHFbin,Form("re%s",hHF[kVar][jcorr1]->GetName()),HFbin);
rehHF[kVar][jcorr2]=(TH1F*)hHF[kVar][jcorr2]->Rebin(nHFbin,Form("re%s",hHF[kVar][jcorr2]->GetName()),HFbin);
rehHF[kVar][0]=(TH1F*)hHF[kVar][0]->Rebin(nHFbin,Form("re%s",hHF[kVar][0]->GetName()),HFbin);
normalizeByBinWidth(rehHF[kVar][jcorr1]);
normalizeByBinWidth(rehHF[kVar][jcorr2]);
normalizeByBinWidth(rehHF[kVar][0]);
}
else{
hFrame->GetXaxis()->SetTitle("track Mult");
rehHF[kVar][jcorr1]=(TH1F*)hHF[kVar][jcorr1]->Clone(Form("re%s",hHF[kVar][jcorr1]->GetName()));
rehHF[kVar][jcorr2]=(TH1F*)hHF[kVar][jcorr2]->Clone(Form("re%s",hHF[kVar][jcorr2]->GetName()));
rehHF[kVar][0]=(TH1F*)hHF[kVar][0]->Clone(Form("re%s",hHF[kVar][0]->GetName()));
}
hFrame->GetXaxis()->SetRangeUser(0.09,50);
hFrame->GetYaxis()->SetTitle("# of events");
if(IsMC)
hFrame->GetYaxis()->SetRangeUser(0.11,1e5);
else
hFrame->GetYaxis()->SetRangeUser(1.1,1e8);
c1->cd(1)->SetLogy();
if(kVar!=3)
c1->cd(1)->SetLogx();
hFrame->DrawCopy();
rehHF[kVar][jcorr1]->SetMarkerColor(1);
rehHF[kVar][jcorr1]->SetMarkerStyle(20);
rehHF[kVar][jcorr1]->SetMarkerSize(1.2);
rehHF[kVar][jcorr2]->SetMarkerColor(2);
rehHF[kVar][jcorr2]->SetMarkerStyle(24);
rehHF[kVar][jcorr2]->SetMarkerSize(1.2);
rehHF[kVar][jcorr1]->Draw("same");
rehHF[kVar][jcorr2]->Draw("same");
leg->AddEntry(rehHF[kVar][jcorr1],Form("%s Eff=%.3f",hHF[kVar][jcorr1]->GetName(),hHF[kVar][0]->Integral()/hHF[kVar][jcorr1]->Integral()),"lp");
leg->AddEntry(rehHF[kVar][jcorr2],Form("%s Eff=%.3f",hHF[kVar][jcorr2]->GetName(),hHF[kVar][0]->Integral()/hHF[kVar][jcorr2]->Integral()),"lp");
leg->Draw("same");
c1->cd(2);
if(kVar!=3)
c1->cd(2)->SetLogx();
hFrame->GetYaxis()->SetRangeUser(0.95,1.049);
hFrame->GetYaxis()->SetTitle(Form("Ratio"));
hFrame->DrawCopy();
TH1F* hRatio=(TH1F*)rehHF[kVar][jcorr1]->Clone("hRatio");
hRatio->Divide(rehHF[kVar][jcorr2]);
hRatio->SetMarkerColor(1);
hRatio->SetMarkerSize(1.2);
hRatio->SetMarkerStyle(20);
hRatio->Draw("same");
if(IsMC)
Tl.DrawLatex(0.4,0.3,Form("Hijing: #frac{%s}{%s}",corr[jcorr1].Data(),corr[jcorr2].Data()));
else
Tl.DrawLatex(0.4,0.3,Form("Data: #frac{%s}{%s}",corr[jcorr1].Data(),corr[jcorr2].Data()));
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
l->Draw("same");
c1->Print("pic/checkDis.png");

}
