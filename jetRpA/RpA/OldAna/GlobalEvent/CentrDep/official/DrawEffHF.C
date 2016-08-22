#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawEffHF(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const int Ncut=5;	const int Nhist=3;	const int NGen=2; const int Ncent=6;
const double HFbin[]={0,0.1,0.2,0.4,0.6,0.8,1,2,4,6,9,12,15,20,30,50,75,100};
const double xerr[Ncent]={0,0,0,0,0,0};
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.92,10.1,19.3,23.7,30.2,100},
{0,1.03,5.00,11.6,15.0,20.3,100},
{0,2.7,12.6,28.6,36.9,49.6,100}};
const double cent[Ncent]={5,15,25,45,75,95};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;
double color[NGen][Nhist]={{2,1,4},{2,1,4}};
double Style[NGen][Nhist]={{20,24,28},{20,24,28}};
//TString cut[Ncut]={"HLT","CutGplus","pPA","pprim","vz",""};
//TString cutName[Ncut]={"HLT","CutG","clean","pprim","vz","All cuts"};
TString cut[Ncut]={"HLT","HLT_CutGplus","HLT_CutGplus_pPA","HLT_CutGplus_pPA_pprim",""};
TString cutName[Ncut]={"HLT","HLT+CutGplus","HLT+CutGplus+clean","HLT+CutG+clean+pprim","All cuts"};
TString Gen[NGen]={"Hijing","Epos"};
TString histoname[Nhist]={"hHFEnergy4","hHFEnergyPlus4","hHFEnergyPlus"};
TString etarange[Nhist]={"HF E_{T} |#eta|>4","HF E_{T} #eta>4","HF E_{T} #eta>0(Pb going side)"};
TString filename[Ncut][NGen];
TFile* f[Ncut][NGen];
TH1F* h_or[Ncut][NGen][Nhist];	TH1F* reh_or[Ncut][NGen][Nhist];
TH1F* h_tr[Ncut][NGen][Nhist];	TH1F* reh_tr[Ncut][NGen][Nhist];
TGraphErrors *g_Eff[Ncut][NGen][Nhist]; double EffValue[Ncut][NGen][Nhist][Ncent];	double EffValueErr[Ncut][NGen][Nhist][Ncent];
TH1F* rehEff[Ncut][NGen][Nhist];	TLatex* T1[Nhist];	TLatex* T2[Ncut];	TLatex* T3[NGen];
c1 = new TCanvas("c1"," ",800,500);
c2 = new TCanvas("c2"," ",800,700);
c4 = new TCanvas("c4"," ",1000,800);
makeMultiPanelCanvas(c1,2,1,0.15,0.15,0.15,0.15,0.03);
makeMultiPanelCanvas(c2,1,1,0.15,0.1,0.15,0.1,0.03);
makeMultiPanelCanvas(c4,3,2,0.15,0.15,0.15,0.15,0.03);

for(int k=0;k<NGen;k++){
T3[k]=new TLatex(0.6,0.9,Gen[k]);
T3[k]->SetNDC();
T3[k]->SetTextAlign(12);
T3[k]->SetTextSize(0.06);
T3[k]->SetTextColor(1);
T3[k]->SetTextFont(42);
}
for(int j=0;j<Nhist;j++){
T1[j]=new TLatex(0.4,0.82,etarange[j]);
T1[j]->SetNDC();
T1[j]->SetTextAlign(12);
T1[j]->SetTextSize(0.06);
T1[j]->SetTextColor(1);
T1[j]->SetTextFont(42);
}
for(int i=0;i<Ncut;i++){
T2[i]=new TLatex(0.7,0.25,cutName[i]);
T2[i]->SetNDC();
T2[i]->SetTextAlign(12);
T2[i]->SetTextSize(0.08);
T2[i]->SetTextColor(1);
T2[i]->SetTextFont(42);
}
for(int k=0;k<NGen;k++){
for(int j=0;j<Nhist;j++){
for(int i=0;i<Ncut;i++){
if(i==Ncut-1)
filename[i][k]=Form("pPbHist_%s%s",Gen[k].Data(),cut[i].Data());
else
filename[i][k]=Form("pPbHist_%s_%s",Gen[k].Data(),cut[i].Data());
f[i][k]=TFile::Open(Form("%s.root",filename[i][k].Data()));
h_or[i][k][j]=(TH1F*)f[i][k]->Get(histoname[j]);
h_tr[i][k][j]=(TH1F*)f[i][k]->Get(Form("%s_tr",histoname[j].Data()));
fixedFontHist(h_or[i][k][j],1.8,2);
fixedFontHist(h_tr[i][k][j],1.8,2);
//reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(Ncent,Form("re%s%s_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin_cent[j]);
//reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(Ncent,Form("re%s%s_tr_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin_cent[j]);
reh_or[i][k][j]=(TH1F*)h_or[i][k][j]->Rebin(nHFbin,Form("re%s%s_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin);
reh_tr[i][k][j]=(TH1F*)h_tr[i][k][j]->Rebin(nHFbin,Form("re%s%s_tr_%d",Gen[k].Data(),histoname[j].Data(),i),HFbin);
normalizeByBinWidth(reh_or[i][k][j]);
normalizeByBinWidth(reh_tr[i][k][j]);
rehEff[i][k][j]=(TH1F*)reh_tr[i][k][j]->Clone(Form("rehEff_%s%s_%d",Gen[k].Data(),histoname[j].Data(),i));
rehEff[i][k][j]->Divide(reh_or[i][k][j]);
//rehEff[i][k][j]->GetXaxis()->SetTitleSize(0.06);
//rehEff[i][k][j]->GetXaxis()->SetLabelSize(0.06);
//rehEff[i][k][j]->GetYaxis()->SetTitleSize(0.06);
//rehEff[i][k][j]->GetYaxis()->SetLabelSize(0.06);
rehEff[i][k][j]->SetMarkerColor(color[k][j]);
rehEff[i][k][j]->SetMarkerStyle(Style[k][j]);
rehEff[i][k][j]->SetMarkerSize(1.5);
rehEff[i][k][j]->GetXaxis()->SetLimits(0.1,99);
rehEff[i][k][j]->SetMaximum(1.09);
rehEff[i][k][j]->SetMinimum(0.4);
rehEff[i][k][j]->SetTitle("");
rehEff[i][k][j]->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
rehEff[i][k][j]->GetYaxis()->SetTitle("Event Selection Efficiency");
for(int icentbin=0;icentbin<Ncent;icentbin++){
//EffValue[i][k][j][icentbin]=h_tr[i][k][j]->Integral(h_tr[i][k][j]->GetXaxis()->FindBin(HFbin_cent[j][Ncent-icentbin-1]),h_tr[i][k][j]->GetXaxis()->FindBin(HFbin_cent[j][Ncent-icentbin]))/h_or[i][k][j]->Integral(h_or[i][k][j]->GetXaxis()->FindBin(HFbin_cent[j][Ncent-icentbin-1]),h_or[i][k][j]->GetXaxis()->FindBin(HFbin_cent[j][Ncent-icentbin]));
EffValue[i][k][j][icentbin]=rehEff[i][k][j]->GetBinContent(Ncent-icentbin);
EffValueErr[i][k][j][icentbin]=rehEff[i][k][j]->GetBinError(Ncent-icentbin);
}
g_Eff[i][k][j]= new TGraphErrors(Ncent,cent,EffValue[i][k][j],xerr,EffValueErr[i][k][j]);
g_Eff[i][k][j]->SetTitle("");
g_Eff[i][k][j]->GetXaxis()->SetTitleSize(0.06);
g_Eff[i][k][j]->GetXaxis()->SetLabelSize(0.06);
g_Eff[i][k][j]->GetYaxis()->SetTitleSize(0.06);
g_Eff[i][k][j]->GetYaxis()->SetLabelSize(0.06);
g_Eff[i][k][j]->SetMarkerSize(0.06);
g_Eff[i][k][j]->GetYaxis()->SetTitle("Event Selection Efficiency");
g_Eff[i][k][j]->GetXaxis()->SetTitle("Centrality");
g_Eff[i][k][j]->GetXaxis()->CenterTitle();
g_Eff[i][k][j]->GetYaxis()->SetTitleOffset(1.2);
g_Eff[i][k][j]->GetXaxis()->SetTitleOffset(0.9);
g_Eff[i][k][j]->GetXaxis()->SetRangeUser(1,99);
g_Eff[i][k][j]->GetYaxis()->CenterTitle();
g_Eff[i][k][j]->SetMarkerSize(1.2);
g_Eff[i][k][j]->SetMarkerColor(color[k][j]);
g_Eff[i][k][j]->SetMarkerStyle(Style[k][j]);
g_Eff[i][k][j]->SetMinimum(0.4);
g_Eff[i][k][j]->SetMaximum(1.09);
}
}
}
/*
c1->cd(i+1)->SetLogy();
reh_or[i]->SetTitle("");
reh_or[i]->GetXaxis()->SetTitle("");
if(i==0)
reh_or[i]->GetYaxis()->SetTitle("Counts");
else
reh_tr[i]->GetYaxis()->SetTitle("");
reh_or[i]->SetMarkerStyle(24);
reh_or[i]->SetMarkerSize(1.2);
reh_or[i]->SetMarkerColor(1);
reh_or[i]->Draw();
reh_tr[i]->SetTitle("");
reh_tr[i]->SetMarkerSize(1.2);
reh_tr[i]->SetMarkerColor(2);
reh_tr[i]->SetMarkerStyle(20);
reh_tr[i]->Draw("same");
TLegend *leg=new TLegend(0.7,0.7,0.8,0.9);
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(reh_or[i],"Before Selection","lp");
leg->AddEntry(reh_tr[i],"After Selection","lp");
leg->Draw("same");

T1[i]->Draw("same");
*//*
i=Ncut-1;	//All cut
for(k=0;k<NGen;k++){
c1->cd(k+1);
if(k<1){
rehEff[i][k][j]->GetYaxis()->SetTitle("");
rehEff[i][k][j]->GetYaxis()->SetLabelSize(0);
}
TLegend *leg=new TLegend(0.2,0.2,0.35,0.4);
for(j=0;j<Nhist;j++){
if(j==0)
rehEff[i][k][j]->Draw();
else
rehEff[i][k][j]->Draw("same");
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(rehEff[i][k][j],etarange[j],"p");
TLine *l=new TLine(cent[0],1,cent[Ncent-1],1);
l->SetLineStyle(2);
}
l->Draw("same");
leg->Draw("same");
T3[k]->Draw("same");
}
c1->Print(Form("pPbHist_%s_HF.png",cut[i].Data()));
*/
i=Ncut-1;j=1;
//c2->cd()->SetLogy();
c2->cd()->SetLogx();
TLegend *leg=new TLegend(0.6,0.2,0.75,0.3);
for(k=0;k<NGen;k++){
rehEff[i][k][j]->GetXaxis()->SetTitleSize(30);
rehEff[i][k][j]->GetXaxis()->SetLabelSize(30);
rehEff[i][k][j]->GetYaxis()->SetTitleSize(30);
rehEff[i][k][j]->GetYaxis()->SetLabelSize(30);
rehEff[i][k][j]->GetYaxis()->SetTitleOffset(0.9);
rehEff[i][k][j]->GetXaxis()->SetTitleOffset(1.1);
rehEff[i][k][j]->GetYaxis()->SetTitleOffset(1);
if(k==0){
rehEff[i][k][j]->Draw();
rehEff[i][k][j]->SetMarkerColor(1);
rehEff[i][k][j]->SetMarkerStyle(20);
}
else{
rehEff[i][k][j]->Draw("same");
rehEff[i][k][j]->SetMarkerColor(2);
rehEff[i][k][j]->SetMarkerStyle(24);
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(rehEff[i][k][j],Gen[k],"p");
TLine *l=new TLine(0,1,cent[Ncent-1],1);
l->SetLineStyle(2);
}
T1[j]->Draw("same");
leg->Draw("same");
l->Draw("same");

//c2->Print(Form("pPbHist_%s%s_HF.png",cut[i].Data(),histoname[j].Data()));

j=1;k=0;
for(i=0;i<Ncut;i++){
//rehEff[i][k][j]->GetYaxis()->SetTitleOffset(2.5);
rehEff[i][k][j]->GetXaxis()->SetTitleOffset(2);
if(i!=0 &&  i!=3){
rehEff[i][k][j]->GetYaxis()->SetTitle("");
rehEff[i][k][j]->GetYaxis()->SetLabelSize(0);
}

if(i<3){
rehEff[i][k][j]->GetXaxis()->SetTitle("");
rehEff[i][k][j]->GetXaxis()->SetLabelSize(0);
}
c4->cd(i+1)->SetLogy();
//c4->cd(i+1)->SetLogx();
rehEff[i][k][j]->Draw();
rehEff[i][k][j]->SetMarkerColor(1);
rehEff[i][k][j]->SetMarkerStyle(20);
T2[i]->Draw("same");
}
T3[k]->Draw("same");
if(Ncut==5){
c4->cd(6);
TH1F *fNull=new TH1F("","",100,0,100);
fNull->GetXaxis()->SetLabelSize(0.06);
fNull->Draw();
}
if(Ncut==5)

//c4->Print(Form("pPbHist_%s_%s_newcut_HF.png",Gen[k].Data(),histoname[j].Data()));
else
//c4->Print(Form("pPbHist_%s_%s_HF.png",Gen[k].Data(),histoname[j].Data()));

}
