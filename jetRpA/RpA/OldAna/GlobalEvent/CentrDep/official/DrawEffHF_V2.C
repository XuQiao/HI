#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void DrawEffHF_V2(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
bool isweight=kFALSE;
//bool isweight=kTRUE;
bool isEvo=kFALSE;
//bool isEvo=kTRUE;
if(isweight)
TString dir="weight";
else
TString dir="noweight";
const int Ncut=6;	const int Nhist=3;	const int NGen=2; const int Ncent=6;
const double HFbin[]={0,0.1,0.2,0.4,0.6,0.8,1,2,4,6,9,12,15,20,30,50,75,100};
const double xerr[Ncent]={0,0,0,0,0,0};
const double HFbin_cent[Nhist][Ncent+1]={
{0,2.92,10.1,19.3,23.7,30.2,100},
{0,1.03,5.00,11.6,15.0,20.3,100},
{0,2.7,12.6,28.6,36.9,49.6,100}};
const double cent[Ncent]={95,75,45,25,15,5};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;
double color[NGen][Nhist]={{2,1,4},{2,1,4}};
double Style[NGen][Nhist]={{20,24,28},{20,24,28}};
if(!isEvo){
TString cut[Ncut]={"HLT","CutGplus","pPA","pprim","vz",""};
TString cutName[Ncut]={"HLT","CutG","pPA","pprim","vz","All cuts"};
}
else{
TString cut[Ncut]={"HLT","HLT_CutGplus","HLT_CutGplus_pPA","HLT_CutGplus_pPA_pprim","",""};
TString cutName[Ncut]={"HLT","HLT+CutGplus","HLT+CutGplus+pPA","#splitline{HLT+CutG}{+pPA+pprim}","#splitline{HLT+CutG+pPA}{+pprim+vz}",""};
}
TString Gen[NGen]={"Hijing","Epos"};
TString histoname[Nhist]={"hHFEnergy4","hHFEnergyPlus4","hHFEnergyPlus"};
TString etarange[Nhist]={" HF E_{T} |#eta|>4"," HF E_{T} #eta>4"," #splitline{HF E_{T}#eta>0}{(Pb going side)}"};
TString filename[Ncut][NGen];
TFile* f[Ncut][NGen];
TH1F* h_or[Ncut][NGen][Nhist];	TH1F* reh_or[Ncut][NGen][Nhist];
TH1F* h_tr[Ncut][NGen][Nhist];	TH1F* reh_tr[Ncut][NGen][Nhist];
TH1F* rehEff[Ncut][NGen][Nhist];	TLatex* T1[Nhist];	TLatex* T2[Ncut];	TLatex* T3[NGen];
c1 = new TCanvas("c1"," ",800,500);
c2 = new TCanvas("c2"," ",800,700);
c3 = new TCanvas("c3"," ",1000,800);
makeMultiPanelCanvas(c1,2,1,0.02,0.02,0.12,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,0.03,0.03,0.03,0.03,0.03);
makeMultiPanelCanvas(c3,3,2,0.02,0.02,0.15,0.15,0.03);

for(int k=0;k<NGen;k++){
T3[k]=new TLatex(0.45,0.9,Gen[k]);
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
T2[i]=new TLatex(0.2,0.25,cutName[i]);
T2[i]->SetNDC();
T2[i]->SetTextAlign(12);
T2[i]->SetTextSize(0.08);
T2[i]->SetTextColor(1);
T2[i]->SetTextFont(42);
}
for(int k=0;k<NGen;k++){
for(int j=0;j<Nhist;j++){
for(int i=0;i<Ncut;i++){
if(cut[i]=="")
filename[i][k]=Form("output/%s/pPbHist_%s%s",dir.Data(),Gen[k].Data(),cut[i].Data());
else
filename[i][k]=Form("output/%s/pPbHist_%s_%s",dir.Data(),Gen[k].Data(),cut[i].Data());
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
rehEff[i][k][j]->SetMarkerColor(color[k][j]);
rehEff[i][k][j]->SetMarkerStyle(Style[k][j]);
rehEff[i][k][j]->SetMarkerSize(1.2);
rehEff[i][k][j]->GetXaxis()->SetLimits(0.1,99);
rehEff[i][k][j]->SetMaximum(1.19);
rehEff[i][k][j]->SetMinimum(0.4);
rehEff[i][k][j]->SetTitle("");
rehEff[i][k][j]->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
rehEff[i][k][j]->GetYaxis()->SetTitle("Event Selection Efficiency");
}
}
}
/*
i=Ncut-1;	//All cut
for(k=0;k<NGen;k++){
c1->cd(k+1)->SetLogx();
TLegend *leg=new TLegend(0.40,0.2,0.55,0.4);
*/
 TH1D * hFrame = new TH1D("","",100,0,100.);
       hFrame-> SetAxisRange(0,100.,"X");
    hFrame->SetAxisRange(0.35, 1.2, "Y") ;
   hFrame->GetXaxis()->SetLimits(0.1,99.9);
   hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
   hFrame->GetYaxis()->SetTitle("Event Selection Efficiency");
   hFrame->GetYaxis()->SetTitleSize(0.05);
	fixedFontHist(hFrame);
  hFrame->GetYaxis()->SetTitleOffset(1.2);/*
hFrame6=(TH1D*)hFrame->Clone("hFrame6");
hFrame6->GetYaxis()->SetTitle("");
hFrame6->GetYaxis()->SetLabelSize(0);
for(j=0;j<Nhist;j++){
if(j==0){
if(k==0)
hFrame->Draw();
else
hFrame6->Draw();
rehEff[i][k][j]->Draw("same");
}
else{
rehEff[i][k][j]->Draw("same");
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(rehEff[i][k][j],etarange[j],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
l->Draw("same");
leg->Draw("same");
T3[k]->Draw("same");
}
c1->Print(Form("pic/pPbHist_%s_%s_HF.png",dir.Data(),cut[i].Data()));
*/
//------------------------------------------------------------------------------------------

i=Ncut-1;j=1;
c2->cd()->SetLogx();
TLegend *leg=new TLegend(0.2,0.6,0.35,0.7);
for(k=0;k<NGen;k++){
rehEff[i][k][j]->GetXaxis()->SetTitleSize(0.05);
rehEff[i][k][j]->GetXaxis()->SetLabelSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetTitleSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetLabelSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetTitleOffset(0.9);
rehEff[i][k][j]->GetXaxis()->SetTitleOffset(0.88);
rehEff[i][k][j]->GetXaxis()->SetLimits(0,100);
hFrame1=(TH1D*)hFrame->Clone("hFrame1");
hFrame1->SetMaximum(1.2);
hFrame1->SetMinimum(0.35);
if(k==0){
hFrame1->Draw();
rehEff[i][k][j]->SetMarkerColor(1);
rehEff[i][k][j]->SetLineColor(1);
rehEff[i][k][j]->SetMarkerStyle(20);
rehEff[i][k][j]->SetMarkerSize(1.6);
rehEff[i][k][j]->Draw("same");
}
else{
rehEff[i][k][j]->SetMarkerColor(2);
rehEff[i][k][j]->SetLineColor(2);
rehEff[i][k][j]->SetMarkerStyle(24);
hFrame1->GetYaxis()->SetTitle("");
hFrame1->GetYaxis()->SetLabelSize(0);
rehEff[i][k][j]->Draw("same");
}
leg->SetTextSize(0.05);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->AddEntry(rehEff[i][k][j],Gen[k],"p");
TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
}
T1[j]->Draw("same");
leg->Draw("same");
l->Draw("same");

//c2->Print(Form("pic/pPbHist_%s_%s_%s_HF.png",dir.Data(),cut[i].Data(),histoname[j].Data()));

//----------------------------------------------------------------------------------

j=1;k=1;
for(i=0;i<Ncut;i++){
rehEff[i][k][j]->GetXaxis()->SetTitleSize(0.05);
rehEff[i][k][j]->GetXaxis()->SetLabelSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetTitleSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetLabelSize(0.05);
rehEff[i][k][j]->GetYaxis()->SetTitleOffset(1.5);
rehEff[i][k][j]->GetXaxis()->SetLimits(0.01,99.99);
hFrame2=(TH1D*)hFrame->Clone("hFrame2");
hFrame2->GetYaxis()->SetTitleOffset(2.2);
hFrame2->GetXaxis()->SetTitleOffset(1.8);
hFrame2->SetMaximum(1.19);
hFrame2->SetMinimum(0.35);
hFrame3=(TH1D*)hFrame2->Clone("hFrame3");
hFrame4=(TH1D*)hFrame2->Clone("hFrame4");
hFrame5=(TH1D*)hFrame2->Clone("hFrame5");

rehEff[i][k][j]->SetMarkerStyle(20);
rehEff[i][k][j]->SetMarkerColor(1);
c3->cd(i+1)->SetLogx();
if(i==0){
hFrame5->Draw();
rehEff[i][k][j]->Draw("same");
}
else if(i==1 || i==2){
hFrame3->Draw();
rehEff[i][k][j]->GetXaxis()->SetTitle("");
hFrame3->GetXaxis()->SetTitle("");
hFrame3->GetYaxis()->SetTitle("");
rehEff[i][k][j]->GetXaxis()->SetLabelSize(0);
hFrame3->GetXaxis()->SetLabelSize(0);
hFrame3->GetYaxis()->SetLabelSize(0);
rehEff[i][k][j]->Draw("same");
}
else if(i==3){
hFrame2->Draw();
rehEff[i][k][j]->Draw("same");
}
else{
hFrame4->Draw();
hFrame4->GetYaxis()->SetLabelSize(0);
hFrame4->GetYaxis()->SetTitle("");
rehEff[i][k][j]->Draw("same");
}

TLine *l=new TLine(0,1,100,1);
l->SetLineStyle(2);
l->Draw("same");
T2[i]->Draw("same");
}
T3[k]->Draw("same");

if(isEvo)
//c3->Print(Form("pic/pPbHist_%s_%s_%s_HF_Evo.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
else{}
//c3->Print(Form("pic/pPbHist_%s_%s_%s_HF.png",dir.Data(),Gen[k].Data(),histoname[j].Data()));
}
