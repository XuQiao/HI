#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

//const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
int ilist0=12, ilist1=6, ilist2=14;
double JetIDcut[2];
double xrange_JetIDcut[2];
const int Neta=8;
const int binbound_pt_pro[]={30,70,100,150,600};
const int Nbin_pt_pro=4;

const TString etabinname[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etabinnameswaped[Neta]={"-20_-15","-15_-10","-10_-5","-5_5","5_10","10_15","15_20","-10_10"};
const double etabin[Neta]={0.5,0.5,0.5,1.0,0.5,0.5,0.5,2.0};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};

TH1D* makehisto(int ilist, int ieta,double cut,int ptmin, int ptmax){
TString JetID;
TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
TString JetIDName = JetIDNameList[ilist];

int ilistw;
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=cut; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}";JetIDcut[0]=0; JetIDcut[1]=cut; ilistw=3;}
else if(ilist== 13) {    JetID = "#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=cut;ilistw=2;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TString histonameID=Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinname[ieta].Data());
TH2F* hdata2F=(TH2F*)fMCPPbxSecJetID->Get(histonameID);
//cout<<Nbin_JetID<<endl;

TH1D* histo=(TH1D*)hdata2F->ProjectionY("histo",hdata2F->GetXaxis()->FindBin(ptmin),hdata2F->GetXaxis()->FindBin(ptmax));
histo=(TH1D*)histo->Rebin(Nbin_JetID,"histo",binbound_JetID);
histo->Draw();
normalizeByBinWidth(histo);
histo->SetMarkerSize(1.2);
histo->SetTitle(JetID);
return histo;
}

void DrawFakeDis(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TH1F* hFrame=new TH1F("","",1000,0,1000);
c1 = new TCanvas("c1","",800,800);
makeMultiPanelCanvas(c1,2,2,0,0,0.15,0.20,0.03);
fixedFontHist(hFrame,1.9,1.4);
hFrame->SetTitle("");
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetRangeUser(0.01,1.95);
TLegend *leg1=new TLegend(0.18,0.60,0.40,0.90);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.04);
leg1->Draw("same");
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
TLine *l = new TLine(0,1,2,1);
l->SetLineStyle(2);
TH1D* histo0[Neta][Nbin_pt_pro];
TH1D* ratio[Neta][Nbin_pt_pro];
for(int i=0;i<Neta;i++){
for(int j=0;j<Nbin_pt_pro;j++){
histo0[i][j] = makehisto(ilist0,i,1.01,binbound_pt_pro[j],binbound_pt_pro[j+1]);
histo0[i][j]->Scale(1./etabin[i]);
}
}
//TH1D* histo0 = makehisto(ilist2,i,0.975);
//TH1D* histo0 = makehisto(ilist1,i,0.05);


for(int j=0;j<Nbin_pt_pro;j++){
c1->cd(j+1)->SetLogy();
//c1->cd(j+1);
if(j==0){
//hFrame->GetYaxis()->SetRangeUser(0.01,10);
hFrame->GetYaxis()->SetRangeUser(6e-10,1);
hFrame->GetYaxis()->SetTitle("d#sigma/dp_{T}d#eta");
}
else if(j==1){
hFrame->GetYaxis()->SetTitle("");
hFrame->GetXaxis()->SetTitle("");
}
else if(j==2){
hFrame->GetYaxis()->SetRangeUser(1e-11,9e-2);
hFrame->GetYaxis()->SetTitle("d#sigma/dp_{T}");
hFrame->GetXaxis()->SetTitle(histo0[0][0]->GetTitle());
}
else if(j==3){
hFrame->GetYaxis()->SetTitle("");
hFrame->GetXaxis()->SetTitle(histo0[0][0]->GetTitle());
}
hFrame->DrawCopy();

//leg1->SetHeader(Form("ratio of %s",etastring[Neta-1].Data()));
for(int i=0;i<Neta;i++){
histo0[i][j]->SetMarkerStyle(markerCode[i]);
histo0[i][j]->SetMarkerColor(colorCode[i]);
histo0[i][j]->SetLineColor(colorCode[i]);
/*
ratio[i][j]=(TH1D*)histo0[i][j]->Clone(Form("ratio_%d_%d",i,j));
if(i!=Neta-1){
ratio[i][j]->Divide(histo0[Neta-1][j]);
ratio[i][j]->Draw("E1same");
if(j==0) leg1->AddEntry(ratio[i][j],Form("%s",etastring[i].Data()),"lp");
}
*/
histo0[i][j]->Draw("E1same");
if(j==0) leg1->AddEntry(histo0[i][j],etastring[i],"lp");

}
T->DrawLatex(0.65,0.9,Form("%d < p_{T} < %d",binbound_pt_pro[j],binbound_pt_pro[j+1]));
if(j==0) T->DrawLatex(0.65,0.85,"Fake jet");
if(j==0) leg1->Draw("same");
l->Draw("same");
}
c1->Print(Form("FakeDis.pdf"));
}
