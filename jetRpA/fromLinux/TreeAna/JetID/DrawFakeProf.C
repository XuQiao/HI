#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

//const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,1000};
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
int ilist=12;//6,12,14
double JetIDcut[2];
double xrange_JetIDcut[2];
const int Neta=8;
const double binbound_pt_pro[]={30,70,100,150,600};
const int Nbin_pt_pro=4;

const TString etabinname[Neta]={"15_20","10_15","5_10","-5_5","-10_-5","-15_-10","-20_-15","-10_10"};
const TString etabinnameswaped[Neta]={"-20_-15","-15_-10","-10_-5","-5_5","5_10","10_15","15_20","-10_10"};
const double etabin[Neta]={0.5,0.5,0.5,1.0,0.5,0.5,0.5,2.0};
const TString etastring[Neta]={"-2.0<#eta_{CM}<-1.5","-1.5<#eta_{CM}<-1.0","-1.0<#eta_{CM}<-0.5","-0.5<#eta_{CM}<0.5","0.5<#eta_{CM}<1.0","1.0<#eta_{CM}<1.5","1.5<#eta_{CM}<2.0","-1.0<#eta_{CM}<1.0"};

TProfile* makehisto(int ilist, int ieta){
TString JetID;
TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
TString JetIDName = JetIDNameList[ilist];
int ilistw;
if(ilist== 6)  {    JetID = "Max p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0]=0; JetIDcut[1]=1; ilistw=1;}
else if(ilist== 7)  {    JetID = "#Sum p_{T}^{charged}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sum p_{T}^{neutral}/p_{T}^{jet}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sum p_{T}^{photons}/p_{T}^{jet}";     JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}";JetIDcut[0]=0; JetIDcut[1]=2; ilistw=3;}
else if(ilist== 13) {    JetID = "#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})";JetIDcut[0]=0;JetIDcut[1]=2;ilistw=2;}

if(ilist==12 || ilist==13 || ilist==14 ){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
//double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
TString histonameID=Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinname[ieta].Data());
TH2F* hdata2F->Scale(1./etabin[ieta]);
TH2F* hdata2F=(TH2F*)fMCPPbJetID->Get(histonameID);
//cout<<Nbin_JetID<<endl;

TProfile* histo=(TProfile*)hdata2F->ProfileX("histo");
histo=(TProfile*)histo->Rebin(Nbin_pt_pro,"histo",binbound_pt_pro);
histo->SetMarkerSize(1.2);
histo->SetTitle(JetID);
return histo;
}

void DrawFakeProf(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);

TH1F* hFrame=new TH1F("","",1000,0,1000);
c1 = new TCanvas("c1","",600,600);
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.17,0.12,0.03);
fixedFontHist(hFrame,1.5,1.8);
hFrame->SetTitle("");
hFrame->GetYaxis()->SetTitleSize(0.035);
hFrame->GetXaxis()->SetRangeUser(29,600);
hFrame->GetYaxis()->SetRangeUser(0.7,2);
TLegend *leg1=new TLegend(0.18,0.70,0.40,0.95);
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
TProfile* histo0[Neta];
for(int i=0;i<Neta;i++){
histo0[i] = makehisto(ilist,i);
}
//TH1D* histo0 = makehisto(ilist2,i,0.975);
//TH1D* histo0 = makehisto(ilist1,i,0.05);


c1->cd(1);
hFrame->GetYaxis()->SetTitle(histo0[0]->GetTitle());
hFrame->GetXaxis()->SetTitle("p_{T} [GeV/c]");
hFrame->DrawCopy();

//leg1->SetHeader(Form("ratio of %s",etastring[Neta-1].Data()));
for(int i=0;i<Neta;i++){
histo0[i]->SetMarkerStyle(markerCode[i]);
histo0[i]->SetMarkerColor(colorCode[i]);
histo0[i]->SetLineColor(colorCode[i]);
/*
ratio[i][j]=(TH1D*)histo0[i][j]->Clone(Form("ratio_%d_%d",i,j));
if(i!=Neta-1){
ratio[i][j]->Divide(histo0[Neta-1][j]);
ratio[i][j]->Draw("E1same");
if(j==0) leg1->AddEntry(ratio[i][j],Form("%s",etastring[i].Data()),"lp");
}
*/
histo0[i]->Draw("same");
leg1->AddEntry(histo0[i],etastring[i],"lp");

}
T->DrawLatex(0.65,0.85,"Fake jet");
leg1->Draw("same");
//l->Draw("same");
c1->Print(Form("FakeProfPPb%s.pdf",JetIDNameList[ilist].Data()));
c1->Print(Form("FakeProfPPb%s.png",JetIDNameList[ilist].Data()));
}
