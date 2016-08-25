#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void P3overlayDataMC(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString listsvarALL[6]={"Max","Sum","Mult","Maxpt","Sumpt","other"};
TString JetIDNameList[18]={"chMax", "chSum", "neuMax", "neuSum", "phoMax", "phoSum", "chMaxpt", "chSumpt", "neuMaxpt", "neuSumpt", "phoMaxpt", "phoSumpt","SumSumpt","SumSumrawpt","neuMaxr","chN","neuN","phoN"};
const double binbound_pt[]={ 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 47, 55, 64,74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 429, 692, 1000};
double binbound_JetID[100];
int Nbin_pt=sizeof(binbound_pt)/sizeof(double)-1;
double xrange_pt[2]={30+1e-4,600-1e-4};
for(int j=0;j<6;j++){
TString svar=listsvarALL[j];
if(svar=="Max"){
int listsvar[]={0,2,4};
double ymin=1e-7,ymax=10;}

else if(svar=="Maxpt"){
int listsvar[]={6,8,10};
double ymin=1e-7,ymax=10;}

else if(svar=="Sum"){
int listsvar[]={1,3,5};
double ymin=1e-7,ymax=10;}

else if(svar=="Sumpt"){
int listsvar[]={7,9,11};
double ymin=1e-7,ymax=10;}

else if(svar=="Mult"){
int listsvar[]={15,16,17};
double ymin=1e-7;ymax=1;}

else if(svar=="other"){
int listsvar[]={12,13,14};
double ymin=1e-8;ymax=10;}
else exit();

int N=sizeof(listsvar)/sizeof(int);
c1 = new TCanvas("c1"," ",1000,400);
makeMultiPanelCanvas(c1,3,1,0,0,0.18,0.16,0.03);

for(int k=0;k<N;k++){
int ilist=listsvar[k];
c1->cd(k+1)->SetLogy();
TString JetID,Unit;  
TString JetIDName=JetIDNameList[ilist];

gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
if(ilist== 0)   JetID = "Max p_{T}^{charged}";
else if(ilist== 1)      JetID = "#Sigma p_{T}^{charged}";
else if(ilist== 2)      JetID = "Max p_{T}^{neutral}";
else if(ilist== 3) JetID = "#Sigma p_{T}^{neutral}";
else if(ilist== 4) JetID = "Max p_{T}^{photon}";
else if(ilist== 5) JetID = "#Sigma p_{T}^{photon}";
else if(ilist== 6)      JetID = "Max p_{T}^{charged}/p_{T}^{jet}";
else if(ilist== 7)      JetID = "#Sigma p_{T}^{charged}/p_{T}^{jet}";
else if(ilist== 8)      JetID = "Max p_{T}^{neutral}/p_{T}^{jet}";
else if(ilist== 9)      JetID = "#Sigma p_{T}^{neutral}/p_{T}^{jet}";
else if(ilist== 10)     JetID = "Max p_{T}^{photon}/p_{T}^{jet}";
else if(ilist== 11)     JetID = "#Sigma p_{T}^{photon}/p_{T}^{jet}";
else if(ilist== 12)     JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photon} +#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{jet}";
else if(ilist== 13)     JetID = "(#Sigma p_{T}^{charged}+#Sigma p_{T}^{photon} +#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e})/p_{T}^{raw}";
else if(ilist== 14)     JetID = "Max p_{T}^{neutral}/Max(#Sigma p_{T}^{charged}, #Sigma p_{T}^{neutral})";
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else{   exit();}
if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.}; int Nbin_JetID=20;}
else if(ilist==6){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.}; int Nbin_JetID=19;}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.}; int Nbin_JetID=20;}
Unit="";
}
else if(JetIDName.Contains("N")){
if(ilist==17){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30}; int Nbin_JetID=17;}
else if(ilist==16){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,10}; int Nbin_JetID=9;}
else if(ilist==15){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,30,40}; int Nbin_JetID=22;}
Unit="";
}
else{
//double binbound_JetID[]={0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300};
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
int Nbin_JetID=20;
Unit="(GeV/c)";
}

TString filename="/scratch/xuq7/RpA/TreeAna/Datacombined.root";
TString PPbfilename="/scratch/xuq7/RpA/TreeAna/MCPPbakPu3PF.root";

TString histoname=Form("jetpt%sCombinedSpectra",JetIDName.Data());
TString PPbhistoname=Form("jetpt%s",JetIDName.Data());

TFile *file=TFile::Open(filename);
TFile *PPbfile=TFile::Open(PPbfilename);

TH2F* hdata2F=(TH2F*)file->Get(histoname);
TH2F* hPPb2F=(TH2F*)PPbfile->Get(PPbhistoname);

TH1D* hdata = (TH1D*)hdata2F->ProjectionY("hdata",hdata2F->GetXaxis()->FindBin(xrange_pt[0]),hdata2F->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb = (TH1D*)hPPb2F->ProjectionY("hPPb",hPPb2F->GetXaxis()->FindBin(xrange_pt[0]),hPPb2F->GetXaxis()->FindBin(xrange_pt[1]));

rehdata=(TH1D*)hdata->Clone("rehdata");
rehPPb=(TH1D*)hPPb->Clone("rehPPb");
rehdata=(TH1D*)hdata->Rebin(Nbin_JetID,"rehdata",binbound_JetID);
rehPPb=(TH1D*)hPPb->Rebin(Nbin_JetID,"rehPPb",binbound_JetID);
normalizeByBinWidth(rehdata);
normalizeByBinWidth(rehPPb);
rehdata->Scale(1./rehdata->Integral(rehdata->FindBin(binbound_JetID[0]),rehdata->FindBin(binbound_JetID[Nbin_JetID])));
rehPPb->Scale(1./rehPPb->Integral(rehPPb->FindBin(binbound_JetID[0]),rehPPb->FindBin(binbound_JetID[Nbin_JetID])));
rehdata->SetMarkerStyle(20);
rehdata->SetMarkerColor(1);
rehdata->SetLineColor(1);
rehdata->SetMarkerSize(1.2);
rehPPb->SetMarkerStyle(0);
rehPPb->SetFillStyle(3004);
rehPPb->SetFillColor(2);
rehPPb->SetLineColor(2);

TH1F* hFrame=new TH1F("","",100000,0,1000);
hFrame->GetYaxis()->SetRangeUser(ymin,ymax);
if(k==0)	hFrame->GetYaxis()->SetTitle("Event Fraction");
else hFrame->GetYaxis()->SetTitle(""); 
hFrame->SetTitle("");
fixedFontHist(hFrame,1.2,1.4);
if(svar=="other")
hFrame->GetXaxis()->SetTitleSize(15);
hFrame->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(),Unit.Data()));
hFrame->GetXaxis()->SetLimits(binbound_JetID[0],binbound_JetID[Nbin_JetID]+1e-3);
hFrame->DrawCopy();
//if(ilist==12 || ilist==13) 
//rehdata->GetXaxis()->SetRangeUser(0.8,1.2);
rehdata->Draw("E1same");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.60,0.80,0.80,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.045);
leg1->AddEntry(rehdata,"Data","lpf");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.045);
T->SetTextColor(1);
T->SetTextFont(42);
if(k==0){
T->DrawLatex(0.20,0.90,Form("p_{T}^{jet}: %.f-%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
T->DrawLatex(0.40,0.85,Form("|#eta_{CM}| < 1"));
leg1->Draw("same");
}
/*c1->cd(2);
TH1D* ratio=(TH1D*)rehdata->Clone("ratio");
ratio->Divide(rehPPb);
ratio->SetMarkerStyle(20);
ratio->SetMarkerSize(1.2);
ratio->SetMarkerColor(1);
ratio->SetMinimum(0);
ratio->SetMaximum(1.9);
ratio->SetTitle("");
ratio->GetYaxis()->SetTitle("Data/MC");
ratio->Draw("E1");
TLine *l =new TLine(binbound_JetID[0],1,binbound_JetID[Nbin_JetID],1);
l->SetLineStyle(2);
l->SetLineColor(1);
l->Draw("same");
*/
}
c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/TreeAna/JetID/pic/overlay_%s_DataMC.png",svar.Data()));
c1->Print(Form("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/TreeAna/JetID/pic/overlay_%s_DataMC.pdf",svar.Data()));
}
}
