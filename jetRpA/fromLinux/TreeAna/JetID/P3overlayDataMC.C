#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"
void P3overlayDataMC(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
int ieta=7;
TString coll="PPb";
TString listsvarALL[6]={"Max","Sum","Mult","Maxpt","Sumpt","other"};
double binbound_JetID[100];
for(int j=0;j<6;j++){
TString svar=listsvarALL[j];
if(svar=="Max"){
int listsvar[]={0,2,4};
double ymin=1e-7,ymax=10;
double margin=8;
}

else if(svar=="Maxpt"){
int listsvar[]={6,8,10};
double ymin=1e-7,ymax=10;
double margin=3e-2;
}
else if(svar=="Sum"){
int listsvar[]={1,3,5};
double ymin=1e-7,ymax=10;
double margin=8;
}

else if(svar=="Sumpt"){
int listsvar[]={7,9,11};
double ymin=1e-7,ymax=10;
double margin=3e-2;
}

else if(svar=="Mult"){
int listsvar[]={15,16,17};
double ymin=1e-7;ymax=10;
double margin=1;
}

else if(svar=="other"){
int listsvar[]={12,13,14};
double ymin=1e-7;ymax=10;
double margin=3e-2;
}
else exit();

int N=sizeof(listsvar)/sizeof(int);
c1 = new TCanvas("c1"," ",1000,400);
makeMultiPanelCanvas(c1,3,1,0,0,0.18,0.16,0.03);

for(int k=0;k<N;k++){
int ilist=listsvar[k];
c1->cd(k+1)->SetLogy();
TString JetID,Unit;  
TString JetIDName=JetIDNameList[ilist];

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

if(coll=="PPb"){
TString histoname=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data());
TH2F* hdata2F=(TH2F*)fDataPPbJetID->Get(histoname);
TH2F* hPPb2F=(TH2F*)fMCPPbJetID->Get(histoname);
}
if(coll=="PbP"){
TString histoname=Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinname[ieta].Data());
TH2F* hdata2F=(TH2F*)fDataPbPJetIDOld->Get(histoname);
TH2F* hPPb2F=(TH2F*)fMCPbPJetID->Get(histoname);
}


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
fixedFontHist(hFrame,1.2,1.5);
hFrame->GetXaxis()->SetTitleSize(0.05);
hFrame->GetXaxis()->SetLabelSize(0.05);
hFrame->GetYaxis()->SetTitleSize(0.05);
hFrame->GetYaxis()->SetLabelSize(0.05);
if(k==0)	hFrame->GetYaxis()->SetTitle("Event Fraction");
else hFrame->GetYaxis()->SetTitle(""); 
hFrame->SetTitle("");
if(svar=="other")
hFrame->GetXaxis()->SetTitleSize(0.04);
hFrame->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(),Unit.Data()));
hFrame->GetXaxis()->SetLimits(binbound_JetID[0]-margin,binbound_JetID[Nbin_JetID]+margin);
hFrame->DrawCopy();
//if(ilist==12 || ilist==13) 
//rehdata->GetXaxis()->SetRangeUser(0.8,1.2);
rehdata->Draw("E1same");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.35,0.80,0.65,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.048);
if(coll=="PPb"){
leg1->AddEntry(rehdata,"Pb going positive Data","lpf");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
}
else if(coll=="PbP"){
leg1->AddEntry(rehdata,"Proton going positive Data","lpf");
leg1->AddEntry(rehPPb,"HIJING+PYTHIA","lpf");
}
TLatex *T=new TLatex();
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.048);
T->SetTextColor(1);
T->SetTextFont(42);
if(k==0)
T->DrawLatex(0.20,0.92,Form("p_{T}^{jet}: %.f-%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
if(k==1)
T->DrawLatex(0.25,0.92,Form("|#eta_{CM}| < 1"));
if(k==2)
leg1->Draw("same");
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
c1->Print(Form("pic/%soverlay_%s_DataMC.png",coll.Data(),svar.Data()));
c1->Print(Form("pic/%soverlay_%s_DataMC.pdf",coll.Data(),svar.Data()));
}
}
