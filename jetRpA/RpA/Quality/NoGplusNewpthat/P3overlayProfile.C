#include "file.h"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void P3overlayProfile(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TString svar="Mult";
if(svar=="Max"){
int listsvar[]={0,2,4,6,8,10};
double ymin1=0.1,ymax1=150;
double ymin2=0,ymax2=0.37;
TString svar1="Max (GeV/c)";
TString svar2="Max / p_{T}";}

else if(svar=="Sum"){
int listsvar[]={1,3,5,7,9,11};
double ymin1=0.1,ymax1=260;
double ymin2=0,ymax2=0.72;
TString svar1="Sum (GeV/c)";
TString svar2="Sum / p_{T}";}

else if(svar=="Mult"){
int listsvar[]={15,16,17};
double ymin1=0;ymax1=25;
TString svar1="Mult (GeV/c)";}
else exit();

int N=sizeof(listsvar)/sizeof(int);
if(N==6){
c1 = new TCanvas("c1"," ",1000,650);
makeMultiPanelCanvas(c1,3,2,0,0,0.22,0.16,0.03);
}
else if(N==3){
c1 = new TCanvas("c1"," ",1000,400);
makeMultiPanelCanvas(c1,3,1,0,0,0.16,0.14,0.03);
}
else exit();

for(int k=0;k<sizeof(listsvar)/sizeof(int);k++){
int ilist=listsvar[k];
c1->cd(k+1);
TString JetID,Unit;     TString JetIDName=JetIDNameList[ilist];
TString histonameID=Form("jetpt%s_0-100%%",JetIDName.Data());
TH2F* hdata2F=(TH2F*)fdata->Get(histonameID);
TH2F* hPPb2F=(TH2F*)fPPb->Get(histonameID);
TH2F* hPP2F=(TH2F*)fPP->Get(histonameID);
if(ilist== 0)   JetID = "chargedMax";
else if(ilist== 1)      JetID = "chargedSum";
else if(ilist== 2)      JetID = "neutralMax";
else if(ilist== 3) JetID = "neutralSum";
else if(ilist== 4) JetID = "photonMax";
else if(ilist== 5) JetID = "photonSum";
else if(ilist== 6)      JetID = "chargedMax/p_{T}";
else if(ilist== 7)      JetID = "chargedSum/p_{T}";
else if(ilist== 8)      JetID = "neutralMax/p_{T}";
else if(ilist== 9)      JetID = "neutralSum/p_{T}";
else if(ilist== 10)     JetID = "photonMax/p_{T}";
else if(ilist== 11)     JetID = "photonSum/p_{T}";
else if(ilist== 12)     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";
else if(ilist== 13)     JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}";
else if(ilist== 14)     JetID = "neutralMax/(neutralMax+chargedMax+photonMax)";
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else{   exit();}

TProfile* hdata=(TProfile*)hdata2F->ProfileX("hdata");
TProfile* hPPb=(TProfile*)hPPb2F->ProfileX("hPPb");
rehdata=(TProfile*)hdata->Clone("rehdata");
rehPPb=(TProfile*)hPPb->Clone("rehPPb");
rehdata=(TProfile*)hdata->Rebin(Nbin_pt,"rehdata",binbound_pt);
rehPPb=(TProfile*)hPPb->Rebin(Nbin_pt,"rehPPb",binbound_pt);
rehPPb->SetMarkerStyle(0);
rehPPb->SetMarkerColor(0);
rehPPb->SetLineColor(2);
rehPPb->SetFillStyle(3004);
rehPPb->SetFillColor(2);
rehdata->SetMarkerStyle(20);
rehdata->SetMarkerSize(1.2);
rehdata->SetMarkerColor(1);
//for(int i=0;i<rehdata->GetNbinsX();i++){
//cout<<rehdata->GetBinContent(i)<<endl;
//}
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(xrange_pt_pre[0],xrange_pt_pre[1]);
if(JetID.Contains("p_{T}")&&svar!="Mult"){
hFrame->GetYaxis()->SetRangeUser(ymin2,ymax2);
hFrame->GetYaxis()->SetTitle(Form("%s",svar2.Data()));}
else{
hFrame->GetYaxis()->SetRangeUser(ymin1,ymax1);
hFrame->GetYaxis()->SetTitle(Form("%s",svar1.Data()));}
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->SetTitle("");
fixedFontHist(hFrame,1.8,2.2);
hFrame->DrawCopy();
TLatex *T1=new TLatex(0.25,0.88,JetID.Data());
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.06);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
rehdata->Draw("E1same");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.60,0.8,0.80,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(rehdata,"Data","lp");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
if(k==0) leg1->Draw("same");
}
c1->Print(Form("pic/JetIDcut/%s/overlay_Profile_%s.gif",algo.Data(),svar.Data()));
c1->Print(Form("pic/JetIDcut/%s/overlay_Profile_%s.pdf",algo.Data(),svar.Data()));
}
