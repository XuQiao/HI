#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
void overlayDataMC(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
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

if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else if(ilist==6){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.};}
Unit="";
}
else if(JetIDName.Contains("N")){
if(ilist==17){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30};}
else if(ilist==16){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,10};}
else if(ilist==15){
double binbound_JetID[]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,30,40};}
Unit="";
}
else{
//double binbound_JetID[]={0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300};
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
Unit="(GeV/c)";
}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TH1D* hdata=(TH1D*)hdata2F->ProjectionY("hdata",hdata2F->GetXaxis()->FindBin(xrange_pt[0]),hdata2F->GetXaxis()->FindBin(xrange_pt[1]));
TH1D* hPPb=(TH1D*)hPPb2F->ProjectionY("hPPb",hPPb2F->GetXaxis()->FindBin(xrange_pt[0]),hPPb2F->GetXaxis()->FindBin(xrange_pt[1]));

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
rehdata->GetYaxis()->SetTitle("Event Fraction");
rehPPb->SetMarkerStyle(0);
rehPPb->SetFillStyle(3004);
rehPPb->SetFillColor(2);
rehPPb->SetLineColor(2);
fixedFontHist(rehdata,1.8,2.6);
fixedFontHist(rehPPb,1.8,2.6);

c1 = new TCanvas("c1","",400,700);
makeMultiPanelCanvas(c1,1,2,-0.12,0,0.18,0.16,0.05);

c1->cd(1)->SetLogy();
rehdata->GetXaxis()->SetTitle(JetID.Data());
rehdata->GetXaxis()->SetRangeUser(binbound_JetID[0],binbound_JetID[Nbin_JetID]);
//if(ilist==12 || ilist==13) 
//rehdata->GetXaxis()->SetRangeUser(0.8,1.2);
rehdata->SetMinimum(TMath::Min(rehdata->GetMinimum(),rehPPb->GetMinimum())/3);
rehdata->SetMaximum(1.5);
rehdata->SetTitle("");
rehdata->Draw("E1");
rehPPb->Draw("HIST E1same");
TLegend *leg1=new TLegend(0.60,0.70,0.80,0.85);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(rehdata,"Data","lpf");
leg1->AddEntry(rehPPb,"PYTHIA+HIJING","lpf");
leg1->Draw("same");
TLatex *T=new TLatex(0.40,0.90,Form("p_{T}^{jet}: %.f-%.f (GeV/c)",xrange_pt[0],xrange_pt[1]));
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
T->Draw("same");
c1->cd(2);
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

//c1->Print(Form("pic/JetIDcut/%s/%s/overlay_DataMC.png",algo.Data(),JetIDName.Data()));

}
