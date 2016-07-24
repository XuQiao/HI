#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include "iostream"
#include "fstream"
#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"

TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
normalizeByBinWidth(rehisto);
}

return rehisto;
}

void JetID2D(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------
int isMatch=1;
int ilist1=12, ilist2=14;
TString JetID1="(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";
TString JetID2="Max h^{0}/Max(#Sigma h^{0},#Sigma h^{#pm})";
/*if(ilist== 0)  { JetID = "chargedMax"; JetIDcut[0]=4; JetIDcut[1]=9999;}
else if(ilist== 1)      JetID = "chargedSum";
else if(ilist== 2)      JetID = "neutralMax";
else if(ilist== 3) JetID = "neutralSum";
else if(ilist== 4) JetID = "photonMax";
else if(ilist== 5) JetID = "photonSum";
else if(ilist== 6)  {    JetID = "Max h^{#pm}/p_{T}"; JetIDcut[0]=0.05; JetIDcut[1]=2;}
else if(ilist== 7)  {    JetID = "#Sigma h^{#pm}/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "Max h^{0}/p_{T}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "#Sigma h^{0}/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "Max #gamma/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "#Sigma #gamma/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{jet}";JetIDcut[0]=0; JetIDcut[1]=1.01;}
else if(ilist== 13) {    JetID = "(#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e)/p_{T}^{raw}"; JetIDcut[0]=0;JetIDcut[1]=1.20;}
else if(ilist== 14)   {  JetID = "Max h^{0}/Max(#Sigma h^{0},#Sigma h^{#pm})";JetIDcut[0]=0;JetIDcut[1]=0.975;}
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else if(ilist== 18)   {  JetID = "PP cut True or False";JetIDcut[0]=1;JetIDcut[1]=2;}
else if(ilist== 19)   {  JetID = "PP cut Tight True or False";JetIDcut[0]=1;JetIDcut[1]=2;}
else{   exit();}
if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.4,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else if(ilist==14){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.84,0.88,0.92,0.96,1.02};}
else if(ilist==8 || ilist==10){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.,1.1,1.5,2};}

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
else if(JetIDName.Contains("PP")){
double binbound_JetID[]={0,1,2};
Unit="";
}
else{
double binbound_JetID[]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,200};
Unit="(GeV/c)";
}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
*/

//--------------------------------------------------------SetUp canvas---------------------------------------------------
c1 = new TCanvas("c1"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.15);
TH2F* hFrame=new TH2F("","",200,0,2,200,0,2);
fixedFontHist(hFrame,2.0,3.0);
hFrame->SetTitle("");
TLegend *leg1=new TLegend(0.36,0.82,0.60,0.95);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.065);
TLatex T1;
T1.SetNDC();
T1.SetTextSize(0.065);
T1.SetTextFont(42);
int i=7;

for(int ipt=0;ipt<Nbin_pt_coarse;ipt++){ 
if(isMatch){
    if(i==7){
	TH2F *hJetID2D= (TH2F*)fMCPPb->Get(Form("jetid12jetid14_%d",ipt));
	TH2F *hJetID2D_real= (TH2F*)fMCPPb->Get(Form("jetid12jetid14_%d_real1",ipt));
	TH2F *hJetID2D_fake= (TH2F*)fMCPPb->Get(Form("jetid12jetid14_%d_fake1",ipt));
    }
    else{
	TH2F *hPPb2D= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[i].Data()));
	TH2F *hPPb2D_real= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s_real1",JetIDName.Data(),etabinnameswap[i].Data()));
	TH2F *hPPb2D_fake= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s_fake1",JetIDName.Data(),etabinnameswap[i].Data()));
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_real1",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_fake1",etabinnameswap[i].Data()));
    }
    }   
else{
    if(i==7){
	TH2F *hPPb2D= (TH2F*)fMCPPb->Get(Form("jetpt%s",JetIDName.Data()));
	TH2F *hPPb2D_real= (TH2F*)fMCPPb->Get(Form("jetpt%s_real",JetIDName.Data()));
	TH2F *hPPb2D_fake= (TH2F*)fMCPPb->Get(Form("jetpt%s_fake",JetIDName.Data()));
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("jetpt"));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPb->Get(Form("jetpt_real"));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPb->Get(Form("jetpt_fake"));
    }
    else{
	TH2F *hPPb2D= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[i].Data()));
	TH2F *hPPb2D_real= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s_real",JetIDName.Data(),etabinnameswap[i].Data()));
	TH2F *hPPb2D_fake= (TH2F*)fMCPPb->Get(Form("jetpt%sEtaBin%s_fake",JetIDName.Data(),etabinnameswap[i].Data()));
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_real= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_real",etabinnameswap[i].Data()));
	TH1F *hPPb_pt_fake= (TH1F*)fMCPPb->Get(Form("jetptEtaBin%s_fake",etabinnameswap[i].Data()));
    }
}

//---------------------------------------PPb JetID f/i&r/i ratio before cut ------------------------------------------ 
c1->cd(ipt+1)->SetLogz();
    if(ipt==0){
        hFrame->GetYaxis()->SetTitle(Form("%s",JetID2.Data()));
        hFrame->SetLabelSize(18,"Y");
    }
    else{
        hFrame->SetLabelSize(0,"Y");
        hFrame->GetYaxis()->SetTitle("");
    }
        hFrame->GetXaxis()->SetTitle("");
hFrame->GetXaxis()->SetNdivisions(510);
fixedFontHist(hFrame,2.0,3.0);
hFrame->SetMinimum(5e-8);
hFrame->SetMaximum(5);
hFrame->DrawCopy();
T1.DrawLatex(0.58,0.78,Form("Mean X= %.2f",hJetID2D_fake->GetMean(1)));
T1.DrawLatex(0.58,0.83,Form("Mean Y= %.2f",hJetID2D_fake->GetMean(2)));
T1.DrawLatex(0.6,0.6,"fake");
T1.DrawLatex(0.25,0.7,Form("%.f < p_{T} < %.f (GeV/c)",binbound_pt_coarse[ipt],binbound_pt_coarse[ipt+1]));
hJetID2D_fake->Draw("colz same");
c1->cd(ipt+1+4)->SetLogz();
    if(ipt==0)
    hFrame->GetYaxis()->SetTitle(Form("%s",JetID2.Data()));
hFrame->GetXaxis()->SetTitle(Form("%s",JetID1.Data()));
hFrame->DrawCopy();
hJetID2D_real->Draw("colz same");
T1.DrawLatex(0.8,0.6,"Real");
T1.DrawLatex(0.58,0.78,Form("Mean X= %.2f",hJetID2D_real->GetMean(1)));
T1.DrawLatex(0.58,0.83,Form("Mean Y= %.2f",hJetID2D_real->GetMean(2)));
if(ipt==0)
T1.DrawLatex(0.4,0.30,etastring[i]);
}
if(isMatch){
c1->Print(Form("pic/jetid12jetid14_2D_cut1.png"));
}
else{
c1->Print(Form("pic/jetid12jetid14_2D.png"));
}
}
