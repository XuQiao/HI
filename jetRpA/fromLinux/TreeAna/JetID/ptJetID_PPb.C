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
#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"

TGraphAsymmErrors* makegraph(TH1* histo1, TH1* histo2){
const int Npoint=1000;
double x[Npoint],y[Npoint],xerrlow[Npoint],xerrhigh[Npoint],yerrlow[Npoint],yerrhigh[Npoint];
TGraphAsymmErrors *g=new TGraphAsymmErrors();
g->BayesDivide(histo1,histo2);
for(int ipoint=0;ipoint<g->GetN();ipoint++){
x[ipoint]=g->GetX()[ipoint];
y[ipoint]=g->GetY()[ipoint];
yerrhigh[ipoint]=g->GetErrorYhigh(ipoint);
yerrlow[ipoint]=g->GetErrorYlow(ipoint);
xerrhigh[ipoint]=0;
xerrlow[ipoint]=0;
}
g0=new TGraphAsymmErrors(g->GetN(),x,y,xerrlow,xerrhigh,yerrlow,yerrhigh);
return g0;
}


TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());

//TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
if(!histoname.Contains("Prof")){
normalizeByBinWidth(rehisto);
}
if(binbound[Nbin]!=1000)
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
else 
rehisto->GetXaxis()->SetRangeUser(47.5,500);

return rehisto;
}

TH1* makejetid(TH2* ptjetid, int ptlow, int pthigh, int Nbin, double binbound[]){
TH1D* hjetid=ptjetid->ProjectionY("hjetid",ptjetid->GetXaxis()->FindBin(ptlow),ptjetid->GetXaxis()->FindBin(pthigh));
hjetid = (TH1D*)myRebin(hjetid, Nbin,binbound);
return hjetid;
}



void ptJetID_PPb(){
int ieta=7;
using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
//gStyle->SetLabelFont(70);
/*for(int iptbin=0;iptbin<=500;iptbin++)
binbound_pt[iptbin]=iptbin*2;
int Nbin_pt=500;*/
//------------------------------------------------------------Get Histograms---------------------------------------------
double xrange_JetIDcut[2]; double JetIDcut[2];
int ilist=0;int Nbin_JetID;
TString JetIDName = JetIDNameList[ilist];
TString JetID,Unit;
if(ilist== 0)  { JetID = "chargedMax"; JetIDcut[0]=4; JetIDcut[1]=300.;}
else if(ilist== 1)   {   JetID = "chargedSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 2)   {   JetID = "neutralMax"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 3)   {   JetID = "neutralSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 4)   {   JetID = "photonMax"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 5)   {   JetID = "photonSum"; JetIDcut[0] =0; JetIDcut[1]=200;}
else if(ilist== 6)  {    JetID = "chargedMax/p_{T}"; JetIDcut[0]=0.05; JetIDcut[1]=1;}
else if(ilist== 7)  {    JetID = "chargedSum/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.6;}
else if(ilist== 8)  {   JetID = "neutralMax/p_{T}"; JetIDcut[0] = 0;  JetIDcut[1]=0.08;}
else if(ilist== 9)  {    JetID = "neutralSum/p_{T}"; JetIDcut[0] = 0; JetIDcut[1]=0.15;}
else if(ilist== 10)  {   JetID = "photonMax/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.50;}
else if(ilist== 11) {    JetID = "photonSum/p_{T}";	JetIDcut[0] = 0; JetIDcut[1]=0.3;}
else if(ilist== 12) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}",JetIDcut[0]=0; JetIDcut[1]=1.01;}
else if(ilist== 13) {    JetID = "#frac{#Sigma p_{T}^{charged}+#Sigma p_{T}^{photons}+#Sigma p_{T}^{neutral}+#Sigma p_{T}^{#mu}+#Sigma p_{T}^{e}}{p_{T}^{jet}}"; JetIDcut[0]=0;JetIDcut[1]=1.13;}
else if(ilist== 14)   {  JetID = "Max p_{T}^{neutral} /Max(#Sigma p_{T}^{charged},#Sigma p_{T}^{neutral})" ;JetIDcut[0]=0;JetIDcut[1]=0.975;}
else if(ilist== 15)     JetID = "charged Multiplicity";
else if(ilist== 16)     JetID = "neutral Multiplicity";
else if(ilist== 17)     JetID = "photon Multiplicity";
else if(ilist== 18)   {  JetID = "#Sigma h^{#pm}+#Sigma #gamma +#Sigma h^{0}+#Sigma #mu+#Sigma e"; JetIDcut[0]=0; JetIDcut[1]=1000;}
else{   exit(0);}

if(JetIDName.Contains("pt") || JetIDName=="neuMaxr"){
if(ilist==12 || ilist==13){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2,1.3,1.4,1.6,1.8,2.};}
else if (ilist==14){
double binbound_JetID[]={0,0.2,0.4,0.6,0.7,0.8,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.1,1.15,1.2};}
else if(ilist==8 || ilist==10){
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5};}
else{
double binbound_JetID[]={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.};}

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
double binbound_JetID[]={1e-1,1,2,3,4,5,6,7,8,10,15,20,30,40,60,75,90,120,140,165,200,240,270,300};
Unit="[GeV/c]";
}
int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;

TH2F *hPPb2D= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
TH2F *hPPb2D_real= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin_real%s",JetIDName.Data(),etabinnameswap[ieta].Data()));
TH2F *hPPb2D_fake= (TH2F*)fMCPPbJetID->Get(Form("jetpt%sEtaBin_fake%s",JetIDName.Data(),etabinnameswap[ieta].Data()));

const double binbound_pt_coarse[]={0,30,70,150,200,600};
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};
xrange_JetIDcut[0]=JetIDcut[0]+1e-4;
xrange_JetIDcut[1]=JetIDcut[1]-1e-4;
TH1D* hPPb_pt=(TH1D*)hPPb2D->ProjectionX("hPPb_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_pt_real=hPPb2D_real->ProjectionX("hPPb_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_pt_fake=hPPb2D_fake->ProjectionX("hPPb_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetID[1]));
TH1D* hPPb_JetIDcut_pt=hPPb2D->ProjectionX("hPPb_JetIDcut_pt",hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_real=hPPb2D_real->ProjectionX("hPPb_JetIDcut_pt_real",hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_real->GetYaxis()->FindBin(xrange_JetIDcut[1]));
TH1D* hPPb_JetIDcut_pt_fake=hPPb2D_fake->ProjectionX("hPPb_JetIDcut_pt_fake",hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[0]),hPPb2D_fake->GetYaxis()->FindBin(xrange_JetIDcut[1]));

TProfile *hProfPPb=(TProfile*)hPPb2D->ProfileX("hProfPPb");
TProfile *hProfPPb_fake=(TProfile*)hPPb2D_fake->ProfileX("hProfPPb_fake");
TProfile *hProfPPb_real=(TProfile*)hPPb2D_real->ProfileX("hProfPPb_real");

//---------------------------------------------------Rebin and Format histograms-----------------------------------------

TH1D* rehisto_hPPb_pt=(TH1D*)myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_real=(TH1D*)myRebin(hPPb_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_pt_fake=(TH1D*)myRebin(hPPb_pt_fake,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt=(TH1D*)myRebin(hPPb_JetIDcut_pt,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt_real=(TH1D*)myRebin(hPPb_JetIDcut_pt_real,Nbin_pt,binbound_pt);
TH1D* rehisto_hPPb_JetIDcut_pt_fake=(TH1D*)myRebin(hPPb_JetIDcut_pt_fake,Nbin_pt,binbound_pt);

TProfile *rehisto_hProfPPb=(TProfile*)myRebin(hProfPPb,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_fake=(TProfile*)myRebin(hProfPPb_fake,Nbin_pt,binbound_pt);
TProfile *rehisto_hProfPPb_real=(TProfile*)myRebin(hProfPPb_real,Nbin_pt,binbound_pt);
/*cout<<"low to high: "<<"PPb real"<<"\t\t"<<"PPb fake"<<"\t\t"<<"PPb Total jets"<<"\t\t"<<"fake fraction"<<endl;
for(int i=1;i<rehisto_hProfPPb_fake->GetXaxis()->GetNbins();i++){
if(rehisto_hPPb_pt->GetBinContent(i)!=0 && rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)>110 && rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)<190){
cout<<rehisto_hPPb_pt_fake->GetXaxis()->GetBinLowEdge(i)<<" to "<<rehisto_hPPb_pt_real->GetXaxis()->GetBinUpEdge(i)<<" : ";
cout<<rehisto_hPPb_pt_real->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt_fake->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt->GetBinContent(i)<<"\t\t"<<rehisto_hPPb_pt_fake->GetBinContent(i)/rehisto_hPPb_pt->GetBinContent(i)<<endl;
}
}
*/
TH1F* hFrame=new TH1F("","",100000,0,1000);
hFrame->SetTitle("");
fixedFontHist(hFrame,1.1,1.4);
hFrame->GetYaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetTitleSize(0.06);
hFrame->GetXaxis()->SetLabelSize(0.06);
hFrame->GetYaxis()->SetLabelSize(0.06);
rehisto_hProfPPb->SetMarkerSize(1.3);
rehisto_hProfPPb_fake->SetMarkerSize(1.3);
rehisto_hProfPPb_real->SetMarkerSize(1.3);
rehisto_hProfPPb->SetMarkerStyle(24);
rehisto_hProfPPb_fake->SetMarkerStyle(20);
rehisto_hProfPPb_real->SetMarkerStyle(34);

rehisto_hProfPPb->SetMarkerColor(1);
rehisto_hProfPPb->SetLineColor(1);
rehisto_hProfPPb_fake->SetMarkerColor(2);
rehisto_hProfPPb_fake->SetLineColor(2);
rehisto_hProfPPb_real->SetMarkerColor(4);
rehisto_hProfPPb_real->SetLineColor(4);

/*
//--------------------------------------------------------pT distribution------------------------------------------------

TCanvas* c1 = new TCanvas("c1"," ",800,600);
makeMultiPanelCanvas(c1,1,1,0,0,0.08,0.1,0.03);
c1->cd(1)->SetLogy();
rehisto_hPP_pt->SetMarkerStyle(20);
rehisto_hPP_pt->Draw("E1");
rehisto_hPPb_pt->Draw("E1same");
TLegend *leg1=new TLegend(0.7,0.7,0.9,0.8);
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->SetLineWidth(0);
leg1->SetTextSize(0.05);
leg1->AddEntry(rehisto_hPP_pt,"PYTHIA","lpf");
leg1->AddEntry(rehisto_hPPb_pt,"PYTHIA+HIJING","lp");
leg1->Draw("same");
c1->Print(Form("pic/%s/jetpt_PPbvsPP.png",JetIDName.Data()));
*/
//-------------------------------------------------PP&PPb JetID distribution---------------------------------------------

TCanvas* c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c2,2,2,0,0,0.20,0.16,0.03);

if(JetIDName.Contains("pt"))
hFrame->GetXaxis()->SetRangeUser(0.01,1.22);
else
hFrame->GetXaxis()->SetRangeUser(0.9e-1,320);
hFrame->GetXaxis()->SetTitle(Form("%s %s",JetID.Data(), Unit.Data()));
const int Npt=4;
int ptbin[5]={30,70,100,150,600};
TH1D* hPPb_JetID[Npt]; TH1D* hPPb_JetID_fake[Npt]; TH1D* hPPb_JetID_real[Npt];
TH1D* ratio_hPPb_JetID_fake[Npt]; TH1D* ratio_hPPb_JetID_real[Npt];
TLatex *T2=new TLatex();
T2->SetTextAlign(12);
T2->SetTextSize(0.055);
T2->SetTextColor(1);
T2->SetTextFont(42);
T2->SetNDC();
for(int ipt=0;ipt<Npt;ipt++){
hPPb_JetID[ipt]=(TH1D*)makejetid(hPPb2D,ptbin[ipt],ptbin[ipt+1],Nbin_JetID,binbound_JetID);
hPPb_JetID_real[ipt]=(TH1D*)makejetid(hPPb2D_real,ptbin[ipt],ptbin[ipt+1],Nbin_JetID,binbound_JetID);
hPPb_JetID_fake[ipt]=(TH1D*)makejetid(hPPb2D_fake,ptbin[ipt],ptbin[ipt+1],Nbin_JetID,binbound_JetID);
hPPb_JetID[ipt]->SetMarkerStyle(24);
hPPb_JetID[ipt]->SetMarkerSize(1.5);
hPPb_JetID[ipt]->SetMarkerColor(1);
hPPb_JetID[ipt]->SetLineColor(1);
hPPb_JetID_fake[ipt]->SetMarkerStyle(20);
hPPb_JetID_fake[ipt]->SetMarkerSize(1.5);
hPPb_JetID_fake[ipt]->SetMarkerColor(2);
hPPb_JetID_fake[ipt]->SetLineColor(2);
hPPb_JetID_real[ipt]->SetMarkerStyle(34);
hPPb_JetID_real[ipt]->SetMarkerSize(1.5);
hPPb_JetID_real[ipt]->SetMarkerColor(4);
hPPb_JetID_real[ipt]->SetLineColor(4);
if(ipt==0||ipt==2)
hFrame->GetYaxis()->SetTitle("d#sigma/dp_{T}");
else
hFrame->GetYaxis()->SetTitle("");
c2->cd(ipt+1)->SetLogy();
if(JetIDName.Contains("pt")){
if(ipt<2)
hFrame->GetYaxis()->SetRangeUser(1e-7,1e-1);
else
hFrame->GetYaxis()->SetRangeUser(1e-10,1.5e-2);
}
else{
c2->cd(ipt+1)->SetLogx();
if(ipt<2)
hFrame->GetYaxis()->SetRangeUser(1.5e-10,1.7e-2);
else
hFrame->GetYaxis()->SetRangeUser(1.5e-12,1.7e-4);
}
hFrame->DrawCopy();
hPPb_JetID[ipt]->Draw("E1same");
hPPb_JetID_fake[ipt]->Draw("E1same");
hPPb_JetID_real[ipt]->Draw("E1same");
if(ipt==1)
T2->DrawLatex(0.10,0.25,Form("%d<p_{T}^{jet}<%d [GeV/c]",ptbin[ipt],ptbin[ipt+1]));
else if(ipt==3)
T2->DrawLatex(0.10,0.35,Form("%d<p_{T}^{jet}<%d [GeV/c]",ptbin[ipt],ptbin[ipt+1]));
else
T2->DrawLatex(0.23,0.40,Form("%d<p_{T}^{jet}<%d [GeV/c]",ptbin[ipt],ptbin[ipt+1]));
}
TLegend *leg2;
if(ilist==7)
leg2=new TLegend(0.68,0.32,0.88,0.45);
else
leg2=new TLegend(0.70,0.78,0.90,0.95);
leg2->SetBorderSize(0);
leg2->SetFillColor(0);
leg2->SetLineWidth(0);
leg2->SetTextSize(0.06);
leg2->AddEntry(hPPb_JetID[0],"Inclusive","lp");
leg2->AddEntry(hPPb_JetID_fake[0],"fake","lp");
leg2->AddEntry(hPPb_JetID_real[0],"real","lp");
leg2->Draw("same");
T2->DrawLatex(0.09,0.85.,"PYTHIA+HIJING");
T2->DrawLatex(0.09,0.75,"|#eta_{CM}|< 1.0");
c2->Print(Form("pic/%s/jetid_PPb.png",JetIDName.Data()));
c2->Print(Form("pic/%s/jetid_PPb.pdf",JetIDName.Data()));

//---------------------------------------- PPb JetID f/i&r/i ratio before cut ------------------------------------------ 
TCanvas* c3 = new TCanvas("c3"," ",600,600);
makeMultiPanelCanvas(c3,2,2,0,0,0.20,0.16,0.03);
TLine *l = new TLine(1e-1,1,300,1);
for(int ipt=0;ipt<4;ipt++){
if(ipt==0||ipt==2)
hFrame->GetYaxis()->SetTitle("Fraction");
else
hFrame->GetYaxis()->SetTitle("");
c3->cd(ipt+1)->SetLogy();
hFrame->GetYaxis()->SetRangeUser(5e-6,8);
if(!JetIDName.Contains("pt"))
c3->cd(ipt+1)->SetLogx();
hFrame->DrawCopy();
ratio_hPPb_JetID_fake[ipt]=(TH1D*)hPPb_JetID_fake[ipt]->Clone("ratio_hPPb_JetID_fake");
ratio_hPPb_JetID_fake[ipt]->Divide(hPPb_JetID[ipt]);
ratio_hPPb_JetID_fake[ipt]->SetTitle("");
ratio_hPPb_JetID_fake[ipt]->Draw("same");
ratio_hPPb_JetID_real[ipt]=(TH1D*)hPPb_JetID_real[ipt]->Clone("ratio_hPPb_JetID_real");
ratio_hPPb_JetID_real[ipt]->Divide(hPPb_JetID[ipt]);
ratio_hPPb_JetID_real[ipt]->Draw("same");
TLegend *leg3=new TLegend(0.45,0.55,0.65,0.75);
leg3->SetBorderSize(0);
leg3->SetFillColor(0);
leg3->SetLineWidth(0);
leg3->SetTextSize(0.065);
leg3->AddEntry(ratio_hPPb_JetID_fake[ipt],"fake/Inclusive","lp");
leg3->AddEntry(ratio_hPPb_JetID_real[ipt],"real/Inclusive","lp");
if(ipt==3)leg3->Draw("same");
T2->Draw("same");
l->Draw("same");
}
l->SetLineStyle(2);
c3->Print(Form("pic/%s/ratio_jetid_PPb.png",JetIDName.Data()));
c3->Print(Form("pic/%s/ratio_jetid_PPb.pdf",JetIDName.Data()));


//--------------------------------------- PPb pt cutvsnocut for f/r/i----------------------------------------------------------
TGraphAsymmErrors* g_fake, *g_real, *g_inc;
TCanvas* c6 = new TCanvas("c6"," ",600,600);
makeMultiPanelCanvas(c6,1,1,-0.10,0,0.12,0.10,0.03);
fixedFontHist(hFrame,1.2,1.6);
c6->cd(1);
hFrame->GetYaxis()->SetRangeUser(0.5,1.1);
hFrame->GetYaxis()->SetTitle("Cut Efficiency");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
hFrame->GetXaxis()->SetRangeUser(47.5,500);
hFrame->DrawCopy();
g_inc=makegraph(rehisto_hPPb_JetIDcut_pt,rehisto_hPPb_pt);
g_inc->SetMarkerSize(1.5);
g_inc->SetMarkerStyle(24);
g_inc->SetMarkerColor(1);
g_inc->SetLineColor(1);
g_inc->Draw("Psame");

/*TH1D* ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcutvsnocut_pt_fake");
TH1D* cl_rehisto_hPPb_pt_fake=(TH1D*)hPPb_pt_fake->Clone("cl_rehisto_hPPb_pt_fake");
cl_rehisto_hPPb_pt_fake=(TH1D*)cl_rehisto_hPPb_pt_fake->Rebin(Nbin_pt_coarse,"cl_rehisto_hPPb_pt_fake",binbound_pt_coarse);
ratio_hPPb_JetIDcutvsnocut_pt_fake=(TH1D*)ratio_hPPb_JetIDcutvsnocut_pt_fake->Rebin(Nbin_pt_coarse,"ratio_hPPb_JetIDcutvsnocut_pt_fake",binbound_pt_coarse);
normalizeByBinWidth(ratio_hPPb_JetIDcutvsnocut_pt_fake);
normalizeByBinWidth(cl_rehisto_hPPb_pt_fake);*/

g_fake=makegraph(rehisto_hPPb_JetIDcut_pt_fake,rehisto_hPPb_pt_fake);
g_fake->SetMarkerSize(1.5);
g_fake->SetMarkerStyle(20);
g_fake->SetMarkerColor(2);
g_fake->SetLineColor(2);
g_fake->Draw("Psame");
g_real=makegraph(rehisto_hPPb_JetIDcut_pt_real,rehisto_hPPb_pt_real);
g_real->SetMarkerSize(1.5);
g_real->SetMarkerStyle(34);
g_real->SetMarkerColor(4);
g_real->SetLineColor(4);
g_real->Draw("Psame");
TLegend *leg6=new TLegend(0.7,0.25,0.9,0.4);
leg6->SetBorderSize(0);
leg6->SetFillColor(0);
leg6->SetLineWidth(0);
leg6->SetTextSize(0.05);
leg6->AddEntry(g_inc,"Inclusive","lp");
leg6->AddEntry(g_fake,"fake","lp");
leg6->AddEntry(g_real,"real","lp");
leg6->Draw("same");
T2->Draw("same");
//T5->Draw("same");
TLine *l = new TLine(47.5,1,500,1);
l->Draw("same");
c6->Print(Form("pic/%s/ratio_jetpt_cutvsnocut_PPb.png",JetIDName.Data()));
c6->Print(Form("pic/%s/ratio_jetpt_cutvsnocut_PPb.pdf",JetIDName.Data()));


//---------------------------------------PPb pT f/i&r/i ratio before and after cut-----------------------------------

TCanvas *c7 = new TCanvas("c7"," ",600,600);
makeMultiPanelCanvas(c7,1,1,-0.10,0,0.12,0.10,0.03);
c7->cd(1)->SetLogy();
c7->cd(1)->SetLogx();
hFrame->GetXaxis()->SetRangeUser(47.5,499.9);
hFrame->GetXaxis()->SetMoreLogLabels();
hFrame->GetXaxis()->SetNoExponent();
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
hFrame->GetYaxis()->SetTitle("Fraction of cross section");
hFrame->GetYaxis()->SetRangeUser(1e-8,5);
hFrame->DrawCopy();

TH1D* ratio_hPPb_pt_fakevsInc=(TH1D*)rehisto_hPPb_pt_fake->Clone("ratio_hPPb_pt_fakevsInc");
ratio_hPPb_pt_fakevsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_fakevsInc->SetMarkerStyle(24);
ratio_hPPb_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_fakevsInc->SetMarkerColor(1);
ratio_hPPb_pt_fakevsInc->SetLineColor(1);
ratio_hPPb_pt_fakevsInc->Draw("same");

/*TH1D* ratio_hPPb_pt_realvsInc=(TH1D*)rehisto_hPPb_pt_real->Clone("ratio_hPPb_pt_realvsInc");
ratio_hPPb_pt_realvsInc->Divide(rehisto_hPPb_pt);
ratio_hPPb_pt_realvsInc->SetMarkerStyle(34);
ratio_hPPb_pt_realvsInc->SetMarkerSize(1.5);
ratio_hPPb_pt_realvsInc->SetMarkerColor(4);
ratio_hPPb_pt_realvsInc->SetLineColor(4);
ratio_hPPb_pt_realvsInc->Draw("same");
*/
TH1D* ratio_hPPb_JetIDcut_pt_fakevsInc=(TH1D*)rehisto_hPPb_JetIDcut_pt_fake->Clone("ratio_hPPb_JetIDcut_pt_fakevsInc");
ratio_hPPb_JetIDcut_pt_fakevsInc->Divide(rehisto_hPPb_JetIDcut_pt);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerStyle(20);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerSize(1.5);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetMarkerColor(2);
ratio_hPPb_JetIDcut_pt_fakevsInc->SetLineColor(2);
ratio_hPPb_JetIDcut_pt_fakevsInc->Draw("same");

TLegend *leg7=new TLegend(0.6,0.80,0.85,0.94);
leg7->SetBorderSize(0);
leg7->SetFillColor(0);
leg7->SetLineWidth(0);
leg7->SetTextSize(0.035);
leg7->AddEntry(ratio_hPPb_pt_fakevsInc,"fake fraction","lp");
leg7->AddEntry(ratio_hPPb_JetIDcut_pt_fakevsInc,"After cut fake fraction","lp");
leg7->Draw("same");
T2->SetTextSize(0.04);
T2->DrawLatex(0.15,0.94,"PYTHIA+HIJING");
T2->DrawLatex(0.15,0.88,"|#eta_{CM}|< 1.0");
T2->DrawLatex(0.15,0.2,Form("Cut: %.2f < %s < %.2f",JetID.Data(),JetIDcut[0],JetIDcut[1]));
c7->Print(Form("pic/%s/ratio_jetpt_%s_PPb.png",JetIDName.Data(),JetIDName.Data()));
c7->Print(Form("pic/%s/ratio_jetpt_%s_PPb.pdf",JetIDName.Data(),JetIDName.Data()));

//------------------------------- PPb pT distribution before and after cut-------------------------------------------
/*TCanvas* c8 = new TCanvas("c8"," ",1200,700);
makeMultiPanelCanvas(c8,2,2,0,0,0.15,0.1,0.03);
c8->cd(1)->SetLogy();
rehisto_hPP_pt->SetMarkerStyle(24);
rehisto_hPP_pt->SetMarkerSize(1.5);
rehisto_hPP_pt->SetMaximum(1e-2);
rehisto_hPP_pt->SetMinimum(1e-13);
rehisto_hPP_pt_real->SetMarkerStyle(34);
rehisto_hPP_pt_real->SetMarkerColor(4);
rehisto_hPP_pt_real->SetLineColor(4);
rehisto_hPP_pt_real->SetMarkerSize(1.5);
rehisto_hPP_pt_fake->SetMarkerStyle(20);
rehisto_hPP_pt_fake->SetMarkerColor(2);
rehisto_hPP_pt_fake->SetLineColor(2);
rehisto_hPP_pt_fake->SetMarkerSize(1.5);
rehisto_hPP_pt->Draw();
rehisto_hPP_pt_real->Draw("E1same");
rehisto_hPP_pt_fake->Draw("E1same");

c8->cd(2)->SetLogy();
rehisto_hPP_JetIDcut_pt->SetMarkerStyle(24);
rehisto_hPP_JetIDcut_pt->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt->SetMaximum(1e-2);
rehisto_hPP_JetIDcut_pt->SetMinimum(1e-13);
rehisto_hPP_JetIDcut_pt->GetYaxis()->SetTitle("");
rehisto_hPP_JetIDcut_pt_real->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt_real->SetMarkerStyle(34);
rehisto_hPP_JetIDcut_pt_real->SetMarkerColor(4);
rehisto_hPP_JetIDcut_pt_real->SetLineColor(4);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerStyle(20);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerColor(2);
rehisto_hPP_JetIDcut_pt_fake->SetLineColor(2);
rehisto_hPP_JetIDcut_pt_fake->SetMarkerSize(1.5);
rehisto_hPP_JetIDcut_pt->Draw();
rehisto_hPP_JetIDcut_pt_real->Draw("E1same");
rehisto_hPP_JetIDcut_pt_fake->Draw("E1same");

c8->cd(3)->SetLogy();
rehisto_hPPb_pt->SetMarkerStyle(24);
rehisto_hPPb_pt->SetMarkerSize(1.5);
rehisto_hPPb_pt->SetMaximum(1e-2);
rehisto_hPPb_pt->SetMinimum(1e-13);
rehisto_hPPb_pt->GetYaxis()->SetTitle("Probability");
rehisto_hPPb_pt_real->SetMarkerStyle(34);
rehisto_hPPb_pt_real->SetMarkerColor(4);
rehisto_hPPb_pt_real->SetLineColor(4);
rehisto_hPPb_pt_real->SetMarkerSize(1.5);
rehisto_hPPb_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_pt_fake->SetMarkerColor(2);
rehisto_hPPb_pt_fake->SetLineColor(2);
rehisto_hPPb_pt_fake->SetMarkerSize(1.5);
rehisto_hPPb_pt->Draw();
rehisto_hPPb_pt_real->Draw("E1same");
rehisto_hPPb_pt_fake->Draw("E1same");

c8->cd(4)->SetLogy();
rehisto_hPPb_JetIDcut_pt->SetMarkerStyle(24);
rehisto_hPPb_JetIDcut_pt->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt->SetMaximum(1e-2);
rehisto_hPPb_JetIDcut_pt->SetMinimum(1e-13);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerStyle(34);
rehisto_hPPb_JetIDcut_pt_real->SetMarkerColor(4);
rehisto_hPPb_JetIDcut_pt_real->SetLineColor(4);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerStyle(20);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerColor(2);
rehisto_hPPb_JetIDcut_pt_fake->SetLineColor(2);
rehisto_hPPb_JetIDcut_pt_fake->SetMarkerSize(1.5);
rehisto_hPPb_JetIDcut_pt->Draw();
rehisto_hPPb_JetIDcut_pt_real->Draw("E1same");
rehisto_hPPb_JetIDcut_pt_fake->Draw("E1same");
TLegend *leg8= new TLegend(0.7,0.75,0.9,0.9);
leg8->SetBorderSize(0);
leg8->SetFillColor(0);
leg8->SetTextSize(0.05);
leg8->SetLineWidth(0);
leg8->AddEntry(rehisto_hPP_pt,"Inclusive","lp");
leg8->AddEntry(rehisto_hPP_pt_real,"real","lp");
leg8->AddEntry(rehisto_hPP_pt_fake,"fake","lp");
leg8->Draw("same");
c8->Print(Form("pic/%s/jetpt_PPb.png",JetIDName.Data()));
*/
//--------------------------------------------------PPb JetID,pT Profile histograms-----------------------------------

TCanvas* c9 = new TCanvas("c9"," ",600,600);
makeMultiPanelCanvas(c9,1,1,-0.10,0,0.12,0.10,0.03);
c9->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
hFrame->GetYaxis()->SetTitle(Form("%s %s",JetID.Data(),Unit.Data()));
if(JetIDName.Contains("pt"))
hFrame->GetYaxis()->SetRangeUser(0.,0.5);
else
hFrame->GetYaxis()->SetRangeUser(0.,150);
fixedFontHist(hFrame,1.2,1.6);
hFrame->DrawCopy();
rehisto_hProfPPb->Draw("same");
rehisto_hProfPPb_fake->Draw("same");
rehisto_hProfPPb_real->Draw("same");
TLegend *leg9=new TLegend(0.65,0.7,0.85,0.85);
leg9->SetBorderSize(0);
leg9->SetFillColor(0);
leg9->SetLineWidth(0);
leg9->SetTextSize(0.05);
leg9->AddEntry(rehisto_hProfPPb,"Inclusive","lp");
leg9->AddEntry(rehisto_hProfPPb_fake,"fake","lp");
leg9->AddEntry(rehisto_hProfPPb_real,"real","lp");
leg9->Draw("same");
T2->DrawLatex(0.68,0.3,"PYTHIA+HIJING");
T2->DrawLatex(0.68,0.2,"|#eta_{CM}|< 1.0");
c9->Print(Form("pic/%s/Profile_PPb.png",JetIDName.Data()));
c9->Print(Form("pic/%s/Profile_PPb.pdf",JetIDName.Data()));

}
