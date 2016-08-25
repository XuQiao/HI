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

bool isMC=kTRUE;



void WriteBins(TH1D *histo){
cout<<"points:"<<endl;
for(int i=1;i<histo->GetXaxis()->GetNbins()+1;i++){
cout<<i<<"th bin, from "<<histo->GetBinLowEdge(i)<<" to "<<histo->GetBinLowEdge(i)+histo->GetBinWidth(i)<<" content: "<<histo->GetBinContent(i)<<endl;
}
}

TH1D* myRebin(TH1D* histo,int Nbin, double binbound[]){

TH1D* rehisto=(TH1D*)histo->Rebin(Nbin,histo->GetName(),binbound);
for(int i=1;i<Nbin+1;i++){
rehisto->SetBinContent(i,rehisto->GetBinContent(i)/rehisto->GetBinWidth(i));
//rehisto->SetBinError(i,rehisto->GetBinError(i)/rehisto->GetBinWidth(i));
rehisto->SetBinError(i,0);

}

//rehisto->GetXaxis()->SetRangeUser(binbound[10],binbound[Nbin-1]);
fixedFontHist(rehisto,1.3,1.1);

return rehisto;
}

void DrawDiffID(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);

//------------------------------------------------------------Get Histograms---------------------------------------------


const int N=3;

//TString JetIDNameList_Simple[N]={"chMaxpt","neuMaxpt","chSumpt","neuSumpt","phoSumpt","SumSumpt"};
TString JetIDNameList_Simple[N]={"chMaxpt","neuMaxr","SumSumpt"};
/*if(algo=="ak3PF"){
double xrange_JetIDcut[N][2]={
{0+1e-4,0.5-1e-4},
{0.0+1e-4,0.15-1e-4},
{0+1e-4,0.8-1e-4},
{0.0+1e-4,0.15-1e-4},
{0.00+1e-4,0.32-1e-4},
{0.00+1e-4,0.99-1e-4}};
}
else{
double xrange_JetIDcut[N][2]={
{0+1e-4,0.3-1e-4},
{0.0+1e-4,0.08-1e-4},
{0+1e-4,0.6-1e-4},
{0.0+1e-4,0.15-1e-4},
{0.00+1e-4,0.30-1e-4},
{0.00+1e-4,1.01-1e-4}};
}
{0.00+1e-4,1.13-1e-4}};

  else if(JetIDNameList_Simple=="phoMaxpt"){
  TString JetID = "photonMax/p_{T}^{jet}";
  TString Unit = "";
  double binbound_JetID[]={0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5};
  double xrange_JetIDcut[2]={0.020+1e-4,0.5-1e-4};}
 */
double xrange_JetIDcut[N][2]={
    {0.05+1e-4,1-1e-4},
    {0+1e-4,0.975-1e-4},
    {0+1e-4,1.01-1e-4}
};

TH2F *hPP2D[N]; TH2F *hPPb2D[N];
TH2F *hPP2D_fake[N]; TH2F *hPPb2D_fake[N];

TH1D* hPPb_pt[N]; TH1D* hPP_pt[N];
TH1D* hPPb_pt_fake[N]; TH1D* hPP_pt_fake[N];
TH1D* hPPb_JetIDcut_pt[N]; TH1D* hPP_JetIDcut_pt[N];
TH1D* hPPb_JetIDcut_pt_fake[N]; TH1D* hPP_JetIDcut_pt_fake[N];

TH1D* rehisto_hPPb_pt[N]; TH1D* rehisto_hPP_pt[N];
TH1D* rehisto_hPPb_pt_fake[N]; TH1D* rehisto_hPP_pt_fake[N];
TH1D* rehisto_hPPb_JetIDcut_pt[N]; TH1D* rehisto_hPP_JetIDcut_pt[N];
TH1D* rehisto_hPPb_JetIDcut_pt_fake[N]; TH1D* rehisto_hPP_JetIDcut_pt_fake[N];


TH1D* ratio_hPP_JetIDcutvsnocut_pt[N];
TH1D* ratio_hPPb_JetIDcutvsnocut_pt[N];
TH1D* ratio_hPP_JetIDcutvsnocut_pt_fake[N];
TH1D* ratio_hPPb_JetIDcutvsnocut_pt_fake[N];
TH1D* cl_rehisto_hPP_pt_fake[N];
TH1D* cl_rehisto_hPPb_pt_fake[N];
//int marker_inc[N]={21,22,29,34,33,20};
int marker_inc[N]={20,21,20};
//int color[N]={2,8,7,4,6,1};
int color[N]={1,4,1};
//int marker_fake[N]={25,26,30,28,27,24};
int marker_fake[N]={24,25,24};

//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//const double binbound_pt[]={0,5,10,15,20,30,45,60,75,90,105,120,140,160,180,200,220,260,300,400,600,1000};
const double binbound_pt_coarse[]={0,30,70,150,200,600,1000};
//const int Nbin_JetID=sizeof(binbound_JetID)/sizeof(double)-1;
const int Nbin_pt_coarse=sizeof(binbound_pt_coarse)/sizeof(double)-1;
//const double xrange_JetID[2]={binbound_JetID[0]+1e-4,binbound_JetID[Nbin_JetID]-1e-4};

c6 = new TCanvas("c6"," ",600,600);
makeMultiPanelCanvas(c6,1,1,-0.12,0,0.13,0.1,0.03);
c6->cd(1)->SetGridx();
TH1F* hFrame=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame,1.2,1.8);
hFrame->GetXaxis()->SetLimits(28,599);
hFrame->SetTitle("");
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetYaxis()->SetTitle(Form("Cut efficiency"));
hFrame->SetMaximum(1.02);
hFrame->SetMinimum(0.92);
hFrame->DrawCopy();
TLegend *leg1=new TLegend(0.40,0.13,0.68,0.43);
TLegend *leg2=new TLegend(0.72,0.13,0.88,0.43);

for(int j=0;j<N;j++){

hPPb2D[j]= (TH2F*)fMCPPb->Get(Form("jetpt%s",JetIDNameList_Simple[j].Data()));
hPPb2D_fake[j]= (TH2F*)fMCPPb->Get(Form("jetpt%s_fake",JetIDNameList_Simple[j].Data()));

hPPb_pt[j]=hPPb2D[j]->ProjectionX(Form("hPPb_pt_%d",j));
hPPb_pt_fake[j]=hPPb2D_fake[j]->ProjectionX(Form("hPPb_pt_fake_%d",j));
hPPb_JetIDcut_pt[j]=hPPb2D[j]->ProjectionX(Form("hPPb_JetIDcut_pt_%d",j),hPPb2D[j]->GetYaxis()->FindBin(xrange_JetIDcut[j][0]),hPPb2D[j]->GetYaxis()->FindBin(xrange_JetIDcut[j][1]));
hPPb_JetIDcut_pt_fake[j]=hPPb2D_fake[j]->ProjectionX(Form("hPPb_JetIDcut_pt_fake_%d",j),hPPb2D_fake[j]->GetYaxis()->FindBin(xrange_JetIDcut[j][0]),hPPb2D_fake[j]->GetYaxis()->FindBin(xrange_JetIDcut[j][1]));

//---------------------------------------------------Rebin and Format histograms-----------------------------------------

rehisto_hPPb_pt[j]=myRebin(hPPb_pt[j],Nbin_pt,binbound_pt);
rehisto_hPPb_pt_fake[j]=myRebin(hPPb_pt_fake[j],Nbin_pt,binbound_pt);
rehisto_hPPb_JetIDcut_pt[j]=myRebin(hPPb_JetIDcut_pt[j],Nbin_pt,binbound_pt);
rehisto_hPPb_JetIDcut_pt_fake[j]=myRebin(hPPb_JetIDcut_pt_fake[j],Nbin_pt,binbound_pt);
}
//----------------------------------------------PP&PPb f&r&i cutvsnocut-------------------------------------------------
/*
c6->cd(1);
ratio_hPP_JetIDcutvsnocut_pt[j]=(TH1D*)rehisto_hPP_JetIDcut_pt[j]->Clone(Form("ratio_hPP_JetIDcut_ptcutvsnocut_pt_%d",j));
ratio_hPP_JetIDcutvsnocut_pt[j]->Divide(rehisto_hPP_pt[j]);
hFrame->GetYaxis()->SetTitle(Form("cut/no cut"));
ratio_hPP_JetIDcutvsnocut_pt[j]->SetMarkerSize(1.5);
ratio_hPP_JetIDcutvsnocut_pt[j]->Draw("same");

ratio_hPP_JetIDcutvsnocut_pt_fake[j]=(TH1D*)rehisto_hPP_JetIDcut_pt_fake[j]->Clone(Form("ratio_hPP_JetIDcutvsnocut_pt_fake_%d",j));
cl_rehisto_hPP_pt_fake[j]=(TH1D*)rehisto_hPP_pt_fake[j]->Clone(Form("cl_rehisto_hPP_pt_fake_%d",j));
cl_rehisto_hPP_pt_fake[j]=(TH1D*)cl_rehisto_hPP_pt_fake[j]->Rebin(Nbin_pt_coarse,Form("cl_rehisto_hPP_pt_fake_%d",j),binbound_pt_coarse);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]=(TH1D*)ratio_hPP_JetIDcutvsnocut_pt_fake[j]->Rebin(Nbin_pt_coarse,Form("ratio_hPP_JetIDcutvsnocut_pt_fake_%d",j),binbound_pt_coarse);
normalizeByBinWidth(ratio_hPP_JetIDcutvsnocut_pt_fake[j]);
normalizeByBinWidth(cl_rehisto_hPP_pt_fake[j]);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->Divide(cl_rehisto_hPP_pt_fake[j]);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->SetMarkerSize(1.5);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->SetMarkerStyle(marker_fake[j]);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->SetMarkerColor(color[j]);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->SetLineColor(color[j]);
ratio_hPP_JetIDcutvsnocut_pt[j]->SetMarkerStyle(marker_inc[j]);
ratio_hPP_JetIDcutvsnocut_pt[j]->SetMarkerColor(color[j]);
ratio_hPP_JetIDcutvsnocut_pt[j]->SetLineColor(color[j]);
ratio_hPP_JetIDcutvsnocut_pt_fake[j]->Draw("same");
leg1->AddEntry(ratio_hPP_JetIDcutvsnocut_pt[j],Form("%s Inc",JetIDNameList_Simple[j].Data()),"lp");
leg2->AddEntry(ratio_hPP_JetIDcutvsnocut_pt_fake[j],Form("%s fake",JetIDNameList_Simple[j].Data()),"lp");
}
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.03);
leg2->SetTextSize(0.03);
leg1->Draw("same");
TLatex *T1=new TLatex(0.5,0.92,"PYTHIA");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
TLine *l3 =new TLine(30,1,600,1);
l3->SetLineStyle(2);
l3->SetLineColor(kBlack);
l3->Draw("same");
*/
for(int j=0;j<2;j++){
ratio_hPPb_JetIDcutvsnocut_pt[j]=(TH1D*)rehisto_hPPb_JetIDcut_pt[j]->Clone(Form("ratio_hPPb_JetIDcut_ptcutvsnocut_pt_%d",j));
ratio_hPPb_JetIDcutvsnocut_pt[j]->Divide(rehisto_hPPb_pt[j]);
ratio_hPPb_JetIDcutvsnocut_pt[j]->SetTitle("");
ratio_hPPb_JetIDcutvsnocut_pt[j]->SetMarkerSize(1.5);
ratio_hPPb_JetIDcutvsnocut_pt[j]->SetMarkerStyle(marker_inc[j]);
ratio_hPPb_JetIDcutvsnocut_pt[j]->SetMarkerColor(color[j]);
ratio_hPPb_JetIDcutvsnocut_pt[j]->SetLineColor(color[j]);

ratio_hPPb_JetIDcutvsnocut_pt_fake[j]=(TH1D*)rehisto_hPPb_JetIDcut_pt_fake[j]->Clone(Form("ratio_hPPb_JetIDcutvsnocut_pt_fake_%d",j));
cl_rehisto_hPPb_pt_fake[j]=(TH1D*)rehisto_hPPb_pt_fake[j]->Clone(Form("cl_rehisto_hPPb_pt_fake_%d",j));
cl_rehisto_hPPb_pt_fake[j]=(TH1D*)cl_rehisto_hPPb_pt_fake[j]->Rebin(Nbin_pt_coarse,Form("cl_rehisto_hPPb_pt_fake_%d",j),binbound_pt_coarse);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]=(TH1D*)ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->Rebin(Nbin_pt_coarse,Form("ratio_hPPb_JetIDcutvsnocut_pt_fake_%d",j),binbound_pt_coarse);
normalizeByBinWidth(ratio_hPPb_JetIDcutvsnocut_pt_fake[j]);
normalizeByBinWidth(cl_rehisto_hPPb_pt_fake[j]);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->Divide(cl_rehisto_hPPb_pt_fake[j]);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->SetMarkerSize(1.5);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->SetMarkerStyle(marker_fake[j]);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->SetMarkerColor(color[j]);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->SetLineColor(color[j]);
ratio_hPPb_JetIDcutvsnocut_pt_fake[j]->Draw("Psame");
ratio_hPPb_JetIDcutvsnocut_pt[j]->Draw("Psame");
leg1->AddEntry(ratio_hPPb_JetIDcutvsnocut_pt[j],Form("%s Inc",JetIDNameList_Simple[j].Data()),"lp");
leg2->AddEntry(ratio_hPPb_JetIDcutvsnocut_pt_fake[j],Form("%s fake",JetIDNameList_Simple[j].Data()),"lp");
}
leg1->SetBorderSize(0);
leg2->SetBorderSize(0);
leg1->SetFillColor(0);
leg2->SetFillColor(0);
leg1->SetTextSize(0.03);
leg2->SetTextSize(0.03);
leg1->Draw("same");

TLatex *T2=new TLatex();
T2->SetNDC();
T2->SetTextAlign(12);
T2->SetTextSize(0.05);
T2->SetTextColor(1);
T2->SetTextFont(42);
T2->DrawLatex(0.4,0.92,"PYTHIA+HIJING");
TLine *l3 =new TLine(30,1,600,1);
l3->SetLineStyle(2);
l3->Draw("same");
T2->SetTextSize(0.045);
T2->DrawLatex(0.15,0.20,etastring[7]);
leg2->Draw("same");

c6->Print(Form("pic/ratio_jetpt_cutvsnocut.png",algo.Data()));


}
