#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawVz(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
const int N=5;
static const int Color[N]={
   1, 2, 4, 46,6, 7,8
};
static const int Style[N] = {
    20, 34, 33, 25,27, 28//,24
};

const double pPbLumi = 20.7 ;

TString TrigName[N]={"Jet20","Jet40","Jet60","Jet80","Jet100"};
TFile *f=TFile::Open("/scratch/xuq7/RpA/GlobalEvent/EventSel");

const int NVz=2;	int Ncut[NVz]={2,3,4,7};	TH1F* hVz[N][NVz];	TH1F* hRatio1[N];
int k0=1;
TH1F* hFrame=new TH1F("","",1000,-50,50);
hFrame->GetXaxis()->SetRangeUser(-15,15);
hFrame->GetYaxis()->SetTitle("Event Fraction");
fixedFontHist(hFrame,1.4,1.6);
TLegend *leg=new TLegend(0.60,0.75,0.8,0.95);
TLegend *leg_=new TLegend(0.60,0.75,0.8,0.95);
//TLine *l=new TLine(-30,1,30,1);
TLine *l=new TLine(-15,1,15,1);
l->SetLineStyle(2);
l->SetLineColor(1);

TLine *l1=new TLine(-15,1.01e-6,-15,2e-1);
TLine *l3=new TLine(-15,0,-15,2);
TLine *l2=new TLine(15,1.01e-6,15,2e-1);
TLine *l4=new TLine(15,0,15,2);
l1->SetLineStyle(2);
l2->SetLineStyle(2);
l3->SetLineStyle(2);
l4->SetLineStyle(2);
l1->SetLineColor(2);
l2->SetLineColor(2);
l3->SetLineColor(2);
l4->SetLineColor(2);
l1->SetLineWidth(4);
l2->SetLineWidth(4);
l3->SetLineWidth(4);
l4->SetLineWidth(4);

leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextFont(42);
leg->SetTextSize(0.035);
/*if(k0==0){
leg->SetHeader("Before Event Selection");
leg_->SetHeader("Before Event Selection");
}
else{
leg->SetHeader("After Event Selection");
leg_->SetHeader("After Event Selection");
}*/
//leg->SetTextAlign(32);
leg_->SetBorderSize(0);
leg_->SetFillColor(0);
leg_->SetTextFont(42);
leg_->SetTextSize(0.035);
for(int i=N-1;i<N;i++){
for(int k=0;k<1;k++){
hVz[i][k]=(TH1F*)f->Get(Form("%svz",TrigName[i].Data()));
hVz[i][k]->Rebin(10);
normalizeByBinWidth(hVz[i][k]);
hVz[i][k]->Scale((double)1.0/hVz[i][k]->Integral());
}
}
c1 = new TCanvas("c1"," ",600,600);
c2 = new TCanvas("c2"," ",600,600);
makeMultiPanelCanvas(c1,1,1,-0.05,0,0.12,0.12,0.03);
makeMultiPanelCanvas(c2,1,1,-0.05,0,0.12,0.12,0.03);

c1->cd(1);//->SetLogy();
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetRangeUser(0,1.2e-1);
hFrame->DrawCopy();
k0=0;
for(i=0;i<N;i++){
hVz[i][k0]->SetMarkerColor(Color[i]);
hVz[i][k0]->SetLineColor(Color[i]);
hVz[i][k0]->SetMarkerStyle(Style[i]);
hVz[i][k0]->SetMarkerSize(1.2);
leg->AddEntry(hVz[i][k0],TrigName[i],"lp");
hVz[i][k0]->Draw("same");
}
leg->Draw("same");
drawCMS(0.2,0.9,pPbLumi);
c1->Print("Vz_zoomDis.pdf");
//l1->Draw("same");
//l2->Draw("same");
c2->cd(1);

hFrame->GetXaxis()->SetTitle("Vz (cm)");
//hFrame->GetYaxis()->SetRangeUser(0.9,1.25);
hFrame->GetYaxis()->SetRangeUser(0.8,1.2);
hFrame->GetYaxis()->SetTitle("Ratio");
hFrame->DrawCopy();
int Selbase=0;
for(i=0;i<N;i++){
if(i!=Selbase){
hRatio1[i]=(TH1F*)hVz[i][k0]->Clone(Form("hRatio1_%d",i));
hRatio1[i]->SetMarkerColor(Color[i]);
hRatio1[i]->SetLineColor(Color[i]);
hRatio1[i]->SetMarkerStyle(Style[i]);
hRatio1[i]->SetMarkerSize(1.2);
leg_->AddEntry(hRatio1[i],Form("%s/%s",TrigName[i].Data(),TrigName[Selbase].Data()),"lp");
hRatio1[i]->Divide(hVz[Selbase][k0]);
hRatio1[i]->Draw("same");
}
}
//l3->Draw("same");
//l4->Draw("same");
leg_->Draw("same");
l->Draw("same");
drawCMS(0.2,0.9,pPbLumi);
c2->Print("AftEvSelVz_zoomRatio.pdf");

/*
int i0=5, k1=0, k2=1;
c2 = new TCanvas("c2"," ",600,800);
makeMultiPanelCanvas(c2,1,2,-0.05,0,0.12,0.12,0.03);
c2->cd(1)->SetLogy();
//hFrame->GetYaxis()->SetRangeUser(1e-5,2e-1);
hFrame->GetYaxis()->SetRangeUser(0.5,1e7);
hFrame->GetXaxis()->SetTitle("");
hFrame->GetYaxis()->SetTitle("Number of Event");
hVz[i0][k1]->SetMarkerColor(1);
hVz[i0][k1]->SetLineColor(1);
hVz[i0][k1]->SetMarkerStyle(20);
hVz[i0][k1]->SetMarkerSize(1.2);
hVz[i0][k2]->SetMarkerColor(2);
hVz[i0][k2]->SetLineColor(2);
hVz[i0][k2]->SetMarkerStyle(24);
hVz[i0][k2]->SetMarkerSize(1.2);
leg->Clear();
leg->SetHeader(TrigName[i0]);
leg->AddEntry(hVz[i0][k1],"Before ES","lp");
leg->AddEntry(hVz[i0][k2],"After ES","lp");
hFrame->DrawCopy();
hVz[i0][k1]->Draw("same");
hVz[i0][k2]->Draw("same");
leg->Draw("same");
c2->cd(2);
hFrame->GetYaxis()->SetRangeUser(0,1.99);
hFrame->GetXaxis()->SetTitle("Vz (cm)");
hFrame->GetYaxis()->SetRangeUser(0.1,1.89);
hFrame->GetYaxis()->SetTitle(Form("After/Before"));
hFrame->DrawCopy();
TH1F* hRatio2=(TH1F*)hVz[i0][k2]->Clone("hRatio2");
hRatio2->Divide(hVz[i0][k1]);
hRatio2->SetMarkerColor(1);
hRatio2->SetLineColor(1);
hRatio2->SetMarkerStyle(20);
hRatio2->Draw("same");
l->Draw("same");
c2->Print(Form("%sEvSel.png",TrigName[i0].Data()));
*/
}
