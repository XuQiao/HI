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

class histogram:public TObject{
public:
TString histoname;
TString xtitle;
TString ytitle;
double xrange[2];
bool islogy;
int Nbin;
double binbound[100];
double dimension;

ClassDef(histogram,1);
} ;

#if !defined(__CINT__)
ClassImp(histogram);
#endif

TH1D* myRebin(TH1D* histo,int Nbin, double binbound[]){
//TH1D* rehisto=(TH1D*)histo->Clone("rehisto");
TH1D* rehisto=(TH1D*)histo->Rebin(Nbin,"rehisto",binbound);
for(int i=1;i<=Nbin;i++){
rehisto->SetBinContent(i,rehisto->GetBinContent(i)/(binbound[i]-binbound[i-1])*histo->GetBinWidth(1));
}
return rehisto;
}


TH1D* extract_2F(TString infilename, TString infileabb, TString histoname,double vmin, double vmax, bool proaxis){ //0: ProX, 1 : ProY
TFile *infile = new TFile(Form("~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/%s",infilename.Data()),"readonly");
TH2F *histo = (TH2F*)infile->Get(histoname.Data());
TString outfilename=histoname+"_"+infileabb;
TH1D *histo_pro; TString str;
if(proaxis){
histo_pro=histo->ProjectionY("histo_pro",histo->GetXaxis()->FindBin(vmin),histo->GetXaxis()->FindBin(vmax));
str="_ProY_";
}
else{
histo_pro=histo->ProjectionX("histo_pro",histo->GetYaxis()->FindBin(vmin),histo->GetYaxis()->FindBin(vmax));
str="_ProX_";
}
histo_pro->GetYaxis()->SetTitle("Counts");
TString outfilename_pro=histoname+str+infileabb;
histo_pro->SetTitle(Form("%.3f -- %.3f",vmin,vmax));
TCanvas *c1 = new TCanvas("c1");
histo->Draw("colz");
c1->Print(Form("pic/%s.png",outfilename.Data()));
histo_pro->Draw();
c1->Print(Form("pic/%s.png",outfilename_pro.Data()));
//delete c1;
return histo_pro;
}

TH1D* extract_1D(TString infilename, TString infileabb, TString histoname){
TFile *infile = new TFile(Form("~/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/%s",infilename.Data()),"readonly");
TH1D *histo = (TH1D*)infile->Get(histoname.Data());
histo->GetYaxis()->SetTitle("Counts");
TString outfilename=histoname+"_"+infileabb;
TCanvas *c2 = new TCanvas("c2");
histo->Draw();
c2->Print(Form("pic/%s.png",outfilename.Data()));
//delete c2;
return histo;
}

void compare(TH1D *histo1, TH1D *histo2, TString filename1, TString filename2,bool isMC, histogram histo){
TCanvas *c3=new TCanvas("c3");
TH1D* histo1_=(TH1D*)histo1->Clone("histo1_");
TH1D* histo2_=(TH1D*)histo2->Clone("histo2_");
if(histo.islogy)
c3->SetLogy();
TH1D* rehisto1=myRebin(histo1_,histo.Nbin-1,histo.binbound);
TH1D* rehisto2=myRebin(histo2_,histo.Nbin-1,histo.binbound);

//if(!histo.islogy){
//histo1_->Scale(1/histo1_->Integral(histo1_->GetXaxis()->FindBin(histo.xrange[0]),histo1_->GetXaxis()->FindBin(histo.xrange[1])));
//histo2_->Scale(1/histo2_->Integral(histo2_->GetXaxis()->FindBin(histo.xrange[0]),histo2_->GetXaxis()->FindBin(histo.xrange[1])));
rehisto1->SetFillStyle(0);
rehisto1->SetMarkerColor(1);
if(!isMC){
rehisto1->SetMarkerStyle(20);
rehisto2->SetFillStyle(3004);
rehisto2->SetMarkerStyle(0);
rehisto2->SetFillColor(2);
}
else{
rehisto1->SetMarkerStyle(29);
rehisto1->SetMarkerSize(1.5);
rehisto2->SetMarkerStyle(24);
rehisto2->SetMarkerSize(1.5);
rehisto2->SetFillStyle(0);
rehisto2->SetLineColor(0);
rehisto2->SetMarkerColor(2);
}
rehisto2->GetXaxis()->SetTitle(histo.xtitle);
rehisto2->GetXaxis()->SetTitleSize(0.04);
rehisto2->GetYaxis()->SetTitleSize(0.04);
rehisto2->GetXaxis()->SetTitleOffset(0.7);
rehisto2->GetYaxis()->SetTitleOffset(0.7);
rehisto2->GetYaxis()->SetTitle("Event Fraction");
rehisto2->GetXaxis()->SetRangeUser(histo.xrange[0],histo.xrange[1]);

if(!isMC)
rehisto2->Draw("Hist");
else
rehisto2->Draw();
rehisto1->Draw("same");
TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.9);
leg1->SetTextSize(0.03);
if(isMC){
leg1->AddEntry(rehisto1,filename1.Data(),"lp");
leg1->AddEntry(rehisto2,filename2.Data(),"lp");
}
else{
leg1->AddEntry(rehisto1,filename1.Data(),"lpf");
leg1->AddEntry(rehisto2,filename2.Data(),"lpf");
}
leg1->Draw("same");
c3->Print(Form("pic/overlay_%.1f_%s_%svs%s.png",histo.dimension,histo.histoname.Data(),filename1.Data(),filename2.Data()));
TCanvas *c4=new TCanvas("c4");
c4->SetGridy();
TH1D *ratio=(TH1D*)rehisto1->Clone("ratio");
ratio->Divide(rehisto2);
ratio->GetYaxis()->SetTitle(Form("#frac{%s}{%s}",filename1.Data(),filename2.Data()));
ratio->SetMarkerSize(0.8);
ratio->SetMarkerStyle(20);
ratio->SetMarkerColor(1);
ratio->GetXaxis()->SetTitle(histo.xtitle);
ratio->GetXaxis()->SetTitleSize(0.04);
ratio->GetYaxis()->SetTitleSize(0.04);
ratio->GetXaxis()->SetRangeUser(histo.xrange[0],histo.xrange[1]);
ratio->GetXaxis()->SetTitleOffset(0.6);
ratio->GetYaxis()->SetTitleOffset(0.6);
ratio->SetMinimum(0);
ratio->SetMaximum(2.0);
ratio->Draw();
TLine *l =new TLine(histo.xrange[0],1,histo.xrange[1],1);
l->SetLineStyle(2);
l->SetLineColor(kBlue);
l->Draw("same");
c4->Print(Form("pic/ratio_%.1f_%s_%svs%s.png",histo.dimension,histo.histoname.Data(),filename1.Data(),filename2.Data()));

//delete c3;
//delete c4;
}

void compare3(TH1D *histo1, TH1D *histo2, TH1D *histo3, TString filename1, TString filename2, TString filename3, histogram histo){
TCanvas *c5=new TCanvas("c5");
TH1D* histo1_=(TH1D*)histo1->Clone("histo1_");
TH1D* histo2_=(TH1D*)histo2->Clone("histo2_");
TH1D* histo3_=(TH1D*)histo3->Clone("histo3_");
if(histo.islogy)
c5->SetLogy();
if(!histo.islogy){
histo1_->Rebin(5);
histo2_->Rebin(5);
histo3_->Rebin(5);
}
//if(!histo.islogy){
//histo1_->Scale(1/histo1_->Integral(histo1_->GetXaxis()->FindBin(histo.xrange[0]),histo1_->GetXaxis()->FindBin(histo.xrange[1])));
//histo2_->Scale(1/histo2_->Integral(histo2_->GetXaxis()->FindBin(histo.xrange[0]),histo2_->GetXaxis()->FindBin(histo.xrange[1])));
//histo3_->Scale(1/histo3_->Integral(histo3_->GetXaxis()->FindBin(histo.xrange[0]),histo3_->GetXaxis()->FindBin(histo.xrange[1])));
histo1_->SetMarkerStyle(29);
histo1_->SetFillStyle(0);
histo1_->SetMarkerColor(1);
histo2_->SetFillStyle(0);
histo2_->SetLineColor(4);
histo2_->SetMarkerStyle(24);
histo2_->SetMarkerSize(1.5);
histo2_->SetMarkerColor(4);
histo3_->SetFillStyle(3004);
histo3_->SetMarkerStyle(0);
histo3_->SetFillColor(2);
double coef;
if(histo.islogy) coef=2;
else coef=1.2;
histo3_->GetXaxis()->SetTitle(histo.xtitle);
histo3_->GetXaxis()->SetLabelFont(0.01);
histo3_->GetXaxis()->SetTitleOffset(0.7);
histo3_->GetYaxis()->SetTitleOffset(0.7);
histo3_->GetYaxis()->SetTitle("Event Fraction");
histo3_->GetXaxis()->SetTitleSize(0.04);
histo3_->GetYaxis()->SetTitleSize(0.04);
histo3_->GetXaxis()->SetRangeUser(histo.xrange[0],histo.xrange[1]);
histo3_->SetMaximum(TMath::Max(histo1_->GetMaximum(),histo2_->GetMaximum())*coef);
histo3_->Draw("Hist");
histo2_->Draw("same");
histo1_->Draw("same");
TLegend *leg1 = new TLegend(0.6,0.7,0.9,0.9);
leg1->SetTextSize(0.03);
leg1->AddEntry(histo1_,filename1.Data(),"lpf");
leg1->AddEntry(histo2_,filename2.Data(),"lpf");
leg1->AddEntry(histo3_,filename3.Data(),"lpf");
leg1->Draw("same");
c5->Print(Form("pic/overlay_%.1f_%s_%svs%svs%s.png",histo.dimension,histo.histoname.Data(),filename1.Data(),filename2.Data(),filename3.Data()));

//delete c5;
}


void Draw(){
#include "parameter.h"

using namespace std;
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);
TH1D *histo[nfile][nhisto];
histogram histolist[nhisto];

for(int i=0;i<nhisto;i++){
histolist[i].histoname=histoname_[i];
histolist[i].xtitle=xtitle_[i];
histolist[i].xrange[0]=xrange_[i][0];
histolist[i].xrange[1]=xrange_[i][1];
histolist[i].islogy=islogy_[i];
histolist[i].dimension=dimension_[i];
histolist[i].Nbin=Nbin_[i];
for(int j=0;j<Nbin_[i];j++){
histolist[i].binbound[j]=binbound_[i][j];
}
}
for(int i=22;i<nfile;i++){
for(int j=0;j<nhisto;j++){
if(histolist[j].dimension>2){
histo[i][j]=extract_2F(infilename[i],infileabb[i],histolist[j].histoname,pro[j][0],pro[j][1],1);
}
else if(histolist[j].dimension<2 && histolist[j].dimension>1){
histo[i][j]=extract_2F(infilename[i],infileabb[i],histolist[j].histoname,pro[j][0],pro[j][1],0);
}
else
histo[i][j]=extract_1D(infilename[i],infileabb[i],histolist[j].histoname);
}
}

for(int j=0;j<noverlayhisto;j++){
for(int i=12;i<noverlayfile;i++){
compare(histo[overlayfile[i][0]][overlayhisto[j]],histo[overlayfile[i][1]][overlayhisto[j]],infileabb[overlayfile[i][0]],infileabb[overlayfile[i][1]],overlayfile[i][2],histolist[overlayhisto[j]]);
}
for(int i=0;i<noverlayfile3;i++){
compare3(histo[overlayfile3[i][0]][overlayhisto[j]],histo[overlayfile3[i][1]][overlayhisto[j]],histo[overlayfile3[i][2]][overlayhisto[j]],infileabb[overlayfile3[i][0]],infileabb[overlayfile3[i][1]],infileabb[overlayfile3[i][2]],histolist[overlayhisto[j]]);
}
}

}
