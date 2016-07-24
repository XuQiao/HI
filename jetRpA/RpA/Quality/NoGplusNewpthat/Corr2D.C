#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
void Corr2D(){
gStyle->SetOptStat(kFALSE);
   TCanvas *c1 = new TCanvas("c1a", "c1",0,0,600,600);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.13);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);
    c1->SetFrameFillStyle(0);
    c1->SetFrameBorderMode(0);

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.04);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.04);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
TCanvas *c2=(TCanvas*)c1->Clone("c2");c2->SetName("c2");
TCanvas *c3=(TCanvas*)c1->Clone("c3");c2->SetName("c3");
TString namelist[3]={"Etaphi","ptEta","ptphi"};
TString histonameID;
TH2D* hFrame=new TH2D("","",20000,-1000,1000,800,-4,4);
for(int i=0;i<3;i++){
TString name=namelist[i];
if(i==0)
c1->cd();
else if(i==1)
c2->cd();
else
c3->cd();

if(name=="Etaphi")
histonameID="jetEtaphi_0-100%";
else if(name=="ptEta")
histonameID="jetptEta_0-100%";
else
histonameID="jetptphi_0-100%";
TH2F* hdata2F=(TH2F*)fdata->Get(histonameID);

if(name.Contains("pt")){
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
hFrame->GetXaxis()->SetRangeUser(50,600);
hdata2F->SetMaximum(1e-5);
hdata2F->SetMinimum(1e-11);
//hFrame->SetAxisRange(0,1e-6,"Z");
}
if(name.Contains("phi")){
hFrame->GetYaxis()->SetTitle("jet #phi");
hFrame->GetYaxis()->SetRangeUser(-3.2,3.2);
}
if(name.Contains("Eta"))
if(name.Contains("pt")){
hFrame->GetYaxis()->SetTitle("jet #eta");
hFrame->GetYaxis()->SetRangeUser(-3.1,3.1);
}
else{
hFrame->GetXaxis()->SetTitle("jet #eta");
hFrame->GetXaxis()->SetRangeUser(-3.2,3.2);
hdata2F->SetMaximum(1e-6);
hdata2F->SetMinimum(1e-9);
}
fixedFontHist(hFrame,1,1);
hFrame->DrawCopy();
c1->SetLogz();
hdata2F->DrawCopy("colz same");
gPad->RedrawAxis();
c1->Print(Form("pic/DataMC/Corr2D_%s.png",name.Data()));
c1->Print(Form("pic/DataMC/Corr2D_%s.pdf",name.Data()));
}
}
