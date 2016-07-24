#include "/home/xuq/Documents/HI/RpA/Quality/root_setting.h"
#include "/home/xuq/Documents/HI/RpA/TreeAna/produceandcheck/file.h"
void Corr2D(){
gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.04);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.04);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

TString namelist[3]={"Etaphi","ptEta","ptphi"};
TString histonameID;
TH2D* hFrame=new TH2D("","",20000,-1000,1000,800,-4,4);	TCanvas *c[3];
for(int i=0;i<3;i++){
    c[i] = new TCanvas(Form("c%da",i), Form("c%d",i),0,0,600,600);
    c[i]->Range(0,0,1,1);
    c[i]->SetFillColor(0);
    c[i]->SetBorderMode(0);
    c[i]->SetBorderSize(2);
    c[i]->SetTickx(1);
    c[i]->SetTicky(1);
    c[i]->SetLeftMargin(0.13);
    c[i]->SetRightMargin(0.13);
    c[i]->SetTopMargin(0.05);
    c[i]->SetBottomMargin(0.1);
    c[i]->SetFrameFillStyle(0);
    c[i]->SetFrameBorderMode(0);
    TString name=namelist[i];
c[i]->cd();
c[i]->SetLogz(1);

if(name=="Etaphi")
histonameID="jetEtaphi";
else if(name=="ptEta")
histonameID="jetptEta";
else
histonameID="jetptphi";
TH2F* hdata2F=(TH2F*)fDataPPb->Get(histonameID);

if(name.Contains("pt")){
hFrame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV/c]");
hFrame->GetXaxis()->SetRangeUser(30,600);
hdata2F->SetMaximum(5e5);
hdata2F->SetMinimum(1);
//hFrame->SetAxisRange(0,1e-6,"Z");
}
if(name.Contains("phi")){
hFrame->GetYaxis()->SetTitle("#phi_{jet}");
hFrame->GetYaxis()->SetRangeUser(-3.2,3.2);
}
if(name.Contains("Eta"))
if(name.Contains("pt")){
hFrame->GetYaxis()->SetTitle("#eta_{lab}^{jet}");
hFrame->GetYaxis()->SetRangeUser(-3.1,3.1);
}
else{
hFrame->GetXaxis()->SetTitle("#eta_{lab}^{jet}");
hFrame->GetXaxis()->SetRangeUser(-3.2,3.2);
}
fixedFontHist(hFrame,1,1);
hFrame->DrawCopy();
hdata2F->DrawCopy("colz same");
gPad->RedrawAxis();
//c[i]->Print(Form("Corr2D_%s.png",name.Data()));
c[i]->Print(Form("Corr2D_%s.pdf",name.Data()));
}
}
