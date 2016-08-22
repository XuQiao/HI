#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
void CompareOldNew(){
gStyle->SetOptStat(kFALSE);
gStyle->SetOptFit(kTRUE);
gStyle->SetErrorX(0);
static const int nColor = 8;
static const int colorCode[nColor] = {
    2, 4, 6, 7, 8, 9, 46,1
};
static const int markerCode[nColor] = {
    33, 34, 29, 24, 30, 28,27,20
};

const int N=8;
const TString etabinname[N]={"12_22","7_12","3_7","-3_3","-7_-3","-12_-7","-22_-12",""};
const double etabin[N]={1.0,0.5,0.4,0.6,0.4,0.5,1,2};
const TString etastring[N]={"-2.2<#eta_{CM}<-1.2","-1.2<#eta_{CM}<-0.7","-0.7<#eta_{CM}<-0.3","-0.3<#eta_{CM}<0.3","0.3<#eta_{CM}<0.7","0.7<#eta_{CM}<1.2","1.2<#eta_{CM}<2.2","-1.0<#eta_{CM}<1.0"};
TH2F* hPPb[N];	 TProfile* hProPPb[N]; 
TH2F* hold[N];	 TProfile* hProold[N]; 
c1 = new TCanvas("c1","c1",600,600);
//c2 = new TCanvas("c2","c2",600,600);
//c2->cd(1)->SetLogx();
makeMultiPanelCanvas(c1,1,1,-0.12,0,0.14,0.12,0.03);
//makeMultiPanelCanvas(c2,1,1,-0.12,0,0.14,0.12,0.03);
TLegend *leg1=new TLegend(0.65,0.7,0.8,0.9);
TLegend *leg2=new TLegend(0.65,0.7,0.8,0.9);
leg1->SetTextSize(0.035);
leg2->SetTextSize(0.04);
TH1F* hFrame=new TH1F("","",1000,0,1000);
hFrame->GetXaxis()->SetLimits(21.9,599.9);
hFrame->SetTitle("");
hFrame->GetYaxis()->SetRangeUser(0.98,1.1);
fixedFontHist(hFrame,1.2,1.6);
TLatex *T2=new TLatex(0.20,0.88,"CMS simulation");
T2->SetNDC();
T2->SetTextAlign(12);
T2->SetTextSize(0.04);
T2->SetTextColor(1);
T2->SetTextFont(42);
TLine *l=new TLine(30.1,1,599.9,1);
l->SetLineStyle(2);
l->SetLineColor(1);
for(int i=0;i<N;i++){
if(i==N-1){
hPPb[i]=(TH2F*)fPPb->Get("rawptJES_0-100%");
hold[i]=(TH2F*)fPPbOld->Get("rawptJES_0-100%");
}
else{
hPPb[i]=(TH2F*)fPPb->Get(Form("rawptJESEtaBin%s_Cen0-100%%",etabinname[i].Data()));
hold[i]=(TH2F*)fPPbOld->Get(Form("rawptJESEtaBin%s_Cen0-100%%",etabinname[i].Data()));
}
hPPb[i]->SetName(Form("hPPb_%d",i));
hold[i]->SetName(Form("hold_%d",i));

hProPPb[i]=(TProfile*)hPPb[i]->ProfileX();
hProold[i]=(TProfile*)hold[i]->ProfileX();
hProPPb[i]=(TProfile*)hProPPb[i]->Rebin(Nbin_pt,hProPPb[i]->GetName(),binbound_pt);
hProold[i]=(TProfile*)hProold[i]->Rebin(Nbin_pt,hProold[i]->GetName(),binbound_pt);
}

c1->cd(1);
hFrame->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
hFrame->GetYaxis()->SetTitle("p_{T}^{jet}/p_{T}^{raw}");
hFrame->DrawCopy();
for(int i=0;i<N;i++){
if(i==2){
TLatex *T1=new TLatex(0.25,0.83,etastring[i]);
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);

hProPPb[i]->SetMarkerStyle(markerCode[i]);
hProold[i]->SetMarkerStyle(markerCode[(i<6)?(6-i):3]);
hProPPb[i]->SetMarkerColor(colorCode[i]);
hProold[i]->SetMarkerColor(colorCode[i]);
hProPPb[i]->SetLineColor(colorCode[i]);
hProold[i]->SetLineColor(colorCode[i]);
hProPPb[i]->SetMarkerSize(1.4);
hProold[i]->SetMarkerSize(1.4);
leg1->AddEntry(hProPPb[i],Form("New %s",etastring[i].Data()),"lp");
leg1->AddEntry(hProold[i],Form("Old %s",etastring[i].Data()),"lp");
hProPPb[i]->DrawCopy("same");
hProold[i]->DrawCopy("same");
}
}
T2->Draw("same");
leg1->SetBorderSize(0);
leg1->SetFillColor(0);
leg1->Draw("same");
c1->Print(Form("OldNewrawptJES_etabin.pdf"));

}
