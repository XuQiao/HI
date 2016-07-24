#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void overlay_dis(){
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
//double HFbin[]={0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,8,10,12,15,20,30,40,50,60,80,100};
//double HFbin[]={0,0.24,0.47,0.67,0.86,1.05,1.24,1.43,1.62,1.83,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.76,4.06,4.38,4.71,5.05,5.42,5.79,6.18,6.58,6.99,7.41,7.84,8.3,8.76,9.24,9.73,10.24,10.77,11.31,11.9,12.49,13.12,13.79,14.5,15.26,16.08,16.97,17.96,19.11,20.44,22.12,24.34,27.91,74.73};
double HFbin[]={0,0.23,0.47,0.67,0.87,1.05,1.24,1.43,1.62,1.82,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.77,4.07,4.38,4.72,5.06,5.41,5.79,6.17,6.57,6.98,7.41,7.85,8.3,8.76,9.24,9.74,10.25,10.79,11.33,11.91,12.51,13.14,13.81,14.51,15.27,16.09,17,17.99,19.14,20.48,22.14,24.35,27.93,90.29};
/*bool isQiao=kTRUE;
if(isQiao){
double HFbin[]={0,0.82,4.62,11.12,14.52,19.80,100};}
else{
double HFbin[]={0,0.92,4.72,11.21,14.60,19.88,100};}
*/
int nHFbin=sizeof(HFbin)/sizeof(double)-1;

TFile *hdata=TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/GlobalEvent/CentrDep/MBCen_weight_RecoHI.root");

TFile *hMC=TFile::Open("pPbHijing_EffCorr_VRw.root");
TString histoname[2]={"HFp4","Mult"};
TString histoname_[2]={"HFPlus4","Ntrk"};
TH1F* fdata[2]; TH1F* fMC[2]; TH1F* fdata_Cl[2];
c1 = new TCanvas("c1"," ",500,700);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.17,0.15,0.03);
for(int i=0;i<2;i++){
fdata[i]=(TH1F*)hdata->Get(Form("pfjet_%s",histoname_[i].Data()));
fMC[i]=(TH1F*)hMC->Get(Form("h%sES",histoname[i].Data()));
if(i!=1){
fdata[i]=(TH1F*)fdata[i]->Rebin(nHFbin,Form("Data_%s",fdata[i]->GetName()),HFbin);
fMC[i]=(TH1F*)fMC[i]->Rebin(nHFbin,Form("MC_%s",fMC[i]->GetName()),HFbin);
}
normalizeByBinWidth(fdata[i]);
normalizeByBinWidth(fMC[i]);
fdata[i]->Scale(1./fdata[i]->Integral());
fMC[i]->Scale(1./fMC[i]->Integral());
fdata_Cl[i]=(TH1F*)fdata[i]->Clone(Form("ratio_%s",fdata[i]->GetName()));
fdata_Cl[i]->Divide(fMC[i]);
fdata[i]->SetMarkerSize(1.2);
fdata[i]->SetMarkerColor(1);
fdata[i]->SetLineColor(1);
fdata[i]->SetMarkerStyle(20);

fdata_Cl[i]->SetMarkerSize(1.2);
fdata_Cl[i]->SetMarkerColor(1);
fdata_Cl[i]->SetLineColor(1);
fdata_Cl[i]->SetMarkerStyle(20);

fMC[i]->SetMarkerSize(1.2);
fMC[i]->SetMarkerColor(2);
fMC[i]->SetLineColor(2);
fMC[i]->SetMarkerStyle(24);
}
i=0;
TH1F* hFrame=new TH1F("hFrame","",10000,0,100);
fixedFontHist(hFrame,2.2,2);
if(i!=1)
hFrame->GetXaxis()->SetRangeUser(0.09,50);
else
hFrame->GetXaxis()->SetRangeUser(0,50);
hFrame->GetYaxis()->SetRangeUser(5e-4,1e-1);
hFrame->GetYaxis()->SetTitle("Events Fraction");

TLegend *leg=new TLegend(0.6,0.2,0.85,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(fdata[i],"Data","lp");
leg->AddEntry(fMC[i],"MC","lp");
c1->cd(1)->SetLogy();
if(i!=1)
c1->cd(1)->SetLogx();
hFrame->GetXaxis()->SetTitle("");
hFrame->DrawCopy();
fMC[i]->Draw("same");
fdata[i]->Draw("same");
leg->Draw("same");
/*if(isQiao)
TLatex *T=new TLatex(0.4,0.9,"Use Qiao's Centrality HFbinning");
else*/
TLatex *T=new TLatex(0.4,0.9,"");
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
T->Draw("same");

c1->cd(2);
if(i!=1)
c1->cd(2)->SetLogx();
if(i!=1)
hFrame->GetXaxis()->SetTitle("HF E_{T} (GeV/c)");
else
hFrame->GetXaxis()->SetTitle("track Mult");
hFrame->GetYaxis()->SetTitle("Data/MC");
hFrame->GetYaxis()->SetRangeUser(0.5,1.59);
hFrame->DrawCopy();
fdata_Cl[i]->Draw("same");
TLine *l=new TLine(HFbin[0],1,HFbin[nHFbin],1);
l->SetLineStyle(2);
l->Draw("same");
c1->Print("CompareWt.png");
}

