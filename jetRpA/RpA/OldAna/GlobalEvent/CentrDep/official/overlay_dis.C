#include "../Quality/root_setting.h"

void overlay_dis(){
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
//double HFbin[]={0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,8,10,12,15,20,30,40,50,60,80,100};
double HFbin[]={0,0.24,0.47,0.67,0.86,1.05,1.24,1.43,1.62,1.83,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.76,4.06,4.38,4.71,5.05,5.42,5.79,6.18,6.58,6.99,7.41,7.84,8.3,8.76,9.24,9.73,10.24,10.77,11.31,11.9,12.49,13.12,13.79,14.5,15.26,16.08,16.97,17.96,19.11,20.44,22.12,24.34,27.91,74.73};
/*bool isQiao=kTRUE;
if(isQiao){
double HFbin[]={0,0.82,4.62,11.12,14.52,19.80,100};}
else{
double HFbin[]={0,0.92,4.72,11.21,14.60,19.88,100};}
*/
int nHFbin=sizeof(HFbin)/sizeof(double)-1;

TFile *hdata=TFile::Open("/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/output/JetTrig/CentrDep/MBCen_weight_bak_1.root");

TFile *hMC=TFile::Open("output/weight/pPbHist_Hijing_DS.root");
TString histoname[4]={"HFEnergy4","HFEnergyPlus","HFEnergyPlus4","Ntrk"};
TString histoname_[4]={"HF4","HFPlus","HFPlus4","Ntrk"};
TH1F* fdata[3]; TH1F* fMC[3]; TH1F* fdata_Cl[3];
c1 = new TCanvas("c1"," ",500,700);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.17,0.15,0.03);
for(int i=0;i<3;i++){
fdata[i]=(TH1F*)hdata->Get(Form("pfjet_%s",histoname_[i].Data()));
fMC[i]=(TH1F*)hMC->Get(Form("h%s_tr",histoname[i].Data()));
fdata[i]=(TH1F*)fdata[i]->Rebin(nHFbin,Form("Data_%s",fdata[i]->GetName()),HFbin);
fMC[i]=(TH1F*)fMC[i]->Rebin(nHFbin,Form("MC_%s",fMC[i]->GetName()),HFbin);
fixedFontHist(fdata[i],2.1,2.5);
fixedFontHist(fMC[i],2.1,2.5);
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
fdata_Cl[i]->GetYaxis()->SetRangeUser(0.8,1.59);

fMC[i]->GetXaxis()->SetTitle("");
fMC[i]->GetYaxis()->SetTitle("Event Fraction");
fdata_Cl[i]->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) #eta>4");
fdata_Cl[i]->GetYaxis()->SetTitle("Data/MC");

fMC[i]->SetMarkerSize(1.2);
fMC[i]->SetMarkerColor(2);
fMC[i]->SetLineColor(2);
fMC[i]->SetMarkerStyle(24);
}

i=2;
TLegend *leg=new TLegend(0.2,0.2,0.45,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(fdata[i],"Data","lp");
leg->AddEntry(fMC[i],"MC","lp");
c1->cd(1)->SetLogy();
//c1->cd(1)->SetLogx();
fMC[i]->Draw();
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
//->SetLogx();
fdata_Cl[i]->Draw();
TLine *l=new TLine(HFbin[0],1,HFbin[nHFbin],1);
l->SetLineStyle(2);
l->Draw("same");
c1->Print("/afs/cern.ch/user/q/qixu/CMSSW_5_3_8_HI/src/jetRpA/RpA/CentrDep/pic/CompareWt.png");
}

