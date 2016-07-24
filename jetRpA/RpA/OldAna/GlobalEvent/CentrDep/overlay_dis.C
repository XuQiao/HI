#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"

void overlay_dis(){
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
double HFbin[]={0,0.1,0.2,0.3,0.5,0.7,1,1.5,2,2.5,3,4,5,6,7,8,10,12,15,20,30};
/*bool isQiao=kTRUE;
if(isQiao){
double HFbin[]={0,0.82,4.62,11.12,14.52,19.80,100};}
else{
double HFbin[]={0,0.92,4.72,11.21,14.60,19.88,100};}
*/
int nHFbin=sizeof(HFbin)/sizeof(double)-1;


TFile *f_E=TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/CentrDep/MBCen_wGplus_wNoise.root");
TFile *f_Q=TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/output/JetTrig/CentrDep/MBCen_Qiaoweight.root");
TString histoname[3]={"HF4","HFPlus","HFPlus4"};
TH1F* h_E[3]; TH1F* h_Q[3]; TH1F* h_E_Cl[3];
c1 = new TCanvas("c1"," ",500,700);
makeMultiPanelCanvas(c1,1,2,-0.05,0,0.17,0.15,0.03);
for(int i=0;i<3;i++){
h_E[i]=(TH1F*)f_E->Get(Form("pfjet_%sWt",histoname[i].Data()));
h_Q[i]=(TH1F*)f_Q->Get(Form("pfjet_%sWt",histoname[i].Data()));
h_E[i]=(TH1F*)h_E[i]->Rebin(nHFbin,Form("Eric_%s",h_E[i]->GetName()),HFbin);
h_Q[i]=(TH1F*)h_Q[i]->Rebin(nHFbin,Form("Qiao_%s",h_Q[i]->GetName()),HFbin);
fixedFontHist(h_E[i],2.1,2.5);
fixedFontHist(h_Q[i],2.1,2.5);
normalizeByBinWidth(h_E[i]);
normalizeByBinWidth(h_Q[i]);
h_E_Cl[i]=(TH1F*)h_E[i]->Clone(Form("ratio_%s",h_E[i]->GetName()));
h_E_Cl[i]->Divide(h_Q[i]);
h_E[i]->SetMarkerSize(1.2);
h_E[i]->SetMarkerColor(1);
h_E[i]->SetLineColor(1);
h_E[i]->SetMarkerStyle(20);

h_E_Cl[i]->SetMarkerSize(1.2);
h_E_Cl[i]->SetMarkerColor(1);
h_E_Cl[i]->SetLineColor(1);
h_E_Cl[i]->SetMarkerStyle(20);
h_E_Cl[i]->GetYaxis()->SetRangeUser(0.8,1.09);

h_Q[i]->GetXaxis()->SetTitle("");
h_Q[i]->GetYaxis()->SetTitle("# of Events");
h_E_Cl[i]->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) #eta>4");
h_E_Cl[i]->GetYaxis()->SetTitle("Eric/Qiao");

h_Q[i]->SetMarkerSize(1.2);
h_Q[i]->SetMarkerColor(2);
h_Q[i]->SetLineColor(2);
h_Q[i]->SetMarkerStyle(24);
}

i=2;
TH1F* h_tr=(TH1F*)f_Q->Get(Form("pfjet_%s",histoname[i].Data()));
h_tr=(TH1F*)h_tr->Rebin(nHFbin,Form("After_sele"),HFbin);
normalizeByBinWidth(h_tr);
fixedFontHist(h_tr);
h_tr->SetMarkerSize(1.2);
h_tr->SetMarkerColor(4);
h_tr->SetLineColor(4);
h_tr->SetMarkerStyle(30);
TLegend *leg=new TLegend(0.2,0.2,0.45,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(h_Q[i],"Weighted by Qiao","lp");
leg->AddEntry(h_E[i],"Weighted by Eric","lp");
leg->AddEntry(h_tr,"Event Selected","lp");
c1->cd(1)->SetLogy();
c1->cd(1)->SetLogx();
h_Q[i]->Draw();
h_E[i]->Draw("same");
h_tr->Draw("same");
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

c1->cd(2)->SetLogx();
h_E_Cl[i]->Draw();

TLine *l=new TLine(HFbin[0],1,HFbin[nHFbin],1);
l->SetLineStyle(2);
l->Draw("same");
/*if(isQiao)
c1->Print("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/pic/CompareWt_Qiao.pdf");
else*/
c1->Print("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/pic/CompareWt.pdf");
}

