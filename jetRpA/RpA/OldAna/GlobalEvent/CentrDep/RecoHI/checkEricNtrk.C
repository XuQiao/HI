#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/TrackCorrector2D.C"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/GlobalEvent/CentrDep/RecoHI/Correction.C"

void checkEricNtrk(){
gStyle->SetOptStat(kFALSE);
const int N=50;
double x[N]; double y[N];	double y1[N],y2[N];	double ratio[N];
//gEff->BayesDivide(h_tr,h_DS);
cout<<"trackMult"<<'\t'<<"Eric"<<'\t'<<"Qiao"<<'\t'<<"ratio Qiao/Eric"<<endl;
TrackCorrector2D corr("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/trackCorrections_HIN12017v5_XSecWeighted.root");
Corrector corr_QiaoReco("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/GlobalEvent/CentrDep/RecoHI/NoVRw/pPbHijing_EffCorr.root");
Corrector corr_QiaoReco_NoGplus("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/GlobalEvent/CentrDep/RecoHI/NoVRw/pPbHijing_NoGplus_EffCorr.root");
TFile *fEff=TFile::Open("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/GlobalEvent/CentrDep/RecoHI/NoVRw/pPbHijing_NoGplus_EffCorr.root");
TH1D* hEff=(TH1D*)fEff->Get("hEffMult");
for(int i=0;i<N;i++){
x[i]=i;
y1[i]=1./corr.getEventWeight(i);
//y1[i]=1./corr.getEventWeightNtrk(i);
y2[i]=1./corr_QiaoReco_NoGplus.getEventWeightNtrk(i);
ratio[i]=(double)y2[i]/y1[i];
cout<<i<<"th: "<<'\t'<<y1[i]<<'\t'<<y2[i]<<'\t'<<ratio[i]<<endl;
}
//TGraph *gEff=new TGraph(N,x,y2);	gEff->SetName("gEff");
//TGraph *gEff_Eric=new TGraph(N,x,y1);	gEff_Eric->SetName("gEff_Eric");
//TGraph *gRatio=new TGraph(N,x,ratio);	gRatio->SetName("gRatio");
TH1D* hEff_Eric=(TH1D*)hEff->Clone("hEff_Eric");
TH1D* hRatio=(TH1D*)hEff->Clone("hRatio");
for(int i=0;i<N;i++){
hEff_Eric->SetBinContent(i,y1[i]);
hEff_Eric->SetBinError(i,0);
hEff->SetBinContent(i,y2[i]);
hRatio->SetBinContent(i,ratio[i]);
}
c1 = new TCanvas("c1"," ",600,800);
makeMultiPanelCanvas(c1,1,2,0,0,0.1,0.12,0.03);
c1->cd(1);
/*
gEff->SetMarkerStyle(20);
gEff->SetMarkerColor(1);
gEff->SetMarkerSize(1.2);
gEff_Eric->SetMarkerStyle(24);
gEff_Eric->SetMarkerColor(2);
gEff_Eric->SetMarkerSize(1.2);
gRatio->SetMarkerStyle(20);
gRatio->SetMarkerColor(1);
gRatio->SetMarkerSize(1.2);
gEff->GetYaxis()->SetTitle("Efficiency");
gEff->GetXaxis()->SetTitleOffset(1.5);
gEff->GetYaxis()->SetTitleOffset(1.5);
gEff->SetTitle("");
gEff->GetYaxis()->SetRangeUser(-0.05,1.1);
gEff->GetXaxis()->SetRangeUser(0,50);
gRatio->SetTitle("");
gRatio->GetXaxis()->SetTitleOffset(1.5);
gRatio->GetYaxis()->SetTitleOffset(1.5);
gRatio->GetYaxis()->SetRangeUser(0.9,1.1);
gRatio->GetXaxis()->SetTitle("Track Mult");
gRatio->GetXaxis()->SetRangeUser(0,50);
gRatio->GetXaxis()->CenterTitle();
gRatio->GetXaxis()->SetTitleSize(18);
gRatio->GetXaxis()->SetTitleFont(44);
gRatio->GetYaxis()->SetTitle("Qiao/Eric");
gRatio->GetYaxis()->CenterTitle();
gRatio->GetYaxis()->SetTitleSize(18);
gRatio->GetYaxis()->SetTitleFont(44);
gEff->GetYaxis()->CenterTitle();
gEff->GetYaxis()->SetTitleSize(18);
gEff->GetYaxis()->SetTitleFont(44);
gEff->Draw("AP");
gEff_Eric->Draw("Psame");
*/
hEff->SetMarkerStyle(20);
hEff->SetMarkerColor(1);
hEff->SetMarkerSize(1.2);
hEff_Eric->SetMarkerStyle(24);
hEff_Eric->SetMarkerColor(2);
hEff_Eric->SetMarkerSize(1.2);
hRatio->SetMarkerStyle(20);
hRatio->SetMarkerColor(1);
hRatio->SetMarkerSize(1.2);
hEff->GetYaxis()->SetTitle("Efficiency");
hEff->GetXaxis()->SetTitleOffset(2);
hEff->GetYaxis()->SetTitleOffset(1.8);
hEff->SetTitle("");
hEff->GetYaxis()->SetRangeUser(-0.05,1.1);
hEff->GetXaxis()->SetRangeUser(0,50);
hRatio->SetTitle("");
hRatio->GetXaxis()->SetTitleOffset(2);
hRatio->GetYaxis()->SetTitleOffset(1.8);
hRatio->GetYaxis()->SetRangeUser(0.9,1.1);
hRatio->GetXaxis()->SetTitle("Track Mult");
hRatio->GetXaxis()->SetRangeUser(0,50);
hRatio->GetXaxis()->CenterTitle();
hRatio->GetXaxis()->SetTitleSize(18);
hRatio->GetXaxis()->SetTitleFont(44);
hRatio->GetYaxis()->SetTitle("Using trackMult/AN-12-377");
hRatio->GetYaxis()->CenterTitle();
hRatio->GetYaxis()->SetTitleSize(18);
hRatio->GetYaxis()->SetTitleFont(44);
hEff->GetYaxis()->CenterTitle();
hEff->GetYaxis()->SetTitleSize(18);
hEff->GetYaxis()->SetTitleFont(44);
hEff->Draw();
hEff_Eric->Draw("same");
TLine *l = new TLine(0,1,50,1);
l->SetLineStyle(2);
l->Draw("same");

TLegend *leg=new TLegend(0.50,0.2,0.75,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
//leg->AddEntry(gEff,"Qiao","lp");
//leg->AddEntry(gEff_Eric,"Eric","lp");
leg->AddEntry(hEff,"Efficiency Using trackMult","lp");
leg->AddEntry(hEff_Eric,"Effciency in AN-12-377","lp");
leg->Draw("same");
c1->cd(2);
//gRatio->Draw("AP");
hRatio->GetYaxis()->SetRangeUser(0.9,1.09);
for(int ibin=0;ibin<hRatio->GetNbinsX();ibin++){
hRatio->SetBinContent(ibin,hRatio->GetBinContent(ibin));
hRatio->SetBinError(ibin,1e-8);
}
hRatio->Draw();
l->Draw("same");
c1->Print("compTrkEff.png");

}

