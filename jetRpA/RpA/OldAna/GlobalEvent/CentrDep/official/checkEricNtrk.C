#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/TrackCorrector2D.C"
#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/Correction_bak.C"

void checkEricNtrk(){
const int N=50;
double x[N]; double y[N];	double y_[N];	double ratio[N];
//gEff->BayesDivide(h_tr,h_DS);
cout<<"trackMult"<<'\t'<<"Eric"<<'\t'<<"Qiao"<<'\t'<<"ratio Qiao/Eric"<<endl;
TrackCorrector2D corr("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Corrections/trackCorrections_HIN12017v5_XSecWeighted.root");
Corrector corr_Qiao("/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/CentrDep/official/output/noweight/pPbHijing_DS_EffCorr.root");
for(int i=0;i<N;i++){
x[i]=i;
y[i]=1./corr.getEventWeight(i);
y_[i]=1./corr_Qiao.getEventWeightNtrk(i);
ratio[i]=(double)y_[i]/y[i];
cout<<i<<"th: "<<'\t'<<y[i]<<'\t'<<y_[i]<<'\t'<<ratio[i]<<endl;
}
TGraph *gEff=new TGraph(N,x,y_);	gEff->SetName("gEff");
TGraph *gEff_Eric=new TGraph(N,x,y);	gEff_Eric->SetName("gEff_Eric");
TGraph *gRatio=new TGraph(N,x,ratio);	gRatio->SetName("gRatio");
c1 = new TCanvas("c1"," ",600,800);
makeMultiPanelCanvas(c1,1,2,0,0,0.08,0.08,0.03);
c1->cd(1);
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
gRatio->SetTitle("");
gRatio->GetXaxis()->SetTitleOffset(1.5);
gRatio->GetYaxis()->SetTitleOffset(1.5);
gRatio->GetYaxis()->SetRangeUser(0.9,1.1);
gRatio->GetXaxis()->SetTitle("Track Mult");
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
TLine *l = new TLine(0,1,100,1);
l->SetLineStyle(2);
l->Draw("same");
gEff_Eric->Draw("Psame");

TLegend *leg=new TLegend(0.55,0.2,0.75,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextSize(0.05);
leg->AddEntry(gEff,"Qiao","lp");
leg->AddEntry(gEff_Eric,"Eric","lp");
leg->Draw("same");
c1->cd(2);
gRatio->Draw("AP");
l->Draw("same");
c1->Print("compTrkEff.png");

}

