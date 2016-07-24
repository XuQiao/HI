#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void effcal(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);


TFile *fileDS = new TFile("../../AnaEff/test/output/NoVRw/DS.root","ReadOnly");
//TFile *fileDS = new TFile("../../AnaEff/test/output/VRw/DS_v2.root","ReadOnly");
TFile *fileES = new TFile("../../ESAnaEff/test/output/NoVRw/ES.root","ReadOnly");
//TFile *fileES = new TFile("../../ESAnaEff/test/output/VRw/ES_v1.root","ReadOnly");

TH1D *hMultDS = (TH1D*) fileDS->Get("demo/hMultDS");	hMultDS->SetName("hMultDS");
TH1D *hMultES = (TH1D*) fileES->Get("demo/hMultES");	hMultES->SetName("hMultES");
TH1D *hMultESDS = (TH1D*) fileES->Get("demo/hMultDS");	hMultESDS->SetName("hMultESDS");
TH1D *hMultESnonDS = (TH1D*) fileES->Get("demo/hMultnonDS");	hMultESnonDS->SetName("hMultESnonDS");

TH1D* hEffDS = (TH1D*)hMultESDS->Clone("hEffDS");	hEffDS->Divide(hMultDS);
TH1D* hContnonDS = (TH1D*)hMultESnonDS->Clone("hContnonDS");	hContnonDS->Divide(hMultES);

TH1D *hHFp4DS = (TH1D*) fileDS->Get("demo/hHFp4DS");    hHFp4DS->SetName("hHFp4DS");
TH1D *hHFp4ES = (TH1D*) fileES->Get("demo/hHFp4ES");    hHFp4ES->SetName("hHFp4ES");

TH1D *hHF4DS = (TH1D*) fileDS->Get("demo/hHF4DS");    hHF4DS->SetName("hHF4DS");
TH1D *hHF4ES = (TH1D*) fileES->Get("demo/hHF4ES");    hHF4ES->SetName("hHF4ES");

/*
const int nHFbin=100;double HFbin[nHFbin+1];
for(int i=0;i<=nHFbin;i++)
HFbin[i]=i;
*/
double HFbin[]={0,0.23,0.47,0.67,0.87,1.05,1.24,1.43,1.62,1.82,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.77,4.07,4.38,4.72,5.06,5.41,5.79,6.17,6.57,6.98,7.41,7.85,8.3,8.76,9.24,9.74,10.25,10.79,11.33,11.91,12.51,13.14,13.81,14.51,15.27,16.09,17,17.99,19.14,20.48,22.14,24.35,27.93,90.29};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;

//TH1D* rehHFp4DS = (TH1D*) hHFp4DS->Rebin(nHFbin,"rehHFp4DS",HFbin);	normalizeByBinWidth(rehHFp4DS);
//TH1D* rehHFp4ES = (TH1D*) hHFp4ES->Rebin(nHFbin,"rehHFp4ES",HFbin);	normalizeByBinWidth(rehHFp4ES);

double hf4bin[]={0,1,2,3,4,6,8,10,13,16,20,25,30,40,55,70,90};
int nhf4bin = 16;
TH1D* rehHFp4DS = (TH1D*) hHFp4DS->Rebin(nhf4bin,"rehHFp4DS",hf4bin);	normalizeByBinWidth(rehHFp4DS);
TH1D* rehHFp4ES = (TH1D*) hHFp4ES->Rebin(nhf4bin,"rehHFp4ES",hf4bin);	normalizeByBinWidth(rehHFp4ES);

TH1D* rehHF4DS = (TH1D*) hHF4DS->Rebin(nhf4bin,"rehHF4DS",hf4bin);     normalizeByBinWidth(rehHF4DS);
TH1D* rehHF4ES = (TH1D*) hHF4ES->Rebin(nhf4bin,"rehHF4ES",hf4bin);     normalizeByBinWidth(rehHF4ES);

cout<<"Zero Multiplicity Fraction= "<<hMultDS->Integral(1,1)/hMultDS->Integral()<<endl;

for(int i=1;i<hMultDS->GetNbinsX();i++){
if(hMultDS->GetBinContent(i)<hMultES->GetBinContent(i)){
//cout<<i<<"th: "<<hMultDS->GetXaxis()->GetBinLowEdge(i)<<"to"<<hMultES->GetXaxis()->GetBinUpEdge(i)<<" "<<hMultDS->GetBinContent(i)<<'\t'<<hMultES->GetBinContent(i)<<endl;
hMultDS->SetBinContent(i,hMultES->GetBinContent(i));
}
}

for(int i=1;i<hHFp4DS->GetNbinsX();i++){
if(hHFp4DS->GetBinContent(i)<hHFp4ES->GetBinContent(i)){
//cout<<i<<"th: "<<hHFp4DS->GetXaxis()->GetBinLowEdge(i)<<"to"<<hHFp4ES->GetXaxis()->GetBinUpEdge(i)<<" "<<hHFp4DS->GetBinContent(i)<<'\t'<<hHFp4ES->GetBinContent(i)<<endl;
hHFp4DS->SetBinContent(i,hHFp4ES->GetBinContent(i));
}
if(rehHFp4DS->GetBinContent(i)<rehHFp4ES->GetBinContent(i)){
//cout<<i<<"th: "<<rehHFp4DS->GetXaxis()->GetBinLowEdge(i)<<"to"<<rehHFp4ES->GetXaxis()->GetBinUpEdge(i)<<" "<<rehHFp4DS->GetBinContent(i)<<'\t'<<rehHFp4ES->GetBinContent(i)<<endl;
rehHFp4DS->SetBinContent(i,hHFp4ES->GetBinContent(i));
}
if(rehHF4DS->GetBinContent(i)<rehHF4ES->GetBinContent(i)){
//cout<<i<<"th: "<<rehHFp4DS->GetXaxis()->GetBinLowEdge(i)<<"to"<<rehHFp4ES->GetXaxis()->GetBinUpEdge(i)<<" "<<rehHFp4DS->GetBinContent(i)<<'\t'<<rehHFp4ES->GetBinContent(i)<<endl;
rehHF4DS->SetBinContent(i,hHF4ES->GetBinContent(i));
}

}

//hMultDS->Draw();
//hMultES->SetLineColor(2);
//hMultES->Draw("same");

TH1D* hFrame=new TH1D("","",1000,0,100);
hFrame->GetXaxis()->CenterTitle();
hFrame->GetYaxis()->CenterTitle();
hFrame->GetYaxis()->SetRangeUser(1e-3,1.15);
hFrame->GetXaxis()->SetRangeUser(1e-1,30);
TLegend *leg=new TLegend(0.3,0.25,0.8,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextFont(42);
leg->SetTextSize(0.04);
leg->SetHeader("Efficiency= Selected/Double Sided");
TLatex *T1=new TLatex(0.23,0.85,"CMS Simulation pPb #sqrt{s}=5.02 TeV");
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.04);
T1->SetTextColor(1);
T1->SetTextFont(42);
  TCanvas *c1 = new TCanvas("c1","Efficiency Fraction",500,500);
  TGraphAsymmErrors *gEffMult = new TGraphAsymmErrors(); 
  gEffMult->SetName("gEffMult");
  gEffMult->BayesDivide(hMultES, hMultDS);
for(int i=0;i<gEffMult->GetN();i++){
gEffMult->SetPointEXhigh(i,0);
gEffMult->SetPointEXlow(i,0);
}
  gEffMult->SetMarkerStyle(20);
  gEffMult->SetLineStyle(1);
  gEffMult->SetLineColor(1);
  gEffMult->SetMarkerColor(1);
  TLine *l=new TLine(0,1,30,1);
  l->SetLineStyle(2);
  hFrame->GetXaxis()->SetTitle("track Multiplicity");
  hFrame->GetYaxis()->SetTitle("Double Sided Efficiency");
  hFrame->DrawCopy();
	leg->AddEntry(gEffMult,"Hijing","lp");
  gEffMult->Draw("Psame");
leg->Draw("same");
  l->Draw("same");
T1->Draw("same");
  c1->Print("Eff_Ntrk.png");
TCanvas *c2 = new TCanvas("c2","Efficiency Fraction",500,500);
  TGraphAsymmErrors *gEffHFp4 = new TGraphAsymmErrors();
  gEffHFp4->SetName("gEffHFp4");
  gEffHFp4->BayesDivide(hHFp4ES, hHFp4DS);
for(int i=0;i<gEffHFp4->GetN();i++){
gEffHFp4->SetPointEXhigh(i,0);
gEffHFp4->SetPointEXlow(i,0);
}
  gEffHFp4->SetMarkerStyle(20);
  gEffHFp4->SetLineStyle(1);
  gEffHFp4->SetLineColor(1);
  gEffHFp4->SetMarkerColor(1);
  hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) #eta>4");
  hFrame->DrawCopy();
  gEffHFp4->Draw("Psame");
leg->Draw("same");
  l->Draw("same");

  TCanvas *c3 = new TCanvas("c3","Efficiency Fraction",500,500);
c3->SetLogx();
  TGraphAsymmErrors *regEffHFp4 = new TGraphAsymmErrors();
  regEffHFp4->SetName("regEffHFp4");
  regEffHFp4->BayesDivide(rehHFp4ES, rehHFp4DS);
for(int i=0;i<regEffHFp4->GetN();i++){
regEffHFp4->SetPointEXhigh(i,0);
regEffHFp4->SetPointEXlow(i,0);
}
  regEffHFp4->SetMarkerStyle(20);
  regEffHFp4->SetLineStyle(1);
  regEffHFp4->SetLineColor(1);
  regEffHFp4->SetMarkerColor(1);
  hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) #eta>4");
  hFrame->DrawCopy();
  regEffHFp4->Draw("Psame");
leg->Draw("same");
  l->Draw("same");
T1->Draw("same");
  c3->Print("Eff_reHF.png");

  TCanvas *c4 = new TCanvas("c4","Efficiency Fraction",500,500);
c4->SetLogx();
  TGraphAsymmErrors *regEffHF4 = new TGraphAsymmErrors();
  regEffHF4->SetName("regEffHF4");
  regEffHF4->BayesDivide(rehHF4ES, rehHF4DS);
for(int i=0;i<regEffHF4->GetN();i++){
regEffHF4->SetPointEXhigh(i,0);
regEffHF4->SetPointEXlow(i,0);
}
  regEffHF4->SetMarkerStyle(20);
  regEffHF4->SetLineStyle(1);
  regEffHF4->SetLineColor(1);
  regEffHF4->SetMarkerColor(1);
  hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} (GeV/c) |#eta|>4");
  hFrame->DrawCopy();
  regEffHF4->Draw("Psame");
leg->Draw("same");
  l->Draw("same");
T1->Draw("same");
  c4->Print("Eff_reHF4.png");


TH1D* hEffMult = (TH1D*)hMultES->Clone("hEffMult");	hEffMult->Divide(hMultDS);
TH1D* hEffHFp4 = (TH1D*)hHFp4ES->Clone("hEffHFp4");	hEffHFp4->Divide(hHFp4DS);
TH1D* rehEffHFp4 = (TH1D*)rehHFp4ES->Clone("rehEffHFp4");	rehEffHFp4->Divide(rehHFp4DS);
TH1D* rehEffHF4 = (TH1D*)rehHF4ES->Clone("rehEffHF4");	rehEffHFp4->Divide(rehHF4DS);

cout<<"----------------TrackMult----------------------------"<<endl;
cout<<"Trigger Efficiency E(M):"<<endl;
for(int i=1;i<=50;i++){
if( i % 10 == 0)  cout<<endl;
cout<<i<<"th: "<<gEffMult->GetXaxis()->GetBinLowEdge(i)<<"to"<<gEffMult->GetXaxis()->GetBinUpEdge(i)<<":"<<gEffMult->GetY()[i]<<endl;
cout<<hEffDS->GetBinContent(i)<<",";
}
cout<<endl;
cout<<"Non-DS contamination F(M):"<<endl;
for(int i=1;i<=50;i++){
if( i % 10 == 0)  cout<<endl;
cout<<hContnonDS->GetBinContent(i)<<",";
}
cout<<endl;

cout<<"DS Efficiency:"<<endl;
for(int i=1;i<=50;i++){
if( i % 10 == 0)  cout<<endl;
cout<<hEffDS->GetBinContent(i)/(1-hContnonDS->GetBinContent(i))<<",";
}
cout<<endl;

cout<<"Integral Efficiency:"<<endl;
cout<<hMultES->Integral()/hMultDS->Integral()<<endl;

cout<<"---------------------------HFp4---------------------------"<<endl;
cout<<"DS Efficiency:"<<endl;
for(int i=1;i<=50;i++){
if( i % 10 == 0)  cout<<endl;
cout<<hEffHFp4->GetBinContent(i)<<",";
}
cout<<endl;

cout<<"Integral Efficiency:"<<endl;
cout<<hHFp4ES->Integral()/hHFp4DS->Integral()<<endl;

cout<<"---------------------------reHFp4---------------------------"<<endl;
cout<<"DS Efficiency:"<<endl;
for(int i=1;i<=50;i++){
if( i % 10 == 0)  cout<<endl;
cout<<rehEffHFp4->GetBinContent(i)<<",";
}
cout<<endl;

cout<<"Integral Efficiency:"<<endl;
cout<<rehHFp4ES->Integral()/rehHFp4DS->Integral()<<endl;


TFile *f = new TFile("pPbHijing_EffCorr_forNBD.root","RECREATE");
f->cd();
hMultES->Write();
hMultDS->Write();
gEffMult->Write();
hEffMult->Write();
hHFp4ES->Write();
hHFp4DS->Write();
gEffHFp4->Write();
hEffHFp4->Write();
rehHFp4ES->Write();
rehHFp4DS->Write();
regEffHFp4->Write();
rehEffHFp4->Write();
rehHF4ES->Write();
rehHF4DS->Write();
regEffHF4->Write();
rehEffHF4->Write();
f->Close();

}
