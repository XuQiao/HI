#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void DrawESDS(){
gStyle->SetOptStat(kFALSE);
gStyle->SetErrorX(0);


//TFile *fileDS = new TFile("../../AnaEff/test/output/NoVRw/DS.root","ReadOnly");
TFile *fileDS = new TFile("../../AnaEff/test/output/VRw/DS_v2.root","ReadOnly");
//TFile *fileES = new TFile("../../ESAnaEff/test/output/NoVRw/ES.root","ReadOnly");
TFile *fileES = new TFile("../../ESAnaEff/test/output/VRw/ES_v1.root","ReadOnly");

TH1D *hMultDS = (TH1D*) fileDS->Get("demo/hMultDS");	hMultDS->SetName("hMultDS");
TH1D *hMultES = (TH1D*) fileES->Get("demo/hMultES");	hMultES->SetName("hMultES");
TH1D *hMultESDS = (TH1D*) fileES->Get("demo/hMultDS");	hMultESDS->SetName("hMultESDS");
TH1D *hMultESnonDS = (TH1D*) fileES->Get("demo/hMultnonDS");	hMultESnonDS->SetName("hMultESnonDS");

TH1D* hEffDS = (TH1D*)hMultESDS->Clone("hEffDS");	hEffDS->Divide(hMultDS);
TH1D* hContnonDS = (TH1D*)hMultESnonDS->Clone("hContnonDS");	hContnonDS->Divide(hMultES);

TH1D *hHFp4DS = (TH1D*) fileDS->Get("demo/hHFp4DS");    hHFp4DS->SetName("hHFp4DS");
TH1D *hHFp4ES = (TH1D*) fileES->Get("demo/hHFp4ES");    hHFp4ES->SetName("hHFp4ES");

/*
const int nHFbin=100;double HFbin[nHFbin+1];
for(int i=0;i<=nHFbin;i++)
HFbin[i]=i;
*/
double HFbin[]={0,0.23,0.47,0.67,0.87,1.05,1.24,1.43,1.62,1.82,2.03,2.24,2.47,2.7,2.95,3.21,3.48,3.77,4.07,4.38,4.72,5.06,5.41,5.79,6.17,6.57,6.98,7.41,7.85,8.3,8.76,9.24,9.74,10.25,10.79,11.33,11.91,12.51,13.14,13.81,14.51,15.27,16.09,17,17.99,19.14,20.48,22.14,24.35,27.93,90.29};
//double HFbin[]={0,0.22,0.46,0.66,0.85,1.04,1.22,1.41,1.6,1.8,2,2.22,2.44,2.68,2.92,3.18,3.45,3.74,4.04,4.35,4.68,5.03,5.38,5.75,6.14,6.54,6.95,7.38,7.82,8.27,8.74,9.22,9.72,10.23,10.77,11.31,11.89,12.49,13.12,13.79,14.5,15.26,16.08,16.99,17.99,19.13,20.48,22.14,24.36,27.94,90.29};
const int nHFbin=sizeof(HFbin)/sizeof(double)-1;

TH1D* rehHFp4DS = (TH1D*) hHFp4DS->Rebin(nHFbin,"rehHFp4DS",HFbin);	normalizeByBinWidth(rehHFp4DS);
TH1D* rehHFp4ES = (TH1D*) hHFp4ES->Rebin(nHFbin,"rehHFp4ES",HFbin);	normalizeByBinWidth(rehHFp4ES);

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
}

//hMultDS->Draw();
//hMultES->SetLineColor(2);
//hMultES->Draw("same");

TH1D* hFrame=new TH1D("","",1000,0,100);
fixedFontHist(hFrame,1.2,1.6);
TLegend *leg=new TLegend(0.3,0.25,0.8,0.4);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetTextFont(42);
leg->SetTextSize(0.04);
c1 = new TCanvas("c1"," ",500,500);
makeMultiPanelCanvas(c1,1,1,-0.06,0,0.14,0.12,0.03);
c1->cd(1)->SetLogy();  
  hMultES->SetMarkerStyle(20);
  hMultES->SetLineStyle(1);
  hMultES->SetLineColor(1);
  hMultES->SetMarkerColor(1);
  hMultDS->SetMarkerStyle(24);
  hMultDS->SetLineStyle(1);
  hMultDS->SetLineColor(2);
  hMultDS->SetMarkerColor(2);
  hFrame->GetXaxis()->SetRangeUser(0,30);
  hFrame->GetYaxis()->SetRangeUser(1,1e6);
  hFrame->GetXaxis()->SetTitle("track Multiplicty");
  hFrame->GetYaxis()->SetTitle("# of Events");
  hFrame->DrawCopy();
	hMultES->Draw("same");
	hMultDS->Draw("same");
	leg->SetHeader("Hijing");
	leg->AddEntry(hMultES,"Event Selected","lp");
	leg->AddEntry(hMultDS,"Double Sided","lp");
	leg->Draw("same");
  c1->Print("DrawESDS_Ntrk.png");
/*
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
*/
c3 = new TCanvas("c3"," ",500,500);
makeMultiPanelCanvas(c3,1,1,-0.06,0,0.14,0.12,0.03);
c3->cd(1)->SetLogy();
c3->cd(1)->SetLogx();
  rehHFp4ES->SetMarkerStyle(20);
  rehHFp4ES->SetLineStyle(1);
  rehHFp4ES->SetLineColor(1);
  rehHFp4ES->SetMarkerColor(1);
  rehHFp4DS->SetMarkerStyle(24);
  rehHFp4DS->SetLineStyle(1);
  rehHFp4DS->SetLineColor(2);
  rehHFp4DS->SetMarkerColor(2);
  hFrame->GetXaxis()->SetRangeUser(1e-1,30);
  hFrame->GetYaxis()->SetRangeUser(5e3,1e6);
  hFrame->GetXaxis()->SetTitle("HF #Sigma E_{T} #eta>4");
  hFrame->GetYaxis()->SetTitle("# of Events");
  hFrame->DrawCopy();
        rehHFp4ES->Draw("same");
        rehHFp4DS->Draw("same");
        leg->Draw("same");
  c3->Print("DrawESDS_reHF.png");
}
