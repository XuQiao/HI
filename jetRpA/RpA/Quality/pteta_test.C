#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TObject.h"
#include <iostream.h>
#include <fstream.h>
#include "root_setting.h"
#include "file.h"

bool isMC=kFALSE;

void printhisto(TH1* histo, ofstream &fout){	
	fout<<histo->GetName()<<endl;
	fout<<"number of bins: "<<histo->GetNbinsX()<<endl;
        for(int i=1;i<=histo->GetNbinsX();i++)
        fout<<i<<" th bin content: "<<histo->GetBinContent(i)<<endl;
        for(int i=1;i<=histo->GetNbinsX();i++)
        fout<<i<<" th bin error: "<<histo->GetBinError(i)<<endl;
}

double* calcerr(TH1* histo1, TH1* histo2, ofstream &fout){
	fout<<histo1->GetName()<<"/"<<histo2->GetName()<<endl;
	fout<<"number of bins: "<<histo1->GetNbinsX()<<endl;
	fout<<"number of bins: "<<histo2->GetNbinsX()<<endl;	
if(histo1->GetNbinsX()==histo2->GetNbinsX()){
	double y[100],yerr[100];
        for(int i=1;i<=histo1->GetNbinsX();i++){
	y[i]=histo1->GetBinContent(i)/histo2->GetBinContent(i);
	fout<<i<<" th bin content :"<<y[i]<<endl;
	}
        for(int i=1;i<=histo1->GetNbinsX();i++){
	yerr[i]=TMath::Sqrt(TMath::Power(histo1->GetBinError(i)/histo1->GetBinContent(i),2)+TMath::Power(histo2->GetBinError(i)/histo2->GetBinContent(i),2))*y[i];
	fout<<i<<" th bin err :"<<yerr[i]<<endl;
	}
return yerr;
}
else{
	cout<<"number of bins not same!"<<endl;
	return 0;
}

}

TH1D* myRebin(TH1D* histo,int Nbin, double binbound[]){

//TH1D* rehisto=(TH1D*) histo->Clone("rehisto");
TH1D* rehisto=(TH1D*)histo->Rebin(Nbin,histo->GetName(),binbound);
for(int i=1;i<Nbin+1;i++){
rehisto->SetBinContent(i,rehisto->GetBinContent(i)/rehisto->GetBinWidth(i));
rehisto->SetBinError(i,rehisto->GetBinError(i)/rehisto->GetBinWidth(i));
}
if(binbound[Nbin]==3.0)
rehisto->GetXaxis()->SetRangeUser(binbound[0],binbound[Nbin]);
else
rehisto->GetXaxis()->SetRangeUser(binbound[10],binbound[Nbin-1]);
rehisto->Scale(1/rehisto->Integral());
fixedFontHist(rehisto);
//cout<<rehisto->Integral()<<endl;
return rehisto;
}

void pteta_test(){

using namespace std;
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetPadBottomMargin(0.12);
gStyle->SetPadTopMargin     (0.025);
gStyle->SetPadLeftMargin   (0.15);
gStyle->SetPadRightMargin (0.025);
gStyle->SetPadTickX           (1);
gStyle->SetPadTickY           (1);

/*
remove("histo_value.txt");
remove("histo_calc.txt");
ofstream foutput,foutcalc;
foutput.open("histo_value.txt",ios::app);
foutcalc.open("histo_calc.txt",ios::app);
*/

//double xrange_eta[2]={-1.465,0.535};
//double binbound_pt[]={30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,260,300,350,500};
//double binbound_pt[]={100,110,120,140,160,180,200,220,260,300,350,500};
//if(!isMC){
double binbound_eta[]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};
/*}
else{
const int Nbin_eta=20;
double binbound_eta[Nbin_eta+1];
for(int i=0;i<Nbin_eta+1;i++)
binbound_eta[i]=xrange_eta[0]+i*(xrange_eta[1]-xrange_eta[0])/Nbin_eta;
}
*/
int Nbin_eta=sizeof(binbound_eta)/sizeof(double)-1;
//for(int i=17;i<18;i++){
double xrange_eta[2]={binbound_eta[0]+1e-4,binbound_eta[Nbin_eta]-1e-4};

if(xrange_pt[0]==30.1 && isMC){
double ratioptmin=0.39, ratioptmax=1.59;
double ratioetamin=0.39, ratioetamax=1.59;
}
else if(xrange_pt[0]!=30.1  && isMC){
double ratioptmin=0.79, ratioptmax=1.19;
double ratioetamin=0.79, ratioetamax=1.19;
}
else{
double ratioptmin=0.49, ratioptmax=1.51;
double ratioetamin=0.49, ratioetamax=1.51;
}

if(!isMC)
TString Titley="Event Fraction";
else TString Titley = "#frac{d^{2}#sigma}{dp_{T}d#eta}";

TH1D* hdata_eta=hdata2F->ProjectionY("hdata_eta",hdata2F->GetXaxis()->FindBin(xrange_pt[0]),hdata2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hdata_eta_all=hdata2F->ProjectionY("hdata_eta_all",hdata2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hdata2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hdata_pt=hdata2F->ProjectionX("hdata_pt",hdata2F->GetYaxis()->FindBin(xrange_eta[0]),hdata2F->GetYaxis()->FindBin(xrange_eta[1]));
//TH1D* hdata_pt_all=hdata2F->ProjectionX("hdata_pt_all",hdata2F->GetYaxis()->FindBin(binbound_eta[0]+1e-4),hdata2F->GetYaxis()->FindBin(binbound_eta[Nbin_eta]-1e-4));

TH1D* hPPb_eta=hPPb2F->ProjectionY("hPPb_eta",hPPb2F->GetXaxis()->FindBin(xrange_pt[0]),hPPb2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPPb_eta_all=hPPb2F->ProjectionY("hPPb_eta_all",hPPb2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPPb2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPPb_pt=hPPb2F->ProjectionX("hPPb_pt",hPPb2F->GetYaxis()->FindBin(xrange_eta[0]),hPPb2F->GetYaxis()->FindBin(xrange_eta[1]));
//TH1D* hPPb_pt_all=hPPb2F->ProjectionX("hPPb_pt_all",hPPb2F->GetYaxis()->FindBin(binbound_eta[0]+1e-4),hPPb2F->GetYaxis()->FindBin(binbound_eta[Nbin_eta]-1e-4));

TH1D* hPP_eta=hPP2F->ProjectionY("hPP_eta",hPP2F->GetXaxis()->FindBin(xrange_pt[0]),hPP2F->GetXaxis()->FindBin(xrange_pt[1]));
//TH1D* hPP_eta_all=hPP2F->ProjectionY("hPP_eta_all",hPP2F->GetXaxis()->FindBin(binbound_pt[10]+1e-4),hPP2F->GetXaxis()->FindBin(binbound_pt[Nbin_pt]-1e-4));

TH1D* hPP_pt=hPP2F->ProjectionX("hPP_pt",hPP2F->GetYaxis()->FindBin(xrange_eta[0]),hPP2F->GetYaxis()->FindBin(xrange_eta[1]));
//TH1D* hPP_pt_all=hPP2F->ProjectionX("hPP_pt_all",hPP2F->GetYaxis()->FindBin(binbound_eta[0]+1e-4),hPP2F->GetYaxis()->FindBin(binbound_eta[Nbin_eta]-1e-4));

TH1D* re_hdata_pt=myRebin(hdata_pt,Nbin_pt,binbound_pt);
TH1D* re_hdata_eta=myRebin(hdata_eta,Nbin_eta,binbound_eta);
TH1D* re_hPPb_pt=myRebin(hPPb_pt,Nbin_pt,binbound_pt);
TH1D* re_hPPb_eta=myRebin(hPPb_eta,Nbin_eta,binbound_eta);
TH1D* re_hPP_pt=myRebin(hPP_pt,Nbin_pt,binbound_pt);
TH1D* re_hPP_eta=myRebin(hPP_eta,Nbin_eta,binbound_eta);
/*
TH1D* re_hdata_pt_all=myRebin(hdata_pt_all,Nbin_pt,binbound_pt);
TH1D* re_hdata_eta_all=myRebin(hdata_eta_all,Nbin_eta,binbound_eta);
TH1D* re_hPPb_pt_all=myRebin(hPPb_pt_all,Nbin_pt,binbound_pt);
TH1D* re_hPPb_eta_all=myRebin(hPPb_eta_all,Nbin_eta,binbound_eta);
TH1D* re_hPP_pt_all=myRebin(hPP_pt_all,Nbin_pt,binbound_pt;
TH1D* re_hPP_eta_all=myRebin(hPP_eta_all,Nbin_eta,binbound_eta);

if(!isMC){
re_hdata_pt->Scale(1/re_hdata_pt->Integral());
re_hdata_eta->Scale(1/re_hdata_eta->Integral());
re_hPPb_pt->Scale(1/re_hPPb_pt->Integral());
re_hPPb_eta->Scale(1/re_hPPb_eta->Integral());
re_hPP_pt->Scale(1/re_hPP_pt->Integral());
re_hPP_eta->Scale(1/re_hPP_eta->Integral());
}
fixedFontHist(re_hdata_pt);
fixedFontHist(re_hPP_pt);
fixedFontHist(re_hPPb_pt);
fixedFontHist(re_hdata_eta);
fixedFontHist(re_hPP_eta);
fixedFontHist(re_hPPb_eta);
*/
//TH1D* re_hPP_eta=hPP_eta->Clone("re_hPP_eta");
//TH1D* re_hPPb_eta=hPPb_eta->Clone("re_hPPb_eta");

re_hdata_pt->SetName("re_hdata_pt");
re_hdata_pt->SetMarkerStyle(20);
re_hdata_pt->SetMarkerSize(1);
re_hdata_pt->SetMarkerColor(1);

re_hPP_pt->SetName("re_hPP_pt");
re_hPP_pt->SetFillColor(2);
re_hPP_pt->SetMarkerColor(2);
re_hPP_pt->SetLineColor(2);

re_hPP_pt->SetFillStyle(3004);
re_hPP_pt->SetMarkerStyle(0);
re_hPP_pt->SetTitle("");

re_hPP_eta->SetFillStyle(3004);
re_hPP_eta->SetMarkerStyle(0);
re_hPP_eta->SetFillColor(2);
re_hPP_eta->SetLineColor(2);
re_hPP_eta->SetTitle("");

if(!isMC){
    re_hPP_pt->SetMaximum(1.5);
    re_hPP_pt->SetMinimum(re_hPP_pt->GetMinimum()/2);
}

re_hPP_pt->GetYaxis()->SetTitle(Titley);
re_hPP_pt->GetXaxis()->SetTitle("p^{jet}_{T} [GeV/c]");
//re_hPP_pt->GetXaxis()->SetTitleOffset(1.2);
//re_hPPb_pt->GetXaxis()->SetTitleOffset(1.2);
re_hPPb_pt->SetName("re_hPPb_pt");
re_hPPb_pt->SetMarkerStyle(28);
re_hPPb_pt->SetMarkerSize(1.2);
re_hPPb_pt->SetMarkerColor(4);
re_hPPb_pt->SetLineColor(4);

re_hdata_eta->SetName("re_hdata_eta");
re_hdata_eta->SetMarkerStyle(20);
re_hdata_eta->SetMarkerSize(1);
re_hdata_eta->SetMarkerColor(1);
if(!isMC){
re_hPP_eta->SetMaximum(re_hPPb_eta->GetMaximum()*1.5);
re_hPP_eta->SetMinimum(re_hPP_eta->GetMinimum()/2);
}
re_hPP_eta->SetName("re_hPP_eta");
re_hPP_eta->GetYaxis()->SetTitle(Titley);
re_hPP_eta->GetXaxis()->SetTitle("#eta");
//re_hPP_eta->GetXaxis()->SetTitleOffset(0.8);
re_hPPb_eta->SetName("re_hPPb_eta");
re_hPPb_eta->SetMarkerStyle(28);
re_hPPb_eta->SetMarkerSize(1.2);
re_hPPb_eta->SetMarkerColor(4);
re_hPPb_eta->SetLineColor(4);
//re_hPPb_eta->GetXaxis()->SetTitleOffset(0.8);

c1 = new TCanvas("c1"," ",700,700);
c2 = new TCanvas("c2"," ",700,700);
makeMultiPanelCanvas(c1,1,2,0.01,0.1,0.14,0.15,0.1);
makeMultiPanelCanvas(c2,1,2,0.01,0.1,0.14,0.15,0.1);

c1->cd(1)->SetLogy();
re_hPP_pt->Draw("HIST E1");
if(!isMC)
re_hdata_pt->Draw("E1same");
re_hPPb_pt->Draw("E1same");
TLatex *T=new TLatex(0.40,0.90,Form("#eta: %.1f--%.1f ",xrange_eta[0],xrange_eta[1]));
T->SetNDC();
T->SetTextAlign(12);
T->SetTextSize(0.05);
T->SetTextColor(1);
T->SetTextFont(42);
T->Draw("same");
TLine *l1 =new TLine(30,1,600,1);
TLegend *leg1=new TLegend(0.65,0.7,0.85,0.9);
leg1->SetBorderSize(0);
leg1->SetFillStyle(0);
//leg1->SetTextFont(63);
leg1->SetTextSize(0.05);
//leg1->SetTextSize(17);
leg1->SetFillColor(0);
//leg1->SetTextSize(0.03);
if(!isMC)
    leg1->AddEntry(re_hdata_pt,"Data","lp");
leg1->AddEntry(re_hPPb_pt,"PYTHIA+HIJING","lp");
leg1->AddEntry(re_hPP_pt,"PYTHIA","lpf");
leg1->Draw("same");
c1->cd(2);
l1->SetLineStyle(2);
if(isMC){ 	// for MC PPb vs Pb check
	TH1D* ratio_pt=(TH1D*)re_hPPb_pt->Clone("ratio_pt");
	ratio_pt->Divide(re_hPP_pt);
	ratio_pt->SetTitle("");
	ratio_pt->GetYaxis()->SetTitle("#frac{PPb}{PP}");
	ratio_pt->GetXaxis()->CenterTitle();
        ratio_pt->GetYaxis()->CenterTitle();
	ratio_pt->SetMarkerStyle(20);
	ratio_pt->SetMarkerSize(1);
	ratio_pt->SetMinimum(ratioptmin);
	ratio_pt->SetMaximum(ratioptmax);
	ratio_pt->Draw("E1");
//	ratio_pt->Fit("pol0");
	l1->Draw("same");
//	c1->Print(Form("pic/test/jetpteta_ProX%.1f-%.1f_%3.f-%3.f_PPbvsPP.png",xrange_eta[0],xrange_eta[1],xrange_pt[0],xrange_pt[1]));
	}

else{	// for Data vs MC study
	TH1D* ratio_DatavsPPb_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPPb_pt");
	ratio_DatavsPPb_pt->Divide(re_hPPb_pt);
	ratio_DatavsPPb_pt->SetTitle("");
	ratio_DatavsPPb_pt->GetYaxis()->SetTitle("#frac{Data}{MC}");
	ratio_DatavsPPb_pt->GetXaxis()->CenterTitle();
	ratio_DatavsPPb_pt->GetYaxis()->CenterTitle();
	ratio_DatavsPPb_pt->SetMarkerStyle(28);
	ratio_DatavsPPb_pt->SetLineColor(4);
	ratio_DatavsPPb_pt->SetMarkerColor(4);
	ratio_DatavsPPb_pt->SetMarkerSize(1);
/*	for(int i=1;i<=ratio_DatavsPPb_pt->GetNbinsX();i++){
	ratio_DatavsPPb_pt->SetBinError(i,calcerr(re_hdata_pt,re_hPPb_pt,foutcalc)[i]);
	}*/
	ratio_DatavsPPb_pt->SetMinimum(ratioptmin);
	ratio_DatavsPPb_pt->SetMaximum(ratioptmax);
	ratio_DatavsPPb_pt->Draw("E1");
	//ratio_DatavsPPb_pt->Fit("pol0");

	TH1D* ratio_DatavsPP_pt=(TH1D*)re_hdata_pt->Clone("ratio_DatavsPP_pt");
	ratio_DatavsPP_pt->Divide(re_hPP_pt);
	ratio_DatavsPP_pt->SetTitle("");
	ratio_DatavsPP_pt->SetMarkerStyle(20);
	ratio_DatavsPP_pt->SetMarkerColor(2);
	ratio_DatavsPP_pt->SetLineColor(2);
	ratio_DatavsPP_pt->SetMarkerSize(1);
	
	ratio_DatavsPP_pt->Draw("E1same");
	//ratio_DatavsPP_pt->Fit("pol0");
	TLegend *leg2=new TLegend(0.65,0.80,0.85,0.95);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
//	leg2->SetTextFont(63);
	leg2->SetTextSize(0.05);
	leg2->SetFillColor(0);
	leg2->AddEntry(ratio_DatavsPPb_pt,"Data/PYTHIA+HIJING","lp");
	leg2->AddEntry(ratio_DatavsPP_pt,"Data/PYTHIA","lp");
	leg2->Draw("same");
	l1->Draw("same");

//	printhisto(ratio_DatavsPPb_pt,foutput);
//	printhisto(ratio_DatavsPP_pt,foutput);

//	c1->Print(Form("pic/test/jetpteta_ProX%.1f-%.1f_%3.f-%3.f_DatavsPPbvsPP.png",xrange_eta[0],xrange_eta[1],binbound_pt[0],binbound_pt[Nbin_pt]));
	}

c2->cd(1);
re_hPP_eta->Draw("HIST E1");
if(!isMC)
re_hdata_eta->Draw("E1same");
	re_hPPb_eta->Draw("E1same");
	TLegend *leg3=new TLegend(0.65,0.7,0.85,0.9);
	leg3->SetTextSize(0.05);
	leg3->SetBorderSize(0);
	leg3->SetFillStyle(0);
//	leg3->SetTextFont(63);
//	leg3->SetTextSize(17);
	leg3->SetFillColor(0);
//	leg3->SetTextSize(0.03);

TLatex *T1=new TLatex(0.40,0.90,Form("p_{T}: %3.f--%3.f [GeV/c]",xrange_pt[0],xrange_pt[1]));
T1->SetNDC();
T1->SetTextAlign(12);
T1->SetTextSize(0.05);
T1->SetTextColor(1);
T1->SetTextFont(42);
T1->Draw("same");
if(!isMC)
	leg3->AddEntry(re_hdata_eta,"Data","lp");
leg3->AddEntry(re_hPPb_eta,"PYTHIA+HIJING","lp");
leg3->AddEntry(re_hPP_eta,"PYTHIA","lpf");
leg3->Draw("same");

c2->cd(2);
//re_hPP_eta->GetXaxis()->SetTitleOffset(0.8);
TLine *l2 =new TLine(30,1,600,1);
l2->SetLineStyle(2);
if(isMC){	//for MC PPb vs PP check
	TH1D* ratio_eta=(TH1D*)re_hPPb_eta->Clone("ratio_eta");
	ratio_eta->Divide(re_hPP_eta);
	ratio_eta->SetTitle("");
	ratio_eta->GetYaxis()->SetTitle("#frac{PPb}{PP}");
	ratio_eta->SetMarkerStyle(20);
	ratio_eta->SetMarkerSize(1);
	ratio_eta->SetMinimum(ratioetamin);
	ratio_eta->SetMaximum(ratioetamax);
	ratio_eta->Draw("E1");
	//ratio_eta->Fit("pol0");
	l2->Draw("same");
//	c2->Print(Form("pic/test/jetpteta_ProY%.1f-%.1f_%3.f-%3.f_PPbvsPP.png",xrange_eta[0],xrange_eta[1],xrange_pt[0],xrange_pt[1]));
	}
else{	// for Data vs MC study
	TH1D* ratio_DatavsPPb_eta=(TH1D*)re_hdata_eta->Clone("ratio_DatavsPPb_eta");
	ratio_DatavsPPb_eta->Divide(re_hPPb_eta);
	ratio_DatavsPPb_eta->SetTitle("");
	ratio_DatavsPPb_eta->GetYaxis()->SetTitle("#frac{Data}{MC}");
	ratio_DatavsPPb_eta->GetXaxis()->CenterTitle();
	ratio_DatavsPPb_eta->GetYaxis()->CenterTitle();
	ratio_DatavsPPb_eta->SetMarkerStyle(28);
	ratio_DatavsPPb_eta->SetMarkerSize(1);
	ratio_DatavsPPb_eta->SetMarkerColor(4);
	ratio_DatavsPPb_eta->SetLineColor(4);
	/*for(int i=1;i<=ratio_DatavsPPb_eta->GetNbinsX();i++){
        ratio_DatavsPPb_eta->SetBinError(i,calcerr(re_hdata_eta,re_hPPb_eta,foutcalc)[i]);
        }
*/
	ratio_DatavsPPb_eta->SetMinimum(ratioetamin);
	ratio_DatavsPPb_eta->SetMaximum(ratioetamax);
	ratio_DatavsPPb_eta->Draw("E1");
	//ratio_DatavsPPb_eta->Fit("pol0");
	
	TH1D* ratio_DatavsPP_eta=(TH1D*)re_hdata_eta->Clone("ratio_DatavsPP_eta");
	ratio_DatavsPP_eta->Divide(re_hPP_eta);
	ratio_DatavsPP_eta->SetTitle("");
	ratio_DatavsPP_eta->SetMarkerStyle(20);
	ratio_DatavsPP_eta->SetMarkerSize(1);
	ratio_DatavsPP_eta->SetMarkerColor(2);
	ratio_DatavsPP_eta->SetLineColor(2);
	ratio_DatavsPP_eta->Draw("E1same");
	//ratio_DatavsPP_eta->Fit("pol0");
	TLegend *leg4=new TLegend(0.65,0.80,0.85,0.95);
	leg4->SetTextSize(0.05);
	leg4->SetBorderSize(0);
	leg4->SetFillStyle(0);
//	leg4->SetTextFont(63);
//	leg4->SetTextSize(17);
	leg4->SetFillColor(0);
	leg4->AddEntry(ratio_DatavsPPb_eta,"Data/PYTHIA+HIJING","lp");
	leg4->AddEntry(ratio_DatavsPP_eta,"Data/PYTHIA","lp");
	leg4->Draw("same");
	l2->Draw("same");

//	c2->Print(Form("pic/test/jetpteta_ProY%.1f-%.1f_%3.f-%3.f_DatavsPPbvsPP.png",xrange_eta[0],xrange_eta[1],xrange_pt[0],xrange_pt[1]));

}
//}
}
