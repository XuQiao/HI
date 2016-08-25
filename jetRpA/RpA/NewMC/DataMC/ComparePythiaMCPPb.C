#include "../../Quality/root_setting.h"
#include "../produceandcheck/file.h"

TH1* myRebin(TH1* histo,int Nbin, double binbound[]){

TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
TString histoname=Form("%s",histo->GetName());
normalizeByBinWidth(rehisto);
//rehisto->Scale(1.0/rehisto->Integral(rehisto->FindBin(xrange_pt[0]),rehisto->FindBin(xrange_pt[1])));
return rehisto;
}

void ComparePythiaMCPPb(){
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptStat(100);
gStyle->SetErrorX(0);
gStyle->SetLabelFont(70);
c1 = new TCanvas("c1"," ",1200,600);
makeMultiPanelCanvas(c1,4,2,0,0,0.25,0.2,0.03);
TH1F* hFrame1=new TH1F("","",1000,0,1000);
fixedFontHist(hFrame1,2.0,3.0);
hFrame1->SetTitle("");
TLatex T1;
T1.SetNDC();
T1.SetTextSize(0.065);
T1.SetTextFont(42);
TString jettype = "refjet";
for(int i=0;i<Neta;i++){
	if(canvas[i]==0){
		double ybase=0.15;	double xbase=0.28;
	}
	else if(canvas[i]>4){
		double ybase=0.30;	double xbase=0.05;
	}
	else{
		double ybase=0.15;	double xbase=0.05;
	}
    if(i==7){
	TH1F *hPPb_pt=  (TH1F*)fMCPPb->Get(Form("%spt",jettype.Data()));
        TH1F *hPP_pt1 = (TH1F*)fMCPP->Get(Form("ak3GenJetSpectrum_QCD10001_%s/JetSpectrum_Fine",etabinnamesym[0].Data()));
        TH1F *hPP_pt2 = (TH1F*)fMCPP->Get(Form("ak3GenJetSpectrum_QCD10001_%s/JetSpectrum_Fine",etabinnamesym[1].Data()));
        TH1F *hPP_pt = (TH1F*)hPP_pt1->Clone();
        hPP_pt->Add(hPP_pt2);
    }
    else{
	TH1F *hPPb_pt= (TH1F*)fMCPPb->Get(Form("%sptEtaBin%s",jettype.Data(),etabinnameswap[i].Data()));
        hPPb_pt->Scale(1.0/etabin[i]);
	TH1F *hPP_pt= (TH1F*)fMCPP->Get(Form("ak3GenJetSpectrum_QCD10001_%s/JetSpectrum_Fine",etabinnamesym[TMath::Abs(i-3)].Data()));
    }
    hPPb_pt->Scale(1/0.87);
    c1->cd(canvas[i]+1)->SetGridx();
    if(canvas[i]==0  || canvas[i]==4){
        hFrame1->GetYaxis()->SetTitle("Pythia+Hijing/Pythia");
        hFrame1->SetLabelSize(18,"Y");
    }
    else{
        hFrame1->SetLabelSize(0,"Y");
        hFrame1->GetYaxis()->SetTitle("");
    }
    if(canvas[i]>=4)
        hFrame1->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
	
   hFrame1->GetXaxis()->SetNdivisions(510);
   hFrame1->GetXaxis()->SetLimits(28,599);
   hFrame1->SetMinimum(0.52);
   hFrame1->SetMaximum(1.49);
   hFrame1->DrawCopy();
   if(canvas[i]!=4){
    TH1F* rehisto_hPPb_pt=(TH1F*)myRebin(hPPb_pt,Nbin_pt,binbound_pt);
    TH1F* rehisto_hPP_pt=(TH1F*)myRebin(hPP_pt,Nbin_pt,binbound_pt);
    TH1F* ratio_pt =(TH1F*)rehisto_hPPb_pt->Clone(Form("ratio_pt_%d",i));
    ratio_pt->Divide(rehisto_hPP_pt);
    ratio_pt->SetMarkerSize(1.4);
    ratio_pt->SetMarkerStyle(20);
    ratio_pt->SetMarkerColor(1);
    ratio_pt->SetLineColor(1);
    ratio_pt->Draw("same");
    TLine *l = new TLine(28,1,599,1);
    l->SetLineStyle(2);
    l->Draw("same");
    T1.DrawLatex(0.4,ybase-0.05,etastring[i]);
   }

    }
}
    

