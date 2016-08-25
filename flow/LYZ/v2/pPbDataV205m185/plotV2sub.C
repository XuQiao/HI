#include "par.h"
void plotV2sub(){
	TFile *f, *ftot;
	if(isSum==0){
	f = TFile::Open("mergedV_Prod_sub.root");
	ftot = TFile::Open("mergedV_Prod.root");
	}
	else{
	f = TFile::Open("mergedV_Sum_sub.root");
	ftot = TFile::Open("mergedV_Sum.root");
	}
        int subsample = (int)(nFileAll/nsamples);
	TH1D* hV2 = new TH1D("hV2","hV2",100,0.05,0.06);
	TH1D* hV2Gau = new TH1D("hV2Gau","hV2Gau",100,-1,1);
	int xbin=0;
	double avgVmean=0, avgdeltaVmean=0;
	double Vmean[nsamples], x[nsamples];
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/Vmean",isample,xbin));
		TVectorD* vecdeltaVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltaVmean",isample,xbin));
		Vmean[isample]=(*vecVmean)[0];
		x[isample]=isample;
		hV2->Fill((*vecVmean)[0]);
		avgVmean+=(*vecVmean)[0];
		avgdeltaVmean+=(*vecdeltaVmean)[0];
	}
	avgVmean/=nsamples;
	avgdeltaVmean/=nsamples;
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/Vmean",isample,xbin));
		TVectorD* vecdeltaVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltaVmean",isample,xbin));
		//hV2Gau->Fill(((*vecVmean)[0]-avgVmean)/(*vecdeltaVmean)[0]);
		hV2Gau->Fill(((*vecVmean)[0]-avgVmean)/avgdeltaVmean);
	}
	TVectorD *vecVmean_ = (TVectorD*)ftot->Get(Form("D_%d/Vmean"));
	TVectorD *vecdeltaVmean_ = (TVectorD*)ftot->Get(Form("D_%d/deltaVmean"));
	TLine *l = new TLine(0,(*vecVmean_)[0],15,(*vecVmean_)[0]);
	TLine *l1 = new TLine(0,(*vecVmean_)[0]+(*vecdeltaVmean_)[0],15,(*vecVmean_)[0]+(*vecdeltaVmean_)[0]);
	l1->SetLineStyle(2);
	TLine *l2 = new TLine(0,(*vecVmean_)[0]-(*vecdeltaVmean_)[0],15,(*vecVmean_)[0]-(*vecdeltaVmean_)[0]);
	l2->SetLineStyle(2);
	TGraph *gV2 = new TGraph(nsamples,x,Vmean);
	gV2->SetMarkerStyle(20);
	gV2->SetMarkerSize(1);
	gV2->GetYaxis()->SetRangeUser((*vecVmean_)[0]-2*(*vecdeltaVmean_)[0],(*vecVmean_)[0]+2*(*vecdeltaVmean_)[0]);
	TCanvas *c1 = new TCanvas();
	hV2->Draw();
	TCanvas *c2 = new TCanvas();
	hV2Gau->Draw();
	TCanvas *c3 = new TCanvas();
	gV2->Draw("AP");
	l->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	c1->SaveAs("hV2.gif");
}
