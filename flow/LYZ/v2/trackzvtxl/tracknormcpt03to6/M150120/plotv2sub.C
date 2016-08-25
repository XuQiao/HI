#include "par.h"
void plotv2sub(){
	TFile *f, *ftot;
	int isSum=0;
	if(isSum==0){
	f = TFile::Open("mergedv_Prod2_sub.root");
	ftot = TFile::Open("mergedv_Prod2.root");
	}
	else{
	f = TFile::Open("mergedv_Prod_sub.root");
	ftot = TFile::Open("mergedv_Prod.root");
	}
        int subsample = (int)(nFileAll/nsamples);
	TH1D* hv2 = new TH1D("hv2","hv2",100,0.0,0.06);
	TH1D* hv2Gau = new TH1D("hv2Gau","hv2Gau",100,-1,1);
	int xbin=0;
	double avgvmean=0, avgdeltavmean=0;
	double vmean[nsamples], x[nsamples];
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("D_%d/s_%d/vmean",xbin,isample));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("D_%d/s_%d/deltavmean",xbin,isample));
		vmean[isample]=(*vecvmean)[0];
		x[isample]=isample;
		hv2->Fill((*vecvmean)[0]);
		avgvmean+=(*vecvmean)[0];
		avgdeltavmean+=(*vecdeltavmean)[0];
	}
	avgvmean/=nsamples;
	avgdeltavmean/=nsamples;
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("D_%d/s_%d//vmean",xbin,isample));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("D_%d/s_%d/deltavmean",xbin,isample));
		//hv2Gau->Fill(((*vecvmean)[0]-avgvmean)/(*vecdeltavmean)[0]);
		hv2Gau->Fill(((*vecvmean)[0]-avgvmean)/avgdeltavmean);
	}
	TVectorD *vecvmean_ = (TVectorD*)ftot->Get(Form("D_%d/vmean"));
	TVectorD *vecdeltavmean_ = (TVectorD*)ftot->Get(Form("D_%d/deltavmean"));
	TLine *l = new TLine(0,(*vecvmean_)[0],nsamples,(*vecvmean_)[0]);
	TLine *l1 = new TLine(0,(*vecvmean_)[0]+(*vecdeltavmean_)[0],nsamples,(*vecvmean_)[0]+(*vecdeltavmean_)[0]);
	l1->SetLineStyle(2);
	TLine *l2 = new TLine(0,(*vecvmean_)[0]-(*vecdeltavmean_)[0],nsamples,(*vecvmean_)[0]-(*vecdeltavmean_)[0]);
	l2->SetLineStyle(2);
	TGraph *gv2 = new TGraph(nsamples,x,vmean);
	gv2->SetMarkerStyle(20);
	gv2->SetMarkerSize(1);
	gv2->GetYaxis()->SetRangeUser((*vecvmean_)[0]-50*(*vecdeltavmean_)[0],(*vecvmean_)[0]+50*(*vecdeltavmean_)[0]);
	TCanvas *c1 = new TCanvas();
	hv2->Draw();
	TCanvas *c2 = new TCanvas();
	hv2Gau->Draw();
	TCanvas *c3 = new TCanvas();
	gv2->Draw("AP");
	l->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	c1->SaveAs("hv2.gif");
}
