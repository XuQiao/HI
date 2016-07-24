#include "par.h"
void plotv2sub_(){
	TFile *f, *ftot;
	if(isSum==0){
	f = TFile::Open("mergedv_Prod2_sub.root");
	ftot = TFile::Open("mergedv_Prod2.root");
	}
	else{
	f = TFile::Open("mergedv_Prod_sub.root");
	ftot = TFile::Open("mergedv_Prod.root");	
	}
        int subsample = (int)(nFileAll/nsamples);
	TH1D* hv2 = new TH1D("hv2","hv2",100,0.053,0.055);
	TH1D* hv2Gau = new TH1D("hv2Gau","hv2Gau",100,-1,1);
	TH1D* hv2Gau_ = new TH1D("hv2Gau_","hv2Gau_",100,-1,1);
	int xbin=0;
	double avgvmean=0, avgdeltavmean=0;
	double vmean[nsamples], deltavmean[nsamples], x[nsamples];
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("s_%d/D_%d/vmean",isample,xbin));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltavmean",isample,xbin));
		vmean[isample]=(*vecvmean)[0];
		deltavmean[isample]=(*vecdeltavmean)[0];
		x[isample]=isample+1;
		hv2->Fill((*vecvmean)[0]);
		avgvmean+=(*vecvmean)[0];
		avgdeltavmean+=(*vecdeltavmean)[0];
	}
	avgvmean/=nsamples;
	avgdeltavmean/=nsamples;
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("s_%d/D_%d/vmean",isample,xbin));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltavmean",isample,xbin));
		hv2Gau->Fill(((*vecvmean)[0]-avgvmean)/(*vecdeltavmean)[0]);
		hv2Gau_->Fill(((*vecvmean)[0]-avgvmean)/avgdeltavmean);
	}
	TVectorD *vecvmean_ = (TVectorD*)ftot->Get(Form("D_%d/vmean"));
	TVectorD *vecdeltavmean_ = (TVectorD*)ftot->Get(Form("D_%d/deltavmean"));
	TLine *l = new TLine(0,(*vecvmean_)[0],15,(*vecvmean_)[0]);
	TLine *l1 = new TLine(0,(*vecvmean_)[0]+(*vecdeltavmean_)[0],15,(*vecvmean_)[0]+(*vecdeltavmean_)[0]);
	l1->SetLineStyle(2);
	TLine *l2 = new TLine(0,(*vecvmean_)[0]-(*vecdeltavmean_)[0],15,(*vecvmean_)[0]-(*vecdeltavmean_)[0]);
	l2->SetLineStyle(2);
	TGraphErrors *gv2 = new TGraphErrors(nsamples,x,vmean,0,deltavmean);
	gv2->SetMarkerStyle(20);
	gv2->SetMarkerSize(1);
	gv2->GetYaxis()->SetRangeUser((*vecvmean_)[0]-10*(*vecdeltavmean_)[0],(*vecvmean_)[0]+10*(*vecdeltavmean_)[0]);
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(2,2);
	c1->cd(1);
	gv2->SetTitle("v2");
	gv2->Draw("AP");
	l->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	c1->cd(2);
	hv2->Draw();
	c1->cd(3);
	hv2Gau->GetXaxis()->SetTitle("(v_{i}-<v_{i}>)/#delta v_{i}");
	hv2Gau->Draw();
	c1->cd(4);
	hv2Gau_->GetXaxis()->SetTitle("(v_{i}-<v_{i}>)/<#delta v_{i}>");
	hv2Gau_->Draw();
	c1->SaveAs("hv2.gif");
}
