#include "par.h"
void plotV2sub(){
	TFile *f, *ftot;
	int isSum=0;
	if(isSum==0){
	f = TFile::Open("mergedV_Prod_sub.root");
	ftot = TFile::Open("mergedV_Prod.root");
	}
	else{
	f = TFile::Open("mergedV_Sum_sub.root");
	ftot = TFile::Open("mergedV_Sum.root");
	}
        int subsample = (int)(nFileAll/nsamples);
	TH1D* hV2 = new TH1D("hV2","hV2",30,0.040,0.070);
	TH1D* hV2Gau = new TH1D("hV2Gau","hV2Gau",1000,-3,3);
	TH1D* hV2Gau_ = new TH1D("hV2Gau_","hV2Gau_",1000,-3,3);
	int xbin=0;
	double avgVmean=0, avgdeltaVmean=0;
	double Vmean[nsamples], deltaVmean[nsamples], x[nsamples];
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/Vmean",isample,xbin));
		TVectorD* vecdeltaVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltaVmean",isample,xbin));
		Vmean[isample]=(*vecVmean)[0];
		deltaVmean[isample]=(*vecdeltaVmean)[0];
		x[isample]=isample+1;
		hV2->Fill((*vecVmean)[0]);
		avgVmean+=(*vecVmean)[0];
		avgdeltaVmean+=(*vecdeltaVmean)[0];
	}
	avgVmean/=nsamples;
	avgdeltaVmean/=nsamples;
	cout<<avgdeltaVmean<<endl;
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/Vmean",isample,xbin));
		TVectorD* vecdeltaVmean = (TVectorD*)f->Get(Form("s_%d/D_%d/deltaVmean",isample,xbin));
		hV2Gau->Fill(((*vecVmean)[0]-avgVmean)/(*vecdeltaVmean)[0]);
		hV2Gau_->Fill(((*vecVmean)[0]-avgVmean)/avgdeltaVmean);
	}
	TVectorD *vecVmean_ = (TVectorD*)ftot->Get(Form("D_%d/Vmean"));
	TVectorD *vecdeltaVmean_ = (TVectorD*)ftot->Get(Form("D_%d/deltaVmean"));
	TLine *l = new TLine(0,avgVmean,nsamples,avgVmean);
	l->SetLineColor(2);
	TLine *l1 = new TLine(0,avgVmean+avgdeltaVmean,nsamples,avgVmean+avgdeltaVmean);
	l1->SetLineColor(2);
	l1->SetLineStyle(2);
	TLine *l2 = new TLine(0,avgVmean-avgdeltaVmean,nsamples,avgVmean-avgdeltaVmean);
	l2->SetLineColor(2);
	l2->SetLineStyle(2);
	TLine *l_ = new TLine(0,(*vecVmean_)[0],nsamples,(*vecVmean_)[0]);
	TLine *l1_ = new TLine(0,(*vecVmean_)[0]+(*vecdeltaVmean_)[0],nsamples,(*vecVmean_)[0]+(*vecdeltaVmean_)[0]);
	l1_->SetLineStyle(2);
	TLine *l2_ = new TLine(0,(*vecVmean_)[0]-(*vecdeltaVmean_)[0],nsamples,(*vecVmean_)[0]-(*vecdeltaVmean_)[0]);
	l2_->SetLineStyle(2);
	TGraphErrors *gV2 = new TGraphErrors(nsamples,x,Vmean,0,deltaVmean);
	gV2->SetMarkerStyle(20);
	gV2->SetMarkerSize(1);
	gV2->GetYaxis()->SetRangeUser((*vecVmean_)[0]-30*(*vecdeltaVmean_)[0],(*vecVmean_)[0]+30*(*vecdeltaVmean_)[0]);
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(2,2);
	c1->cd(1);
	gV2->SetTitle("V2");
	gV2->Draw("AP");
	l_->Draw("same");
	l1_->Draw("same");
	l2_->Draw("same");
	l->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	c1->cd(2);
	hV2->Draw();
	c1->cd(3);
	hV2Gau->GetXaxis()->SetTitle("(V_{i}-<V_{i}>)/#delta V_{i}");
	hV2Gau->Draw();
	c1->cd(4);
	hV2Gau_->GetXaxis()->SetTitle("(V_{i}-<V_{i}>)/<#delta V_{i}>");
	hV2Gau_->Draw();
	c1->SaveAs("hV2.gif");
}
