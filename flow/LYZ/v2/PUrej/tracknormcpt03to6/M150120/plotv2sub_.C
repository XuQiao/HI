#include "par.h"
void plotv2sub_(){
	TFile *f, *ftot;
	int isSum==0;
	if(isSum==0){
	f = TFile::Open("mergedv_Prod2_sub.root");
	ftot = TFile::Open("mergedv_Prod2.root");
	}
	else{
	f = TFile::Open("mergedv_Prod_sub.root");
	ftot = TFile::Open("mergedv_Prod.root");	
	}
        int subsample = (int)(nFileAll/nsamples);
	TH1D* hv2 = new TH1D("hv2","hv2",100,0.0,0.10);
	TH1D* hv2Gau = new TH1D("hv2Gau","hv2Gau",100,-5,5);
	TH1D* hv2Gau_ = new TH1D("hv2Gau_","hv2Gau_",100,-5,5);
	int xbin=0;
	int xptbin=2;
	double avgvmean=0, avgdeltavmean=0;
	double vmean[nsamples], deltavmean[nsamples], x[nsamples];
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("D_%d/s_%d/vmean",xbin,isample));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("D_%d/s_%d/deltavmean",xbin,isample));
		vmean[isample]=(*vecvmean)[xptbin];
		deltavmean[isample]=(*vecdeltavmean)[xptbin];
		x[isample]=isample+1;
		hv2->Fill((*vecvmean)[xptbin]);
		avgvmean+=(*vecvmean)[xptbin];
		avgdeltavmean+=(*vecdeltavmean)[xptbin];
	}
	avgvmean/=nsamples;
	avgdeltavmean/=nsamples;
        for(int isample=0;isample<nsamples;isample++){
		TVectorD* vecvmean = (TVectorD*)f->Get(Form("D_%d/s_%d/vmean",xbin,isample));
		TVectorD* vecdeltavmean = (TVectorD*)f->Get(Form("D_%d/s_%d/deltavmean",xbin,isample));
		hv2Gau->Fill(((*vecvmean)[xptbin]-avgvmean)/(*vecdeltavmean)[xptbin]);
		hv2Gau_->Fill(((*vecvmean)[xptbin]-avgvmean)/avgdeltavmean);
	}
	TVectorD *vecvmean_ = (TVectorD*)ftot->Get(Form("D_%d/vmean"));
	TVectorD *vecdeltavmean_ = (TVectorD*)ftot->Get(Form("D_%d/deltavmean"));
	TLine *l = new TLine(0,(*vecvmean_)[xptbin],nsamples,(*vecvmean_)[xptbin]);
	TLine *l1 = new TLine(0,(*vecvmean_)[xptbin]+(*vecdeltavmean_)[xptbin],nsamples,(*vecvmean_)[xptbin]+(*vecdeltavmean_)[xptbin]);
	l1->SetLineStyle(2);
	TLine *l2 = new TLine(0,(*vecvmean_)[xptbin]-(*vecdeltavmean_)[xptbin],nsamples,(*vecvmean_)[xptbin]-(*vecdeltavmean_)[xptbin]);
	l2->SetLineStyle(2);
	TLine *l3 = new TLine(0,avgvmean-avgdeltavmean,nsamples,avgvmean-avgdeltavmean);
	l3->SetLineColor(2);
	l3->SetLineStyle(2);
	TLine *l4 = new TLine(0,avgvmean,nsamples,avgvmean);
	l4->SetLineColor(2);
	TLine *l5 = new TLine(0,avgvmean+avgdeltavmean,nsamples,avgvmean+avgdeltavmean);
	l5->SetLineColor(2);
	l5->SetLineStyle(2);
	TGraphErrors *gv2 = new TGraphErrors(nsamples,x,vmean,0,deltavmean);
	gv2->SetMarkerStyle(20);
	gv2->SetMarkerSize(1);
	gv2->GetYaxis()->SetRangeUser((*vecvmean_)[0]-50*(*vecdeltavmean_)[0],(*vecvmean_)[0]+50*(*vecdeltavmean_)[0]);
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(2,2);
	c1->cd(1);
	gv2->SetTitle("v2");
	gv2->Draw("AP");
	l->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
	l5->Draw("same");
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
