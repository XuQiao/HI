#include <iostream>
#include <iomanip>
#include "par.h"
void plotgvsr(){
	ifstream outG1("G2mod2.dat");
	ifstream outR("R.dat");
	TFile *f = TFile::Open("mergedV_Sum.root");
	int xbin=0;
	int xpt=0;
	int xtheta=0;
	TVectorD *vecDr = f->Get(Form("D_%d/r",xbin));
	TVectorD *vecDg2 = f->Get(Form("D_%d/D_%d/D_%d/G2",xbin,xpt,xtheta));
	TVectorD *vecDsigma2 = f->Get(Form("D_%d/D_%d/sigma2",xbin,xpt));
	TVectorD *vecDV = f->Get(Form("D_%d/D_%d/V",xbin,xpt));
	TVectorD *vecDavgmult = f->Get(Form("D_%d/avgmult",xbin));

	double *r = vecDr->GetMatrixArray();
	double *g2 = vecDg2->GetMatrixArray();
	double *sigma2 = vecDsigma2->GetMatrixArray();
	double *V = vecDV->GetMatrixArray();
	double *avgmult = vecDavgmult->GetMatrixArray();

	double r_[nstepr], g2_[nstepr];
	double temp;
	for(int ir=0;ir<nstepr;ir++){
		outG1>>temp;
		g2_[ir]=temp;
		outR>>temp;
		r_[ir]=temp;
	}
	for(int ir=0;ir<10;ir++){
	cout<<r[ir]<<"\t"<<r_[ir]<<"\t\t"<<g2[ir]<<"\t"<<g2_[ir]<<endl;
	}
	TCanvas *c1 = new TCanvas;
	c1->SetLogy();
	TGraph *gr1=new TGraph(nstepr,r_,g2_);
	TGraph *gr2=new TGraph(nstepr,r,g2);
	gr1->SetMarkerSize(1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(1);
	gr1->GetYaxis()->SetRangeUser(1e-7,2);
	gr2->SetLineColor(2);
	gr1->Draw("AP");
	gr2->Draw("Csame");
}
