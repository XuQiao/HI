#include "fstream.h"
#include <vector>
#include "/phenix/plhf/xuq/phenix/flow/PP/sort.h"

void findpt()
{
	const int color[] = {2,3,4,5,6,7,8,4,3,2,6,7,8};
//	const double ptbin[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0};
	//const double ptbin[] = {0.2,1.0,2.0,3.0,5.0};
	const double ptbin[] = {2.5,4.0};
	const double erf = 1e-3;
	int npt = (int)sizeof(ptbin)/sizeof(double)-1;
	TFile *f = TFile::Open("merged.root");
	TH1F* ptdis = (TH1F*)f->Get(Form("ptforedis_%d",0));
	for(int xcent=1; xcent < 6; xcent++){
		TH1F* ptdis_t = (TH1F*)f->Get(Form("ptforedis_%d",xcent));
		ptdis->Add(ptdis_t);
	}
	ptdis->SetTitle("track pt distribution");
	ptdis->GetXaxis()->SetTitle("p_{T}^{track}");	
	ptdis->GetYaxis()->SetTitle("# of tracks");	

	ofstream fstr("ptbin.txt");
	fstr<<"ptbin"<<"\t"<<"pt ratio"<<"\t"<<"pt mean"<<endl;
	
	for(int ipt = 0; ipt < npt; ipt++){
		TH1F* pt_t = (TH1F*)ptdis->Clone(Form("pt_%d",ipt));
		pt_t->GetXaxis()->SetRangeUser(ptbin[ipt]-erf,ptbin[ipt+1]+erf);
		double ptmean = pt_t->GetMean(1);
		fstr<<ptbin[ipt]<<"\t"<<ptdis->Integral(ptdis->FindBin(ptbin[0]+erf),ptdis->FindBin(ptbin[ipt]-erf))/ptdis->Integral(ptdis->FindBin(ptbin[0]-erf),ptdis->FindBin(ptbin[npt-1]+erf))<<"\t"<<ptmean<<endl;
	cout<<ptmean<<", ";
	}
cout<<endl;

	TCanvas *c1 = new TCanvas();
	c1->SetLogy();
	ptdis->Draw();
	ptdis->GetXaxis()->SetRangeUser(0,6);
	for(int ipt = 0; ipt < npt; ipt++){
		TH1F* pt_t = (TH1F*)ptdis->Clone(Form("pt_%d",ipt));
		pt_t->GetXaxis()->SetRangeUser(ptbin[ipt]-erf,ptbin[ipt+1]+erf);
		pt_t->SetFillColor(color[ipt]);
		pt_t->Draw("HIST same");
	}
	c1->Print("ptdis.png");
	
}
