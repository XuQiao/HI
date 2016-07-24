#include "fstream.h"
#include <vector>
#include "/home/xuq7/HI/utilities/sort.h"

void findcent()
{
	const int color[10] = {2,3,4,5,6,7,8,4,3,2};
	//const double centbin[] = {0,0.05,0.1,0.2,0.4,0.6,1.0};
	//const double centbin[] = {0,0.001,0.005,0.01,0.05,0.1,0.2,0.4,0.6,1.0};
        double centbin[101];
        for(int i=0;i<101;i++){
            centbin[i]=i*0.01;
        }
	const double erf = 1e-5;
	int ncent = (int)sizeof(centbin)/sizeof(double);
	//TFile *f = TFile::Open("merged_AnapAlfvtxor.root");
	//TFile *f = TFile::Open("merged_AnapAlfvtxand.root");
	//TFile *f = TFile::Open("merged_AnapAlfvtxsouth.root");
	TFile *f = TFile::Open("merged_AnapAlmb.root");
//	TH2F* pc1hitsbbcdis = (TH2F*)f->Get("pc1hitsbbcdis"); 
//	TH1F* pc1hitsdis = (TH1F*)f->Get("pc1hitsdis"); 
	TH2F* hbbcntrk = (TH2F*)f->Get("hbbcsbbcn");
//	hbbcntrk->RebinX(3);
//	TH1F* bbcsouthdis = (TH1F*)pc1hitsbbcdis->ProjectionX(0,-1);
	TH1F* bbcsouthdis = (TH1F*)hbbcntrk->ProjectionX(0,-1);
//	TH2F* bbcsouthdis = (TH2F*)hbbcntrk->Clone("bbcsouthdis");
	bbcsouthdis->SetTitle("BBC charge distribution");	
	bbcsouthdis->GetXaxis()->SetTitle("BBC charge");	
	bbcsouthdis->GetYaxis()->SetTitle("# of Events");	
//	pc1hitsdis->SetTitle("PC1 hits distribution");	
//	pc1hitsdis->GetXaxis()->SetTitle("PC1 hits");	
//	pc1hitsdis->GetYaxis()->SetTitle("# of Events");	


	ofstream fstr("centbin_bbcmb.txt");
//	fstr<<"centbin"<<"\t"<<"bbcsouth"<<"\t"<<"pc1hits"<<"\t"<<"bbc ratio"<<"\t"<<"pc1hits ratio"<<endl;
	fstr<<"centbin"<<"\t"<<"bbc"<<"\t"<<"bbc ratio"<<endl;
	
	vector<double> bbcs;
//	vector<double> pc1hs;

	for(int icent = 0; icent < ncent; icent++){
		double bbcskp = findpoint(bbcsouthdis,centbin[icent]);
//		double pc1hskp = findpoint(pc1hitsdis,centbin[icent]);
		if(bbcskp<0) bbcskp = 0;
//		if(pc1hskp<0) pc1hskp = 0;
		bbcs.push_back(bbcskp);
cout<<bbcskp<<",";
//		pc1hs.push_back(pc1hskp);
	}
	
	for(int icent = 0; icent < ncent; icent++){
	//	fstr<<centbin[icent]<<"\t"<<bbcs[icent]<<"\t"<<pc1hs[icent]<<"\t"<<bbcsouthdis->Integral(bbcsouthdis->FindBin(bbcs[icent]+erf),bbcsouthdis->FindBin(bbcs[0]-erf))/bbcsouthdis->Integral(bbcsouthdis->FindBin(bbcs[ncent-1]-erf),bbcsouthdis->FindBin(bbcs[0]+erf))<<"\t"<<pc1hitsdis->Integral(pc1hitsdis->FindBin(pc1hs[icent]+erf),pc1hitsdis->FindBin(pc1hs[0]-erf))/pc1hitsdis->Integral(pc1hitsdis->FindBin(pc1hs[ncent-1]-erf),pc1hitsdis->FindBin(pc1hs[0]+erf))<<endl;
		fstr<<centbin[icent]<<"\t"<<bbcs[icent]<<"\t"<<bbcsouthdis->Integral(bbcsouthdis->FindBin(bbcs[icent]+erf),bbcsouthdis->FindBin(bbcs[0]-erf))/bbcsouthdis->Integral(bbcsouthdis->FindBin(bbcs[ncent-1]-erf),bbcsouthdis->FindBin(bbcs[0]+erf))<<endl;
		//fstr<<centbin[icent]<<"\t"<<bbcs[icent]<<"\t"<<bbcsouthdis->Integral(bbcsouthdis->GetXaxis()->FindBin(bbcs[icent]+erf),bbcsouthdis->GetXaxis()->FindBin(bbcs[0]-erf),bbcsouthdis->GetYaxis()->FindBin(bbcs[icent]+erf),bbcsouthdis->GetYaxis()->FindBin(bbcs[0]-erf))/bbcsouthdis->Integral(bbcsouthdis->FindBin(bbcs[ncent-1]-erf,bbcs[ncent-1]-erf),bbcsouthdis->FindBin(bbcs[0]+erf,bbcs[0]+erf))<<endl;
	}

	TCanvas *c1 = new TCanvas();
	//c1->SetLogx();
	c1->SetLogy();
//	c1->SetLogz();
	bbcsouthdis->GetXaxis()->SetRangeUser(0,500);
//	bbcsouthdis->GetYaxis()->SetRangeUser(0,100);
	bbcsouthdis->Draw();
	for(int icent = 0; icent < ncent; icent++){
		TH1F* bbcs_t = (TH1F*)bbcsouthdis->Clone(Form("bbcs_%d",icent));
		bbcs_t->GetXaxis()->SetRangeUser(bbcs[icent+1]-erf,bbcs[icent]+erf);
		//bbcs_t->GetYaxis()->SetRangeUser(bbcs[icent+1]-erf,bbcs[icent]+erf);
		bbcs_t->SetFillColor(color[icent]);
		bbcs_t->Draw("HIST same");
	}
	c1->Print("bbcaddcent.png");
	/*
 	TCanvas *c2 = new TCanvas();
	c2->SetLogy();
	pc1hitsdis->Draw();
	for(int icent = 0; icent < ncent; icent++){
		TH1F* pc1hs_t = (TH1F*)pc1hitsdis->Clone(Form("pc1hs_%d",icent));
		pc1hs_t->GetXaxis()->SetRangeUser(pc1hs[icent+1]-erf,pc1hs[icent]+erf);
		pc1hs_t->SetFillColor(color[icent]);
		pc1hs_t->Draw("HIST same");
	}
	c2->Print("pc1hits.png");
*/
}
