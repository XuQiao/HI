#include "../parameter.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include "par.h"
void DrawTHStack(){
        static const int nColor = 9;
        static const int colorCode[nColor] = {
        1,2,3,6,4,8,1,9,5};
	gStyle->SetOptStat(kFALSE);
	int sth=0; int Gth=0;
	TFile *f = TFile::Open(outG);
	TObjArray *hs[3];
        std::vector<TH1D*> h[3];
	TCanvas *c1 = new TCanvas("c1","stacked hists",10,10,600,600);
	c1->SetLogy();
	TString dirname;
	TString name;	
        TH1D* hFrame = new TH1D("","",4000,0,4000);
	hFrame->GetXaxis()->SetRangeUser(0,500);
	hFrame->GetYaxis()->SetRangeUser(1,1e5);
	hFrame->GetXaxis()->SetTitle("N_{part}");
	hFrame->GetYaxis()->SetTitle("# of Events");
        hFrame->DrawCopy();
       	if(Gth==0)
        	name = "G0";
        else if(Gth<nGlau && Gth>0)
        	name = Form("Glau_%d",Gth);
        else
        	name = Form("bin_%d",Gth-nGlau+1);
        TVectorD *centbin = (TVectorD*)f->Get(Form("std/G0/centbin"));
        TVectorD *kpoint = (TVectorD*)f->Get(Form("std/G0/kpoint"));

	for(int sth=0; sth<1; sth++){
        	if(sth==0){dirname = "std";}
        	else if(sth==1){dirname ="Gri055";}
        	else {dirname ="Gri101";}
		hs[sth] = (TObjArray*)f->Get(Form("%s/%s/Npartdis",dirname.Data(),name.Data()));
                h[sth].resize(N);
                for(int i=0;i<N;i++){
//		h[sth][i] = (TH1D*)hs[sth]->GetStack()->At(N-1);
		h[sth][i] = (TH1D*)hs[sth]->At(i);
		h[sth][i]->SetLineColor(colorCode[i]);
		h[sth][i]->SetMarkerSize(1.0);
		h[sth][i]->SetMarkerStyle(20);
		h[sth][i]->SetMarkerColor(colorCode[i]);
		h[sth][i]->Draw("Psame");
//                hs[sth]->Draw("Psame");
                }
	}
    std::vector<TString> label(N);
	TLegend *leg0 = new TLegend(0.60,0.60,0.90,0.88);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.030);
    for(int i=0;i<N-1;i++){
	if(method==0)label[i] = Form("%.2f-%.2f%%",(*centbin)[i]*100,(*centbin)[i+1]*100);
	else label[i] = Form("%.2f-%.2f",(*kpoint)[i],(*kpoint)[i+1]);
//    leg0->AddEntry(h[0][i],"standard","p");
    leg0->AddEntry(h[0][i],label[i],"p");
    }
 //   leg0->AddEntry(h[1],"Gribov #Omega=0.55","p");
 //   leg0->AddEntry(h[2],"Gribov #Omega=1.01","p");
        leg0->Draw();
c1->Print("Npartdis.pdf");
c1->Print("Npartdis.png");
}
	
