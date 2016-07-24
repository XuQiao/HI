#include "../parameter.h"
#include "par.h"
void DrawTHStack(){
	int color1[3]={1,4,2};
	gStyle->SetOptStat(kFALSE);
	int sth=0; int Gth=0;
	TFile *f = TFile::Open("G.root");
	THStack *hs[3];
	TH1D* h[3];
	TCanvas *c1 = new TCanvas("c1","stacked hists",10,10,700,900);
	c1->SetLogy();
	TString dirname;
	TString name;	
       	if(Gth==0)
        	name = "G0";
        else if(Gth<nGlau && Gth>0)
        	name = Form("Glau_%d",Gth);
        else
        	name = Form("bin_%d",Gth-nGlau+1);

	for(int sth=0; sth<3; sth++){
        	if(sth==0){dirname = "std";}
        	else if(sth==1){dirname ="Gri055";}
        	else {dirname ="Gri101";}
		hs[sth] = (THStack*)f->Get(Form("%s/%s/Npartdis",dirname.Data(),name.Data()));
		h[sth] = (TH1D*)hs[sth]->GetStack()->At(N-1);
		h[sth]->GetXaxis()->SetRangeUser(0,60);
		h[sth]->GetXaxis()->SetTitle("N_{part}");
		h[sth]->GetYaxis()->SetTitle("Event Fraction P(N_{part})");
		h[sth]->SetLineColor(color1[sth]);
		h[sth]->SetMarkerSize(1.0);
		h[sth]->SetMarkerStyle(20);
		h[sth]->SetMarkerColor(color1[sth]);
		h[sth]->Draw("Psame");
	}
	TLegend *leg0 = new TLegend(0.60,0.68,0.90,0.88);
    leg0->SetFillColor(10);
    leg0->SetBorderSize(0.035);
    leg0->SetTextFont(42);
    leg0->SetTextSize(0.030);
    leg0->AddEntry(h[0],"standard","p");
    leg0->AddEntry(h[1],"Gribov #Omega=0.55","p");
    leg0->AddEntry(h[2],"Gribov #Omega=1.01","p");
        leg0->Draw();
c1->Print("Npartdis_Gri.gif");
}
	
