#include "parameter.h"
void Npartdistr(){
	gStyle->SetOptStat(kFALSE);
	int color1[3]={1,4,2};
	TFile *fout = new TFile("GlauNpartdis.root","Recreate");
	TH1D* hNpart[3];
	TFile *fGlauber;
	TTree *t;
	TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
        c1->SetLogy();
	for(int i=0;i<3;i++){
	if(i==0)    fGlauber = TFile::Open(stdGlaulist[0]);    
 	else if(i==1)  fGlauber = TFile::Open(Gri055Glaulist[0]);    
	else    fGlauber = TFile::Open(Gri101Glaulist[0]);    
	hNpart[i] = new TH1D("","",100,0,100);
	hNpart[i]->Sumw2();
	t = (TTree*) fGlauber ->Get("nt_p_Pb");
	Float_t Npart;
	t->SetBranchAddress("Npart",&Npart);
	for(int ievt=0;ievt<t->GetEntries();ievt++){
		t->GetEntry(ievt);
		hNpart[i]->Fill(Npart);
		}
	hNpart[i]->Scale(1.0/hNpart[i]->Integral());
	fout->cd();
	if(i==0)	hNpart[i]->Write("stdNpart");
	else if(i==1)	hNpart[i]->Write("Gri055Npart");
	else	hNpart[i]->Write("Gri101Npart");
	c1->cd();
        hNpart[i]->GetXaxis()->SetTitle("N_{part}");
        hNpart[i]->GetYaxis()->SetTitle("Event Fraction P(N_{part})");
        hNpart[i]->GetXaxis()->SetRangeUser(0,60);
        hNpart[i]->SetLineColor(color1[i]);
        hNpart[i]->SetMarkerSize(1.0);
        hNpart[i]->SetMarkerStyle(20);
        hNpart[i]->SetMarkerColor(color1[i]);
	if(i==0)
        hNpart[i]->Draw();
	else
        hNpart[i]->Draw("Psame");
        }
        TLegend *leg0 = new TLegend(0.60,0.68,0.90,0.88);
	leg0->SetFillColor(10);
	leg0->SetBorderSize(0.035);
    	leg0->SetTextFont(42);
    	leg0->SetTextSize(0.030);
    	leg0->AddEntry(hNpart[0],"standard","p");
    	leg0->AddEntry(hNpart[1],"Gribov #Omega=0.55","p");
    	leg0->AddEntry(hNpart[2],"Gribov #Omega=1.01","p");
    	leg0->Draw("same");
	fout->Close();
}
