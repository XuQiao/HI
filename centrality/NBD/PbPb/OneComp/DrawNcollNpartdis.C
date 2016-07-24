#include "parameter.h"
void DrawNcollNpartdis(){
        int Gth=0;
	TFile *fGlauber = TFile::Open(stdGlaulist[Gth]);
	TTree *t = (TTree*)fGlauber->Get("nt_Pb_Pb");
        Float_t Ncoll, Npart, B;	Long_t Nevent;
	t->SetBranchAddress("Ncoll",&Ncoll);
	t->SetBranchAddress("Npart",&Npart);
	t->SetBranchAddress("B",&B);
        Nevent = (Long_t) t->GetEntries();
    	Long_t Ev;	Int_t Bino;	Double_t Para, Bi_Para;
        TH1D* hNcoll = new TH1D("","",4000,0,4000);
        TH1D* hNpart = new TH1D("","",4000,0,4000);
	for (Ev=0; Ev<Nevent; Ev++){
	//if(Ev%100000==0)	cout<<"\t"<<"Have run "<<Ev<<" events"<<endl;
        	t->GetEntry(Ev);
                hNcoll->Fill(Ncoll);
                hNpart->Fill(Npart);
        }
	TString name;
        if(Gth==0)
                name = "G0";
        else if(Gth<nGlau)
                name = Form("Glau_%d",Gth);
        else
                name = Form("bin_%d",Gth-nGlau+1);
        TCanvas *c1 = new TCanvas();
        c1->SetLogy();
        hNcoll->Scale(1./hNcoll->Integral());
        hNcoll->GetXaxis()->SetRangeUser(0,3000);
        hNcoll->GetXaxis()->SetTitle("Ncoll");
        hNcoll->GetYaxis()->SetTitle("Event Fraction");
        hNcoll->Draw();
        c1->Print(Form("%s_Ncoll.png",name.Data()));
        TCanvas *c2 = new TCanvas();
        c2->SetLogy();
        hNpart->Scale(1./hNpart->Integral());
        hNpart->GetXaxis()->SetRangeUser(0,500);
        hNpart->GetXaxis()->SetTitle("Npart");
        hNpart->GetYaxis()->SetTitle("Event Fraction");
        hNpart->Draw();
        c2->Print(Form("%s_Npart.png",name.Data()));
}
