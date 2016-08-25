#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void pthatCentOrder(){
	gStyle->SetOptStat(kFALSE);
	h1 = DrawTH1("hiBin",15);
	h1->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitle("Centrality");
	h1->GetYaxis()->SetTitle("Event Fraction");
	h1->GetXaxis()->SetRangeUser(0,100);
	h1->GetYaxis()->SetRangeUser(0,0.07);

	h2 = DrawTH1("hiBin",120);
	h3 = DrawTH1("hiBin",460);
	h1->SetMarkerStyle(20);
	h1->SetMarkerColor(1);
	h1->SetLineColor(1);
	h2->SetMarkerStyle(29);
	h2->SetMarkerColor(2);
	h2->SetLineColor(2);
	h3->SetMarkerStyle(24);
	h3->SetMarkerColor(4);
	h3->SetLineColor(4);
	TCanvas *c1 = new TCanvas();
	c1->cd();
	h1->Draw("P");
	h2->Draw("Psame");
	h3->Draw("Psame");
	TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
	leg->AddEntry(h1,"pthat = 15","p");
	leg->AddEntry(h2,"pthat = 120","p");
	leg->AddEntry(h3,"pthat = 460","p");
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.04);
	leg->Draw("same");
	c1->Print("hiBinpthatOrder.png");

	h1_ = DrawTH1("hiHFplusEta4",15);
	h1_->SetTitle("");
        h1_->GetXaxis()->SetTitleSize(0.05);
        h1_->GetYaxis()->SetTitleSize(0.05);
        h1_->GetXaxis()->SetTitle("Centrality");
        h1_->GetYaxis()->SetTitle("Event Fraction");
        h1_->GetXaxis()->SetRangeUser(0,100);
	h2_ = DrawTH1("hiHFplusEta4",120);
	h3_ = DrawTH1("hiHFplusEta4",460);

	TCanvas *c2 = new TCanvas();
	c2->cd();
	c2->SetLogy();
	h1_->GetXaxis()->SetTitle("HF Energy #eta>4");
	h1_->GetYaxis()->SetTitle("Event Fraction");
	h1_->GetYaxis()->SetRangeUser(1e-6,1e-1);
	h1_->SetMarkerStyle(20);
        h1_->SetMarkerColor(1);
        h1_->SetLineColor(1);
        h2_->SetMarkerStyle(29);
        h2_->SetMarkerColor(2);
        h2_->SetLineColor(2);
        h3_->SetMarkerStyle(24);
        h3_->SetMarkerColor(4);
        h3_->SetLineColor(4);
        h1_->Draw("P");
        h2_->Draw("Psame");
        h3_->Draw("Psame");
        leg->Draw("same");
	c2->Print("hiHFpthatOrder.png");
}

TH1D* DrawTH1(TString var, int pthat){
	TFile *f = TFile::Open(Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/QCD_Pt_%d_TuneZ2_Mix_5p02TeV_pythia.root",pthat));
	TTree *t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
	Float_t hiHFplusEta4, vz;
	Int_t hiBin;
	TH1D* h = new TH1D("h","h",100,0,100);
	t->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
	t->SetBranchAddress("vz",&vz);
	t->SetBranchAddress("hiBin",&hiBin);
	for(int Ev=0;Ev<t->GetEntries();Ev++){
		t->GetEntry(Ev);
		if(TMath::Abs(vz)>=15) continue;
		if(var == "hiBin")
		h->Fill(hiBin);
		else if(var == "hiHFplusEta4")
		h->Fill(hiHFplusEta4);
	}
//	t->Draw(Form("%s>>h",var.Data()),"TMath::Abs(vz)<15");
//	h = (TH1D*)gDirectory->Get("h");
	h->Scale(1.0/h->Integral());
	return h;
}
