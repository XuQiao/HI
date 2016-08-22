void plotdis(){
	TFile *f = TFile::Open("../jobs/pPbAnaTree_merged.root");
	TH1D *hvz = (TH1D*)f->Get("demo/vz");
	TH1D *hpt = (TH1D*)f->Get("demo/pt");
	TH1D *hphi = (TH1D*)f->Get("demo/phi");
	
	TH1D *hNtrk = (TH1D*)f->Get("demo/hntrk");
	
	TCanvas *c1 = new TCanvas();
	hvz->Draw();
	c1->SaveAs("vz.gif");
	
	TCanvas *c2 = new TCanvas();
	c2->SetLogy();
	hpt->Draw();
	c2->SaveAs("pt.gif");	

	TCanvas *c3 = new TCanvas();
	hphi->Draw();
	c3->SaveAs("phi.gif");
	
	TCanvas *c4 = new TCanvas();
	c4->SetLogy();
	hNtrk->Draw();
	c4->SaveAs("hNtrk.gif");
}

	
