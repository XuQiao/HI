void Draw(){
	gStyle->SetOptStat(kFALSE);
	TFile *f = TFile::Open("histo.root");
	TH1D* hpt = (TH1D*)f->Get("hpt");
	TH1D* hpfpt_id1 = (TH1D*)f->Get("hpfpt_id1");
	TH1D* heta = (TH1D*)f->Get("heta");
	TH1D* hpfeta_id1 = (TH1D*)f->Get("hpfeta_id1");
	TH1D* hphi = (TH1D*)f->Get("hphi");
	TH1D* hpfphi = (TH1D*)f->Get("hpfphi_id1");
	cout<<"pt:"<<hpt->Integral()<<"\t"<<hpfpt_id1->Integral()<<endl;
	cout<<"eta:"<<heta->Integral()<<"\t"<<hpfeta_id1->Integral()<<endl;


	TCanvas *c1 = new TCanvas("c1","c1",1200,400);
	c1->Divide(3,1);
	c1->cd(1)->SetLogy();
	hpt->SetLineColor(1);
	hpt->SetMarkerColor(1);
	hpt->SetMarkerStyle(20);
	hpt->SetMarkerSize(1.2);
	hpt->Draw("P");
	hpfpt_id1->SetLineColor(2);
	hpfpt_id1->SetMarkerColor(2);
	hpfpt_id1->SetMarkerStyle(24);
	hpfpt_id1->SetMarkerSize(1.2);
	hpfpt_id1->Draw("Psame");

        c1->cd(2);
        heta->SetLineColor(1);
        heta->SetMarkerColor(1);
        heta->SetMarkerStyle(20);
        heta->SetMarkerSize(1.2);
        heta->Draw("P");
        hpfeta_id1->SetLineColor(2);
        hpfeta_id1->SetMarkerColor(2);
        hpfeta_id1->SetMarkerStyle(24);
        hpfeta_id1->SetMarkerSize(1.2);
        hpfeta_id1->Draw("Psame");

        c1->cd(3);
        hphi->SetLineColor(1);
        hphi->SetMarkerColor(1);
        hphi->SetMarkerStyle(20);
        hphi->SetMarkerSize(1.2);
        hphi->Draw("P");
        hpfphi_id1->SetLineColor(2);
        hpfphi_id1->SetMarkerColor(2);
        hpfphi_id1->SetMarkerStyle(24);
        hpfphi_id1->SetMarkerSize(1.2);
        hpfphi_id1->Draw("Psame");
	c1->SaveAs("plot.pdf");

	TCanvas *c2 = new TCanvas("c2","c2",1200,400);
	c2->Divide(3,1);
	c2->cd(1);
	TH1D* hr_pt = (TH1D*)hpt->Clone("ratio_pt");
	hr_pt->Divide(hpfpt_id1);
	hr_pt->SetMarkerSize(1);
	hr_pt->Draw("P");
	TLine *l1 = new TLine(0,1,20,1);
	l1->SetLineStyle(2);
	l1->Draw("same");
	c2->cd(2);
	TH1D* hr_eta = (TH1D*)heta->Clone("ratio_eta");
	hr_eta->Divide(hpfeta_id1);
	hr_eta->SetMarkerSize(1);
	hr_eta->Draw("P");
	TLine *l2 = new TLine(-5,1,5,1);
	l2->SetLineStyle(2);
	l2->Draw("same");
	c2->cd(3);
	TH1D* hr_phi = (TH1D*)hphi->Clone("ratio_phi");
	hr_phi->Divide(hpfphi_id1);
	hr_phi->SetMarkerSize(1);
	hr_phi->Draw("P");
	TLine *l3 = new TLine(-4,1,4,1);
	l3->SetLineStyle(2);
	l3->Draw("same");
}
