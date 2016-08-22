#include "/phenix/u/xuq/util/SimplifyLife.C"
void Drawnvtxntrk(){
	TFile *f = new TFile("output.root");
	TCanvas *c1;
	TLegend *leg = new TLegend(0.5,0.7,0.7,0.85);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.048);
	TH2F* nbbcntrk = (TH2F*)f->Get(Form("hnvtxntrk_0"));
	c1= new TCanvas();
	c1->cd();
	//c1->SetLogz();
	//SetTitle(*nbbcntrk,"nvtx hits","ntracks","");
	//SetStyle(*nbbcntrk,0,0,0,0,0);
	SetRange(*nbbcntrk,0,0,100,100);
	TProfile *hProf = (TProfile*)nbbcntrk->ProfileX();
	SetTitle(*hProf,"nvtx hits","ntracks","");
//	nbbcntrk->Draw("colz");	
	hProf->SetMarkerSize(1.4);
	hProf->Draw("P");

	c1->Print(Form("nvtxntrkProfX.png"));
}
