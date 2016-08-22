#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
void DrawQuality(){
	TString dir1 = "PP2010trackpt03to6";
	TString dir2 = "tracknormcpt03to6";
	TString mult = "M150120";
	TH1D* hpt_PFon = plot(dir1,"","pt",1,20);
	TH1D* heta_PFon = plot(dir1,"","eta",1,20);
	TH1D* hphi_PFon = plot(dir1,"","phi",1,20);
        TH1D* hpt_t03to6 = plot(dir2,mult,"pt",2,20);
        TH1D* heta_t03to6 = plot(dir2,mult,"eta",2,20);
        TH1D* hphi_t03to6 = plot(dir2,mult,"phi",2,20);

	c1 = new TCanvas("c1"," ",1200,400);
        c1->Divide(3,1);
	c1->cd(1)->SetLogy();
//	makeMultiPanelCanvas(c1,3,1,0,0,0.25,0.2,0.03);
        gStyle->SetOptFit(1);
        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);
        gStyle->SetErrorX(0);
	TH1D* hFrame_pt = new TH1D("hFrame_pt","",100,0,10);
        hFrame_pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hFrame_pt->GetXaxis()->SetRangeUser(0,10);
        hFrame_pt->GetYaxis()->SetRangeUser(1,1e9);
        hFrame_pt->GetYaxis()->SetTitle("# of particles");
        hFrame_pt->GetXaxis()->SetTitleOffset(1);
	hFrame_pt->Draw();
        hpt_PFon->Draw("Psame");
        hpt_t03to6->Draw("Psame");
	TH1D* ratio = (TH1D*)hpt_PFon->Clone();
	ratio->Divide(hpt_t03to6);
	TCanvas *c2 = new TCanvas();
	ratio->Draw();
        c1->cd(2)->SetLogy();
	TH1D* hFrame_eta = new TH1D("hFrame_eta","",100,-5,5);
        hFrame_eta->GetXaxis()->SetTitle("#eta");
        hFrame_eta->GetXaxis()->SetRangeUser(-5,5);
        hFrame_eta->GetYaxis()->SetRangeUser(1,5e7);
        hFrame_eta->GetYaxis()->SetTitle("# of particles");
        hFrame_eta->GetXaxis()->SetTitleOffset(1);
        hFrame_eta->Draw();
        heta_PFon->Draw("Psame");
        heta_t03to6->Draw("Psame");
        c1->cd(3)->SetLogy();
	TH1D* hFrame_phi = new TH1D("hFrame_phi","",100,-4,4);
	hFrame_phi->GetXaxis()->SetTitle("#phi");
        hFrame_phi->GetXaxis()->SetRangeUser(-4,4);
        hFrame_phi->GetYaxis()->SetRangeUser(1e1,1e7);
        hFrame_phi->GetYaxis()->SetTitle("# of particles");
        hFrame_phi->GetXaxis()->SetTitleOffset(1);
        hFrame_phi->Draw();
        hphi_PFon->Draw("Psame");
        hphi_t03to6->Draw("Psame");
	TLegend *tl = new TLegend(0.3,0.3,0.6,0.5);
	tl->SetTextSize(0.045);
	tl->SetBorderSize(0);
	tl->SetFillColor(0);
	//tl->AddEntry(hpt_PFon,"PFcandidate charged","p");
	tl->AddEntry(hpt_PFon,"PFcandidate","p");
	tl->AddEntry(hpt_t03to6,"track loose cut 0.1<p_{T}<10.0","p");
	tl->Draw();
	c1->SaveAs("PFQualify.png");

}
	
TH1D* plot(TString dir, TString Tmult, TString var,int color, int marker){
	TString pf;
	if(dir.Contains("PF"))	pf="pf";
	else pf="";
	int xbin=0;
	if(Tmult==""){
	TFile *fhisto = TFile::Open(Form("%s/histomerged.root",dir.Data()));
	}
	else
	TFile *fhisto = TFile::Open(Form("%s/%s/histomerged.root",dir.Data(),Tmult.Data()));
	TH1D* histo = (TH1D*)fhisto->Get(Form("D_%d/h%s%s",xbin,pf.Data(),var.Data()));
	int rebin=histo->GetNbinsX()/100;
	histo->Rebin(rebin);
        histo->SetMarkerColor(color);
        histo->SetLineColor(color);
        histo->SetMarkerStyle(marker);
        histo->SetMarkerSize(1);
	cout<<histo->GetEntries()<<endl;
	return histo;
}
