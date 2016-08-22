#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
TString coll = "PPb";
void OldvsNewCentpthat(){
	gStyle->SetOptStat(kFALSE);
	int pthat = 0;
	int hibinmin = 0;
	int hibinmax = 10;
/*	h1_inc = DrawTH1(0,"hiHFEta4",pthat,0,100,1,30);
	h1 = DrawTH1(0,"hiHFEta4",pthat,hibinmin,hibinmax,1,20);
	h1_old = DrawTH1(1,"hiHFEta4",pthat,hibinmin,hibinmax,2,20);
	h2 = DrawTH1(0,"hiHFplusEta4",pthat,hibinmin,hibinmax,1,24);
	h2_inc_old = DrawTH1(1,"hiHFplusEta4",pthat,0,100,2,29);
	h2_old = DrawTH1(1,"hiHFplusEta4",pthat,hibinmin,hibinmax,2,24);
*/	
	h1_data = DrawTH1(2,"hiHFEta4",pthat,hibinmin,hibinmax,1,20);
	h1_inc_data = DrawTH1(2,"hiHFEta4",pthat,0,100,1,30);
	h2_data = DrawTH1(2,"hiHFplusEta4",pthat,hibinmin,hibinmax,2,24);
	h2_inc_data = DrawTH1(2,"hiHFplusEta4",pthat,0,100,2,29);
/*	h1->SetTitle("");
	h1->GetXaxis()->SetTitleSize(0.05);
	h1->GetYaxis()->SetTitleSize(0.05);
	h1->GetXaxis()->SetTitle("HF Energy");
	h1->GetYaxis()->SetTitle("# of Events");
	h1->GetXaxis()->SetRangeUser(0,100);
	h1->GetYaxis()->SetRangeUser(1,1e7);*/
        h1_data->SetTitle("");
        h1_data->GetXaxis()->SetTitleSize(0.05);
        h1_data->GetYaxis()->SetTitleSize(0.05);
        h1_data->GetXaxis()->SetTitle("HF Energy");
        h1_data->GetYaxis()->SetTitle("# of Events");
        h1_data->GetXaxis()->SetRangeUser(0,100);
        h1_data->GetYaxis()->SetRangeUser(1,1e9);

	TCanvas *c1 = new TCanvas();
	c1->cd();
	c1->SetLogy();
	//h1->Draw();
	//h1_inc->Draw("same");
	//h2_inc_old->Draw("same");
	//h1_old->Draw("same");
	//h2_old->Draw("same");
	//h2->Draw("same");
	h1_data->Draw();
	h1_inc_data->Draw("same");
	h2_data->Draw("same");
	h2_inc_data->Draw("same");
	TLegend *leg = new TLegend(0.45,0.70,0.80,0.9);
//	leg->AddEntry(h1,Form("HF |#eta| > 4, %d%% < centrality < %d%%, New",hibinmin,hibinmax),"p");
//	leg->AddEntry(h1_inc,Form("HF |#eta| > 4, inclusive, New"),"p");
//	leg->AddEntry(h1_old,Form("HF |#eta| > 4, %d%% < centrality < %d%%, Old",hibinmin,hibinmax),"p");
//	leg->AddEntry(h2,Form("HF #eta > 4, %d%% < centrality < %d%%, New",hibinmin,hibinmax),"p");
//	leg->AddEntry(h2_inc_old,Form("HF #eta > 4, inclusive, Old"),"p");
//	leg->AddEntry(h2_old,Form("HF #eta > 4, %d%% < centrality < %d%%, Old",hibinmin,hibinmax),"p");
	leg->AddEntry(h1_data,Form("HF |#eta| > 4, %d%% < centrality < %d%%, Data",hibinmin,hibinmax),"p");
	leg->AddEntry(h1_inc_data,Form("HF |#eta| > 4, inclusive, Data",hibinmin,hibinmax),"p");
	leg->AddEntry(h2_data,Form("HF #eta > 4, %d%% < centrality < %d%%, Data",hibinmin,hibinmax),"p");
	leg->AddEntry(h2_inc_data,Form("HF #eta > 4, inclusive, Data",hibinmin,hibinmax),"p");
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetTextSize(0.035);
	leg->Draw("same");
	TLatex T;
	T.SetNDC();
	T.SetTextAlign(12);
	T.SetTextSize(0.05);
	T.SetTextColor(1);
	T.SetTextFont(42);
/*	if(pthat!=0)
		T.DrawLatex(0.35,0.85,Form("pthat = %d, %d%% < centrality < %d%%",pthat,hibinmin,hibinmax));
	else
		T.DrawLatex(0.35,0.85,Form("merged, %d%% < centrality < %d%%",hibinmin,hibinmax));
*/
/*	TH1D* ratio = (TH1D*)h1->Clone("ratio");
	ratio->Divide(h1_old);
	TLegend *leg = new TLegend(0.65,0.6,0.8,0.8);
	leg->AddEntry(h1,Form("pthat = %d",pthat),"p");
	TLine *l = new TLine(0,1,100,1);
	l->SetLineStyle(2);
	ratio->Draw();
	leg->Draw("same");
	l->Draw("same");
	c1->Print(Form("pthat%dratio.png",pthat));
*/
	c1->Print(Form("Data%sDiffvar%d_%d.png",coll.Data(),hibinmin,hibinmax));
}

TH1D* DrawTH1(int isOld, TString var, int pthat, int hibinmin, int hibinmax, int color, int marker){
	if(pthat!=0){
		if(isOld==0)
			TFile *f = TFile::Open(Form("/cms/store/user/ymao/pA5TEV/Mixing/STARTHI53V27/merged/QCD_Pt_%d_TuneZ2_Mix_5p02TeV_pythia.root",pthat));
		else if(isOld==1)
			TFile *f = TFile::Open(Form("/cms/store/user/qixu/jetRpA/OldForest/pPb_Full_BTagForest%d_Fix3output.root",pthat));
		TTree *t = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
	}
	else{
		if(isOld==0)
			TFile *f = TFile::Open(Form("/store/user/qixu/jetRpA/skimTree/MC%sakPu3PFskimfile0_10.root",coll.Data()));
		else if(isOld==1)
			TFile *f = TFile::Open(Form("/store/user/qixu/jetRpA/skimTree/MC%sakPu3PFskimfile0_8.root",coll.Data()));
		else if(isOld==2)
			TFile *f = TFile::Open(Form("/store/user/qixu/jetRpA/skimTree/DATA%sakPu3PFskimfile.root",coll.Data()));
		TTree *t = (TTree*)f->Get("nt");
		
	}
	Float_t hiHFplusEta4, hiHFminusEta4, vz;
	Int_t hiBin;
	TH1D* h = new TH1D("h","h",100,0,100);
	h->Sumw2();
	t->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
	t->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
	t->SetBranchAddress("vz",&vz);
	t->SetBranchAddress("hiBin",&hiBin);
	cout<<"Total Events="<<t->GetEntries()<<endl;
	for(int Ev=0;Ev<t->GetEntries();Ev++){
		if(Ev%1000000==0) cout<<"has processed "<<Ev<<endl;
		t->GetEntry(Ev);
		if(hiBin > hibinmax || hiBin < hibinmin) continue;
		if(TMath::Abs(vz)>=15) continue;
		if(var == "hiBin")
		h->Fill(hiBin);
		else if(var == "hiHFplusEta4")
		h->Fill(hiHFplusEta4);
		else if(var == "hiHFEta4")
		h->Fill(hiHFplusEta4+hiHFminusEta4);
	}
//	t->Draw(Form("%s>>h",var.Data()),"TMath::Abs(vz)<15");
//	h = (TH1D*)gDirectory->Get("h");
//	h->Scale(1.0/h->Integral());
	h->SetMarkerColor(color);
	h->SetMarkerStyle(marker);
	h->SetLineColor(color);
	h->SetMarkerSize(1.3);
	return h;
}
