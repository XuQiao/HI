#include "/home/xuq7/HI/jetRpA/RpA/Quality/root_setting.h"
#include "file.h"
TH1* myRebin(TH1* histo,int Nbin, double binbound[]){
//rehisto->Scale(1/rehisto->Integral());
TH1* rehisto=histo->Clone("rehisto");
TH1* rehisto=(TH1*)histo->Rebin(Nbin,histo->GetName(),binbound);
fixedFontHist(rehisto,1.2,1.3);
return rehisto;
}

void Drawjtptvsjtpu(){
	gStyle->SetOptStat(kFALSE);

	TH2F* jtptvsjtpu = new TH2F("jtptvsjtpu","",1000,0,1000,1000,0,1000);
	TH1F* hFrame=new TH1F("","",1000,0,1000);
	hFrame->GetXaxis()->SetLimits(30.1,599.9);
	hFrame->SetTitle("");
	hFrame->GetYaxis()->SetLabelSize(20);
	hFrame->SetMaximum(100);
	hFrame->SetMinimum(0);
	TFile *f = TFile::Open("/cms/store/user/qixu/jetRpA/skimTree/MCPPbakPu3PFskimfile0_9.root");
	TTree *nt = (TTree*)f->Get("nt");
	Float_t jtpt[100],jteta[100],jtpu[100],rawpt[100];
	Int_t nref;
  Float_t weight;
  Float_t vz;

  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("nref",&nref);
  nt->SetBranchAddress("jtpt",jtpt);
  nt->SetBranchAddress("jtpu",jtpu);
  nt->SetBranchAddress("jteta",jteta);
  nt->SetBranchAddress("rawpt",rawpt);
  nt->SetBranchAddress("weight",&weight);

  Int_t nentries = nt->GetEntries();

  for(int i=0; i<nentries; i++){
	 nt->GetEntry(i);
	if(TMath::Abs(vz)>15) continue;
    if(i%100000==1)    cout<<"analyzing "<< i <<" th event"<<endl;  
	for(int j4i = 0; j4i < nref; j4i++){
    double jet_pt = jtpt[j4i];
    double jet_eta = jteta[j4i];
    double raw_pt = rawpt[j4i];
    double jet_pu = jtpu[j4i];
        if(raw_pt<20.0 || fabs(jet_eta)>5.0) continue;
	if(TMath::Abs(jet_eta+0.465)<=1.) {
		jtptvsjtpu->Fill(jet_pt,jet_pu,weight);
	}
	}
	}
	
	c1 = new TCanvas("c1"," ",600,600);
//	makeMultiPanelCanvas(c1,2,1,0,0,0.1,0.15,0.03);		
//	c1->Divide(1,1,0.01,0.01);
	c1->SetRightMargin(0.15);
//	t->Draw("jtpt:jtpu>>jtptvsjtpu","weight*(rawpt>=20 && jteta<=1-0.465 && jteta>=-1-0.465)");
	jtptvsjtpu->SetTitle("");
	jtptvsjtpu->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
	jtptvsjtpu->GetXaxis()->SetRangeUser(0,130);
	jtptvsjtpu->GetYaxis()->SetTitle("jtpu (GeV/c)");
	jtptvsjtpu->GetYaxis()->SetRangeUser(0,25);
	//TProfile *hProf = (TProfile*)jtptvsjtpu->ProfileX("hProfX",1,-1);
	//TProfile *rehisto_hProf=(TProfile*)myRebin(hProf,Nbin_pt,binbound_pt);
	hFrame->GetXaxis()->SetTitle("jet p_{T} (GeV/c)");
	hFrame->GetYaxis()->SetTitle("average jtpu (GeV/c)");
	c1->cd(1)->SetLogz();
	jtptvsjtpu->Draw("colz");
//	c1->cd(2);
//	hFrame->DrawCopy();
//	rehisto_hProf->Draw("same");
	c1->Print("jtptvsjtpu.png");
}
	
