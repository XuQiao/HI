void NjetNeventratio(){
	const TString instring = "/store/user/tuos/mergedForestPythiaMiniBias.root";
	TFile *fin = TFile::Open(instring);
	TFile *fout = new TFile("PythiaMBhistos.root","Recreate");
	Float_t hiHFplusEta4, hiHFminusEta4, hiHFplus, hiHF;
	Float_t genpt[1000],geneta[1000];
	int ngen;
	//int Nevt, Njet;
	TH1D::SetDefaultSumw2();
	TH1D* hHFplus4Njet = new TH1D("hHFplus4Njet","",50,0,50);
	TH1D* hHFplus4Njet_ = new TH1D("hHFplus4Njet_","",50,0,50);
	TH1D* hHFplus4 = new TH1D("hHFplus4","",50,0,50);
        TH1D* hHFplusNjet = new TH1D("hHFplusNjet","",100,0,100);
        TH1D* hHFplusNjet_ = new TH1D("hHFplusNjet_","",100,0,100);
        TH1D* hHFplus = new TH1D("hHFplus","",100,0,100);
        TH1D* hHF4 = new TH1D("hHF4","",1000,0,100);
        TH1D* hHF = new TH1D("hHF","",1000,0,100);
	TTree *tree = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
	tree->AddFriend("akPu3PFJetAnalyzer/t");
	tree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
	tree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
	tree->SetBranchAddress("hiHFplus",&hiHFplus);
	tree->SetBranchAddress("hiHF",&hiHF);
	tree->SetBranchAddress("ngen",&ngen);
	tree->SetBranchAddress("genpt",genpt);
	tree->SetBranchAddress("geneta",geneta);
	for(int evi = 0; evi < tree->GetEntries(); evi++) {
		if(evi % 100000 == 0) cout<<"has processed "<<evi<<" events"<<endl;
	        tree->GetEntry(evi);
		hHFplus4->Fill(hiHFplusEta4);
		hHFplus->Fill(hiHFplus);
		hHF4->Fill(hiHFplusEta4+hiHFminusEta4);
		hHF->Fill(hiHF);
		for(int j4i = 0; j4i<ngen; j4i++){
			if(genpt[j4i]>20 && TMath::Abs(geneta[j4i])<2.8){
			hHFplus4Njet->Fill(hiHFplusEta4);
			hHFplusNjet->Fill(hiHFplus);
			}
		}
		if(ngen!=0){
			hHFplus4Njet_->Fill(hiHFplusEta4);
			hHFplusNjet_->Fill(hiHFplus);
		}
	}
	TH1D* hratioplus4 = (TH1D*)hHFplus4Njet->Clone("NjetratioPlus4");
	TH1D* hratioplus4_ = (TH1D*)hHFplus4Njet_->Clone("NjetratioPlus4_");
	hratioplus4->Divide(hHFplus4);
	hratioplus4_->Divide(hHFplus4);
	TH1D* hratioplus = (TH1D*)hHFplusNjet->Clone("NjetratioPlus");
        TH1D* hratioplus_ = (TH1D*)hHFplusNjet_->Clone("NjetratioPlus_");
        hratioplus->Divide(hHFplus);
        hratioplus_->Divide(hHFplus);
	fout->cd();
	hHFplus4Njet->Write();
	hHFplus4Njet_->Write();
	hHFplus4->Write();
	hratioplus4->Write();
	hratioplus4_->Write();
        hHFplusNjet->Write();
        hHFplusNjet_->Write();
        hHFplus->Write();
        hHF4->Write();
        hHF->Write();
	hratioplus->Write();
	hratioplus_->Write();
	fin->Close();
	fout->Close();
}
