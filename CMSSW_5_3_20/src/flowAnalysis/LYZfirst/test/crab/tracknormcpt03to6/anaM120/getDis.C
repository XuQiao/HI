void getDis(){

	system("hadd merged.root /lio/lfs/cms/store/user/qixu/flow/PACorrHM/tracktcpt03to6/multiM120/*/pPbAnaTree*.root");
	TFile *fmerged = TFile::Open("merged.root");
	TH1D* hNtrkvsMult = (TH1D*)fmerged->Get("demo/hNtrkvsMult");
	TH1D* hpt = (TH1D*)fmerged->Get("demo/pt");
	TH1D* heta = (TH1D*)fmerged->Get("demo/eta");
	TH1D* hphi = (TH1D*)fmerged->Get("demo/phi");
	TFile *fout = new TFile("merged_Dis.root","Recreate");
	fout->cd();
	hNtrkvsMult->Write();	
	hpt->Write();
	heta->Write();	
	hphi->Write();	
	remove("merged.root");
}
