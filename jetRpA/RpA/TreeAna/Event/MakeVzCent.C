void MakeVzCent(){
TFile *f=new TFile("/scratch/xuq7/RpA/TreeAna/PbPskimfile0_6815.root");
const int Ntrig=5;
TString trigName[Ntrig]={"Jet20","Jet40","Jet60","Jet80","Jet100"};

TTree *nt=f->Get("nt");
TFile *fout=new TFile("PbPJetTrigvzandcent.root","RECREATE");

TH1F *vzhisto[Ntrig];	TH1F *centhisto[Ntrig];
TH1F *vzh;	TH1F *centh;

 for(int itrig=0; itrig<Ntrig; itrig++){
vzhisto[itrig]=new TH1F(Form("%svz",trigName[itrig].Data()),Form("%svz",trigName[itrig].Data()),600,-30,30); vzhisto[itrig]->Sumw2();
centhisto[itrig]=new TH1F(Form("%scent",trigName[itrig].Data()),Form("%scent",trigName[itrig].Data()),100,0,100); centhisto[itrig]->Sumw2();
}
vzh=new TH1F(Form("vz"),Form("vz"),600,-30,30); vzh->Sumw2();
centh=new TH1F(Form("cent"),Form("cent"),100,0,100); centh->Sumw2();

  Int_t jetTr2[Ntrig];
  Int_t HLT_PAJet20_noJetID_v1;
  Int_t HLT_PAJet40_noJetID_v1;
  Int_t HLT_PAJet60_noJetID_v1;
  Int_t HLT_PAJet80_noJetID_v1;
  Int_t HLT_PAJet100_noJetID_v1;
  Double_t vz;
  Int_t hiBin;
  Int_t pVertexFilterCutGplus;
  Int_t pPAcollisionEventSelectionPA;
  Int_t pBeamScrapingFilter;
  Int_t pprimaryvertexFilter;
  Int_t phfPosFilter1;
  Int_t phfNegFilter1;
  Int_t pHBHENoiseFilter;

  nt->SetBranchAddress("HLT_PAJet20_noJetID_v1",&HLT_PAJet20_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet40_noJetID_v1",&HLT_PAJet40_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet60_noJetID_v1",&HLT_PAJet60_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet80_noJetID_v1",&HLT_PAJet80_noJetID_v1);
  nt->SetBranchAddress("HLT_PAJet100_noJetID_v1",&HLT_PAJet100_noJetID_v1);
  nt->SetBranchAddress("vz",&vz);
  nt->SetBranchAddress("hiBin",&hiBin);
  nt->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  nt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  nt->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);

for(Long_t i; i<nt->GetEntries(); i++){
if(i%100000==0)
cout<<"has processed "<<i<<" th Events"<<endl; 
nt->GetEntry(i);
int  flag=0;
   if(TMath::Abs(vz)>15) continue;
   if(!pHBHENoiseFilter || !pprimaryvertexFilter || !pPAcollisionEventSelectionPA ) continue;
      jetTr2[0] = HLT_PAJet20_noJetID_v1;
      jetTr2[1] = HLT_PAJet40_noJetID_v1;
      jetTr2[2] = HLT_PAJet60_noJetID_v1;
      jetTr2[3] = HLT_PAJet80_noJetID_v1;
      jetTr2[4] = HLT_PAJet100_noJetID_v1;
   for(int itrig=0; itrig<Ntrig; itrig++){
	if(jetTr2[itrig]==1){
	vzhisto[itrig]->Fill(vz);
	centhisto[itrig]->Fill(hiBin);
	flag =1;
	}
	}
	
	if(flag==1){
	 vzh->Fill(vz);
 	centh->Fill(hiBin);}
		
}	//event loop

fout->cd();
for(int itrig=0; itrig<Ntrig; itrig++){
	vzhisto[itrig]->Write();
	centhisto[itrig]->Write();
}
vzh->Write();
centh->Write();
fout->Close();
}
