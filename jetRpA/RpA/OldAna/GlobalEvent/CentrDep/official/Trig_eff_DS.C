void Calc(TString Generatorfile,TString outfile,bool Isweight){

TFile *fin = TFile::Open(Generatorfile);
TFile *fout = TFile::Open(outfile,"RECREATE"); 

TH1F *hHFEnergy4 = new TH1F("hHFEnergy4","",10000,0,100);	hHFEnergy4->Sumw2();
TH1F *hHFEnergy4_tr = new TH1F("hHFEnergy4_tr","",10000,0,100); hHFEnergy4_tr->Sumw2();
TH1F *hHFEnergy4_DS = new TH1F("hHFEnergy4_DS","",10000,0,100); hHFEnergy4_DS->Sumw2();
TH1F *hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","",10000,0,100);	hHFEnergyPlus4->Sumw2();
TH1F *hHFEnergyPlus4_tr = new TH1F("hHFEnergyPlus4_tr","",10000,0,100); 	hHFEnergyPlus4_tr->Sumw2();
TH1F *hHFEnergyPlus4_DS = new TH1F("hHFEnergyPlus4_DS","",10000,0,100); 	hHFEnergyPlus4_DS->Sumw2();
TH1F *hHFEnergyPlus = new TH1F("hHFEnergyPlus","",10000,0,100);	hHFEnergyPlus->Sumw2();
TH1F *hHFEnergyPlus_tr = new TH1F("hHFEnergyPlus_tr","",10000,0,100); 	hHFEnergyPlus_tr->Sumw2();
TH1F *hHFEnergyPlus_DS = new TH1F("hHFEnergyPlus_DS","",10000,0,100); 	hHFEnergyPlus_DS->Sumw2();

TH1F *hNtrk = new TH1F("hNtrk","",500,0,500);   hNtrk->Sumw2();
TH1F *hNtrk_tr = new TH1F("hNtrk_tr","",500,0,500); hNtrk_tr->Sumw2();
TH1F *hNtrk_DS = new TH1F("hNtrk_DS","",500,0,500); hNtrk_DS->Sumw2();

fin->cd();
TTree* HiTree=(TTree*)hiEvtAnalyzer->Get("HiTree");
/*tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");
//tree->AddFriend("genpartTree=genparticles/hi");
tree->AddFriend("genparticle=HiGenParticleAna/hi");
//tree->AddFriend("pptracks/trackTree");
tree->AddFriend("tracktree=ppTrack/trackTree");
*/
TTree *skimTree=(TTree*)skimanalysis->Get("HltTree");
TTree *hltTree=(TTree*)hltanalysis->Get("HltTree");
TTree *genTree=(TTree*)genparticles->Get("hi");
TTree *ppTree=(TTree*)pptracks->Get("trackTree");
/*double HFbin[100];
for(int i=0;i<100;i++)
HFbin[i]=i;
int nHFbin=sizeof(HFbin)/sizeof(double)-1;*/

const int Maxmult=10000;const int MaxTrk=1000;
Float_t hiHF, vz, hiHFplus,eta[Maxmult];
Float_t hiHFminusEta4, hiHFplusEta4;
Int_t nTrk;
Float_t trkDz1[MaxTrk], trkDxy1[MaxTrk], trkDzError1[MaxTrk], trkDxyError1[MaxTrk],trkPtError[MaxTrk], trkPt[MaxTrk];
Bool_t highPurity[MaxTrk];
Long_t Ev,Nevent, Nevent_counter, Nleft;
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter,phfPosFilter1, phfNegFilter1,pBeamScrapingFilter;
Int_t Nskim=0;	double weight;	int mult, nDS=0, nPlus35, nMinus35;	int sta[Maxmult],pdg[Maxmult];
Float_t phi[Maxmult],pt[Maxmult];;
Nevent_counter=0; Nleft=0;
if(HiTree->GetEntries()==skimTree->GetEntries())
Nevent=HiTree->GetEntries();

HiTree->SetBranchAddress("hiHF",&hiHF);//HF energy
HiTree->SetBranchAddress("hiHFplus",&hiHFplus);//HF energy positive
HiTree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);//HF energy positive 4 to 5
HiTree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);//HF energy negative -5 to -4
HiTree->SetBranchAddress("vz",&vz);
hltTree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
skimTree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
skimTree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
//skimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
skimTree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
skimTree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
skimTree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
skimTree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
genTree->SetBranchAddress("eta",&eta);
genTree->SetBranchAddress("sta",&sta);
genTree->SetBranchAddress("pt",&pt);
genTree->SetBranchAddress("pdg",&pdg);
genTree->SetBranchAddress("phi",&phi);
genTree->SetBranchAddress("mult",&mult);
ppTree->SetBranchAddress("nTrk",&nTrk);
ppTree->SetBranchAddress("highPurity",&highPurity);
ppTree->SetBranchAddress("trkDz1",&trkDz1);
ppTree->SetBranchAddress("trkDxy1",&trkDxy1);
ppTree->SetBranchAddress("trkDzError1",&trkDzError1);
ppTree->SetBranchAddress("trkDxyError1",&trkDxyError1);
ppTree->SetBranchAddress("trkPtError",&trkPtError);
ppTree->SetBranchAddress("trkPt",&trkPt);

TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);
TLorentzVector * v = new TLorentzVector();

for(Ev=0;Ev<Nevent;Ev++){

HiTree->GetEntry(Ev);
hltTree->GetEntry(Ev);
skimTree->GetEntry(Ev);
genTree->GetEntry(Ev);
ppTree->GetEntry(Ev);
weight=1.;
if(Isweight)
weight*=fVz->Eval(vz);
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;
//cout<<"weight="<<weight<<endl;

hHFEnergy4->Fill(hiHFplusEta4+hiHFminusEta4,weight);
hHFEnergyPlus4->Fill(hiHFplusEta4,weight);
hHFEnergyPlus->Fill(hiHFplus,weight);
Nevent_counter++;

 int trackMult = 0;
      for(int j=0;j<nTrk;j++) {
         if(!((highPurity[j])
             && (fabs(trkDz1[j]/trkDzError1[j])<3)
             && (fabs(trkDxy1[j]/trkDxyError1[j])<3)
             && (trkPtError[j]/trkPt[j]<0.1)
            ))
            continue;
         trackMult++;
      }
hNtrk->Fill(trackMult,weight);

//if(trackMult==0)	continue;

nPlus35=0; nMinus35=0;

for(int j=0; j<mult;j++){
float m = -1.;
      if(fabs(pdg[j]) == 111) m = 0.13498; //pion 0
      else if(fabs(pdg[j]) == 211) m = 0.13957; //pion +-
      else if(fabs(pdg[j]) == 321) m = 0.49368; //kaon +-
      else if(fabs(pdg[j]) == 130) m = 0.49761; //kaon 0 long
      else if(fabs(pdg[j]) == 2212) m = 0.93827; //proton
      else if(fabs(pdg[j]) == 2112) m = 0.93957; //neutron
      else if(fabs(pdg[j]) == 22) m = 0.; //photon
      else if(fabs(pdg[j]) == 11) m = 0.000511; //electron
      else if(fabs(pdg[j]) == 13) m = 0.10566; //muon
	      else continue;

v->SetPtEtaPhiM(pt[j],eta[j],phi[j],m);
if(eta[j]>3 &&eta[j]<5 && v->E()>3. && (sta[j]==1 || fabs(pdg[j]) == 111)) nPlus35++;
if(eta[j]>-5 && eta[j]<-3 && v->E()>3. && (sta[j]==1 || fabs(pdg[j]) == 111)) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<=15){
hHFEnergy4_DS->Fill(hiHFplusEta4+hiHFminusEta4,weight);
hHFEnergyPlus4_DS->Fill(hiHFplusEta4,weight);
hHFEnergyPlus_DS->Fill(hiHFplus,weight);
hNtrk_DS->Fill(trackMult,weight);
nDS++;
}

//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pHBHENoiseFilter && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<=15))
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))
{Nskim++; continue;}

hHFEnergy4_tr->Fill(hiHFplusEta4+hiHFminusEta4,weight);
hHFEnergyPlus4_tr->Fill(hiHFplusEta4,weight);
hHFEnergyPlus_tr->Fill(hiHFplus,weight);
hNtrk_tr->Fill(trackMult,weight);
Nleft++;

}
cout<<Nevent_counter<<endl;
cout<<nDS<<endl;
cout<<Nleft<<endl;


fout->cd();
hHFEnergy4->Write();
hHFEnergy4_tr->Write();
hHFEnergy4_DS->Write();

hHFEnergyPlus4->Write();
hHFEnergyPlus4_tr->Write();
hHFEnergyPlus4_DS->Write();

hHFEnergyPlus->Write();
hHFEnergyPlus_tr->Write();
hHFEnergyPlus_DS->Write();

hNtrk->Write();
hNtrk_tr->Write();
hNtrk_DS->Write();

fin->Close();
fout->Close();

gROOT->Reset();

}

void Trig_eff_DS(){
bool Isweight=kFALSE;

//Char_t *Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root";
//Char_t *Epos_file = "/afs/cern.ch/work/k/kjung/public/EposForest_CandyCane_Merged1.root";
TString Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_EposLHC_official_500k.root";

//TString Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pPb_Hijing_MB02/prod12/HiForest_v68_merged01/H02_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod12/HiForest_v68_merged02/H01_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod09/HiForest_v55_merged02/merged_forest_0.root";
TString Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_Hijing_official_500k.root";
if(!Isweight)
TString dir="output/noweight";
else
TString dir="output/weight";

TString Epos_outfile = Form("%s/pPbHist_Epos_DS.root",dir.Data());
TString Hijing_outfile = Form("%s/pPbHist_Hijing_DS.root",dir.Data());

//Calc(Epos_file,Epos_outfile,Isweight);
Calc(Hijing_file,Hijing_outfile,Isweight);

}
