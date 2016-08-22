void Calc(TString Generatorfile,TString outfile){

TFile *fin = TFile::Open(Generatorfile);
TFile *fout = TFile::Open(outfile,"RECREATE"); 

TH1F *hHFEnergy4 = new TH1F("hHFEnergy4","HFEnergy Distribution before Trigger selection",30000,0,300);	hHFEnergy4->Sumw2();
TH1F *hHFEnergy4_tr = new TH1F("hHFEnergy4_tr","HFEnergy Distribution after Trigger selection",30000,0,300); hHFEnergy4_tr->Sumw2();
TH1F *hHFEnergy4_DS = new TH1F("hHFEnergy4_DS","HFEnergy Distribution after Doudle Side event selection",30000,0,300); hHFEnergy4_DS->Sumw2();
TH1F *hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","HFEnergy Distribution before Trigger selection",30000,0,300);	hHFEnergyPlus4->Sumw2();
TH1F *hHFEnergyPlus4_tr = new TH1F("hHFEnergyPlus4_tr","HFEnergy Distribution after Trigger selection",30000,0,300); 	hHFEnergyPlus4_tr->Sumw2();
TH1F *hHFEnergyPlus4_DS = new TH1F("hHFEnergyPlus4_DS","HFEnergy Distribution after Double Side event selection",30000,0,300); 	hHFEnergyPlus4_DS->Sumw2();
TH1F *hHFEnergyPlus = new TH1F("hHFEnergyPlus","HFEnergy Distribution before Trigger selection",30000,0,300);	hHFEnergyPlus->Sumw2();
TH1F *hHFEnergyPlus_tr = new TH1F("hHFEnergyPlus_tr","HFEnergy Distribution after Trigger selection",30000,0,300); 	hHFEnergyPlus_tr->Sumw2();
TH1F *hHFEnergyPlus_DS = new TH1F("hHFEnergyPlus_DS","HFEnergy Distribution after Double Side event selection",30000,0,300); 	hHFEnergyPlus_DS->Sumw2();

TTree *tree;
fin->cd();
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");
//tree->AddFriend("genpartTree=genparticles/hi");

const int Maxmult=10000;
Float_t hiHF, vz, hiHFplus,eta[Maxmult];
Float_t hiHFminusEta4, hiHFplusEta4;
Long_t Ev,Nevent, Nevent_counter, Nleft;
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter,phfPosFilter1, phfNegFilter1,pBeamScrapingFilter,pPApileupVertexFilterCutGplus;
Int_t Nskim=0;	double weight;	int mult, nDS=0, nPlus35, nMinus35;	int sta[Maxmult];
Nevent_counter=0; Nleft=0;
Nevent=tree->GetEntries();

tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplus",&hiHFplus);//HF energy positive
tree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pPApileupVertexFilterCutGplus",&pPApileupVertexFilterCutGplus);
//tree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
tree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
//tree->SetBranchAddress("eta",&eta);
//tree->SetBranchAddress("sta",&sta);
//tree->SetBranchAddress("mult",&mult);

TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);

for(Ev=0;Ev<Nevent;Ev++){

tree->GetEntry(Ev);
weight=1.;
//weight*=fVz->Eval(vz);
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;
//cout<<"weight="<<weight<<endl;

hHFEnergy4->Fill(hiHFplusEta4+hiHFminusEta4);
hHFEnergyPlus4->Fill(hiHFplusEta4);
hHFEnergyPlus->Fill(hiHFplus);
Nevent_counter++;
/*
nPlus35=0; nMinus35=0;
for(int j=0; j<mult;j++){
if(eta[j]>3&&eta[j]<5&&sta[j]==1) nPlus35++;
if(eta[j]>-5&&eta[j]<-3&&sta[j]==1) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<15){
hHFEnergy4_DS->Fill(hiHFplusEta4+hiHFminusEta4);
hHFEnergyPlus4_DS->Fill(hiHFplusEta4);
hHFEnergyPlus_DS->Fill(hiHFplus);
nDS++;
}
*/
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pPApileupVertexFilterCutGplus && pprimaryvertexFilter && TMath::Abs(vz)<15))
{Nskim++; continue;}

hHFEnergy4_tr->Fill(hiHFplusEta4+hiHFminusEta4);
hHFEnergyPlus4_tr->Fill(hiHFplusEta4);
hHFEnergyPlus_tr->Fill(hiHFplus);
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

fin->Close();
fout->Close();

gROOT->Reset();

}

void Trig_eff_DS(){

//Char_t *Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root";
//Char_t *Epos_file = "/afs/cern.ch/work/k/kjung/public/EposForest_CandyCane_Merged1.root";
//TString Epos_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";

//Char_t *Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pPb_Hijing_MB02/prod12/HiForest_v68_merged01/H02_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod12/HiForest_v68_merged02/H01_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod09/HiForest_v55_merged02/merged_forest_0.root";

TString Hijing_file="/store/user/tuos/MinBiasTree_v3_210614.root";

TString dir=".";
TString Hijing_outfile = Form("%s/pPbHist_Hijing_withGplus_DS.root",dir.Data());

Calc(Hijing_file,Hijing_outfile);

}
