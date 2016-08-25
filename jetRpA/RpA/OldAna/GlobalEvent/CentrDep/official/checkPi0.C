void Calc(TString Generatorfile,bool Isweight){

TFile *fin = TFile::Open(Generatorfile);
//TFile *fout = TFile::Open(outfile,"RECREATE"); 

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
int nPi0=0;
int nPi=0;
int nK=0;
int nK0L=0;
int nPro=0;
int nNeu=0;
int nPho=0;
int ne=0;
int nmu=0;
int nPi0raw=0;
int nPi0Minus99=0;
int nPi01=0;
int nPi02=0;
int nPi03=0;
int nPiPraw=0;
int nPiNraw=0;
int nKPraw=0;
int nKNraw=0;
int nK0Lraw=0;
int nProraw=0;
int nNeuraw=0;
int nPhoraw=0;
int neraw=0;
int nmuraw=0;
int ntau=0;
int nSum=0;
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
double tau[Maxmult];
Nevent_counter=0; Nleft=0;
if(HiTree->GetEntries()==skimTree->GetEntries())
Nevent=HiTree->GetEntries();
nPi0;
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
//genTree->SetBranchAddress("tau",&tau);
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

for(Ev=0;Ev<10000;Ev++){

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
for(int j=0; j<mult;j++){
nSum++;
if(pdg[j]==111 && sta[j]==-99) nPi0Minus99++;
if(pdg[j]==111 && sta[j]==1) nPi01++;
if(pdg[j]==111 && sta[j]==2) nPi02++;
if(pdg[j]==111 && sta[j]==3) nPi03++;
if(pdg[j]==111) nPi0raw++;
if(pdg[j]==211) nPiPraw++;
if(pdg[j]==-211) nPiNraw++;
if(pdg[j]==321) nKPraw++;
if(pdg[j]==-321) nKNraw++;
if(pdg[j]==130) nK0Lraw++;
if(pdg[j]==2212) nProraw++;
if(pdg[j]==2112) nNeuraw++;
if(pdg[j]==22) nPhoraw++;
if(pdg[j]==11) neraw++; 
if(pdg[j]==13) nmuraw++;
}
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

//if(trackMult==0)	continue;

nPlus35=0; nMinus35=0;

for(int j=0; j<mult;j++){
float m = -1.;
      if(fabs(pdg[j]) == 111) {m = 0.13498; } //pion 0 
      else if(fabs(pdg[j]) == 211) {m = 0.13957;}//pion +-
      else if(fabs(pdg[j]) == 321) {m = 0.49368;}//kaon +-
      else if(fabs(pdg[j]) == 130) {m = 0.49761;}//kaon 0 long
      else if(fabs(pdg[j]) == 2212) {m = 0.93827; }//proton
      else if(fabs(pdg[j]) == 2112) {m = 0.93957; }//neutron
      else if(fabs(pdg[j]) == 22) {m = 0.; }//photon
      else if(fabs(pdg[j]) == 11) {m = 0.000511;}//electron
      else if(fabs(pdg[j]) == 13) {m = 0.10566; }//muon
	      else continue;

v->SetPtEtaPhiM(pt[j],eta[j],phi[j],m);
if(eta[j]>3 &&eta[j]<5 && v->E()>3. && (sta[j]==1 || fabs(pdg[j])==111)) {nPlus35++;

if(fabs(pdg[j]) == 111) {nPi0++;} //pion 0 
      else if(fabs(pdg[j]) == 211) {nPi++;}//pion +-
      else if(fabs(pdg[j]) == 321) {nK++;}//kaon +-
      else if(fabs(pdg[j]) == 130) {nK0L++;}//kaon 0 long
      else if(fabs(pdg[j]) == 2212) {nPro++;}//proton
      else if(fabs(pdg[j]) == 2112) {nNeu++;}//neutron
      else if(fabs(pdg[j]) == 22) {nPho++;}//photon
      else if(fabs(pdg[j]) == 11) {ne++;}//electron
      else if(fabs(pdg[j]) == 13) {nmu++;}//muon
}
if(eta[j]>-5 && eta[j]<-3 && v->E()>3. && (sta[j]==1 || fabs(pdg[j])==111)){ nMinus35++;
if(fabs(pdg[j]) == 111) {nPi0++;} //pion 0 
      else if(fabs(pdg[j]) == 211) {nPi++;}//pion +-
      else if(fabs(pdg[j]) == 321) {nK++;}//kaon +-
      else if(fabs(pdg[j]) == 130) {nK0L++;}//kaon 0 long
      else if(fabs(pdg[j]) == 2212) {nPro++;}//proton
      else if(fabs(pdg[j]) == 2112) {nNeu++;}//neutron
      else if(fabs(pdg[j]) == 22) {nPho++;}//photon
      else if(fabs(pdg[j]) == 11) {ne++;}//electron
      else if(fabs(pdg[j]) == 13) {nmu++;}//muon
}
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<=15){
nDS++;
}


if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<=15))
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))
{Nskim++; continue;}
Nleft++;
}
/*
cout<<nPi0<<endl;
cout<<nPi<<endl;
cout<<nK<<endl;
cout<<nK0L<<endl;
cout<<nPro<<endl;
cout<<nNeu<<endl;
cout<<nPho<<endl;
cout<<ne<<endl;
cout<<nmu<<endl;*/
cout<<"ntau:"<<ntau<<endl;
cout<<"nSum:"<<nSum<<endl;
cout<<"nPi0 sta=-99:"<<nPi0Minus99<<endl;
cout<<"nPi0 sta=1:"<<nPi01<<endl;
cout<<"nPi0 sta=2:"<<nPi02<<endl;
cout<<"nPi0 sta=3:"<<nPi03<<endl;
cout<<"pi0:"<<nPi0raw<<endl;
cout<<"pi+:"<<nPiPraw<<endl;
cout<<"pi-:"<<nPiNraw<<endl;
cout<<"K+:"<<nKPraw<<endl;
cout<<"K-:"<<nKNraw<<endl;
cout<<"K0L:"<<nK0Lraw<<endl;
cout<<"Proton:"<<nProraw<<endl;
cout<<"Neutral:"<<nNeuraw<<endl;
cout<<"Photon:"<<nPhoraw<<endl;
cout<<"electron:"<<neraw<<endl;
cout<<"muon:"<<nmuraw<<endl;
cout<<Nevent_counter<<endl;
cout<<nDS<<endl;
cout<<nleft<<endl;

fin->Close();

}

void checkPi0(){
bool Isweight=kFALSE;

//Char_t *Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root";
//Char_t *Epos_file = "/afs/cern.ch/work/k/kjung/public/EposForest_CandyCane_Merged1.root";
TString Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_EposLHC_official_500k.root";

//TString Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pPb_Hijing_MB02/prod12/HiForest_v68_merged01/H02_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod12/HiForest_v68_merged02/H01_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod09/HiForest_v55_merged02/merged_forest_0.root";
//TString Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_Hijing_official_500k.root";
TString Hijing_file="/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";

//TString Epos_outfile = Form("%s/pPbHist_Epos_DS.root",dir.Data());
//TString Hijing_outfile = Form("%s/pPbHist_Hijing_DS.root",dir.Data());

Calc(Epos_file,Isweight);
//Calc(Hijing_file,Isweight);

}
