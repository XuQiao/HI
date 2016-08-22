void Calc(TString Generatorfile,TString outfile){

TFile *fin = TFile::Open(Generatorfile);
TFile *fout = TFile::Open(outfile,"RECREATE"); 

TH1F *hNtrk = new TH1F("hNtrk","",100,0,100);	hNtrk->Sumw2();
TH1F *hNtrk_tr = new TH1F("hNtrk_tr","",100,0,100); hNtrk_tr->Sumw2();
TH1F *hNtrk_DS = new TH1F("hNtrk_DS","",100,0,100); hNtrk_DS->Sumw2();

TTree *tree;
fin->cd();
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");
tree->AddFriend("genpartTree=genparticles/hi");
tree->AddFriend("pptracks/trackTree");

const int Maxmult=10000; const int MaxTrk=1000;
Float_t hiHF, vz, hiHFplus,eta[Maxmult];
Float_t hiHFminusEta4, hiHFplusEta4;
Long_t Ev,Nevent, Nevent_counter, Nleft;
Int_t nTrk;
Float_t trkDz1[MaxTrk], trkDxy1[MaxTrk], trkDzError1[MaxTrk], trkDxyError1[MaxTrk],trkPtError[MaxTrk], trkPt[MaxTrk];
Bool_t highPurity[MaxTrk];
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter,phfPosFilter1, phfNegFilter1,pBeamScrapingFilter;
Int_t Nskim=0;	double weight;	int mult, nDS=0, nPlus35, nMinus35;	int sta[Maxmult];
Nevent_counter=0; Nleft=0;
Nevent=tree->GetEntries();

tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
//tree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
tree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
tree->SetBranchAddress("eta",&eta);
tree->SetBranchAddress("sta",&sta);
tree->SetBranchAddress("mult",&mult);
tree->SetBranchAddress("nTrk",&nTrk);
tree->SetBranchAddress("highPurity",&highPurity);
tree->SetBranchAddress("trkDz1",&trkDz1);
tree->SetBranchAddress("trkDxy1",&trkDxy1);
tree->SetBranchAddress("trkDzError1",&trkDzError1);
tree->SetBranchAddress("trkDxyError1",&trkDxyError1);
tree->SetBranchAddress("trkPtError",&trkPtError);
tree->SetBranchAddress("trkPt",&trkPt);


TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);

for(Ev=0;Ev<Nevent;Ev++){

tree->GetEntry(Ev);
weight=1.;
weight*=fVz->Eval(vz);
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;
//cout<<"weight="<<weight<<endl;
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
//if(trackMult==0) continue;
hNtrk->Fill(trackMult,weight);
Nevent_counter++;

nPlus35=0; nMinus35=0;
for(int j=0; j<mult;j++){
if(eta[j]>3&&eta[j]<5&&sta[j]==1) nPlus35++;
if(eta[j]>-5&&eta[j]<-3&&sta[j]==1) nMinus35++;
}
if(nPlus35>=1 && nMinus35>=1 && TMath::Abs(vz)<15){
hNtrk_DS->Fill(trackMult,weight);
nDS++;
}
if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
//if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))
{Nskim++; continue;}

hNtrk_tr->Fill(trackMult,weight);
Nleft++;
}
cout<<Nevent_counter<<endl;
cout<<nDS<<endl;
cout<<Nleft<<endl;
fout->cd();
hNtrk->Write();
hNtrk_tr->Write();
hNtrk_DS->Write();

fin->Close();
fout->Close();

gROOT->Reset();

}

void Trig_eff_Ntrk_DS(){

//Char_t *Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root";
//Char_t *Epos_file = "/afs/cern.ch/work/k/kjung/public/EposForest_CandyCane_Merged1.root";
TString Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_EposLHC_official_500k.root";

//Char_t *Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pPb_Hijing_MB02/prod12/HiForest_v68_merged01/H02_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod12/HiForest_v68_merged02/H01_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod09/HiForest_v55_merged02/merged_forest_0.root";
TString Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/pPb_MC_MinBiasTree/pPb_MC_MinBiasTree_v6_Hijing_official_500k.root";


TString dir="output/weight";
TString Epos_outfile = Form("%s/pPbHist_Epos_Ntrk_DS.root",dir.Data());
TString Hijing_outfile = Form("%s/pPbHist_Hijing_Ntrk_DS.root",dir.Data());

Calc(Epos_file,Epos_outfile);
Calc(Hijing_file,Hijing_outfile);

}
