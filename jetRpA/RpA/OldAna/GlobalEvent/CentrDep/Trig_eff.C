void Calc(TString Generatorfile,TString outfile,TString cutName){

TFile *fin = TFile::Open(Generatorfile);
TFile *fout = TFile::Open(outfile,"RECREATE"); 

TH1F *hHFEnergy4 = new TH1F("hHFEnergy","HFEnergy Distribution before Trigger selection",10000,0,100);	hHFEnergy4->Sumw2();
TH1F *hHFEnergy4_tr = new TH1F("hHFEnergy_tr","HFEnergy Distribution after Trigger selection",10000,0,100); hHFEnergy4_tr->Sumw2();
TH1F *hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","HFEnergy Distribution before Trigger selection",10000,0,100);	hHFEnergyPlus4->Sumw2();
TH1F *hHFEnergyPlus4_tr = new TH1F("hHFEnergyPlus4_tr","HFEnergy Distribution after Trigger selection",10000,0,100); 	hHFEnergyPlus4_tr->Sumw2();
TH1F *hHFEnergyPlus = new TH1F("hHFEnergyPlus","HFEnergy Distribution before Trigger selection",10000,0,100);	hHFEnergyPlus->Sumw2();
TH1F *hHFEnergyPlus_tr = new TH1F("hHFEnergyPlus_tr","HFEnergy Distribution after Trigger selection",1000,0,100); 	hHFEnergyPlus_tr->Sumw2();

TTree *tree;
fin->cd();
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");

Float_t hiHF, vz, hiHFplus;
Float_t hiHFminusEta4, hiHFplusEta4;
Long_t Ev,Nevent, Nevent_counter, Nleft;
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter;
Int_t Nskim=0;	double weight;
Nevent_counter=0; Nleft=0;
Nevent=tree->GetEntries();

tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplus",&hiHFplus);//HF energy positive
tree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
//tree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
tree->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);

TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
fVz->SetParameters(1.60182e+00,1.08425e-03,-1.29156e-02,-7.24899e-06,2.80750e-05);

for(Ev=0;Ev<Nevent;Ev++){

tree->GetEntry(Ev);
weight=1.;
weight*=fVz->Eval(vz);
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<'\t';
cout<<"weight="<<weight<<endl;

hHFEnergy4->Fill(hiHFplusEta4+hiHFminusEta4,weight);
hHFEnergyPlus4->Fill(hiHFplusEta4,weight);
hHFEnergyPlus->Fill(hiHFplus,weight);
Nevent_counter++;
//(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15){
if(cutName.Contains("HLT") && !(HLT_PAZeroBiasPixel_SingleTrack_v1))
{Nskim++; continue;}
if(cutName.Contains("CutGplus") && !(pVertexFilterCutGplus))
{Nskim++; continue;}
if(cutName.Contains("pPA") && !(pPAcollisionEventSelectionPA))
{Nskim++; continue;}
if(cutName.Contains("pprim") && !(pprimaryvertexFilter))
{Nskim++; continue;}
if(cutName.Contains("vz") && !(TMath::Abs(vz)<15))
{Nskim++; continue;}
if(cutName=="" && !(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15))
{Nskim++; continue;}

hHFEnergy4_tr->Fill(hiHFplusEta4+hiHFminusEta4,weight);
hHFEnergyPlus4_tr->Fill(hiHFplusEta4,weight);
hHFEnergyPlus_tr->Fill(hiHFplus,weight);
Nleft++;

}

fout->cd();
hHFEnergy4->Write("hHFEnergy4",TObject::kOverwrite);
hHFEnergy4_tr->Write("hHFEnergy4_tr",TObject::kOverwrite);

hHFEnergyPlus4->Write("hHFEnergyPlus4",TObject::kOverwrite);
hHFEnergyPlus4_tr->Write("hHFEnergyPlus4_tr",TObject::kOverwrite);

hHFEnergyPlus->Write("hHFEnergyPlus",TObject::kOverwrite);
hHFEnergyPlus_tr->Write("hHFEnergyPlus_tr",TObject::kOverwrite);

fin->Close();
fout->Close();

gROOT->Reset();

}

void Trig_eff(){

//Char_t *Epos_file = "root://eoscms//eos/cms/store/group/phys_heavyions/azsigmon/HiForest_pPb_Epos_336800.root";
//Char_t *Epos_file = "/afs/cern.ch/work/k/kjung/public/EposForest_CandyCane_Merged1.root";
TString Epos_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";

//Char_t *Hijing_file="root://eoscms//eos/cms/store/group/phys_heavyions/dgulhan/pPb_Hijing_MB02/prod12/HiForest_v68_merged01/H02_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod12/HiForest_v68_merged02/H01_merged_hiforest_v68_prod12.root";
//Char_t *Hijing_file="root://eoscms//eos/cms/store/caf/user/dgulhan/pPb_Hijing_MB01/prod09/HiForest_v55_merged02/merged_forest_0.root";
TString Hijing_file="/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";
//bool isEvo=kTRUE;	
bool isEvo=kFALSE;


if(!isEvo){
const int Ncut=6;
TString cut[Ncut]={"HLT","CutGplus","pPA","pprim","vz",""};
}
else{
const int Ncut=5;
TString cut[Ncut]={"HLT","HLT_CutGplus","HLT_CutGplus_pPA","HLT_CutGplus_pPA_pprim",""};
}
TString Epos_outfile[Ncut];
TString Hijing_outfile[Ncut];
TString dir="output/weight";
for(int i=0;i<Ncut;i++){
if(i!=Ncut-1){
Epos_outfile[i] = Form("%s/pPbHist_Epos_%s.root",dir.Data(),cut[i].Data());
Hijing_outfile[i] = Form("%s/pPbHist_Hijing_%s.root",dir.Data(),cut[i].Data());
}
else{
Epos_outfile[i] = Form("%s/pPbHist_Epos%s.root",dir.Data(),cut[i].Data());
Hijing_outfile[i] = Form("%s/pPbHist_Hijing%s.root",dir.Data(),cut[i].Data());
}

Calc(Epos_file,Epos_outfile[i],cut[i]);
Calc(Hijing_file,Hijing_outfile[i],cut[i]);

}
}
