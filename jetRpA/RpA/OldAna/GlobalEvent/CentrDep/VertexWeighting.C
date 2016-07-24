#include "/home/xuq7/CMSSW_6_2_3_patch1/src/jetRpA/RpA/Quality/root_setting.h"
void VertexWeighting(){
TString Gen="Hijing";
//TString Gen="Epos";
if(Gen=="Hijing"){
TString Hijing_file="/store/user/tuos/pPb_MC_MinBiasTree_v4_Hijing_180k.root";
TFile *fin=TFile::Open(Hijing_file);
TFile *fout=TFile::Open("output/Hijing_vz.root","RECREATE");
}
else if(Gen=="Epos"){
TString Epos_file = "/store/user/tuos/pPb_MC_MinBiasTree_v4_Epos_189k.root";
TFile *fin=TFile::Open(Epos_file);
TFile *fout=TFile::Open("output/Epos_vz.root","RECREATE");
}

TH1F *hVz_MC=new TH1F("hVz_MC","Vz distribution",400, -20., 20.);	hVz_MC->Sumw2();
TH2F *hVzHF4=new TH2F("hVzHF4","Vz HF distribution",400,-20.,20, 500,0.,100.);	hVzHF4->Sumw2();
TH2F *hVzHFPlus4=new TH2F("hVzHFPlus4","Vz HF distribution",400,-20.,20, 500,0.,100.);	hVzHFPlus4->Sumw2();
TH2F *hVzHFPlus=new TH2F("hVzHFPlus","Vz HF distribution",400,-20.,20, 500,0.,100.);	hVzHFPlus->Sumw2();

TTree *tree;
fin->cd();
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
tree->AddFriend("skimTree=skimanalysis/HltTree");
tree->AddFriend("hltTree=hltanalysis/HltTree");

Float_t hiHF, vz, hiHFplus;
Float_t hiHFminusEta4, hiHFplusEta4;
Long_t Ev,Nevent, Nevent_counter, Nleft;
Int_t n, HLT_PAZeroBiasPixel_SingleTrack_v1, pprimaryvertexFilter, pVertexFilterCutGplus, pPAcollisionEventSelectionPA, pHBHENoiseFilter;
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
if(Ev%5000==0)
cout<<"Have run "<<Ev<<" Events"<<endl;

tree->GetEntry(Ev);
if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pPAcollisionEventSelectionPA && pprimaryvertexFilter && TMath::Abs(vz)<15)) continue;
hVz_MC->Fill(vz);
hVzHF4->Fill(vz,hiHFplusEta4+hiHFminusEta4);
hVzHFPlus4->Fill(vz,hiHFplusEta4);
hVzHFPlus->Fill(vz,hiHFplus);
}

fout->cd();
hVz_MC->Write();
hVzHF4->Write();
hVzHFPlus4->Write();
hVzHFPlus->Write();

fout->Close();
}
