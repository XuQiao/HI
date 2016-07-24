{

//Char_t *HiForestFile = "/store/user/tuos/pPb_MinBiasTree_v5_211256_json.root";
Char_t *HiForestFile = "root://eoscms.cern.ch//eos/cms/store/group/phys_heavyions/velicanu/forest/HIRun2015/HIMinimumBias2/Merged/HiForestPromptReco_262726.root";
Char_t *outfile="./PbPbHistPlusEta4.root";

TFile *foutput = TFile::Open(outfile,"update");//outfile;
TFile *hHiForestFile = TFile::Open(HiForestFile);

TTree *tree;
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
//Char_t *tr2 = "pptracks/trackTree";
Char_t *tr3 = "skimanalysis/HltTree";
Char_t *tr4 = "hltanalysis/HltTree";
//Char_t *tr5 = "hfrechits/hfTree";

//tree.AddFriend(tr2);
tree.AddFriend(tr3);
tree.AddFriend(tr4);
//tree.AddFriend(tr5);

Float_t hiHF, vz;
Float_t hiHFminus4, hiHFplus4;
Long_t Nevent;
Int_t Ntrack, n, pcollisionEventSelection, HLT_PAZeroBiasPixel_SingleTrack_v1, phltPixelClusterShapeFilter, phfPosFilter1, phfNegFilter1, pprimaryvertexFilter, pBeamScrapingFilter, pVertexFilterCutGplus;
Int_t Nskim=0;

Nevent=tree->GetEntries();

tree->SetBranchAddress("hiNtracks",&Ntrack);//number of tracks
tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplusEta4",&hiHFplus4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminus4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
//tree->SetBranchAddress("n",&n);//HF hits
//tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
//tree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
tree->SetBranchAddress("phfPosFilter3",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter3",&phfNegFilter1);
tree->SetBranchAddress("pprimaryVertexFilter",&pprimaryvertexFilter);
//tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
//tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);


//-------------------------------------------------------------
//Event Selection

TH1D* hNtrack = new TH1D("hNtrack","hNtrack;track hits;# of events",300,0,300);
TH1F* hHFEnergy = new TH1F("hHFEnergy","hHF Deposit Energy;Energy;# of events",5000,0,5000);
TH1F* hHFEnergy4 = new TH1F("hHFEnergy4","hHF Deposit Energy;Energy;# of events",5000,0,5000);
TH1F* hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","hHF Deposit Energy eta range;Energy;# of events",5000,0,5000);
TH1F* hHFEnergyPlus4_Rebin = new TH1F("hHFEnergyPlus4_Rebin","hHF Deposit Energy eta range;Energy;# of events",1000,0,5000);
TH1F* hHFEnergyMinus4 = new TH1F("hHFEnergyMinus4","hHF Deposit Energy eta range;Energy;# of events",5000,0,5000);
TH1D* hHFHit = new TH1D("hHFHit","hHFHit;HF Hits;# of events",1000,0,1000);

Long_t Ev;

for(Ev=0; Ev<Nevent; Ev++){
if(Ev%5000==0) cout<<"Ev = "<<Ev<<endl;
tree->GetEntry(Ev);
if(TMath::Abs(vz)>=15) continue;
//if(!(phltPixelClusterShapeFilter && HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter))	// && TMath::Abs(vz)<15
if(!(phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && pcollisionEventSelection))	// && TMath::Abs(vz)<15
{
Nskim++;
continue;
}

else{
hNtrack->Fill(Ntrack);
hHFEnergy->Fill(hiHF);
hHFEnergy4->Fill(hiHFplus4+hiHFminus4);
hHFEnergyPlus4->Fill(hiHFplus4);
hHFEnergyPlus4_Rebin->Fill(hiHFplus4);
hHFEnergyMinus4->Fill(hiHFminus4);
//hHFHit->Fill(n);
}

}

cout<<"Nevent="<<Nevent<<endl<<"Nskim="<<Nskim<<endl<<"Percentage="<<(double)Nskim/Nevent<<endl;

foutput->cd();
hNtrack->Write("hNtrack",TObject::kOverwrite);
hHFEnergy->Write("hHFEnergy",TObject::kOverwrite);
hHFEnergy4->Write("hHFEnergy4",TObject::kOverwrite);
hHFEnergyPlus4->Write("hHFEnergyPlus4",TObject::kOverwrite);
hHFEnergyPlus4_Rebin->Write("hHFEnergyPlus4_Rebin",TObject::kOverwrite);
hHFEnergyMinus4->Write("hHFEnergyMinus4",TObject::kOverwrite);
//hHFHit->Write("hHFHit",TObject::kOverwrite);

foutput->Close();
}
