{
const Int_t nTrkMax = 50000;

Char_t *HiForestFile = "/store/user/tuos/pPb_MinBiasTree_v5_211256_json.root";
Char_t *outfile="./pPbHistNtrackcut.root";

TFile *foutput = TFile::Open(outfile,"update");//outfile;
TFile *hHiForestFile = TFile::Open(HiForestFile);

TTree *tree;
tree=(TTree*)hiEvtAnalyzer->Get("HiTree");
Char_t *tr2 = "pptracks/trackTree";
Char_t *tr3 = "skimanalysis/HltTree";
Char_t *tr4 = "hltanalysis/HltTree";
Char_t *tr5 = "hfrechits/hfTree";

tree.AddFriend(tr2);
tree.AddFriend(tr3);
tree.AddFriend(tr4);
tree.AddFriend(tr5);

Float_t hiHF, vz;
Float_t hiHFminus4, hiHFplus4;
Bool_t highPurity[nTrkMax];
Float_t trkPt[nTrkMax], trkPtError[nTrkMax], trkEta[nTrkMax], trkDxy1[nTrkMax], trkDxyError1[nTrkMax], trkDz1[nTrkMax], trkDzError1[nTrkMax];
Long_t Nevent;
Int_t Ntrack, n, HLT_PAZeroBiasPixel_SingleTrack_v1, phltPixelClusterShapeFilter, phfPosFilter1, phfNegFilter1, pprimaryvertexFilter, pBeamScrapingFilter, pVertexFilterCutGplus;
Int_t nTrk, iTrk, nTrk_lc;

Int_t Nskim=0;

Nevent=tree->GetEntries();

tree->SetBranchAddress("hiNtracks",&Ntrack);//number of tracks
tree->SetBranchAddress("hiHF",&hiHF);//HF energy
tree->SetBranchAddress("hiHFplusEta4",&hiHFplus4);//HF energy positive 4 to 5
tree->SetBranchAddress("hiHFminusEta4",&hiHFminus4);//HF energy negative -5 to -4
tree->SetBranchAddress("vz",&vz);
tree->SetBranchAddress("n",&n);//HF hits
tree->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
tree->SetBranchAddress("phltPixelClusterShapeFilter",&phltPixelClusterShapeFilter);
tree->SetBranchAddress("phfPosFilter1",&phfPosFilter1);
tree->SetBranchAddress("phfNegFilter1",&phfNegFilter1);
tree->SetBranchAddress("pprimaryvertexFilter",&pprimaryvertexFilter);
tree->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
tree->SetBranchAddress("nTrk",&nTrk);
tree->SetBranchAddress("trkPt",&trkPt);
tree->SetBranchAddress("trkPtError",&trkPtError);
tree->SetBranchAddress("trkEta",&trkEta);
tree->SetBranchAddress("highPurity",&highPurity);
tree->SetBranchAddress("trkDxy1",&trkDxy1);
tree->SetBranchAddress("trkDxyError1",&trkDxyError1);
tree->SetBranchAddress("trkDz1",&trkDz1);
tree->SetBranchAddress("trkDzError1",&trkDzError1);

//-------------------------------------------------------------
//Event Selection

TH1D* hNtrack = new TH1D("hNtrack","hNtrack;track hits;# of events",300,0,300);
TH1F* hHFEnergy = new TH1F("hHFEnergy","hHF Deposit Energy;Energy;# of events",300,0,300);
TH1F* hHFEnergy4 = new TH1F("hHFEnergy4","hHF Deposit Energy;Energy;# of events",300,0,300);
TH1F* hHFEnergyPlus4 = new TH1F("hHFEnergyPlus4","hHF Deposit Energy eta range;Energy;# of events",300,0,100);
TH1F* hHFEnergyMinus4 = new TH1F("hHFEnergyMinus4","hHF Deposit Energy eta range;Energy;# of events",300,0,300);
TH1D* hHFHit = new TH1D("hHFHit","hHFHit;HF Hits;# of events",800,0,800);
TH1D* hNtrack_lc = new TH1D("hNtrack_lc","hNtrack loose cut;track hits;# of events",300,0,300);

Long_t Ev;

for(Ev=0; Ev<Nevent; Ev++){
if(Ev%5000==0) cout<<"Ev = "<<Ev<<endl;
tree->GetEntry(Ev);

if(!(phltPixelClusterShapeFilter && HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))

{
Nskim++;
continue;
}

else{

nTrk_lc=0;
for(iTrk=0;iTrk<nTrk;iTrk++){
if(highPurity[iTrk] && trkPt[iTrk] > 0.1 && fabs(trkEta[iTrk])<2.4 && trkPtError[iTrk]/trkPt[iTrk]<0.1 && fabs(trkDz1[iTrk]/trkDzError1[iTrk])<3.0 && fabs(trkDxy1[iTrk]/trkDxyError1[iTrk])<3.0)
nTrk_lc++;
}

hNtrack->Fill(Ntrack);
hHFEnergy->Fill(hiHF);
hHFEnergy4->Fill(hiHFplus4+hiHFminus4);
hHFEnergyPlus4->Fill(hiHFplus4);
hHFEnergyMinus4->Fill(hiHFminus4);
hHFHit->Fill(n);
hNtrack_lc->Fill(nTrk_lc);

}

}

cout<<"Nevent="<<Nevent<<endl<<"Nskim="<<Nskim<<endl<<"Percentage="<<(double)Nskim/Nevent<<endl;

foutput->cd();
hNtrack->Write("hNtrack",TObject::kOverwrite);
hHFEnergy->Write("hHFEnergy",TObject::kOverwrite);
hHFEnergy4->Write("hHFEnergy4",TObject::kOverwrite);
hHFEnergyPlus4->Write("hHFEnergyPlus4",TObject::kOverwrite);
hHFEnergyMinus4->Write("hHFEnergyMinus4",TObject::kOverwrite);
hHFHit->Write("hHFHit",TObject::kOverwrite);
hNtrack_lc->Write("hNtrack_lc",TObject::kOverwrite);

foutput->Close();
}
