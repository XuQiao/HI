{


Char_t *HiForestFile1 = "/store/user/tuos/MinBiasTree_v4_211607.root";
Char_t *HiForestFile2 = "/store/user/tuos/pPb_MinBiasTree_v5_211256_json.root";
Char_t *outfile="./pPbHistEta4.root";
Char_t *outfile_origin="./pPbHist_origin.root";

TFile *foutput = TFile::Open(outfile,"update");//outfile;

Char_t *tr1 = "hiEvtAnalyzer/HiTree";
Char_t *tr2 = "pptracks/trackTree";
Char_t *tr3 = "skimanalysis/HltTree";
Char_t *tr4 = "hltanalysis/HltTree";
Char_t *tr5 = "hfrechits/hfTree";

tree = new TChain(tr1);
t2 = new TChain(tr2);
t3 = new TChain(tr3);
t4 = new TChain(tr4);
t5 = new TChain(tr5);

tree.Add(HiForestFile1);
tree.Add(HiForestFile2);
t2.Add(HiForestFile1);
t2.Add(HiForestFile2);
t3.Add(HiForestFile1);
t3.Add(HiForestFile2);
t4.Add(HiForestFile1);
t4.Add(HiForestFile2);
t5.Add(HiForestFile1);
t5.Add(HiForestFile2);

tree.AddFriend(t2,tr2);
tree.AddFriend(t3,tr3);
tree.AddFriend(t4,tr4);
tree.AddFriend(t5,tr5);

Float_t hiHF, vz;
Float_t hiHFminus4, hiHFplus4;
Long_t Nevent;
Int_t Ntrack, n, HLT_PAZeroBiasPixel_SingleTrack_v1, phltPixelClusterShapeFilter, phfPosFilter1, phfNegFilter1, pprimaryvertexFilter, pBeamScrapingFilter, pVertexFilterCutGplus;
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


//-------------------------------------------------------------
//Event Selection

TH1D* hNtrack = new TH1D("hNtrack","hNtrack;track hits;# of events",300,0,300);	hNtrack->Sumw2();
TH1F* hHFEnergy = new TH1F("hHFEnergy","hHF Deposit Energy;Energy;# of events",300,0,300);	hHFEnergy->Sumw2();
TH1F* hHFEnergy4 = new TH1F("hHFEnergy4","hHF Deposit Energy eta range;Energy;# of events",300,0,300);	hHFEnergy4->Sumw2();
TH1D* hHFHit = new TH1D("hHFHit","hHFHit;HF Hits;# of events",800,0,800);	hHFHit->Sumw2();

Long_t Ev;
cout<<Nevent<<endl;
for(Ev=0; Ev<Nevent; Ev++){
if(Ev%5000==0) cout<<"Ev = "<<Ev<<endl;
tree->GetEntry(Ev);

if(!(HLT_PAZeroBiasPixel_SingleTrack_v1 && pVertexFilterCutGplus && pBeamScrapingFilter && phfPosFilter1 && phfNegFilter1 && pprimaryvertexFilter && TMath::Abs(vz)<15))
{
Nskim++;
continue;
}

else{
hNtrack->Fill(Ntrack);
hHFEnergy->Fill(hiHF);
hHFEnergy4->Fill(hiHFplus4+hiHFminus4);
hHFHit->Fill(n);
}

}

cout<<"Nevent="<<Nevent<<endl<<"Nskim="<<Nskim<<endl<<"Percentage="<<(double)Nskim/Nevent<<endl;

foutput->cd();
hNtrack->Write("hNtrack",TObject::kOverwrite);
hHFEnergy->Write("hHFEnergy",TObject::kOverwrite);
hHFEnergy4->Write("hHFEnergy4",TObject::kOverwrite);
hHFHit->Write("hHFHit",TObject::kOverwrite);

foutput->Close();
}
